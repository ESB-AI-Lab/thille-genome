#!/usr/bin/env python3
"""One stacked coverage figure per gene, driven by one BED file per assembly.

For every gene name in the BED files, draws a fixed-width window centred on the
gene for each coverage track. The panels against the B reference (Thille) are on
top and the panels against the A reference (Hass) are below, each using its own
BED coordinates. Every position is coloured by its depth relative to its own
chromosome median, with the same classes as graph_coverage.py:

    low <= 0.25M < half <= 0.75M < medium <= 1.25M < high <= 1.5M < outlier

Coverage files are read from <coverage-dir>/<track>.<contig>.coverage.gz
(3 columns: contig, 1-based position, depth; every position present). Each
(track, contig) pair is scanned once and cached as an int32 .npy plus a .json
holding the median, so later runs take seconds.

Examples:
    # Thille only
    python3 plot_gene_stacks.py -b genes_thille.bed
    # both references, joined by gene name (needs a contig map for CM names)
    python3 plot_gene_stacks.py -a genes_hass.bed -b genes_thille.bed \\
        --contig-names-a hass_names.tsv
    # a third reference C (e.g. the West Indian T2T), drawn below A
    python3 plot_gene_stacks.py -b genes_thille.bed -a genes_hass.bed -c genes_wi.bed \\
        --tracks-c gwen.clr.wi thille.hifi.wi --labels-c "Gwen → WI" "Thille → WI"

Every figure shows the same 100 kb window at the same scale (bp per inch), so
genes can be compared directly.
"""

import argparse
import dataclasses
import json
import multiprocessing
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

__version__ = "0.10.0"

WINDOW_SIZE = 100_000  # fixed so every gene is drawn at the same resolution
# tandem-array figures: the smallest tier that holds the whole array (Edwin, 2026-09-26),
# so figures within a tier share one resolution
TANDEM_TIERS = [100_000, 250_000, 500_000, 1_000_000, 2_500_000, 10_000_000]
TIER_MARGIN = 1_000  # covers the +-500 bp shift from rounding the window start to kb
RASTERIZE_ABOVE = 100_000  # wider windows embed the scatter as a 600 dpi image; 100 kb figures stay vector
ROUNDING = 3  # window start rounded to the nearest 10**ROUNDING bp, as in graph_coverage.py
BIN_FACTORS = [0, 0.25, 0.75, 1.25, 1.5, np.inf]
CLASS_LABELS = ["low", "half", "medium", "high", "outlier"]
# display scale shared by every coverage panel (Edwin, 2026-09-26): depths above the cap are drawn
# CAP_OFFSET above it so they read as one band; classes are still set against the median before capping
Y_LIMITS = (-5, 125)
DISPLAY_CAP = 120
CAP_OFFSET = 2
BOUNDARY_ALPHA = 0.4
SHADE_ALPHA = 0.08
MEDIAN_COLOR = "#4d4d4d"
OPACITY = 0.05
# --no-text (Edwin, 2026-10-09: no figure text inside the PDFs/PNGs; titles, subtitles and explanatory numbers go to the
# legend md): only axis/tick labels, the color key and short row/genome identifiers are drawn. Set once in main().
NO_TEXT = False
# axes are placed in inches, not with tight_layout, so the data area (and so bp per inch)
# is identical in every figure regardless of tick label or title widths
FIG_WIDTH = 20
AXES_LEFT = 1.0
AXES_WIDTH = 18.5
AXES_HEIGHT = 3.0
PANEL_GAP = 1.0
TOP_MARGIN = 1.0
BOTTOM_MARGIN = 0.6
CHUNK_ROWS = 5_000_000
# annotation tracks sit under a group's coverage panels; heights are fixed so figures stay comparable
GENE_TRACK_HEIGHT = 1.2
TE_TRACK_HEIGHT = 0.8
ANNOTATION_GAP = 0.7   # room for the coverage panel's x tick labels above the first track
TRACK_GAP = 0.15
PRIMARY_FILL = sns.color_palette("colorblind")[6]
GENE_FILL = "#b8bdb5"
MAX_TE_CLASSES = 20
COPY_TRACK_HEIGHT = 0.6
CLASS_ROW_HEIGHT = 0.3  # one row per coverage track in the class strip
# meanings of the itemRgb colours, from the tandem_copies_*.bed track description
TANDEM_COLOR_MEANINGS = {(90, 110, 160): "copy unresolved", (230, 120, 20): "flowering array copy",
                         (200, 0, 0): "curated copy (or closest paralog, ortholog absent)"}  # more classes than this are drawn as "other" so colours stay distinguishable


def log(message):
    print(message, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- BED input

@dataclass(frozen=True)
class Gene:
    name: str
    contig: str
    start: int  # 1-based, inclusive (BED start + 1)
    stop: int   # 1-based, inclusive (BED end)
    strand: str
    label: str = ""  # display name; the key (name) when empty

    @property
    def display(self):
        return self.label or self.name


def read_bed(path, name_fields=None):
    """Return {key: Gene}. BED3+ with name in column 4, optional strand in column 6.

    name_fields=N keys each row by the first N '|'-separated fields of the name, so rows from
    different assemblies (e.g. FAMILY|ver1_id|model|...) share a key.
    """
    genes = {}
    duplicates = []
    with open(path) as bed:
        for line_number, line in enumerate(bed, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: need at least 4 tab-separated columns (contig, start, end, name)")
            contig, start, end, name = fields[0], int(fields[1]), int(fields[2]), fields[3]
            if name_fields:
                name = "|".join(name.split("|")[:name_fields])
            strand = fields[5] if len(fields) >= 6 else "+"
            if strand not in ("+", "-"):
                raise ValueError(f"{path}:{line_number}: strand must be '+' or '-', got {strand!r}")
            if end <= start:
                raise ValueError(f"{path}:{line_number}: end {end} <= start {start}")
            if name in genes:
                duplicates.append(name)
            genes[name] = Gene(name, contig, start + 1, end, strand)
    if duplicates:
        raise ValueError(f"{path}: duplicate gene names {sorted(set(duplicates))}; names must be unique per BED")
    return genes


def display_names(keys, name_fields):
    """key -> display name: FAMILY, or FAMILY (ver1_id) where a family has more than one key; the key itself otherwise."""
    if not name_fields:
        return {key: key for key in keys}
    families = {}
    for key in keys:
        families.setdefault(key.split("|")[0], []).append(key)
    names = {}
    for family, members in families.items():
        for key in members:
            parts = key.split("|")
            names[key] = family if len(members) == 1 or len(parts) < 2 else f"{family} ({parts[1]})"
    return names


class ContigNames:
    """Map a BED contig to the contig used in the coverage files, and to a panel label."""

    def __init__(self, path=None):
        self.table = {}
        if path is None:
            return
        with open(path) as names:
            for line_number, line in enumerate(names, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 3:
                    raise ValueError(f"{path}:{line_number}: need 3 columns: bed_contig, coverage_contig, label")
                self.table[fields[0]] = (fields[1], fields[2])

    def coverage_contig(self, contig):
        return self.table.get(contig, (contig, None))[0]

    def label(self, contig):
        if contig in self.table:
            return self.table[contig][1]
        match = re.fullmatch(r"chr(\d+)", contig)
        return f"Chromosome {int(match.group(1))}" if match else contig


# ---------------------------------------------------------------- coverage cache

def coverage_path(coverage_dir, track, contig):
    return Path(coverage_dir) / f"{track}.{contig}.coverage.gz"


def cache_paths(cache_dir, track, contig):
    # plain concatenation: Path.with_suffix would treat ".chr10" as the suffix and replace it
    base = f"{track}.{contig}"
    return Path(cache_dir) / f"{base}.depth.npy", Path(cache_dir) / f"{base}.json"


def source_stamp(path):
    stat = os.stat(path)
    return {"source": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def cache_is_valid(source, meta_path):
    if not meta_path.is_file():
        return False
    meta = json.loads(meta_path.read_text())
    return {k: meta.get(k) for k in ("source", "size", "mtime_ns")} == source_stamp(source)


def scan_contig(source, contig, npy_path, meta_path):
    """Extract one contig's depths from a whole-genome coverage file into the cache.

    grep keeps only the contig's rows so pandas parses ~50 M rows instead of ~860 M;
    the 1..L position check below catches anything the filter got wrong.
    """
    started = time.time()
    log(f"[scan] {source.name}: {contig} start")
    decompress = subprocess.Popen(["pigz", "-dc", str(source)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    select = subprocess.Popen(["grep", "-F", f"{contig}\t"], stdin=decompress.stdout,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, env={**os.environ, "LC_ALL": "C"})
    decompress.stdout.close()  # grep owns the pipe now, so pigz gets SIGPIPE if grep dies
    positions, depths = [], []
    rows = 0
    try:
        reader = pd.read_csv(select.stdout, sep="\t", header=None, usecols=[0, 1, 2],
                             names=["contig", "position", "depth"],
                             dtype={"contig": "category", "position": "int32", "depth": "int32"},
                             chunksize=CHUNK_ROWS)
        for chunk in reader:
            chunk = chunk[chunk["contig"] == contig]
            positions.append(chunk["position"].to_numpy())
            depths.append(chunk["depth"].to_numpy())
            rows += len(chunk)
    except pd.errors.EmptyDataError:
        pass  # grep found nothing; reported below with the pipeline status
    select_status = select.wait()
    decompress_status = decompress.wait()
    if decompress_status != 0:
        raise RuntimeError(f"pigz failed on {source} (exit {decompress_status}): {decompress.stderr.read().decode().strip()}")
    if select_status not in (0, 1):
        raise RuntimeError(f"grep failed on {source} (exit {select_status}): {select.stderr.read().decode().strip()}")
    if rows == 0:
        raise RuntimeError(f"contig {contig} not found in {source}")

    position = np.concatenate(positions)
    depth = np.concatenate(depths)
    expected = np.arange(1, rows + 1, dtype=np.int32)
    if not np.array_equal(position, expected):
        first_bad = int(np.flatnonzero(position != expected)[0])
        raise RuntimeError(f"{source}: {contig} positions are not 1..{rows} in order "
                           f"(row {first_bad} has position {position[first_bad]}, expected {first_bad + 1})")
    median = float(np.median(depth))

    # write to .tmp then rename, so an interrupted scan never leaves a cache that looks valid
    npy_tmp = npy_path.with_name(npy_path.name + ".tmp")
    with open(npy_tmp, "wb") as handle:
        np.save(handle, depth)
    os.replace(npy_tmp, npy_path)
    meta = {**source_stamp(source), "contig": contig, "length": rows, "median": median,
            "max_depth": int(depth.max()), "version": __version__}
    meta_tmp = meta_path.with_name(meta_path.name + ".tmp")
    meta_tmp.write_text(json.dumps(meta, indent=1))
    os.replace(meta_tmp, meta_path)
    log(f"[scan] {source.name}: {contig} done, {rows:,} positions, median {median}, {time.time() - started:.0f}s")
    return meta


def ensure_cache(pairs, coverage_dir, cache_dir, jobs):
    """Scan every (track, contig) pair whose cache is missing or stale, in parallel."""
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    todo = []
    for track, contig in sorted(pairs):
        source = coverage_path(coverage_dir, track, contig)
        npy_path, meta_path = cache_paths(cache_dir, track, contig)
        if npy_path.is_file() and cache_is_valid(source, meta_path):
            continue
        todo.append((source, contig, npy_path, meta_path))
    if not todo:
        log(f"[cache] all {len(pairs)} track/contig pairs cached")
        return
    log(f"[cache] scanning {len(todo)} of {len(pairs)} track/contig pairs with {min(jobs, len(todo))} workers")
    failures = []
    with ProcessPoolExecutor(max_workers=min(jobs, len(todo))) as pool:
        futures = {pool.submit(scan_contig, *item): item for item in todo}
        for future in as_completed(futures):
            source, contig = futures[future][:2]
            try:
                future.result()
            except Exception as error:  # collect every failed scan before stopping
                failures.append(f"{source} {contig}: {error}")
    if failures:
        raise RuntimeError("coverage scan failed:\n  " + "\n  ".join(failures))


def load_track(cache_dir, track, contig):
    npy_path, meta_path = cache_paths(cache_dir, track, contig)
    meta = json.loads(meta_path.read_text())
    return np.load(npy_path, mmap_mode="r"), meta


# ---------------------------------------------------------------- window and classes

def window(gene, size):
    """Window [start, stop] (1-based, inclusive) of `size` bp centred on the gene, start rounded to kb.

    Same arithmetic as graph_start_stop in graph_coverage.py; Python round() sends ties to the even kb.
    """
    mid = (gene.stop + gene.start) // 2
    start = round(mid - size // 2, -ROUNDING)
    stop = start + size - 1
    if not (start <= gene.start and gene.stop <= stop):
        raise ValueError(f"{gene.name}: gene {gene.start}-{gene.stop} does not fit in window {start}-{stop}")
    return start, stop


def classify(depth, median):
    """Bin by the median (same bins as graph_coverage.py), then cap the drawn value. Returns (drawn depth, categories).

    The category comes from the true depth; only the plotted height is capped at DISPLAY_CAP.
    """
    depth = pd.Series(depth, dtype="int32")
    categories = pd.cut(depth, bins=np.array(BIN_FACTORS) * median, labels=CLASS_LABELS, include_lowest=True)
    drawn = depth.mask(depth.gt(DISPLAY_CAP), other=DISPLAY_CAP + CAP_OFFSET)
    return drawn, categories


# ---------------------------------------------------------------- annotation

@dataclass
class GeneModel:
    gene_id: str
    contig: str
    start: int
    stop: int
    strand: str
    exons: list  # [(start, stop)], 1-based inclusive; empty when the file has gene spans only
    cds: list
    transcript_id: str = ""  # the transcript the exons/CDS come from (longest span); empty for span-only genes


def gtf_attribute(attributes, key):
    match = re.search(rf'{key} "([^"]*)"', attributes)
    return match.group(1) if match else ""


def read_gene_models(path, contigs):
    """GTF -> {contig: [GeneModel]} for the given contigs, one transcript per gene (the longest span).

    Files with only `gene` lines (e.g. the GenBank Hass GTF) give span-only models.
    """
    frame = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 2, 3, 4, 6, 8],
                        names=["contig", "feature", "start", "stop", "strand", "attributes"],
                        dtype={"contig": str, "feature": str, "start": "int64", "stop": "int64", "strand": str, "attributes": str})
    frame = frame[frame["contig"].isin(contigs)]
    frame["gene_id"] = frame["attributes"].map(lambda a: gtf_attribute(a, "gene_id"))
    frame["transcript_id"] = frame["attributes"].map(lambda a: gtf_attribute(a, "transcript_id"))

    genes = {row.gene_id: row for row in frame[frame["feature"] == "gene"].itertuples()}
    transcripts = frame[frame["feature"].isin(["transcript", "mRNA"])]
    orphans = sorted(set(transcripts["gene_id"]) - set(genes))
    if orphans and genes:
        raise ValueError(f"{path}: transcripts whose gene_id has no gene line: {orphans[:10]}")
    longest = {}
    for row in transcripts.itertuples():
        best = longest.get(row.gene_id)
        if best is None or row.stop - row.start > best.stop - best.start:
            longest[row.gene_id] = row
    parts = frame[frame["feature"].isin(["exon", "CDS"])]
    blocks = {}
    for row in parts.itertuples():
        blocks.setdefault((row.transcript_id, row.feature), []).append((row.start, row.stop))

    models = {}
    spans = genes if genes else {row.gene_id: row for row in longest.values()}
    for gene_id, row in spans.items():
        transcript = longest.get(gene_id)
        tid = transcript.transcript_id if transcript is not None else None
        models.setdefault(row.contig, []).append(GeneModel(
            gene_id, row.contig, int(row.start), int(row.stop), row.strand,
            sorted(blocks.get((tid, "exon"), [])), sorted(blocks.get((tid, "CDS"), [])), tid or ""))
    for contig in models:
        models[contig].sort(key=lambda model: model.start)
    return models


def te_class_from_name(name):
    """Class implied by a RepeatModeler/RepeatMasker family name when no classified library is given."""
    if re.fullmatch(r"\([ACGTN]+\)n", name):
        return "Simple_repeat"
    if name.endswith("-rich"):
        return "Low_complexity"
    if name.startswith("ltr-"):
        return "LTR (unclassified)"
    if name.startswith("rnd-"):
        return "Unclassified"
    return "Other"


def read_te_library(path):
    """RepeatModeler families FASTA -> {family name: class}, from headers like >rnd-5_family-4001#LTR/Gypsy."""
    classes = {}
    with open(path) as library:
        for line in library:
            if line.startswith(">"):
                name, _, rest = line[1:].split()[0].partition("#")
                classes[name] = rest or "Unclassified"
    return classes


def read_te(path, contigs, library=None):
    """RepeatMasker GFF -> ({contig: DataFrame start/stop/cls}, class counts over the whole file)."""
    frame = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 3, 4, 8],
                        names=["contig", "start", "stop", "attributes"],
                        dtype={"contig": str, "start": "int64", "stop": "int64", "attributes": str})
    frame["family"] = frame["attributes"].str.extract(r'Motif:([^"]+)"', expand=False)
    missing = frame["family"].isna()
    if missing.any():
        raise ValueError(f"{path}: {int(missing.sum())} rows without a Target \"Motif:...\" family, e.g. {frame.loc[missing, 'attributes'].iloc[0]!r}")
    families = frame["family"].unique()
    lookup = {family: (library.get(family) if library else None) or te_class_from_name(family) for family in families}
    if library:
        unmatched = sum(family not in library for family in families)
        if unmatched:
            log(f"[te] {unmatched} of {len(families)} families in {path} are not in the library; classed by name")
    frame["cls"] = frame["family"].map(lookup)
    counts = frame["cls"].value_counts()
    frame = frame[frame["contig"].isin(contigs)]
    return {contig: part[["start", "stop", "cls"]].reset_index(drop=True) for contig, part in frame.groupby("contig")}, counts


def te_colors(counts):
    """Fixed class -> colour table for the run: most common classes first, the tail pooled as 'other'."""
    classes = list(counts.index)
    palette_ = sns.color_palette("tab20", MAX_TE_CLASSES)
    colors = {cls: palette_[i] for i, cls in enumerate(classes[:MAX_TE_CLASSES])}
    for cls in classes[MAX_TE_CLASSES:]:
        colors[cls] = (0.6, 0.6, 0.6)
    return colors


def pack_rows(intervals):
    """Greedy packing: each (start, stop) goes in the first row whose last stop is before it. Returns row per interval."""
    row_ends = []
    rows = []
    for start, stop in intervals:
        for row, end in enumerate(row_ends):
            if start > end:
                row_ends[row] = stop
                rows.append(row)
                break
        else:
            row_ends.append(stop)
            rows.append(len(row_ends) - 1)
    return rows


@dataclass
class Annotation:
    """What to draw under one panel group."""
    genes: list          # GeneModel overlapping the window (or None when no gene file was given)
    tes: object          # DataFrame of TEs overlapping the window (or None when no TE file was given)
    te_colors: dict
    te_source: str       # legend note on where TE classes came from
    primary: Gene
    bed_genes: list      # every Gene of this reference's BED
    copies: list = None  # TandemCopy overlapping the window (None when no tandem BED was given or none overlap)
    members: list = None  # GeneModel of the array members (--member-models), drawn in their own track
    member_style: dict = None  # gene_id -> (color, legend label)


@dataclass(frozen=True)
class TandemCopy:
    group: str       # first two '|' fields of the name, e.g. KCS11|ChrU.ver1.g440620 or jg101413|jg101413
    name: str
    contig: str
    start: int       # 1-based inclusive
    stop: int
    strand: str
    rgb: tuple
    curated: bool


def read_tandem(path):
    """BED9 with itemRgb -> {contig: [TandemCopy]} sorted by start."""
    copies = {}
    with open(path) as bed:
        for line_number, line in enumerate(bed, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(f"{path}:{line_number}: need BED9 with itemRgb, got {len(fields)} columns")
            rgb = tuple(int(v) for v in fields[8].split(","))
            parts = fields[3].split("|")
            copy = TandemCopy("|".join(parts[:2]), fields[3], fields[0], int(fields[1]) + 1, int(fields[2]),
                              fields[5], rgb, "CURATED" in parts)
            copies.setdefault(copy.contig, []).append(copy)
    unknown = sorted({c.rgb for v in copies.values() for c in v} - set(TANDEM_COLOR_MEANINGS))
    if unknown:
        log(f"[tandem] {path}: colours without a known meaning {unknown}; shown as 'rgb r,g,b'")
    for contig in copies:
        copies[contig].sort(key=lambda c: c.start)
    return copies


def tandem_label(rgb):
    return TANDEM_COLOR_MEANINGS.get(rgb, "rgb " + ",".join(map(str, rgb)))


# ---------------------------------------------------------------- plotting

@dataclass
class Panel:
    track_label: str
    gene: Gene
    chrom_label: str
    depth: np.ndarray  # whole chromosome, 1-based position i at index i-1
    median: int
    shades: list = None  # [(start, stop)] shaded yellow; the gene span when None


def palette():
    colors = sns.color_palette("colorblind")
    return {"outlier": colors[3], "high": colors[4], "medium": colors[0], "half": colors[2], "low": colors[1]}


def draw_panel(ax, panel, start, stop, colors, rasterize=False):
    first = max(start, 1)
    last = min(stop, len(panel.depth))
    positions = np.arange(first, last + 1)
    clipped, categories = classify(np.asarray(panel.depth[first - 1:last]), panel.median)
    frame = pd.DataFrame({"Position": positions, "Coverage": clipped, "Category": categories})
    sns.scatterplot(ax=ax, data=frame, x="Position", y="Coverage", hue="Category",
                    palette=colors, s=15, linewidth=0, alpha=OPACITY, legend=False, rasterized=rasterize)
    ax.axhline(y=panel.median, color=MEDIAN_COLOR, linestyle="solid")
    for y in [0] + [factor * panel.median for factor in BIN_FACTORS[1:-1]] + [DISPLAY_CAP]:
        ax.axhline(y=y, color=MEDIAN_COLOR, linestyle="dashed", linewidth=0.8, alpha=BOUNDARY_ALPHA)
    ax.set_ylim(*Y_LIMITS)
    for shade_start, shade_stop in panel.shades or [(panel.gene.start, panel.gene.stop)]:
        ax.axvspan(shade_start, shade_stop, facecolor="yellow", alpha=SHADE_ALPHA, edgecolor="none")
        for x in (shade_start, shade_stop):
            ax.axvline(x=x, color="goldenrod", linewidth=0.6)
    if NO_TEXT:
        ax.set_title(panel.track_label, loc="left")
    else:
        ax.set_title(f"{panel.track_label}: {panel.gene.display} on {panel.chrom_label} "
                     f"({len(panel.depth) // 1000}kb), median {panel.median}", loc="left")
    ax.set_xlabel("Position")
    ax.set_ylabel("Coverage")


def draw_gene_track(ax, annotation, lo, hi):
    """Gene models packed into rows. Primary gene filled, other BED genes bold-bordered and labelled."""
    ax.set_yticks([])
    ax.set_ylabel("Genes", rotation=0, ha="right", va="center")
    ax.grid(False)
    primary = annotation.primary  # None in a tandem figure whose group is not a BED gene
    others = [gene for gene in annotation.bed_genes if primary is None or gene.name != primary.name]
    models = annotation.genes
    pad = (hi - lo) * 0.01
    rows = pack_rows([(model.start, model.stop + pad) for model in models])
    n_rows = max(rows, default=0) + 1
    ax.set_ylim(n_rows, 0)
    primary_found = False
    for model, row in zip(models, rows):
        is_primary = primary is not None and model.start <= primary.stop and model.stop >= primary.start
        labels = [gene.display for gene in others if gene.contig == model.contig
                  and model.start <= gene.stop and model.stop >= gene.start]
        primary_found |= is_primary
        if is_primary:
            fill, edge, width, label = PRIMARY_FILL, "black", 0.8, annotation.primary.display
        elif labels:
            fill, edge, width, label = GENE_FILL, "black", 2.0, "/".join(labels)
        else:
            fill, edge, width, label = GENE_FILL, "none", 0.0, None
        # box sits in the lower part of its row so the label above it stays inside the axis
        y = row + 0.65
        ax.plot([model.start, model.stop], [y, y], color="black", lw=0.8, zorder=1)
        # strand chevrons along the gene line so the direction reads at any zoom
        step = (hi - lo) / 60
        marks = np.arange(model.start + step / 2, model.stop, step)
        ax.plot(marks, np.full(len(marks), y), linestyle="none", marker=">" if model.strand == "+" else "<",
                markersize=6, markeredgewidth=0.8, color="black", zorder=1)
        blocks = model.exons or [(model.start, model.stop)]
        ax.broken_barh([(a, b - a + 1) for a, b in blocks], (y - 0.15, 0.3),
                       facecolors=fill, edgecolors=edge, linewidths=width, zorder=2)
        if model.cds:
            ax.broken_barh([(a, b - a + 1) for a, b in model.cds], (y - 0.3, 0.6),
                           facecolors=fill, edgecolors=edge, linewidths=width, zorder=3)
        elif not model.exons:
            ax.broken_barh([(model.start, model.stop - model.start + 1)], (y - 0.3, 0.6),
                           facecolors=fill, edgecolors=edge, linewidths=width, zorder=3)
        if label:
            x = (max(model.start, lo) + min(model.stop, hi)) / 2
            ax.text(x, y - 0.33, label, ha="center", va="bottom", fontsize=9, zorder=4, clip_on=True)
    if primary is not None and not primary_found:
        # the BED gene has no annotated model: show its interval so it is visibly missing, not silently absent
        gene = annotation.primary
        ax.broken_barh([(gene.start, gene.stop - gene.start + 1)], (0.2, n_rows - 0.4),
                       facecolors="none", edgecolors="black", linestyles="dashed", linewidths=1.0, zorder=4)
        ax.text((max(gene.start, lo) + min(gene.stop, hi)) / 2, 0.3,
                gene.display if NO_TEXT else f"{gene.display} (no annotated model)",
                ha="center", va="top", fontsize=9, zorder=4, clip_on=True)
    return primary is None or primary_found


def draw_members_track(ax, annotation, lo, hi):
    """Array copies: the member gene models (exon and CDS structure) filled with their --member-models colors."""
    ax.set_yticks([])
    ax.set_ylabel("Array\ncopies", rotation=0, ha="right", va="center")
    ax.grid(False)
    models = annotation.members
    pad = (hi - lo) * 0.01
    rows = pack_rows([(model.start, model.stop + pad) for model in models])
    n_rows = max(rows, default=0) + 1
    ax.set_ylim(n_rows, 0)
    for model, row in zip(models, rows):
        fill = annotation.member_style[model.gene_id][0]
        y = row + 0.5
        ax.plot([model.start, model.stop], [y, y], color="black", lw=0.8, zorder=1)
        step = (hi - lo) / 60
        marks = np.arange(model.start + step / 2, model.stop, step)
        ax.plot(marks, np.full(len(marks), y), linestyle="none", marker=">" if model.strand == "+" else "<",
                markersize=6, markeredgewidth=0.8, color="black", zorder=1)
        blocks = model.exons or [(model.start, model.stop)]
        ax.broken_barh([(a, b - a + 1) for a, b in blocks], (y - 0.15, 0.3), facecolors=fill, edgecolors="black",
                       linewidths=0.6, zorder=2)
        if model.cds:
            ax.broken_barh([(a, b - a + 1) for a, b in model.cds], (y - 0.3, 0.6), facecolors=fill,
                           edgecolors="black", linewidths=0.6, zorder=3)
    return {annotation.member_style[m.gene_id] for m in models}


def draw_te_track(ax, annotation):
    ax.set_yticks([])
    ax.set_ylabel("TEs", rotation=0, ha="right", va="center")
    ax.grid(False)
    tes = annotation.tes
    rows = pack_rows(list(zip(tes["start"], tes["stop"])))
    n_rows = max(rows, default=0) + 1
    ax.set_ylim(n_rows, 0)
    for row in range(n_rows):
        part = tes[np.array(rows) == row]
        ax.broken_barh(list(zip(part["start"], part["stop"] - part["start"] + 1)), (row + 0.1, 0.8),
                       facecolors=[annotation.te_colors[cls] for cls in part["cls"]], linewidths=0)
    return set(tes["cls"])


def draw_copies_track(ax, copies):
    """Tandem copies in the BED's own colours; copies flagged CURATED in the name get a bold edge."""
    ax.set_yticks([])
    ax.set_ylabel("Tandem\ncopies", rotation=0, ha="right", va="center")
    ax.grid(False)
    rows = pack_rows([(c.start, c.stop) for c in copies])
    n_rows = max(rows, default=0) + 1
    ax.set_ylim(n_rows, 0)
    for copy, row in zip(copies, rows):
        ax.broken_barh([(copy.start, copy.stop - copy.start + 1)], (row + 0.1, 0.8),
                       facecolors=[tuple(v / 255 for v in copy.rgb)],
                       edgecolors="black" if copy.curated else "none", linewidths=1.5 if copy.curated else 0)
    return {c.rgb for c in copies}


def draw_class_strip(ax, group, start, stop, colors):
    """One row per coverage panel of the group: every position colored by its class (same classes and colors as the
    panels, from the true depth against the chromosome median)."""
    first = max(start, 1)
    last = min(stop, min(len(p.depth) for p in group))
    order = list(CLASS_LABELS)
    codes = np.vstack([np.asarray(classify(np.asarray(p.depth[first - 1:last]), p.median)[1].cat.codes) for p in group])
    cmap = matplotlib.colors.ListedColormap([colors[c] for c in order])
    ax.imshow(codes, aspect="auto", interpolation="nearest", cmap=cmap, vmin=-0.5, vmax=len(order) - 0.5,
              extent=(first - 0.5, last + 0.5, len(group) - 0.5, -0.5), rasterized=True)
    ax.set_yticks(range(len(group)), [p.track_label for p in group], fontsize=8)
    ax.set_ylim(len(group) - 0.5, -0.5)
    ax.grid(False)
    if not NO_TEXT:
        ax.set_title("Coverage class per position (rows as the panels above)", loc="left", fontsize=9)


def plot_stack(name, groups, size, out_path, annotations=None, title=None, shade_label=None,
               pdf=True, png_dpi=None, class_strip=False):
    """groups: list of panel lists, top to bottom (B reference first). Panels in a group share an x axis.

    annotations: optional list parallel to groups; an Annotation adds gene/TE tracks under that group.
    Without annotations the layout is exactly the coverage-only layout.
    """
    annotations = annotations or [None] * len(groups)
    colors = palette()
    sns.set_style("whitegrid")
    # (kind, height, gap above) for every axis, top to bottom
    rows = []
    for group, annotation in zip(groups, annotations):
        for _panel in group:
            rows.append(("coverage", AXES_HEIGHT, PANEL_GAP))
        if class_strip:
            rows.append(("classes", CLASS_ROW_HEIGHT * len(group), ANNOTATION_GAP))
        if annotation is not None and annotation.members:
            rows.append(("members", GENE_TRACK_HEIGHT, ANNOTATION_GAP))
        if annotation is not None and annotation.genes is not None:
            rows.append(("genes", GENE_TRACK_HEIGHT, TRACK_GAP if rows[-1][0] == "members" else ANNOTATION_GAP))
        if annotation is not None and annotation.tes is not None:
            gap = TRACK_GAP if rows[-1][0] == "genes" else ANNOTATION_GAP
            rows.append(("tes", TE_TRACK_HEIGHT, gap))
        if annotation is not None and annotation.copies:
            gap = TRACK_GAP if rows[-1][0] in ("genes", "tes") else ANNOTATION_GAP
            rows.append(("copies", COPY_TRACK_HEIGHT, gap))
    height = TOP_MARGIN + BOTTOM_MARGIN + sum(h for _k, h, _g in rows) + sum(g for _k, _h, g in rows[1:])
    fig = plt.figure(figsize=(FIG_WIDTH, height))
    axes = []
    top = height - TOP_MARGIN
    for i, (_kind, axis_height, gap) in enumerate(rows):
        if i:
            top -= gap
        bottom = top - axis_height
        axes.append(fig.add_axes([AXES_LEFT / FIG_WIDTH, bottom / height, AXES_WIDTH / FIG_WIDTH, axis_height / height]))
        top = bottom
    kb = FuncFormatter(lambda x, _pos: f"{int(x / 1000)}kb")
    index = 0
    te_classes_shown = set()
    copy_colors_shown = set()
    member_styles_shown = set()
    missing_models = []
    for group, annotation in zip(groups, annotations):
        gene = group[0].gene
        start, stop = window(gene, size)
        first_ax = axes[index]
        for panel in group:
            ax = axes[index]
            if ax is not first_ax:
                ax.sharex(first_ax)
            draw_panel(ax, panel, start, stop, colors, rasterize=size > RASTERIZE_ABOVE)
            ax.xaxis.set_major_formatter(kb)
            index += 1
        if class_strip:
            ax = axes[index]
            ax.sharex(first_ax)
            draw_class_strip(ax, group, start, stop, colors)
            ax.xaxis.set_major_formatter(kb)
            index += 1
        if annotation is not None:
            if annotation.members:
                ax = axes[index]
                ax.sharex(first_ax)
                member_styles_shown |= draw_members_track(ax, annotation, start, stop)
                ax.xaxis.set_major_formatter(kb)
                index += 1
            if annotation.genes is not None:
                ax = axes[index]
                ax.sharex(first_ax)
                if not draw_gene_track(ax, annotation, start, stop):
                    missing_models.append(annotation.primary.display)
                if annotation.members:
                    ax.set_ylabel("Other\ngenes", rotation=0, ha="right", va="center")
                ax.xaxis.set_major_formatter(kb)
                index += 1
            if annotation.tes is not None:
                ax = axes[index]
                ax.sharex(first_ax)
                te_classes_shown |= draw_te_track(ax, annotation)
                ax.xaxis.set_major_formatter(kb)
                index += 1
            if annotation.copies:
                ax = axes[index]
                ax.sharex(first_ax)
                copy_colors_shown |= draw_copies_track(ax, annotation.copies)
                ax.xaxis.set_major_formatter(kb)
                index += 1
            # only the group's bottom track keeps its x tick labels
            n_tracks = (bool(annotation.members) + (annotation.genes is not None) + (annotation.tes is not None)
                        + bool(annotation.copies))
            for track_ax in axes[index - n_tracks:index - 1]:
                track_ax.tick_params(labelbottom=False)
        # a - strand gene is drawn right to left so it reads like a + strand gene
        first_ax.set_xlim((stop, start) if gene.strand == "-" else (start, stop))

    legend = [
        Line2D([0], [0], color=MEDIAN_COLOR, lw=2, label="Chromosome median"),
        Line2D([0], [0], color=MEDIAN_COLOR, lw=0.8, linestyle="dashed", alpha=BOUNDARY_ALPHA,
               label=f"0, 0.25/0.75/1.25/1.5x Median, {DISPLAY_CAP}"),
        Line2D([0], [0], marker="s", color="w", label=shade_label or f"{name} gene", markerfacecolor=sns.color_palette("colorblind")[8],
               alpha=0.6, markersize=10),
        Line2D([0], [0], marker="o", color="w", label=f"Outliers > 1.5x Median (> {DISPLAY_CAP} drawn at {DISPLAY_CAP + CAP_OFFSET})",
               markerfacecolor=colors["outlier"], markersize=10),
        Line2D([0], [0], marker="o", color="w", label="High <= 1.5x Median", markerfacecolor=colors["high"], markersize=10),
        Line2D([0], [0], marker="o", color="w", label="Medium <= 1.25x Median", markerfacecolor=colors["medium"], markersize=10),
        Line2D([0], [0], marker="o", color="w", label="Half <= 0.75x Median", markerfacecolor=colors["half"], markersize=10),
        Line2D([0], [0], marker="o", color="w", label="Low <= 0.25x Median", markerfacecolor=colors["low"], markersize=10),
    ]
    gene_annotations = [a for a in annotations if a is not None and a.genes is not None]
    if gene_annotations:
        if any(a.primary is not None for a in gene_annotations):
            legend.append(Line2D([0], [0], marker="s", color="w", label=f"{name} gene model", markerfacecolor=PRIMARY_FILL,
                                 markeredgecolor="black", markersize=10))
        legend += [
            Line2D([0], [0], marker="s", color="w", label="Other BED gene", markerfacecolor=GENE_FILL,
                   markeredgecolor="black", markeredgewidth=2, markersize=10),
            Line2D([0], [0], marker="s", color="w", label="Other annotated gene", markerfacecolor=GENE_FILL, markersize=10),
        ]
    te_annotation = next((a for a in annotations if a is not None and a.tes is not None), None)
    if te_annotation is not None and te_classes_shown:
        legend.append(Line2D([0], [0], color="none", label=f"TE class ({te_annotation.te_source})"))
        for cls in [c for c in te_annotation.te_colors if c in te_classes_shown]:
            legend.append(Line2D([0], [0], marker="s", color="w", label=cls,
                                 markerfacecolor=te_annotation.te_colors[cls], markersize=10))
    if member_styles_shown:
        legend.append(Line2D([0], [0], color="none", label="Array copies (gene models)"))
        for color, label in sorted(member_styles_shown, key=lambda s: s[1]):
            legend.append(Line2D([0], [0], marker="s", color="w", label=label, markerfacecolor=color,
                                 markeredgecolor="black", markersize=10))
    if copy_colors_shown:
        legend.append(Line2D([0], [0], color="none", label="Tandem copies"))
        for rgb in [c for c in TANDEM_COLOR_MEANINGS if c in copy_colors_shown] + sorted(copy_colors_shown - set(TANDEM_COLOR_MEANINGS)):
            legend.append(Line2D([0], [0], marker="s", color="w", label=tandem_label(rgb),
                                 markerfacecolor=tuple(v / 255 for v in rgb), markersize=10))
    if NO_TEXT:
        # what the title and panel titles said goes to a sidecar table for the figure legend
        with open(Path(out_path).with_name(Path(out_path).stem + "_notes.tsv"), "w") as fh:
            fh.write("title\ttrack\tchrom\tchrom_kb\tchrom_median\n")
            for group in groups:
                for p in group:
                    fh.write(f"{title or name}\t{p.track_label}\t{p.chrom_label}\t{len(p.depth) // 1000}\t{p.median}\n")
    else:
        fig.suptitle(title or f"{size // 10**ROUNDING}kb surrounding {name}", fontsize=16, y=1 - 0.3 / height)
    fig.legend(handles=legend, loc="center left", bbox_to_anchor=((AXES_LEFT + AXES_WIDTH + 0.2) / FIG_WIDTH, 0.5))
    # out_path names the PDF; the PNG goes next to it with the same stem
    if pdf:
        fig.savefig(out_path, dpi=600, bbox_inches="tight")
    if png_dpi:
        fig.savefig(Path(out_path).with_suffix(".png"), dpi=png_dpi, bbox_inches="tight")
    plt.close(fig)
    return missing_models


# ---------------------------------------------------------------- driver

@dataclass
class Reference:
    key: str        # "b", "a" or "c"
    genes: dict
    tracks: list
    labels: list
    names: ContigNames
    gene_models: dict = None   # {contig: [GeneModel]} when --genes-gtf is given
    tes: dict = None           # {contig: DataFrame} when --te-gff is given
    te_colors: dict = None
    te_source: str = ""
    tandem: dict = None        # {contig: [TandemCopy]} when --tandem is given


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", "--bed-a", help="BED of genes on the A reference (Hass); drawn in the bottom panels")
    parser.add_argument("-b", "--bed-b", help="BED of genes on the B reference (Thille); drawn in the top panels")
    parser.add_argument("--tracks-a", nargs="+", default=["gwen.clr.hass", "thille.hifi.hass"],
                        help="coverage prefixes on A, top to bottom (default: %(default)s)")
    parser.add_argument("--labels-a", nargs="+", default=["Gwen → Hass", "Thille → Hass"],
                        help="panel labels for --tracks-a (default: %(default)s)")
    parser.add_argument("--tracks-b", nargs="+", default=["gwen.clr.thille", "thille.hifi.thille"],
                        help="coverage prefixes on B, top to bottom (default: %(default)s)")
    parser.add_argument("--labels-b", nargs="+", default=["Gwen → Thille", "Thille → Thille"],
                        help="panel labels for --tracks-b (default: %(default)s)")
    parser.add_argument("-c", "--bed-c", help="BED of genes on an optional third reference C; drawn below the A panels")
    parser.add_argument("--tracks-c", nargs="+", default=[], help="coverage prefixes on C, top to bottom (needed with --bed-c)")
    parser.add_argument("--labels-c", nargs="+", default=[], help="panel labels for --tracks-c")
    parser.add_argument("--contig-names-a", help="TSV bed_contig, coverage_contig, label for A (e.g. CM056818.1 chr10 'Chromosome 10')")
    parser.add_argument("--contig-names-b", help="TSV bed_contig, coverage_contig, label for B")
    parser.add_argument("--contig-names-c", help="TSV bed_contig, coverage_contig, label for C")
    for key, ref_name in (("a", "Hass"), ("b", "Thille"), ("c", "third reference")):
        parser.add_argument(f"--genes-gtf-{key}", help=f"GTF gene annotation on {key.upper()} ({ref_name}); adds a gene-model track under the {key.upper()} panels")
        parser.add_argument(f"--te-gff-{key}", help=f"RepeatMasker GFF on {key.upper()}; adds a TE track under the {key.upper()} panels")
        parser.add_argument(f"--te-library-{key}", help=f"RepeatModeler families FASTA (>name#Class/Family) giving TE classes for --te-gff-{key}; "
                                                     "without it classes come from the family name prefix")
    parser.add_argument("--tandem-b", help="tandem-copies BED9 (itemRgb) on B; adds a copies track where copies fall in the window")
    parser.add_argument("--tandem-a", help="tandem-copies BED9 (itemRgb) on A")
    parser.add_argument("--tandem-c", help="tandem-copies BED9 (itemRgb) on C")
    parser.add_argument("--coverage-dir", default="covfiles", help="directory of <track>.<contig>.coverage.gz (default: %(default)s)")
    parser.add_argument("--cache-dir", default="cache", help="per track/contig depth cache (default: %(default)s)")
    parser.add_argument("--genes", nargs="+", help="only plot these genes, by key or display name (default: every gene in the BED files)")
    parser.add_argument("--title-map", help="TSV key<TAB>file_stem<TAB>title[<TAB>label]: per gene key (BED name fields), the output "
                        "file stem, the figure title and optionally the gene label used in panel titles, legend and gene "
                        "track; {chrom} in the title becomes the chromosome label of the first reference. Genes not "
                        "listed keep the default name and title (default: none)")
    parser.add_argument("--name-fields", type=int, help="key genes by the first N '|'-separated fields of the BED name "
                        "(use 2 for the flowering BEDs: FAMILY|ver1_id); default: the whole name")
    parser.add_argument("--member-models", help="TSV gene_id<TAB>color<TAB>legend label: these gene models (array copies) "
                        "are drawn in their own 'Array copies' track in the given colors, the gene track shows the other "
                        "genes, and tandem copy boxes are not drawn (default: none; output unchanged)")
    parser.add_argument("--tandem-figures", action="store_true",
                        help="draw one figure per tandem group (from --tandem-a/-b) instead of one per BED gene; the window is "
                             f"the smallest tier in {[t // 1000 for t in TANDEM_TIERS]} kb that holds the whole array")
    parser.add_argument("--class-strip", action="store_true",
                        help="add a class strip under each reference's coverage panels: one row per panel, every position "
                             "colored by its coverage class (default off; output without it is unchanged)")
    parser.add_argument("--no-text", action="store_true",
                        help="no titles or explanatory numbers in the figure (axis labels, color key and track identifiers "
                             "only); the title, chromosome, window and chromosome medians go to <stem>_notes.tsv "
                             "(default off; output without it is unchanged)")
    parser.add_argument("--no-pdf", action="store_true", help="write only the PNGs (default: PDF and PNG for every figure)")
    parser.add_argument("--png-dpi", type=int, default=600, help="PNG resolution (default: %(default)s, same as the PDFs)")
    parser.add_argument("--outdir", help="output directory (default: Results, or Results_tandem with --tandem-figures)")
    parser.add_argument("--jobs", type=int, default=8, help="parallel coverage scans (default: %(default)s)")
    parser.add_argument("--plot-jobs", type=int, default=1,
                        help="draw figures in this many processes (default: %(default)s); each 10 Mb figure uses ~3 GB RAM")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args(argv)
    if args.tandem_figures:
        if not args.tandem_a and not args.tandem_b:
            parser.error("--tandem-figures needs --tandem-a and/or --tandem-b")
    elif not args.bed_a and not args.bed_b:
        parser.error("give at least one of -a/--bed-a or -b/--bed-b")
    args.member_style = None
    if args.member_models:
        args.member_style = {}
        for line in open(args.member_models):
            if line.strip() and not line.startswith("#"):
                gene_id, color, label = line.rstrip("\n").split("\t")[:3]
                args.member_style[gene_id] = (color, label)
    if args.plot_jobs < 1:
        parser.error("--plot-jobs must be at least 1")
    if args.png_dpi <= 0:
        parser.error("--png-dpi must be positive")
    if args.outdir is None:
        args.outdir = "Results_tandem" if args.tandem_figures else "Results"
    if args.bed_c and not args.tracks_c:
        parser.error("--bed-c needs --tracks-c")
    for key in ("a", "b", "c"):
        if len(getattr(args, f"tracks_{key}")) != len(getattr(args, f"labels_{key}")):
            parser.error(f"--tracks-{key} and --labels-{key} need the same number of entries")
        if getattr(args, f"te_library_{key}") and not getattr(args, f"te_gff_{key}"):
            parser.error(f"--te-library-{key} needs --te-gff-{key}")
    return args


@dataclass
class Figure:
    """One output figure: per reference, the region to centre on and what to shade."""
    key: str
    display: str
    size: int
    out_name: str
    regions: dict        # ref key -> Gene (the BED gene, or the tandem array span)
    shades: dict         # ref key -> [(start, stop)] or None (shade the gene span)
    primaries: dict      # ref key -> Gene or None (gene track highlight)
    title: str = None
    shade_label: str = None


def tier_for(span):
    """Smallest tandem tier that holds an array of `span` bp after kb rounding; None when none does."""
    return next((tier for tier in TANDEM_TIERS if span <= tier - TIER_MARGIN), None)


def safe_file_name(text):
    return re.sub(r"_+", "_", re.sub(r"[^\w.+-]", "_", text)).strip("_")


def gene_figures(references, names, display):
    figures = []
    for name in names:
        regions = {ref.key: ref.genes[name] for ref in references if name in ref.genes}
        figures.append(Figure(name, display[name], WINDOW_SIZE,
                              f"{safe_file_name(display[name])}_{WINDOW_SIZE // 10**ROUNDING}kb_stacked.pdf",
                              regions, {k: None for k in regions}, dict(regions)))
    return figures


def tandem_figures(references, display, skipped):
    """One figure per tandem group; both references use the tier that fits the larger array."""
    keys = sorted({c.group for ref in references if ref.tandem for v in ref.tandem.values() for c in v})
    figures = []
    for key in keys:
        parts = key.split("|")
        name = display.get(key) or (parts[0] if len(parts) < 2 or parts[0] == parts[1] else f"{parts[0]} ({parts[1]})")
        regions, shades, primaries, spans = {}, {}, {}, {}
        for ref in references:
            copies = [c for v in (ref.tandem or {}).values() for c in v if c.group == key]
            if not copies:
                continue
            contigs = {c.contig for c in copies}
            if len(contigs) > 1:
                raise ValueError(f"tandem group {key} spans several contigs on {ref.key.upper()}: {sorted(contigs)}")
            start, stop = min(c.start for c in copies), max(c.stop for c in copies)
            regions[ref.key] = Gene(key, contigs.pop(), start, stop, "+", name)
            shades[ref.key] = [(c.start, c.stop) for c in copies]
            primaries[ref.key] = ref.genes.get(key)
            spans[ref.key] = stop - start + 1
        tier = tier_for(max(spans.values()))
        if tier is None:
            skipped.append(f"{name}: array span {max(spans.values()):,} bp exceeds the largest tier {TANDEM_TIERS[-1]:,}")
            continue
        n_copies = sum(len(v) for v in shades.values())
        span_text = ", ".join(f"{k.upper()} {v:,} bp" for k, v in spans.items())
        figures.append(Figure(key, name, tier, f"{safe_file_name(name)}_{tier // 10**ROUNDING}kb_tandem.pdf",
                              regions, shades, primaries,
                              title=f"{tier // 10**ROUNDING}kb around tandem array {name} ({n_copies} copies; span {span_text})",
                              shade_label="Tandem copies (shaded)"))
    return figures


_RENDER_CONTEXT = None  # (figures, references, args), set by main before drawing; inherited by forked workers


def render_figure_by_index(index):
    figures, references, args = _RENDER_CONTEXT
    return render_figure(figures[index], references, args)


def render_figure(figure, references, args):
    """Build the panels and tracks of one figure and save it. Returns a result dict; never raises for a bad figure."""
    result = {"out_path": None, "failed": None, "no_model": [], "skipped": []}
    size = figure.size
    groups = []
    annotations = []
    too_long = []
    for ref in references:
        gene = figure.regions.get(ref.key)
        if gene is None:
            continue
        if gene.stop - gene.start > size:
            too_long.append(f"{ref.key.upper()} {gene.stop - gene.start + 1}bp")
            continue
        contig = ref.names.coverage_contig(gene.contig)
        group = []
        for track, label in zip(ref.tracks, ref.labels):
            depth, meta = load_track(args.cache_dir, track, contig)
            group.append(Panel(label, gene, ref.names.label(gene.contig), depth, int(meta["median"]),
                               figure.shades.get(ref.key)))
        groups.append(group)
        try:
            lo, hi = window(gene, size)
        except ValueError:
            annotations.append(None)  # plot_stack raises the same error inside the try below
            continue
        copies_in = None
        if ref.tandem is not None:
            copies_in = [c for c in ref.tandem.get(gene.contig, []) if c.start <= hi and c.stop >= lo] or None
        if ref.gene_models is None and ref.tes is None and copies_in is None:
            annotations.append(None)
            continue
        genes_in = None
        if ref.gene_models is not None:
            genes_in = [m for m in ref.gene_models.get(gene.contig, []) if m.start <= hi and m.stop >= lo]
        tes_in = None
        if ref.tes is not None:
            part = ref.tes.get(gene.contig)
            tes_in = (part[(part["start"] <= hi) & (part["stop"] >= lo)] if part is not None
                      else pd.DataFrame({"start": [], "stop": [], "cls": []}))
        members_in = None
        if args.member_style and genes_in is not None:
            # array members get their own track; the copies BED then only sizes and shades the window
            members_in = [m for m in genes_in if m.gene_id in args.member_style] or None
            genes_in = [m for m in genes_in if m.gene_id not in args.member_style]
            copies_in = None
        annotations.append(Annotation(genes_in, tes_in, ref.te_colors, ref.te_source, figure.primaries.get(ref.key),
                                      list(ref.genes.values()), copies_in, members_in, args.member_style))
    if too_long:
        result["skipped"].append(f"{figure.display} at {size}bp: longer than window ({', '.join(too_long)})")
    if not groups:
        return result
    out_path = Path(args.outdir) / figure.out_name
    started = time.time()
    try:
        zero = [f"{p.track_label} on {p.chrom_label}" for group in groups for p in group if p.median == 0]
        if zero:
            raise ValueError(f"median coverage is 0 for {zero}")
        result["no_model"] = plot_stack(figure.display, groups, size, out_path, annotations, pdf=not args.no_pdf,
                                        png_dpi=args.png_dpi, title=figure.title, shade_label=figure.shade_label,
                                        class_strip=args.class_strip)
    except ValueError as error:  # a bad figure is reported, the rest still get drawn
        result["failed"] = f"{figure.display} at {size}bp: {error}"
        return result
    result["out_path"] = out_path
    outputs = ([out_path.name] if not args.no_pdf else []) + [out_path.with_suffix(".png").name]
    log(f"[plot] {args.outdir}/{' + '.join(outputs)}  {time.time() - started:.0f}s")
    return result


def main(argv=None):
    args = parse_args(argv)
    global NO_TEXT
    NO_TEXT = args.no_text   # module global so forked plot workers inherit it
    started = time.time()
    mode = "tandem figures" if args.tandem_figures else f"gene figures, window={WINDOW_SIZE}"
    log(f"plot_gene_stacks {__version__}  {mode}  bed-b={args.bed_b}  bed-a={args.bed_a}  bed-c={args.bed_c}  tandem-b={args.tandem_b}  "
        f"tandem-a={args.tandem_a}  coverage-dir={args.coverage_dir}  cache-dir={args.cache_dir}  outdir={args.outdir}  jobs={args.jobs}")

    references = []  # top to bottom: B above A, then C when given
    for key in ("b", "a", "c"):
        bed = getattr(args, f"bed_{key}")
        tandem_path = getattr(args, f"tandem_{key}")
        if bed or (args.tandem_figures and tandem_path):
            references.append(Reference(key, read_bed(bed, args.name_fields) if bed else {}, getattr(args, f"tracks_{key}"),
                                        getattr(args, f"labels_{key}"), ContigNames(getattr(args, f"contig_names_{key}"))))
    for ref in references:
        tandem_path = getattr(args, f"tandem_{ref.key}")
        if tandem_path:
            ref.tandem = read_tandem(tandem_path)
            log(f"[tandem] {ref.key.upper()}: {sum(len(v) for v in ref.tandem.values())} copies in "
                f"{len({c.group for v in ref.tandem.values() for c in v})} groups from {tandem_path}")

    names = sorted(set().union(*(ref.genes for ref in references)))
    display = display_names(names, args.name_fields)
    tmap = {}
    if args.title_map:
        for line in open(args.title_map):
            if line.strip() and not line.startswith("#"):
                f = line.rstrip("\n").split("\t")
                tmap[f[0]] = (f[1], f[2])
                if len(f) > 3 and f[3] and f[0] in display:
                    display[f[0]] = f[3]
    if len(set(display.values())) != len(display):
        sys.exit("display names are not unique; use a larger --name-fields")
    for ref in references:
        ref.genes = {key: dataclasses.replace(gene, label=display[key]) for key, gene in ref.genes.items()}

    skipped = []
    figures = tandem_figures(references, display, skipped) if args.tandem_figures else gene_figures(references, names, display)
    if args.genes:
        known = {f.key for f in figures} | {f.display for f in figures}
        unknown = sorted(set(args.genes) - known)
        if unknown:
            sys.exit(f"--genes not found: {unknown}")
        figures = [f for f in figures if f.key in args.genes or f.display in args.genes]
    if tmap:  # e.g. match the variant-matrix figures' gene labels and headers
        renamed = []
        for f in figures:
            if f.key in tmap:
                stem, title = tmap[f.key]
                ref = next(r for r in references if r.key in f.regions)
                title = title.replace("{chrom}", ref.names.label(f.regions[ref.key].contig))
                f = dataclasses.replace(f, out_name=f"{safe_file_name(stem)}_{WINDOW_SIZE // 10**ROUNDING}kb_stacked.pdf",
                                        title=title)
            renamed.append(f)
        figures = renamed
        if len({f.out_name for f in figures}) != len(figures):
            sys.exit("--title-map gives two figures the same file name")
    if not figures:
        sys.exit("nothing to plot")
    unmatched = [f.display for f in figures if len(references) > 1 and len(f.regions) == 1]

    pairs = set()
    missing = []
    for figure in figures:
        for ref in references:
            if ref.key not in figure.regions:
                continue
            contig = ref.names.coverage_contig(figure.regions[ref.key].contig)
            for track in ref.tracks:
                pairs.add((track, contig))
                if not coverage_path(args.coverage_dir, track, contig).is_file():
                    missing.append(f"{coverage_path(args.coverage_dir, track, contig)} ({figure.display})")
    if missing:
        sys.exit("missing coverage files:\n  " + "\n  ".join(sorted(set(missing))))

    ensure_cache(pairs, args.coverage_dir, args.cache_dir, args.jobs)

    for ref in references:
        contigs = {f.regions[ref.key].contig for f in figures if ref.key in f.regions}
        gtf = getattr(args, f"genes_gtf_{ref.key}")
        if gtf:
            ref.gene_models = read_gene_models(gtf, contigs)
            absent = sorted(contigs - set(ref.gene_models))
            if absent:
                log(f"[genes] {ref.key.upper()}: no annotated genes on contigs {absent} (check contig names in {gtf})")
            log(f"[genes] {ref.key.upper()}: {sum(len(v) for v in ref.gene_models.values()):,} gene models on {len(ref.gene_models)} contigs from {gtf}")
        te_gff = getattr(args, f"te_gff_{ref.key}")
        if te_gff:
            library_path = getattr(args, f"te_library_{ref.key}")
            library = read_te_library(library_path) if library_path else None
            ref.tes, counts = read_te(te_gff, contigs, library)
            ref.te_colors = te_colors(counts)
            ref.te_source = "from library" if library else "from family name"
            log(f"[te] {ref.key.upper()}: {int(counts.sum()):,} TEs, {len(counts)} classes ({ref.te_source}) from {te_gff}")

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    global _RENDER_CONTEXT
    _RENDER_CONTEXT = (figures, references, args)
    results = []
    if args.plot_jobs > 1 and len(figures) > 1:
        # fork: workers inherit the loaded annotations and coverage memmaps instead of pickling them
        workers = min(args.plot_jobs, len(figures))
        log(f"[plot] drawing {len(figures)} figures with {workers} workers")
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("fork")) as pool:
            futures = [pool.submit(render_figure_by_index, i) for i in range(len(figures))]
            results = [future.result() for future in futures]
    else:
        results = [render_figure_by_index(i) for i in range(len(figures))]
    written = [r["out_path"] for r in results if r["out_path"] is not None]
    failed = [r["failed"] for r in results if r["failed"]]
    no_model = [name for r in results for name in r["no_model"]]
    skipped += [line for r in results for line in r["skipped"]]

    log(f"\nsummary: {len(written)} figures in {args.outdir}/, {len(skipped)} skipped, {len(failed)} failed, "
        f"{time.time() - started:.0f}s")
    if args.tandem_figures:
        tiers = {}
        for figure in figures:
            tiers[figure.size] = tiers.get(figure.size, 0) + 1
        log("  tiers: " + ", ".join(f"{t // 1000} kb: {n}" for t, n in sorted(tiers.items())))
    for line in skipped:
        log(f"  skipped {line}")
    for name in unmatched:
        log(f"  only in one reference: {name}")
    for line in sorted(set(no_model)):
        log(f"  no annotated gene model overlaps BED gene: {line}")
    for line in failed:
        log(f"  FAILED {line}")
    return 1 if failed else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError) as error:
        log(f"error: {error}")
        sys.exit(1)
