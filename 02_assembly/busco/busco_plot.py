#!/usr/bin/env python3
"""
busco_plot.py - Publication-ready BUSCO plotting tool in 100% pure Python.

Drop-in replacement for the legacy R/ggplot2-dependent generate_plot.py script.
Generates both high-resolution PNG and vector PDF figures without needing R,
Rscript, or ggplot2 installed.

Uses colorblind-safe palettes (Okabe-Ito / Color Universal Design by default)
ensuring full accessibility across protanopia, deuteranopia, and tritanopia.

Supports both BUSCO short_summary.txt and short_summary.json files across all
BUSCO versions (v2 - v6).
"""

import os
import sys
import re
import json
import glob
import argparse
from pathlib import Path
from dataclasses import dataclass
from typing import List, Optional, Tuple, Dict, Any

import matplotlib
matplotlib.use("Agg")  # Headless backend suitable for server environments
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np


# Available palettes:
# Default: Okabe-Ito (Color Universal Design - Bang Wong, Nature Methods 2011).
# S and F were already Okabe-Ito Sky Blue and Yellow in legacy BUSCO.
# D and M in legacy BUSCO were non-colorblind-safe (#3492C7 and #F04442);
# they are replaced by the closest colorblind-safe colors in the palette:
# D -> #0072B2 (Blue) and M -> #D55E00 (Vermilion).
PALETTES = {
    "colorblind": {
        "S": "#56B4E9",  # Okabe-Ito Sky Blue (Complete single-copy)
        "D": "#0072B2",  # Okabe-Ito Blue (Complete duplicated)
        "F": "#F0E442",  # Okabe-Ito Yellow (Fragmented)
        "M": "#D55E00",  # Okabe-Ito Vermilion (Missing)
    },
    "okabe-ito": {
        "S": "#56B4E9",
        "D": "#0072B2",
        "F": "#F0E442",
        "M": "#D55E00",
    },
    "tol": {
        "S": "#66CCEE",  # Paul Tol Cyan
        "D": "#4477AA",  # Paul Tol Blue
        "F": "#CCBB44",  # Paul Tol Yellow
        "M": "#EE6677",  # Paul Tol Red
    },
    "legacy": {
        "S": "#56B4E9",  # Legacy BUSCO Sky Blue
        "D": "#3492C7",  # Legacy BUSCO Steel Blue (non-standard)
        "F": "#F0E442",  # Legacy BUSCO Yellow
        "M": "#F04442",  # Legacy BUSCO Red (not colorblind safe)
    },
}

CATEGORY_LABELS = {
    "S": "Complete (C) and single-copy (S)",
    "D": "Complete (C) and duplicated (D)",
    "F": "Fragmented (F)",
    "M": "Missing (M)",
}


@dataclass
class BuscoResult:
    label: str
    file_path: str
    dataset: str
    version: str
    single: int
    duplicated: int
    fragmented: int
    missing: int
    total: int
    single_pct: float
    duplicated_pct: float
    fragmented_pct: float
    missing_pct: float
    complete_pct: float
    summary_string: str = ""

    @property
    def complete(self) -> int:
        return self.single + self.duplicated


def parse_busco_txt(file_path: str) -> BuscoResult:
    """Parse a BUSCO short_summary.txt file."""
    path = Path(file_path)
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()

    # Extract dataset name
    dataset_match = re.search(r"The lineage dataset is:\s*([^\s,(]+)", text)
    dataset = dataset_match.group(1).strip() if dataset_match else ""

    # Extract BUSCO version
    version_match = re.search(r"BUSCO version is:\s*([^\s\n]+)", text)
    version = version_match.group(1).strip() if version_match else "unknown"

    # Extract one line summary if present
    # Format: C:97.3%[S:93.3%,D:4.0%],F:0.1%,M:2.7%,n:1658[,E:10.5%]
    summary_match = re.search(
        r"C:([0-9.]+)%\[S:([0-9.]+)%,D:([0-9.]+)%\],F:([0-9.]+)%,M:([0-9.]+)%,n:(\d+)",
        text,
    )

    # Extract raw counts
    s_match = re.search(r"(\d+)\s+Complete and single-copy BUSCOs \(S\)", text)
    d_match = re.search(r"(\d+)\s+Complete and duplicated BUSCOs \(D\)", text)
    f_match = re.search(r"(\d+)\s+Fragmented BUSCOs \(F\)", text)
    m_match = re.search(r"(\d+)\s+Missing BUSCOs \(M\)", text)
    n_match = re.search(r"(\d+)\s+Total BUSCO groups searched", text)

    if not (s_match and d_match and f_match and m_match and n_match):
        # Fallback search if formatting differs slightly
        s_match = s_match or re.search(r"(\d+)\s+Complete and single-copy", text)
        d_match = d_match or re.search(r"(\d+)\s+Complete and duplicated", text)
        f_match = f_match or re.search(r"(\d+)\s+Fragmented", text)
        m_match = m_match or re.search(r"(\d+)\s+Missing", text)
        n_match = n_match or re.search(r"(\d+)\s+Total BUSCO", text)

    if not (s_match and d_match and f_match and m_match and n_match):
        raise ValueError(f"Could not parse BUSCO counts from file: {file_path}")

    s = int(s_match.group(1))
    d = int(d_match.group(1))
    f = int(f_match.group(1))
    m = int(m_match.group(1))
    total = int(n_match.group(1))

    if total == 0:
        total = s + d + f + m

    if summary_match:
        c_pct = float(summary_match.group(1))
        s_pct = float(summary_match.group(2))
        d_pct = float(summary_match.group(3))
        f_pct = float(summary_match.group(4))
        m_pct = float(summary_match.group(5))
        summary_str = summary_match.group(0)
    else:
        s_pct = round(s / total * 100, 1)
        d_pct = round(d / total * 100, 1)
        f_pct = round(f / total * 100, 1)
        m_pct = round(100.0 - s_pct - d_pct - f_pct, 1)
        c_pct = round(s_pct + d_pct, 1)
        summary_str = f"C:{c_pct}%[S:{s_pct}%,D:{d_pct}%],F:{f_pct}%,M:{m_pct}%,n:{total}"

    label = derive_label(path, dataset)

    return BuscoResult(
        label=label,
        file_path=str(path.resolve()),
        dataset=dataset,
        version=version,
        single=s,
        duplicated=d,
        fragmented=f,
        missing=m,
        total=total,
        single_pct=s_pct,
        duplicated_pct=d_pct,
        fragmented_pct=f_pct,
        missing_pct=m_pct,
        complete_pct=c_pct,
        summary_string=summary_str,
    )


def parse_busco_json(file_path: str) -> BuscoResult:
    """Parse a BUSCO short_summary.json file."""
    path = Path(file_path)
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    res = data.get("results", {})
    params = data.get("parameters", {})
    lineage = data.get("lineage_dataset", {})
    versions = data.get("versions", {})

    dataset = lineage.get("name", "")
    version = str(versions.get("busco", "unknown"))

    s = int(res.get("Single copy BUSCOs", 0))
    d = int(res.get("Multi copy BUSCOs", 0))
    f = int(res.get("Fragmented BUSCOs", 0))
    m = int(res.get("Missing BUSCOs", 0))
    total = int(res.get("n_markers", s + d + f + m))

    s_pct = float(res.get("Single copy percentage", round(s / total * 100, 1) if total else 0.0))
    d_pct = float(res.get("Multi copy percentage", round(d / total * 100, 1) if total else 0.0))
    f_pct = float(res.get("Fragmented percentage", round(f / total * 100, 1) if total else 0.0))
    m_pct = float(res.get("Missing percentage", round(m / total * 100, 1) if total else 0.0))
    c_pct = float(res.get("Complete percentage", round(s_pct + d_pct, 1)))

    summary_str = res.get(
        "one_line_summary",
        f"C:{c_pct}%[S:{s_pct}%,D:{d_pct}%],F:{f_pct}%,M:{m_pct}%,n:{total}",
    )

    label = params.get("plot_label") or params.get("out") or derive_label(path, dataset)

    return BuscoResult(
        label=label,
        file_path=str(path.resolve()),
        dataset=dataset,
        version=version,
        single=s,
        duplicated=d,
        fragmented=f,
        missing=m,
        total=total,
        single_pct=s_pct,
        duplicated_pct=d_pct,
        fragmented_pct=f_pct,
        missing_pct=m_pct,
        complete_pct=c_pct,
        summary_string=summary_str,
    )


def derive_label(path: Path, dataset: str = "") -> str:
    """
    Intelligently derive a clean sample label from the filename.
    Standard BUSCO format: short_summary.[generic|specific].<dataset>.<label>.txt
    """
    name = path.stem
    if name.endswith(".txt") or name.endswith(".json"):
        name = Path(name).stem

    # Remove short_summary prefix
    name_clean = re.sub(r"^short_summary\.(?:specific|generic)\.", "", name)

    # If dataset is known, remove it if it prefixes the rest
    if dataset and name_clean.startswith(f"{dataset}."):
        name_clean = name_clean[len(dataset) + 1 :]

    # Fallback to stem if pattern didn't isolate label
    if not name_clean:
        name_clean = path.stem

    return name_clean


def load_busco_file(file_path: str) -> BuscoResult:
    """Load a BUSCO file based on extension (.txt or .json)."""
    p = Path(file_path)
    if not p.is_file():
        raise FileNotFoundError(f"File not found: {file_path}")

    if p.suffix.lower() == ".json":
        return parse_busco_json(file_path)
    else:
        return parse_busco_txt(file_path)


def find_summary_files(
    directory: str = ".",
    run_type: str = "both",
) -> List[Path]:
    """Find all BUSCO short_summary files in a directory."""
    dir_path = Path(directory)
    if not dir_path.is_dir():
        raise NotADirectoryError(f"Directory not found: {directory}")

    patterns = []
    if run_type in ("specific", "both"):
        patterns.extend([
            "short_summary.specific.*.txt",
            "short_summary.specific.*.json",
        ])
    if run_type in ("generic", "both"):
        patterns.extend([
            "short_summary.generic.*.txt",
            "short_summary.generic.*.json",
        ])
    if run_type == "all":
        patterns.extend([
            "short_summary*.txt",
            "short_summary*.json",
        ])

    found_files = []
    for pat in patterns:
        found_files.extend(list(dir_path.glob(pat)))

    # Deduplicate: if both .txt and .json exist for the same base run, keep .txt
    unique_files: Dict[str, Path] = {}
    for f in found_files:
        stem = f.stem
        if stem.endswith(".txt") or stem.endswith(".json"):
            stem = Path(stem).stem
        if stem not in unique_files or f.suffix == ".txt":
            unique_files[stem] = f

    return sorted(list(unique_files.values()))


def plot_busco(
    results: List[BuscoResult],
    output_prefix: str = "busco_figure",
    title: Optional[str] = "BUSCO Assessment Results",
    annotation_type: str = "counts",
    width: Optional[float] = None,
    height: Optional[float] = None,
    bar_height: float = 0.55,
    dpi: int = 300,
    colors: Optional[Dict[str, str]] = None,
    legend_pos: str = "top",
    font_family: str = "sans-serif",
    save_formats: Optional[List[str]] = None,
    quiet: bool = False,
) -> Tuple[str, ...]:
    """
    Generate publication-quality stacked bar plot for BUSCO assessments.
    Saves figure in both PNG and PDF format.
    """
    if not results:
        raise ValueError("No BUSCO results provided to plot.")

    if colors is None:
        colors = PALETTES["colorblind"]

    if save_formats is None:
        save_formats = ["png", "pdf"]

    n_samples = len(results)

    # Dynamic sizing if not explicitly given
    if width is None:
        width = 10.5
    if height is None:
        # Single sample ~ 3.2 inches, each additional adds ~0.55 inches
        height = max(3.0, 1.2 + 0.55 * n_samples)

    # Reverse results so the first sample appears at the top of the horizontal bar chart
    rev_results = list(reversed(results))
    y_positions = np.arange(n_samples)

    # Convert values into percentages array normalized to exactly 100%
    s_arr = np.array([r.single / r.total * 100.0 for r in rev_results])
    d_arr = np.array([r.duplicated / r.total * 100.0 for r in rev_results])
    f_arr = np.array([r.fragmented / r.total * 100.0 for r in rev_results])
    m_arr = np.array([r.missing / r.total * 100.0 for r in rev_results])

    # Y-axis labels: if multiple different datasets are present, append dataset name
    distinct_datasets = {r.dataset for r in results if r.dataset}
    if len(distinct_datasets) > 1:
        y_labels = [f"{r.label} ({r.dataset})" for r in rev_results]
    else:
        y_labels = [r.label for r in rev_results]

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)

    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Helvetica", "Arial"]
    plt.rcParams["font.family"] = font_family

    # Plot stacked horizontal bars
    # S: Single copy
    ax.barh(
        y_positions,
        s_arr,
        height=bar_height,
        color=colors["S"],
        edgecolor="none",
        label=CATEGORY_LABELS["S"],
        zorder=3,
    )
    # D: Duplicated
    ax.barh(
        y_positions,
        d_arr,
        left=s_arr,
        height=bar_height,
        color=colors["D"],
        edgecolor="none",
        label=CATEGORY_LABELS["D"],
        zorder=3,
    )
    # F: Fragmented
    ax.barh(
        y_positions,
        f_arr,
        left=s_arr + d_arr,
        height=bar_height,
        color=colors["F"],
        edgecolor="none",
        label=CATEGORY_LABELS["F"],
        zorder=3,
    )
    # M: Missing
    ax.barh(
        y_positions,
        m_arr,
        left=s_arr + d_arr + f_arr,
        height=bar_height,
        color=colors["M"],
        edgecolor="none",
        label=CATEGORY_LABELS["M"],
        zorder=3,
    )

    # Bar annotation text (responsive font size)
    base_font_size = max(7.5, min(10.0, 9.5 * (height / (n_samples + 2))))

    for y_pos, r in zip(y_positions, rev_results):
        if annotation_type == "none":
            break

        c_total = r.single + r.duplicated
        if annotation_type == "percentages":
            text_str = (
                f"C:{r.complete_pct}% [S:{r.single_pct}%, D:{r.duplicated_pct}%], "
                f"F:{r.fragmented_pct}%, M:{r.missing_pct}%, n:{r.total}"
            )
        elif annotation_type == "both":
            text_str = (
                f"C:{r.complete_pct}% ({c_total}) [S:{r.single_pct}% ({r.single}), "
                f"D:{r.duplicated_pct}% ({r.duplicated})], F:{r.fragmented_pct}% ({r.fragmented}), "
                f"M:{r.missing_pct}% ({r.missing}), n:{r.total}"
            )
        else:  # counts (classic BUSCO standard)
            text_str = (
                f"C:{c_total} [S:{r.single}, D:{r.duplicated}], "
                f"F:{r.fragmented}, M:{r.missing}, n:{r.total}"
            )

        text_obj = ax.text(
            2.5,
            y_pos,
            text_str,
            ha="left",
            va="center",
            fontsize=base_font_size,
            fontweight="bold",
            color="#111111",
            zorder=4,
        )
        text_obj.set_path_effects([
            pe.withStroke(linewidth=2.5, foreground="white", alpha=0.90)
        ])

    # X-axis
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xticklabels(["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=10)
    ax.set_xlabel("% BUSCOs", fontsize=11, fontweight="bold", labelpad=8)

    # Y-axis
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=11, fontweight="bold")
    ax.set_ylim(-0.6, n_samples - 0.4)

    # Gridlines
    ax.xaxis.grid(True, linestyle="--", linewidth=0.7, color="#CCCCCC", alpha=0.7, zorder=0)
    ax.yaxis.grid(False)

    # Spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_color("#222222")
    ax.tick_params(axis="x", which="major", length=5, width=1.0, color="#222222")
    ax.tick_params(axis="y", which="major", length=0)

    # Legend
    if legend_pos != "none":
        handles, leg_labels = ax.get_legend_handles_labels()
        if legend_pos == "top":
            ncol = 4 if width >= 9.0 else 2
            ax.legend(
                handles,
                leg_labels,
                loc="lower center",
                bbox_to_anchor=(0.5, 1.02),
                ncol=ncol,
                frameon=False,
                fontsize=8.5,
                columnspacing=1.0,
                handletextpad=0.4,
            )
            if title:
                ax.set_title(title, fontsize=13, fontweight="bold", pad=32)
        elif legend_pos == "bottom":
            ncol = 4 if width >= 9.0 else 2
            ax.legend(
                handles,
                leg_labels,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.18),
                ncol=ncol,
                frameon=False,
                fontsize=8.5,
                columnspacing=1.0,
                handletextpad=0.4,
            )
            if title:
                ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
        elif legend_pos == "right":
            ax.legend(
                handles,
                leg_labels,
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                frameon=False,
                fontsize=8.5,
            )
            if title:
                ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    else:
        if title:
            ax.set_title(title, fontsize=13, fontweight="bold", pad=12)

    plt.tight_layout()

    # Save output files
    saved_files = []
    out_base = Path(output_prefix)
    out_base.parent.mkdir(parents=True, exist_ok=True)

    if out_base.suffix.lower() in (".png", ".pdf", ".svg"):
        out_base = out_base.with_suffix("")

    for fmt in save_formats:
        fmt_clean = fmt.lower().lstrip(".")
        target_file = out_base.with_suffix(f".{fmt_clean}")
        fig.savefig(
            target_file,
            dpi=dpi,
            bbox_inches="tight",
            format=fmt_clean,
        )
        saved_files.append(str(target_file.resolve()))
        if not quiet:
            print(f"[BUSCO Plot] Saved figure: {target_file}")

    plt.close(fig)
    return tuple(saved_files)


def main():
    parser = argparse.ArgumentParser(
        description="BUSCO Assessment Plotting Tool (Pure Python - No R required)\n"
                    "Generates publication-ready figures in both PNG and PDF format\n"
                    "using colorblind-safe palettes (Okabe-Ito / CUD by default).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Input options
    input_group = parser.add_argument_group("Input options")
    input_group.add_argument(
        "-wd",
        "--working_directory",
        metavar="PATH",
        default=None,
        help="Directory containing BUSCO short_summary files (default: current directory).",
    )
    input_group.add_argument(
        "-i",
        "--input",
        nargs="+",
        metavar="FILE",
        help="One or more specific BUSCO short_summary files (txt or json).",
    )
    input_group.add_argument(
        "-rt",
        "--run_type",
        choices=["generic", "specific", "both", "all"],
        default="both",
        help="Type of summary to search for when scanning directory (default: both).",
    )
    input_group.add_argument(
        "--labels",
        nargs="+",
        metavar="NAME",
        help="Custom label(s) for the samples (must match number of input files).",
    )

    # Output options
    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        "-o",
        "--output",
        metavar="PREFIX",
        default="busco_figure",
        help="Output file prefix/path (default: busco_figure). Creates <PREFIX>.png and <PREFIX>.pdf.",
    )
    output_group.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Resolution for raster PNG output in dots per inch (default: 300).",
    )
    output_group.add_argument(
        "--svg",
        action="store_true",
        help="Also export vector SVG figure.",
    )

    # Palette and plot customization
    style_group = parser.add_argument_group("Plot customization")
    style_group.add_argument(
        "-p",
        "--palette",
        choices=["colorblind", "okabe-ito", "tol", "legacy"],
        default="colorblind",
        help="Color palette to use (default: 'colorblind' [Okabe-Ito CUD]).\n"
             "  - colorblind / okabe-ito: S=#56B4E9, D=#0072B2, F=#F0E442, M=#D55E00\n"
             "  - tol: Paul Tol qualitative colorblind palette\n"
             "  - legacy: Original BUSCO palette (#56B4E9, #3492C7, #F0E442, #F04442)",
    )
    style_group.add_argument(
        "-t",
        "--title",
        metavar="TEXT",
        default="BUSCO Assessment Results",
        help="Plot title (default: 'BUSCO Assessment Results').",
    )
    style_group.add_argument(
        "--no-title",
        action="store_true",
        help="Omit title from plot (useful for manuscript figures).",
    )
    style_group.add_argument(
        "--annotation",
        choices=["counts", "percentages", "both", "none"],
        default="counts",
        help="Summary text displayed on bars (default: counts, classic BUSCO format).",
    )
    style_group.add_argument(
        "--bar-height",
        type=float,
        default=0.55,
        help="Relative height/thickness of each bar (0.1 to 1.0, default: 0.55).",
    )
    style_group.add_argument(
        "--width",
        type=float,
        default=None,
        help="Figure width in inches (default: 10.5).",
    )
    style_group.add_argument(
        "--height",
        type=float,
        default=None,
        help="Figure height in inches (default: auto-calculated based on sample count).",
    )
    style_group.add_argument(
        "--legend",
        choices=["top", "bottom", "right", "none"],
        default="top",
        help="Legend position (default: top).",
    )
    style_group.add_argument(
        "--sort-by",
        choices=["input", "name", "complete", "single", "missing"],
        default="input",
        help="Order samples on plot: 'input', 'name', 'complete', 'single', or 'missing'.",
    )
    style_group.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress informational console output.",
    )

    args = parser.parse_args()

    # Collect summary files
    summary_files: List[Path] = []
    if args.input:
        for item in args.input:
            matched = glob.glob(item)
            if matched:
                summary_files.extend([Path(m) for m in matched])
            else:
                summary_files.append(Path(item))
    else:
        scan_dir = args.working_directory if args.working_directory else "."
        summary_files = find_summary_files(scan_dir, run_type=args.run_type)

    if not summary_files:
        print("Error: No BUSCO summary files found.", file=sys.stderr)
        if args.working_directory:
            print(f"Checked directory: {args.working_directory}", file=sys.stderr)
        sys.exit(1)

    # Parse summary files
    results: List[BuscoResult] = []
    for f in summary_files:
        try:
            res = load_busco_file(str(f))
            results.append(res)
            if not args.quiet:
                print(f"[BUSCO Plot] Loaded: {f.name} -> {res.label} ({res.summary_string})")
        except Exception as e:
            print(f"Warning: Failed to parse {f}: {e}", file=sys.stderr)

    if not results:
        print("Error: Could not parse any valid BUSCO results from files.", file=sys.stderr)
        sys.exit(1)

    # Custom labels
    if args.labels:
        if len(args.labels) != len(results):
            print(
                f"Error: Number of custom labels ({len(args.labels)}) does not match "
                f"number of samples ({len(results)}).",
                file=sys.stderr,
            )
            sys.exit(1)
        for r, lab in zip(results, args.labels):
            r.label = lab

    # Sorting
    if args.sort_by == "name":
        results.sort(key=lambda r: r.label.lower())
    elif args.sort_by == "complete":
        results.sort(key=lambda r: (r.complete_pct, r.single_pct), reverse=True)
    elif args.sort_by == "single":
        results.sort(key=lambda r: r.single_pct, reverse=True)
    elif args.sort_by == "missing":
        results.sort(key=lambda r: r.missing_pct)

    plot_title = None if args.no_title else args.title

    formats = ["png", "pdf"]
    if args.svg:
        formats.append("svg")

    out_prefix = args.output
    if args.working_directory and not os.path.isabs(out_prefix) and not os.path.dirname(out_prefix):
        out_prefix = os.path.join(args.working_directory, out_prefix)

    # Selected palette
    selected_colors = PALETTES.get(args.palette, PALETTES["colorblind"])

    # Generate plot
    plot_busco(
        results=results,
        output_prefix=out_prefix,
        title=plot_title,
        annotation_type=args.annotation,
        width=args.width,
        height=args.height,
        bar_height=args.bar_height,
        dpi=args.dpi,
        colors=selected_colors,
        legend_pos=args.legend,
        save_formats=formats,
        quiet=args.quiet,
    )

    if not args.quiet:
        print("\nPlot generation complete!")
        print(f"  Palette: {args.palette} (S: {selected_colors['S']}, D: {selected_colors['D']}, F: {selected_colors['F']}, M: {selected_colors['M']})")
        print(f"  PNG: {out_prefix}.png")
        print(f"  PDF: {out_prefix}.pdf")
        if args.svg:
            print(f"  SVG: {out_prefix}.svg")


if __name__ == "__main__":
    main()
