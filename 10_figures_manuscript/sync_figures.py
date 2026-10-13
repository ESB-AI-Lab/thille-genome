#!/usr/bin/env python3
"""Keep /gpool/edwin/avocados/thille_paper/figures/ holding only the current version of every figure (Edwin, 2026-10-05).

One subfolder per figure (Figure_1 ... Figure_S12; S12 = inversion breakpoints); each panel is copied as PDF and PNG with a
name giving its order and panel (e.g. Figure_6/01_SDMYB.pdf). Sources are the rows of manuscript/figure_paths_v4.tsv
(scripts/figure_paths_v4.py) plus the inversion breakpoint figures. Files are replaced when the source differs, and any
PDF/PNG in figures/ that is not a current figure is removed, so the folder never holds deprecated versions. Each figure's
legend md (manuscript/figure_legends/Figure_<id>.md; Edwin, 2026-10-09: figure text lives in the legend) is copied into
its folder as Figure_<id>/Figure_<id>.md.
Usage (from the repo root): python3 scripts/figure_paths_v4.py && python3 scripts/sync_figures.py
"""
import csv
import filecmp
import re
import shutil
import sys
from pathlib import Path

REPO = Path("/gpool/edwin/avocados/thille_paper")
OUT = REPO / "figures"
INV = Path("/epool/edwin/pamericana/avoflow/assemblies/thille.hic/compare/v4/inversions/figures_v2_notext")  # no-text renders (2026-10-09; figures_v2 kept)
INV_NAMES = {"hass": "Hass_UoQ", "westindian": "Unnamed_Lowland", "puebla": "Puebla"}


def clean(s):
    return re.sub(r"[^A-Za-z0-9.-]+", "_", s).strip("_")


def main():
    targets = {}  # destination -> source
    order = {}
    for r in csv.DictReader(open(REPO / "manuscript/figure_paths_v4.tsv"), delimiter="\t"):
        if not r["item"].startswith("Figure"):
            continue
        folder = r["item"].replace("Figure ", "Figure_")
        order[folder] = order.get(folder, 0) + 1
        panel = clean(r["panel"].replace("(crop pending)", "uncropped")) or folder
        stem = f"{order[folder]:02d}_{panel}"
        src = Path(r["path"])
        for ext in (".pdf", ".png"):
            targets[OUT / folder / (stem + ext)] = src.with_suffix(ext)
    for d, name in INV_NAMES.items():
        for src in sorted((INV / d).glob("*.p[dn][fg]")):
            targets[OUT / "Figure_S12" / f"{name}_{src.name.split('_', 1)[1]}"] = src
    for src in sorted(INV.glob("overview_split_fraction.p[dn][fg]")):
        targets[OUT / "Figure_S12" / ("00_" + src.name)] = src
    for md in sorted((REPO / "manuscript/figure_legends").glob("Figure_*.md")):
        targets[OUT / md.stem / md.name] = md
    missing = [s for s in targets.values() if not s.is_file()]
    if missing:
        sys.exit("missing sources:\n" + "\n".join(map(str, missing)))
    copied = same = 0
    for dst, src in targets.items():
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.is_file() and filecmp.cmp(src, dst, shallow=False):
            same += 1
            continue
        shutil.copy2(src, dst)
        copied += 1
    removed = []
    for f in OUT.rglob("*"):
        if f.is_file() and f.suffix in (".pdf", ".png", ".md") and f not in targets:
            f.unlink()  # Edwin: no deprecated figures in figures/ (the sources stay where they are)
            removed.append(f)
    for d in sorted(OUT.iterdir()):
        if d.is_dir() and not any(d.iterdir()):
            d.rmdir()
    print(f"{OUT}: {len(targets)} files ({copied} copied or replaced, {same} unchanged, {len(removed)} removed)")
    for f in removed:
        print(f"  removed {f}")


if __name__ == "__main__":
    main()
