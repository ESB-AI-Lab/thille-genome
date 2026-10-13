#!/usr/bin/env python3
"""Dot input for render_karyotype.py --dots: one line per Thille gene locus and set.

Sets: fst = rows of the Fst flowering table plus rows of the Gwen-paper table in a top 1% Fst set (contig or
scaffold); coverage = rows of the coverage-selected flowering table; florid = table rows whose Gwen gene has a
FLOR-ID match in the candidate screen (T1_florid). Rows sharing a Thille locus within a set give one dot.

Usage: make_karyotype_dots.py > fst_check/karyotype_dots.tsv
"""
import csv
import re

T = "fst_check/tables"
CANDS = "/epool/edwin/pamericana/gwen/annotation/fstpeaksearch/results/final_candidates/flowering_candidates_all_v4.csv"
LOC = re.compile(r"(chr\d+):([\d,]+)-([\d,]+)")


def rows(name):
    return list(csv.DictReader(open(f"{T}/{name}.tsv"), delimiter="\t"))


def main():
    florid = {r["ver1_id"] for r in csv.DictReader(open(CANDS)) if r["T1_florid"].strip()}
    sel = [("fst", r) for r in rows("fst_flowering_v2")]
    sel += [("fst", r) for r in rows("coverage_gwen_v2") if r["Top 1% Fst set"] != "neither"]
    sel += [("coverage", r) for r in rows("no_fst_flowering_v2")]
    allrows = rows("fst_flowering_v2") + rows("no_fst_flowering_v2") + rows("coverage_gwen_v2")
    sel += [("florid", r) for r in allrows if set(re.findall(r"Chr\w+\.ver1\.g\d+", r["Gwen gene (ver1 / jg)"])) & florid]
    seen = set()
    print("chrom\tstart\tend\tgene\tset")
    for s, r in sel:
        m = LOC.match(r["Thille location"])
        c, a, b = m.group(1), int(m.group(2).replace(",", "")) - 1, int(m.group(3).replace(",", ""))
        if (s, c, a, b) in seen:
            continue
        seen.add((s, c, a, b))
        gene = re.sub(r"^(T2|hemi) ", "", r["SwissProt gene"]).split(":")[0].split(" (")[0]
        print(f"{c}\t{a}\t{b}\t{gene}\t{s}")


if __name__ == "__main__":
    main()
