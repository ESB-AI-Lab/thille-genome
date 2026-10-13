"""Convert the manuscript docx into a bioRxiv LaTeX source: main.tex, citations.bib, and figures/ (source PDFs).

The text, captions, and tables come from the docx; citations become natbib \\citep/\\citet keys built from the
manuscript's own reference list; figure panels are the source PDFs of the images in the docx (matched by image
content, see FIGURES), copied into figures/. Panels with no PDF are listed in MISSING and drawn as a framed
placeholder until a file is given.

Usage: python3 build_latex.py <manuscript.docx> <out_dir>
"""
import csv, html, re, shutil, sys, unicodedata, zipfile
from pathlib import Path

DOCX, OUT = sys.argv[1], Path(sys.argv[2])
REPO = Path("/gpool/edwin/avocados/thille_paper")
TH = Path("/epool/edwin/pamericana/avoflow/assemblies/thille.hic")
INV_OVERVIEW = TH / "compare/v4/inversions/figures_v2_notext/overview_split_fraction.pdf"  # no-text render (2026-10-09)
MS = REPO / "manuscript"  # Figures 1, 2 and 3B: PDFs provided by Edwin (2026-10-02), stored only in figures/

xml = zipfile.ZipFile(DOCX).read("word/document.xml").decode("utf-8")
body = xml[xml.index("<w:body>") + 8:]
BLOCKS = re.findall(r"<w:tbl>.*?</w:tbl>|<w:p[ >].*?</w:p>|<w:p/>", body, flags=re.S)


def plain(b):
    return html.unescape("".join(re.findall(r"<w:t(?: [^>]*)?>([^<]*)</w:t>", b)))


def find_block(prefix, start=0):
    for i in range(start, len(BLOCKS)):
        if not BLOCKS[i].startswith("<w:tbl") and plain(BLOCKS[i]).strip().startswith(prefix):
            return i
    raise SystemExit(f"block not found: {prefix}")


# ------------------------------------------------------------------ text -> LaTeX
CHARS = {"×": "x",  # ASCII x (Edwin, 2026-10-02)
         "\u2009": r"\,",  # thin space between initials (E.\,A.\,S.)
         "′": r"$'$", "≤": r"$\leq$", "≥": r"$\geq$", "°": r"\textdegree{}", "β": r"$\beta$",
         "σ": r"$\sigma$", "’": "'", "‘": "`", "“": "``", "”": "''", "\u00a0": "~", "µ": r"$\mu$", "±": r"$\pm$"}
SPECIAL = {"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "$": r"\$", "&": r"\&", "#": r"\#", "%": r"\%",
           "_": r"\_", "^": r"\^{}", "~": r"\textasciitilde{}", "<": r"\textless{}", ">": r"\textgreater{}"}
unknown_chars = set()


def esc(text):
    out = []
    for ch in text:
        if ch in SPECIAL:
            out.append(SPECIAL[ch])
        elif ch in CHARS:
            out.append(CHARS[ch])
        elif ord(ch) > 127 and unicodedata.category(ch)[0] != "L":
            unknown_chars.add(ch); out.append(ch)
        else:
            out.append(ch)
    return "".join(out)


def on(rpr, tag):
    m = re.search(rf"<w:{tag}(?: w:val=\"([^\"]*)\")?/>", rpr)
    return bool(m) and m.group(1) not in ("0", "false")


def segments(p):
    """[(text, bold, italic, vert)] for the runs of a paragraph; '\n' marks a line break."""
    segs = []
    for r in re.findall(r"<w:r(?: [^>]*)?>.*?</w:r>", p, flags=re.S):
        rpr = (re.search(r"<w:rPr>.*?</w:rPr>", r, flags=re.S) or [""])[0]
        if isinstance(rpr, re.Match):
            rpr = rpr.group(0)
        b, i = on(rpr, "b"), on(rpr, "i")
        v = (re.search(r'<w:vertAlign w:val="(\w+)"', rpr) or [None, None])[1]
        for tok in re.findall(r"<w:t(?: [^>]*)?>[^<]*</w:t>|<w:br[^>]*/>|<w:tab/>", r):
            if tok.startswith("<w:br"):
                segs.append(("\n", b, i, v))
            elif tok == "<w:tab/>":
                segs.append((" ", b, i, v))
            else:
                segs.append((html.unescape(re.sub(r"<[^>]+>", "", tok)), b, i, v))
    merged = []
    for s in segs:
        if merged and merged[-1][1:] == s[1:] and s[0] != "\n" and merged[-1][0] != "\n":
            merged[-1] = (merged[-1][0] + s[0],) + s[1:]
        else:
            merged.append(s)
    return merged


def latex(segs, bold=True):
    out = []
    for text, b, i, v in segs:
        if text == "\n":
            out.append("\n"); continue
        t = esc(text)
        if not t.strip():
            out.append(t); continue
        if v == "superscript":
            t = rf"\textsuperscript{{{t}}}"
        elif v == "subscript":
            t = rf"\textsubscript{{{t}}}"
        if i:
            t = rf"\textit{{{t}}}"
        if b and bold:
            t = rf"\textbf{{{t}}}"
        out.append(t)
    s = "".join(out)
    s = re.sub(r"(https?://[^\s,;)]+[^\s,;).])", lambda m: r"\url{" + m.group(1).replace("\\_", "_").replace("\\%", "%").replace("\\#", "#").replace("\\&", "&") + "}", s)
    return s


# ------------------------------------------------------------------ references -> BibTeX
def ascii_fold(s):
    return unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode()


def bib_author(name):
    name = name.strip()
    if name == "et al.":
        return "others"
    if name == "UniProt Consortium":
        return "{UniProt Consortium}"
    parts = name.split()
    if len(parts) >= 2 and re.fullmatch(r"[A-Z]{1,4}", parts[-1]):
        initials = " ".join(f"{c}." for c in parts[-1])
        return f"{' '.join(parts[:-1])}, {initials}"
    return "{" + name + "}"


def bib_escape(s):
    return s.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")


refs = []
r0 = find_block("References")
for b in BLOCKS[r0 + 1:find_block("Supplementary Figures", r0)]:
    t = plain(b).strip()
    if not t:
        continue
    if "Cabanettes F" in t:   # D-GENIES was not used in this study (Edwin, 2026-10-02)
        continue
    m = re.match(r"^(.*?)\.? ((?:19|20)\d\d(?:-(?:19|20)\d\d)?[a-z]?)\. (.*)$", t)   # 2025a, 2025b (v37)
    if not m:
        raise SystemExit(f"cannot parse reference: {t}")
    auth_s, year, rest = m.groups()
    parts = [a.strip() for a in re.split(r", | et al\.?$", auth_s) if a.strip()]
    authors = [a for a in parts if a.rstrip(".") != "et al"]   # "Smith NJ, et al." must not become an author
    if re.search(r"(?:^|,| ) ?et al\.?$", auth_s) or len(authors) < len(parts):
        authors.append("et al.")
    surname = authors[0].rsplit(" ", 1)[0] if authors[0] != "UniProt Consortium" else "UniProt Consortium"
    n = len([a for a in authors if a != "et al."])
    if surname == "UniProt Consortium":
        label = "UniProt Consortium"
    elif "et al." in authors or n >= 3:
        label = f"{surname} et al."
    elif n == 2:
        label = f"{surname} and {authors[1].rsplit(' ', 1)[0]}"
    else:
        label = surname
    refs.append(dict(text=t, authors=authors, surname=surname, year=year, rest=rest, label=f"{label} {year}"))

keys = {}
for r in refs:
    base = re.sub(r"[^a-z]", "", ascii_fold(r["surname"]).lower()) + r["year"][:4]
    k = base
    n = 1
    while k in keys.values():
        n += 1; k = f"{base}{chr(96 + n)}"
    keys[r["label"]] = k
    r["key"] = k


def bib_entry(r):
    rest = r["rest"]
    doi = None
    m = re.search(r"\s*doi:(\S+?)\.?$", rest)
    if m:
        doi = m.group(1); rest = rest[:m.start()].rstrip()
    url = None
    m = re.search(r"\s*(https?://\S+?)\.?$", rest)   # references without a DOI that end in a URL (e.g. Delphey 1960)
    if m and not doi:
        url = m.group(1); rest = rest[:m.start()].rstrip()
    authors = " and ".join(bib_author(a) for a in r["authors"])
    f = {"author": authors, "year": r["year"].replace("-", "--")}
    if doi:
        f["doi"] = doi
    if url:
        f["url"] = url
    k = r["key"]
    if rest.startswith("Aligning sequence reads"):          # arXiv preprint
        title, tail = rest.split(". ", 1)
        f.update(title="{" + title + "}", eprint="1303.3997", archivePrefix="arXiv", primaryClass="q-bio.GN",
                 howpublished="arXiv:1303.3997 [q-bio.GN]")
        typ = "misc"
    elif ". In: " in rest:                                   # conference papers
        title, tail = rest.split(". In: ", 1)
        book, pages = re.match(r"(.*?)\. p\. (\S+?)\.?$", tail).groups()
        f.update(title="{" + title + "}", booktitle=book, pages=pages.replace("-", "--"))
        typ = "inproceedings"
    elif rest.startswith("Manual of tropical"):
        f.update(title="{Manual of tropical and subtropical fruits}", publisher="Macmillan", address="New York (NY)",
                 note="Facsimile edition, New York (NY): Hafner Press, 1974")
        typ = "book"
    elif re.search(r"\. bioRxiv \d+\.?$", rest):              # bioRxiv preprints (e.g. Dudchenko et al. 2018)
        title, num = re.match(r"^(.*?)\. bioRxiv (\d+)\.?$", rest).groups()
        f.update(title="{" + title.replace("$", r"\$") + "}", howpublished=f"bioRxiv {num}", note="Preprint")
        typ = "misc"
    elif rest.startswith("RepeatMasker"):
        f.update(title="{RepeatMasker Open-4.0}", howpublished=r"\url{http://www.repeatmasker.org}",
                 note="Software")
        typ = "misc"
    else:
        m = re.match(r"^(.*?)\. (.*?)\.? (\d+)(?:\(([^)]+)\))?:(\S+?)\.?$", rest)
        if not m:
            raise SystemExit(f"cannot parse journal part of {k}: {rest}")
        title, journal, vol, num, pages = m.groups()
        f.update(title="{" + title + "}", journal=journal, volume=vol, pages=pages.replace("-", "--"))
        if num:
            f["number"] = num
        typ = "article"
    lines = [f"@{typ}{{{k},"] + [f"  {a} = {{{bib_escape(v) if a not in ('howpublished', 'url') else v}}}," for a, v in f.items()]
    return "\n".join(lines) + "\n}\n"


# ------------------------------------------------------------------ citations in text
LABELS = sorted(keys, key=len, reverse=True)
unresolved = set()
NAME = r"[A-Z][A-Za-zÀ-ÿ'\-]+(?: [A-Z][A-Za-zÀ-ÿ'\-]+)?"


def label_key(s):
    s = s.strip()
    return keys.get(s)


def cite(s):
    # narrative: "Groh et al. (2026)"
    def narr(m):
        k = label_key(f"{m.group(1)} {m.group(2)}")
        if k:
            return rf"\citet{{{k}}}"
        unresolved.add(m.group(0)); return m.group(0)
    s = re.sub(rf"((?:{NAME}) et al\.|{NAME} and {NAME}|{NAME}) \(((?:19|20)\d\d)\)", narr, s)

    def paren(m):
        inner = m.group(1)
        parts = inner.split("; ")
        out, all_cites = [], True
        for p in parts:
            k = label_key(p)
            if k:
                out.append(("c", k))
            else:
                all_cites = False
                if re.search(r"(?:et al\.|[A-Z][a-z]+) (?:19|20)\d\d", p):
                    unresolved.add(p)
                out.append(("t", p))
        if all_cites:
            return r"\citep{" + ",".join(k for _t, k in out) + "}"
        if not any(t == "c" for t, _ in out):
            return m.group(0)
        return "(" + "; ".join(rf"\citealt{{{v}}}" if t == "c" else v for t, v in out) + ")"
    s = re.sub(r"\(([^()]*)\)", paren, s)
    for lab in LABELS:                                   # bare labels left in running text
        s = re.sub(rf"(?<![\w{{]){re.escape(lab)}(?!\d)", rf"\\citealt{{{keys[lab]}}}", s)
    return s


# ------------------------------------------------------------------ body
SECTIONS = {"Abstract", "Introduction", "Methods", "Results", "Discussion", "Data Availability", "Funding", "Acknowledgments", "Summary:"}


def paragraph(b):
    segs = segments(b)
    text = "".join(s[0] for s in segs).strip()
    if not text:
        return ""
    if text.rstrip(":") in {x.rstrip(":") for x in SECTIONS} and len(text) < 30:
        return f"\\section*{{{esc(text.rstrip(':'))}}}\n"
    # run-in bold heading followed by a line break
    if segs and segs[0][1] and any(s[0] == "\n" for s in segs):
        n = next(i for i, s in enumerate(segs) if s[0] == "\n")
        if all(s[1] for s in segs[:n]):
            head = "".join(s[0] for s in segs[:n]).strip()
            rest = latex([s for s in segs[n + 1:]], bold=False).replace("\n", " ").strip()
            return f"\\subsection*{{{esc(head)}}}\n{cite(rest)}\n"
    if all(s[1] for s in segs if s[0].strip()) and len(text) < 120:
        return f"\\subsection*{{{esc(text)}}}\n"
    return cite(latex(segs).replace("\n", " \\\\\n").strip()) + "\n"


def table(b, landscape=False):
    widths = [float(w) for w in re.findall(r'<w:gridCol w:w="([^"]+)"', b)]
    total = sum(widths)
    avail = 0.98 * (23 if landscape else 16.5)  # cm of text width
    cols = "".join(rf">{{\raggedright\arraybackslash}}p{{{avail * w / total - 0.2:.2f}cm}}" for w in widths)
    rows = []
    for tr in re.findall(r"<w:tr[ >].*?</w:tr>", b, flags=re.S):
        cells = []
        for tc in re.findall(r"<w:tc>.*?</w:tc>", tr, flags=re.S):
            span = int((re.search(r'<w:gridSpan w:val="(\d+)"', tc) or [0, 1])[1])
            content = r" \newline ".join(cite(latex(segments(p), bold=True)).replace("\n", " ").strip()
                                         for p in re.findall(r"<w:p[ >].*?</w:p>", tc, flags=re.S)
                                         if plain(p).strip())
            cells.append(content if span == 1 else rf"\multicolumn{{{span}}}{{l}}{{{content}}}")
        rows.append(" & ".join(cells) + r" \\")
    head = rows[0]
    body_rows = "\n".join(rows[1:])
    return (f"\\begin{{longtable}}{{{cols}}}\n\\toprule\n{head}\n\\midrule\n\\endfirsthead\n\\toprule\n{head}\n"
            f"\\midrule\n\\endhead\n{body_rows}\n\\bottomrule\n\\end{{longtable}}\n")


# ------------------------------------------------------------------ figures: the v4 sources listed in figure_paths_v4.tsv
# (scripts/figure_paths_v4.py; the same files as /gpool/edwin/avocados/thille_paper/figures/). Gene figures are laid out
# two panels per row, labelled A, B, ... per row as in the docx; Figure 5 keeps its three SyRI overlays in row A.
MISSING = {}  # name -> description, for panels without a source PDF
FP = {}
for r in csv.DictReader(open(REPO / "manuscript/figure_paths_v4.tsv"), delimiter="\t"):
    if r["item"].startswith("Figure"):
        FP.setdefault(r["item"], []).append((r["panel"], Path(r["path"])))


def fname(item, k, panel, path):
    tag = re.sub(r"[^A-Za-z0-9.-]+", "_", f"{item.replace('Figure ', 'fig')}_{k:02d}_{panel}").strip("_")
    return f"{tag}{path.suffix}"


def pairs(item, labelled=True):
    ps = [(p, fname(item, k, lab, p)) for k, (lab, p) in enumerate(FP[item], start=1)]
    rows = [ps[i:i + 2] for i in range(0, len(ps), 2)]
    return [[(chr(65 + n) if labelled and c == 0 else "", ps_) for c, ps_ in enumerate(row)] for n, row in enumerate(rows)]


def single(item, labels):
    ps = [(p, fname(item, k, lab, p)) for k, (lab, p) in enumerate(FP[item], start=1)]
    return [[(lab, ps[i])] for i, lab in enumerate(labels)]


f5 = [(p, fname("Figure 5", k, lab, p)) for k, (lab, p) in enumerate(FP["Figure 5"], start=1)]
f6 = [(p, fname("Figure 6", k, lab, p)) for k, (lab, p) in enumerate(FP["Figure 6"], start=1)]
FIGURES = {
    "Figure 1": single("Figure 1", [""]),
    "Figure 2": single("Figure 2", [""]),
    "Figure 3": [[("A", (FP["Figure 3"][0][1], fname("Figure 3", 1, "A", FP["Figure 3"][0][1]))),
                  ("B", (FP["Figure 3"][1][1], fname("Figure 3", 2, "B", FP["Figure 3"][1][1])))]],
    "Figure 4": single("Figure 4", [""]),
    "Figure 5": [[("A", f5[0]), ("", f5[1]), ("", f5[2])], [("B", f5[3]), ("C", f5[4])]],
    "Figure 6": [[("A", f6[0])], [("B", f6[1]), ("", f6[2])], [("C", f6[3]), ("", f6[4])], [("D", f6[5]), ("", f6[6])]],
    "Figure 7": pairs("Figure 7"),
    "Supplementary Figure S2": single("Figure S2", [""]),
    "Supplementary Figures S5 and S6": pairs("Figure S5", False) + pairs("Figure S6", False),
    "Supplementary Figure S7": single("Figure S7", [""]),
    "Supplementary Figure S8": single("Figure S8", [""]),
    "Supplementary Figure S9": single("Figure S9", ["A", "B"]),
    "Supplementary Figure S10": single("Figure S10", ["A", "B"]),
    "Supplementary Figure S11": pairs("Figure S11", False),
    # S12: the overview only; the 152 breakpoint figures are supplementary files (figures/Figure_S12/)
    "Supplementary Figure S12": [[("", (INV_OVERVIEW, "figS12_overview_split_fraction.pdf"))]],
    "Supplementary Figure S13": pairs("Figure S13", False),
    "Supplementary Figure S14": [[("A", None), ("B", None)], [("C", None), ("D", None)]],
    "Supplementary Figure S15": single("Figure S15", ["A", "B", "C"]),
}
s14 = [(p, fname("Figure S14", k, lab, p)) for k, (lab, p) in enumerate(FP["Figure S14"], start=1)]
FIGURES["Supplementary Figure S14"] = [[("A", s14[0]), ("B", s14[1])], [("C", s14[2]), ("D", s14[3])]]

figdir = OUT / "figures"
figdir.mkdir(parents=True, exist_ok=True)
copied = []
for rows in FIGURES.values():
    for row in rows:
        for _lab, (path, name) in row:
            if path is None:
                continue
            if not path.exists():
                raise SystemExit(f"figure source missing: {path}")
            dst = figdir / name
            if dst.exists() and dst.read_bytes() != path.read_bytes():
                raise SystemExit(f"refusing to overwrite a different {dst}")
            if not dst.exists():
                shutil.copy2(path, dst)
            copied.append((str(path), name))


def figure_env(key, caption_tex, max_rows=4, ncol=None, legend_page=True):
    """One float per page; figures with more than max_rows panel rows continue on the next page. A legend longer than
    about half a page goes on its own page before the panels, which then carry the figure name only."""
    rows = FIGURES[key]
    if legend_page and len(re.sub(r"\\[a-zA-Z]+\*?|[{}]", "", caption_tex)) > 2400:
        head = re.match(r"\\textbf\{([^}]*):\}", caption_tex)
        name = head.group(1) if head else key
        return ("\\begin{figure}[p]\n" + rf"\caption*{{{caption_tex}}}" + "\n\\end{figure}\n"
                + figure_env(key, rf"\textbf{{{name}}} (legend on the previous page).", max_rows, ncol, legend_page=False))
    ncol = ncol or max(len(r) for r in rows)  # a short last row keeps the panel size and leaves its neighbor empty
    if len(rows) > max_rows:
        chunks = [rows[i:i + max_rows] for i in range(0, len(rows), max_rows)]
        head = re.match(r"\\textbf\{([^}]*):\}", caption_tex)
        out = []
        for n, chunk in enumerate(chunks):
            FIGURES[f"{key} part {n}"] = chunk
            cap = caption_tex if n == 0 else rf"\textbf{{{head.group(1) if head else key} (continued).}}"
            out.append(figure_env(f"{key} part {n}", cap, max_rows, ncol, legend_page=False))
        return "\n".join(out)
    # panels share the page with the legend (v38 legends carry the text moved out of the figures): estimate the caption
    # height (small font, ~90 characters per 12.5 pt line, text height ~650 pt) and give the panels the rest
    chars = len(re.sub(r"\\[a-zA-Z]+\*?|[{}]", "", caption_tex))
    cap_frac = (chars / 90 + 1) * 12.5 / 650
    row_h = max(0.06, min(0.42, (0.93 - cap_frac - len(rows) * 0.03) / len(rows)))
    out = ["\\begin{figure}[p]", "\\centering"]
    for row in rows:
        w = 0.98 / ncol
        cells = []
        for lab, (path, name) in row:
            g = (rf"\includegraphics[width=\linewidth,height={row_h:.3f}\textheight,keepaspectratio]{{figures/{name}}}"
                 if path is not None else rf"\missingfig{{{name}}}")
            lab_tex = rf"\textbf{{{lab})}}\\" if lab else r"\mbox{}\\"
            cells.append(rf"\begin{{minipage}}[t]{{{w:.3f}\linewidth}}\raggedright{lab_tex}{g}\end{{minipage}}")
        cells += [rf"\begin{{minipage}}[t]{{{w:.3f}\linewidth}}\mbox{{}}\end{{minipage}}"] * (ncol - len(row))
        out.append("\\hfill".join(cells) + r"\par\medskip")
    out += [rf"\caption*{{{caption_tex}}}", "\\end{figure}", ""]
    return "\n".join(out)


def caption(title_block, legend_block=None):
    t = latex(segments(BLOCKS[title_block]), bold=False).strip()
    m = re.match(r"^((?:Supplementary )?(?:Figures?|Table) S?\d+[A-Z]?(?: and S\d+)?):\s*(.*)$", t, flags=re.S)
    head, rest = (m.group(1), m.group(2)) if m else ("", t)
    cap = rf"\textbf{{{head}:}} {rest}" if head else rest
    if legend_block is not None:
        cap = cap.rstrip()
        cap += ("" if cap[-1] in ".:;)" else ".") + " " + latex(segments(BLOCKS[legend_block]), bold=False).strip()
    return cite(cap)


# ------------------------------------------------------------------ tables from TSV summaries
SUM = REPO / "summaries"
SUM4 = REPO / "v4work/summaries"  # tables rebuilt on the Thille v4 assembly (2026-10-05)
LOWLAND = {"WI": "Unnamed Lowland", "Lowland (WI) HiFi": "Unnamed Lowland HiFi"}


def cell(text):
    """Escape a TSV cell; GENE^x becomes a superscript allele letter; long ;-lists may break after each ;."""
    parts = re.split(r"\^([A-Za-z?]+)", text)
    out = "".join(esc(p) if i % 2 == 0 else rf"\textsuperscript{{{esc(p)}}}" for i, p in enumerate(parts))
    return out.replace(";", r";\allowbreak{}")


def tsv_rows(path):
    lines = [l.rstrip("\n") for l in open(path) if l.strip()]
    notes = [l.lstrip("# ").strip() for l in lines if l.startswith("#")]
    rows = [l.split("\t") for l in lines if not l.startswith("#")]
    return rows[0], rows[1:], notes


def longtable(header, rows, widths, landscape=False, size=r"\scriptsize"):
    avail = 0.98 * (23 if landscape else 16.5)
    total = sum(widths)
    cols = "".join(rf">{{\raggedright\arraybackslash}}p{{{avail * w / total - 0.2:.2f}cm}}" for w in widths)
    head = " & ".join(rf"\textbf{{{h}}}" for h in header) + r" \\"
    body = "\n".join(" & ".join(r) + r" \\" for r in rows)
    return (f"{{{size}\n\\begin{{longtable}}{{{cols}}}\n\\toprule\n{head}\n\\midrule\n\\endfirsthead\n\\toprule\n{head}\n"
            f"\\midrule\n\\endhead\n{body}\n\\bottomrule\n\\end{{longtable}}}}\n")


def table2():
    _h, rows, _n = tsv_rows(SUM / "tables/annotation.tsv")
    out, last = [], None
    for section, metric, value in rows:
        if section != last:
            out.append([rf"\multicolumn{{2}}{{l}}{{\textit{{{esc(section)}}}}}"]); last = section
        indent = r"\qquad " if metric.startswith("  ") else r"\quad "
        out.append([indent + cell(metric.strip()), cell(value)])
    return longtable(["Metric", "Value"], [r if len(r) == 2 else [r[0]] for r in out], [3, 2], size=r"\small")


def table7():
    h, rows, _n = tsv_rows(SUM4 / "tables/tandem_arrays.tsv")
    keep = [i for i, c in enumerate(h) if not c.startswith("_")]
    return longtable([esc(h[i]) for i in keep], [[cell(r[i]) for i in keep] for r in rows], [3, 2.2, 2.6, 2.2, 4, 0.9, 1.1], landscape=True)


def table_s3():
    h, rows, _n = tsv_rows(SUM4 / "genome_wide_coverage_thille.tsv")
    ix = {c: i for i, c in enumerate(h)}
    classes = ["low", "half", "medium", "high", "outlier"]
    hdr = ["Read set", "Chromosome", "Median depth", "Length (bp)"] + [f"{c.capitalize()} bp (\\%)" for c in classes]
    body = [[cell(LOWLAND.get(r[ix["read_set"]], r[ix["read_set"]])), cell(r[ix["chromosome"]]), cell(r[ix["median_depth"]]), f"{int(r[ix['length_bp']]):,}"]
            + [f"{int(r[ix[c + '_bp']]):,} ({r[ix[c + '_pct']]})" for c in classes] for r in rows]
    return longtable(hdr, body, [1.4, 1.2, 0.9, 1.5] + [1.7] * 5, landscape=True)


def table_s4():
    h, rows, _n = tsv_rows(SUM4 / "fst_alleles/a_shared_alleles.tsv")
    hdr = ["Gene", "A-shared allele", "Amino-acid changes", "In Unnamed Lowland", "Carriers", "Contig Fst rank (mean Fst)",
           "Scaffold Fst rank (mean Fst)", "Top 1\\%", "Note"]
    return longtable(hdr, [[cell(c) for c in r] + [""] * (len(h) - len(r)) for r in rows], [1.3, 1.5, 4, 1.2, 2.2, 1.6, 1.6, 1.1, 1.2], landscape=True)


def table_s5():
    h, rows, _n = tsv_rows(SUM4 / "fst_alleles/af_skew_check.tsv")
    rows = [[LOWLAND.get(c, c) for c in r] for r in rows]
    hdr = [esc(c.replace("_", " ")) for c in h]
    w = [1.1, 1.1, 0.8] + [0.55] * 2 + [0.6] * 9 + [6, 1.4]
    return longtable(hdr, [[cell(c) for c in r] for r in rows], w, landscape=True, size=r"\tiny")


def table_s6():
    h, rows, notes = tsv_rows(SUM4 / "fst_alleles/pigment_fst.tsv")
    hdr = ["Gene", "Thille model", "Gwen copy", "A vs B contig rank (mean Fst)", "A vs B scaffold rank (mean Fst)",
           "Gua vs Mex contig rank (mean Fst)", "Gua vs Mex scaffold rank (mean Fst)", "Gua vs Mex top 1\\%", "Note"]
    tab = longtable(hdr, [[cell(c) for c in r] + [""] * (len(h) - len(r)) for r in rows], [1, 1, 2.2, 1.4, 1.4, 1.4, 1.4, 1.5, 4], landscape=True)
    keep = [n for n in notes if not n.startswith(("Source:", "report "))]
    return tab + "{\\footnotesize " + " ".join(cell(n) for n in keep) + "}\n"


# ------------------------------------------------------------------ assemble
ti = find_block("Title:")
title = esc(plain(BLOCKS[ti]).split(":", 1)[1].strip())
authors = []
for b in BLOCKS[find_block("Co-Authors") + 1:find_block("Affiliations")]:
    segs = [s for s in segments(b) if s[0].strip()]
    if not segs:
        continue
    name = "".join(s[0] for s in segs if s[3] != "superscript").strip()
    sup = "".join(s[0] for s in segs if s[3] == "superscript").strip()
    if "," in name:
        last, first = [x.strip() for x in name.split(",", 1)]
        name = f"{first} {last}"
    star = name.endswith("*")
    name = name.rstrip("*").strip()
    if sup.endswith("*") or star:
        sup = sup.rstrip("*")
        star = True
    sup = ",".join(sup.replace(",", "")) if sup else ""
    authors.append(esc(name) + (rf"\textsuperscript{{{sup}{',*' if star else ''}}}" if sup or star else ""))
affils = []
for b in BLOCKS[find_block("Affiliations") + 1:find_block("Potential Keywords")]:
    t = plain(b).strip()
    if t:
        m = re.match(r"^([A-Z*])\.\s*(.*)$", t)
        affils.append(rf"\textsuperscript{{{m.group(1)}}}{esc(m.group(2))}" if m else esc(t))

preamble = r"""\documentclass[11pt]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{textcomp}
\usepackage[margin=2.5cm]{geometry}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{booktabs,longtable,array}
\usepackage{pdflscape}
\usepackage[authoryear,round,semicolon]{natbib}
\usepackage{xurl}
\usepackage[hidelinks]{hyperref}
\usepackage{setspace}
\onehalfspacing
\captionsetup{font=small,labelformat=empty,singlelinecheck=false}
\setlength{\LTcapwidth}{\linewidth}
\newcommand{\missingfig}[1]{\fbox{\parbox[c][5cm][c]{0.9\linewidth}{\centering Figure file not yet provided:\\\texttt{\detokenize{#1}}}}}
"""
byline = ", ".join(authors) + "\\par\n\\medskip\n" + "{\\small " + "\\\\\n".join(affils) + "\\par}\n\\end{center}\n\n"
parts = [preamble, "\\begin{document}\n", "\\begin{center}\n{\\LARGE\\bfseries " + title + "\\par}\n\\bigskip\n" + byline]
supp_title = ("\\begin{center}\n{\\LARGE\\bfseries Supplementary Material\\par}\n\\medskip\n{\\large " + title
              + "\\par}\n\\bigskip\n" + byline)

start = find_block("Abstract")
end = find_block("Figures and Tables")
for b in BLOCKS[start:end]:
    parts.append(table(b) if b.startswith("<w:tbl") else paragraph(b))
parts.append("\n\\bibliographystyle{plainnat}\n\\bibliography{citations}\n\n\\clearpage\n\\section*{Figures}\n\n")

# main figures
# every main figure has a legend paragraph right after its title (v38: legends built from manuscript/figure_legends/)
for key, tb, lb in (("Figure 1", "Figure 1:", 1), ("Figure 2", "Figure 2:", 1), ("Figure 3", "Figure 3:", 1),
                    ("Figure 4", "Figure 4:", 1), ("Figure 5", "Figure 5:", 1),
                    ("Figure 6", "Figure 6:", 1), ("Figure 7", "Figure 7:", 1)):
    i = find_block(tb)
    parts.append(figure_env(key, caption(i, i + lb if lb else None)))

# tables
parts.append("\\clearpage\n\\section*{Tables}\n\n")
i1 = find_block("Table 1:")
parts.append(rf"\noindent {caption(i1)}" + "\n\n" + table(BLOCKS[i1 + 1]) + "\n" + cite(latex(segments(BLOCKS[i1 + 2]), bold=False)) + "\n\n" + table(BLOCKS[i1 + 3]))
i2 = find_block("Table 2:")
parts.append("\n" + rf"\noindent {caption(i2)}" + "\n\n" + table(next(b for b in BLOCKS[i2:i2 + 4] if b.startswith("<w:tbl"))))
i3 = find_block("Table 3:")
parts.append(rf"\noindent {caption(i3)}" + "\n\n" + table(BLOCKS[i3 + 1]))
for lab in ("Table 4:", "Table 5:", "Table 6:"):
    i = find_block(lab)
    tb = next(j for j in range(i, i + 5) if BLOCKS[j].startswith("<w:tbl"))
    parts.append("\\begin{landscape}\n" + rf"\noindent {caption(i, 1)}" + "\n\n" + table(BLOCKS[tb], landscape=True) + "\\end{landscape}\n")
i7 = find_block("Table 7:")
parts.append("\\begin{landscape}\n" + rf"\noindent {caption(i7, 1)}" + "\n\n" + table7() + "\\end{landscape}\n")

parts.append("\\end{document}\n")
main_parts, parts = parts, [preamble, "\\begin{document}\n", supp_title]

# supplementary (supplemental.tex; numbering as cited in the main text)
iS2 = find_block("Supplementary Table S2:")
tbS2 = next(j for j in range(iS2, iS2 + 5) if BLOCKS[j].startswith("<w:tbl"))
parts.append("\\section*{Supplementary Tables}\n\n\\begin{landscape}\n" + rf"\noindent {caption(iS2, 1)}" + "\n\n"
             + table(BLOCKS[tbS2], landscape=True) + "\\end{landscape}\n")
# Supplementary Text S1 (v37: carotenoid and chlorophyll section), between its heading and "Supplementary Figures"
iT = find_block("Supplementary Text S1")
iF = find_block("Supplementary Figures", iT)
parts.append("\\clearpage\n\\section*{" + esc(plain(BLOCKS[iT]).strip()) + "}\n\n")
for b in BLOCKS[iT + 1:iF]:
    parts.append(table(b) if b.startswith("<w:tbl") else paragraph(b))
parts.append("\\clearpage\n\\section*{Supplementary Figures}\n\n")
for key, tb in (("Supplementary Figure S2", "Supplementary Figure S2:"), ("Supplementary Figures S5 and S6", "Supplementary Figures S5 and S6:"),
                ("Supplementary Figure S7", "Supplementary Figure S7:"), ("Supplementary Figure S8", "Supplementary Figure S8:"),
                ("Supplementary Figure S9", "Supplementary Figure S9:"), ("Supplementary Figure S10", "Supplementary Figure S10:"),
                ("Supplementary Figure S11", "Supplementary Figure S11:"), ("Supplementary Figure S12", "Supplementary Figure S12:"),
                ("Supplementary Figure S13", "Supplementary Figure S13:"), ("Supplementary Figure S14", "Supplementary Figure S14:"),
                ("Supplementary Figure S15", "Supplementary Figure S15:")):
    i = find_block(tb)
    parts.append(figure_env(key, caption(i), max_rows=4))
for lab, fn in (("Supplementary Table S3:", table_s3), ("Supplementary Table S4:", table_s4),
                ("Supplementary Table S5:", table_s5), ("Supplementary Table S6:", table_s6)):
    cap = re.sub(r"\s*\((?:summaries/[^)]*)\)", "", caption(find_block(lab)))   # drop the internal file path
    parts.append("\\begin{landscape}\n" + rf"\noindent {cap}" + "\n\n" + fn() + "\\end{landscape}\n")
parts.append("\\end{document}\n")

supp = "".join(parts)
if re.search(r"\\cite", supp):
    supp = supp.replace("\\end{document}\n", "\n\\bibliographystyle{plainnat}\n\\bibliography{citations}\n\\end{document}\n")
(OUT / "supplemental.tex").write_text(supp, encoding="utf-8")
tex = "".join(main_parts)
TEXT_FIXES = []  # corrections made for v20 now live in docx v21 (build_v21.py)
for old, new in TEXT_FIXES:
    if tex.count(old) != 1:
        raise SystemExit(f"text fix matched {tex.count(old)} times: {old[:70]}")
    tex = tex.replace(old, new)
(OUT / "main.tex").write_text(tex, encoding="utf-8")
cited = set(re.findall(r"\\cite(?:p|t|alt)\{([^}]*)\}", tex + supp))
cited = {k for g in cited for k in g.split(",")}
(OUT / "citations.bib").write_text("".join(bib_entry(r) for r in refs), encoding="utf-8")

print(f"main.tex: {len(tex):,} chars; references {len(refs)}; cited keys {len(cited)}")
print("uncited references:", sorted(set(keys.values()) - cited))
print("unresolved citation text:", sorted(unresolved))
print("unmapped characters:", sorted(unknown_chars))
print(f"figures copied: {len(copied)}; missing:", *[f"\n  {n}: {d}" for n, d in MISSING.items()])
