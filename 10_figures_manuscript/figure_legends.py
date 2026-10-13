"""Figure legends from manuscript/figure_legends/Figure_<id>.md into the docx (Edwin, 2026-10-09: no figure text inside
the PDFs/PNGs; each figure's legend lives in its own md file, split by panel, and the manuscript legends are built from
those files). build_latex.py takes captions from the docx, so the LaTeX build follows.

md format (agreed with the thille-hic-bf session):
  # Figure 6                      or  # Supplementary Figure S7
  **<title>.**
  [shared text]                   optional, before the panels
  (A-C) <text shared by panels A to C>   shared details first (Edwin, 2026-10-10)
  (A) <panel A text>              then one line per panel ID; a figure without panels has one body line instead
Plain text, *italic* for gene and species names, and X^y for a one-letter allele superscript (SGRL^d).

Docx result: main figures keep the title paragraph ("Figure 6: <title>", bold) and get one legend paragraph after it
("(A) ... (B) ..."; inserted when the figure had none). Supplementary figures are one paragraph: the bold
"Supplementary Figure S7: <title>" followed by the legend; S5 and S6 share the paragraph "Supplementary Figures S5 and
S6:". Usage: apply_legends(document_xml, md_dir, ids) -> (document_xml, report).
"""
import html
import re
from pathlib import Path

T_RE = re.compile(r"(<w:t(?: [^>]*)?>)([^<]*)(</w:t>)")
P_RE = re.compile(r"<w:p[ >].*?</w:p>", flags=re.S)
txt = lambda p: "".join(html.unescape(m.group(2)) for m in T_RE.finditer(p))
PANEL = re.compile(r"\([A-Z](?:-[A-Z])?(?:, [A-Z](?:-[A-Z])?)*\) ")   # (A), (A-C), (A, C), (B-D, F)
INLINE = re.compile(r"\*([^*]+)\*|([A-Za-z0-9-]+)\^([a-z?])")


def parse(path):
    lines = [l.rstrip() for l in Path(path).read_text().splitlines() if l.strip()]
    head = re.fullmatch(r"# (Figure \w+|Supplementary Figure S\d+)", lines[0])
    title = re.fullmatch(r"\*\*(.+)\*\*", lines[1])
    if not head or not title:
        raise SystemExit(f"{path}: expected '# Figure <id>' and '**<title>.**' on the first two lines")
    body = lines[2:]
    is_panel = [bool(PANEL.match(l)) for l in body]
    if any(is_panel) and not all(is_panel[is_panel.index(True):]):
        raise SystemExit(f"{path}: a non-panel line after the panel lines")   # shared text may precede the panels
    return head.group(1), title.group(1), " ".join(body)


def runs(text, rpr, bold=False):
    """w:r elements for text with *italic* and X^y superscripts, using the paragraph's run properties."""
    base = re.sub(r"<w:(?:b|bCs|i|iCs)(?: [^>]*)?/>|<w:vertAlign [^>]*/>", "", rpr or "<w:rPr></w:rPr>")
    def r(t, b=False, i=False, sup=False):
        extra = ("<w:b/><w:bCs/>" if b else "") + ("<w:i/><w:iCs/>" if i else "") + ('<w:vertAlign w:val="superscript"/>' if sup else "")
        pr = base.replace("<w:rPr>", "<w:rPr>" + extra, 1) if extra else base
        return f'<w:r>{pr}<w:t xml:space="preserve">{html.escape(t, quote=False)}</w:t></w:r>'
    out, pos = [], 0
    for m in INLINE.finditer(text):
        if m.start() > pos:
            out.append(r(text[pos:m.start()], bold))
        if m.group(1):
            out.append(r(m.group(1), bold, i=True))
        else:
            out.append(r(m.group(2), bold)); out.append(r(m.group(3), bold, sup=True))
        pos = m.end()
    if pos < len(text):
        out.append(r(text[pos:], bold))
    return "".join(out)


def rebuild(p, pieces):
    """Paragraph p with its runs replaced by pieces [(text, bold)], keeping pPr and the first run's rPr."""
    ppr = re.search(r"<w:pPr>.*?</w:pPr>", p, flags=re.S)
    rpr = re.search(r"<w:r>(?:\s*)(<w:rPr>.*?</w:rPr>)", p, flags=re.S) or re.search(r"<w:r [^>]*>(?:\s*)(<w:rPr>.*?</w:rPr>)", p, flags=re.S)
    open_tag = re.match(r"<w:p(?: [^>]*)?>", p).group(0)
    # image runs stay (Figure 1's image sits in its title paragraph, after the text)
    drawings = "".join(m.group(0) for m in re.finditer(r"<w:r(?: [^>]*)?>(?:(?!</w:r>).)*?<w:drawing>.*?</w:r>", p, flags=re.S))
    return (open_tag + (ppr.group(0) if ppr else "") + "".join(runs(t, rpr.group(1) if rpr else None, b) for t, b in pieces)
            + drawings + "</w:p>")


def apply_legends(d, md_dir, ids, legend_template_prefix):
    """ids: figure ids to apply (e.g. ["1", "2", "6", "S7"]); legend_template_prefix: start of an existing legend
    paragraph whose formatting is cloned for main figures without one."""
    paras = lambda: list(P_RE.finditer(d))
    tmpl = next(m.group(0) for m in paras() if txt(m.group(0)).startswith(legend_template_prefix))
    report = []
    md = {i: parse(Path(md_dir) / f"Figure_{i}.md") for i in ids}
    for i, (head, title, body) in md.items():
        if i.startswith("S") and i in ("S5", "S6"):
            continue
        if not i.startswith("S"):
            ps = paras()
            k = next((n for n, m in enumerate(ps) if txt(m.group(0)).startswith(f"Figure {i}:")), None)
            if k is None:
                raise SystemExit(f"no title paragraph 'Figure {i}:' in the document")
            t = ps[k]
            new_t = rebuild(t.group(0), [(f"Figure {i}: {title}", True)])
            nxt = ps[k + 1]
            nt = txt(nxt.group(0)).strip()
            is_legend = (bool(nt) and not re.fullmatch(r"[A-Z]\)", nt) and "<w:drawing" not in nxt.group(0)
                         and not re.match(r"(Supplementary )?Figures? ", nt))   # Figures 1 and 2 have consecutive titles
            if is_legend:
                d = d[:t.start()] + new_t + d[t.end():nxt.start()] + rebuild(nxt.group(0), [(body, False)]) + d[nxt.end():]
                report.append(f"Figure {i}: title and legend replaced")
            else:
                d = d[:t.start()] + new_t + rebuild(tmpl, [(body, False)]) + d[t.end():]
                report.append(f"Figure {i}: title replaced, legend inserted")
        else:
            ps = paras()
            m = next(m for m in ps if txt(m.group(0)).startswith(f"Supplementary Figure {i}:"))
            d = d[:m.start()] + rebuild(m.group(0), [(f"Supplementary Figure {i}: {title}", True), (" " + body, False)]) + d[m.end():]
            report.append(f"Supplementary Figure {i}: replaced")
    if "S5" in md and "S6" in md:
        ps = paras()
        m = next(m for m in ps if txt(m.group(0)).startswith("Supplementary Figures S5 and S6:"))
        pieces = [("Supplementary Figures S5 and S6:", True)]
        for i in ("S5", "S6"):
            _, title, body = md[i]
            pieces += [(f" {i}: {title}", True), (" " + body, False)]
        d = d[:m.start()] + rebuild(m.group(0), pieces) + d[m.end():]
        report.append("Supplementary Figures S5 and S6: replaced")
    return d, report
