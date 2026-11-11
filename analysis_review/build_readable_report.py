"""Render a Markdown report as a self-contained HTML document and printable PDF.

Default: python analysis_review/build_readable_report.py
One section: python analysis_review/build_readable_report.py --source reports/04_classifier_sequences.md
"""
from pathlib import Path
import argparse
import base64
import html
import os
import re
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "venomics-mpl"))
import markdown
from bs4 import BeautifulSoup
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.utils import ImageReader
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, KeepTogether, Preformatted
from matplotlib.font_manager import findfont, FontProperties

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports/exports"
OUT.mkdir(parents=True, exist_ok=True)


# A single catalogue drives batch builds and links between the rendered reports.
REPORTS = {
    "reports/FINAL_ANALYSIS.md": ("FINAL_ANALYSIS", "Detailed results", 9.7),
    "reports/01_shortlists.md": ("01_shortlists", "Section 1: shortlists", 10.5),
    "reports/02_concentration_rules.md": ("02_concentration_rules", "Section 2: concentration rules", 10.5),
    "reports/03_matched_targets.md": ("03_matched_targets", "Section 3: matched targets", 10.5),
    "reports/04_classifier_sequences.md": ("04_classifier_sequences", "Section 4: classifier sequences", 10.5),
    "reports/05_sequence_responses.md": ("05_sequence_responses", "Section 5: sequence responses", 11),
    "reports/06_venomics.md": ("06_venomics", "Section 6: Venomics", 10.5),
    "reports/07_interpretation.md": ("07_interpretation", "Section 7: interpretation", 10.5),
    "reports/SOURCE_COMPARISONS.md": ("SOURCE_COMPARISONS", "Source comparisons", 10.5),
    "reports/PRELIMINARY_ANALYSES.md": ("PRELIMINARY_ANALYSES", "Other preliminary analyses", 10.5),
    "reports/REPRODUCIBILITY.md": ("REPRODUCIBILITY", "Reproducibility", 10.5),
    "reports/METHODS.md": ("METHODS", "Methods and analysis map", 9.7),
}
REPORT_PAGES = {(ROOT / source).resolve(): stem + ".html" for source, (stem, _, _) in REPORTS.items()}


def render_report(source_path, output_stem, review_date):
    source_path = Path(source_path).resolve()
    if Path(output_stem).name != output_stem:
        raise ValueError("output_stem must be a filename stem without directories")
    source_relative = os.path.relpath(source_path, ROOT)
    spec = REPORTS.get(source_relative, (output_stem, "Detailed results", 9.7))
    _, report_label, body_font = spec
    rendered = markdown.markdown(source_path.read_text(), extensions=["tables", "fenced_code", "toc", "md_in_html"])
    soup = BeautifulSoup(rendered, "html.parser")
    title = soup.find("h1").get_text() if soup.find("h1") else "DAMP results"
    for img in soup.find_all("img"):
        path = source_path.parent / img["src"]
        img["src"] = "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()
    for link in soup.find_all("a", href=True):
        if not re.match(r"^[a-z]+:|^#", link["href"]):
            target, separator, fragment = link["href"].partition("#")
            target_path = (source_path.parent / target).resolve()
            destination = REPORT_PAGES.get(target_path, os.path.relpath(target_path, OUT))
            link["href"] = destination + (separator + fragment if separator else "")
    style = """
    :root { color-scheme:light; --ink:#20313c; --accent:#225b78; }
    * {box-sizing:border-box} body {margin:0;background:#f4f1eb;color:var(--ink);
    font:18px/1.65 -apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;}
    main {max-width:1050px;margin:40px auto;padding:48px 64px;background:white;
    border-top:7px solid var(--accent);box-shadow:0 8px 35px #20313c0c}
    h1 {font-size:2.1rem;line-height:1.2;letter-spacing:-.025em;margin-top:0}
    h2 {font-size:1.35rem;line-height:1.35;color:var(--accent);margin-top:2.2em}
    a {color:#165d83;text-underline-offset:3px} p {margin:1em 0}
    table {border-collapse:collapse;width:100%;font-size:.92rem;margin:1.3em 0}
    th,td {text-align:left;padding:12px;border-bottom:1px solid #d9e1e6;vertical-align:top}
    th {background:#eaf1f4} tbody tr:nth-child(even) {background:#f8fafb}
    img {display:block;width:100%;height:auto;margin:1.5em 0 .6em}
    details {margin:1.3em 0;padding:.7em 1em;border:1px solid #d9e1e6;border-radius:6px}
    summary {cursor:pointer;color:var(--accent);font-weight:600}
    blockquote {border-left:4px solid #287a65;margin:1.4em 0;padding:.2em 1.3em;
    background:#f2f7f5;font-size:.98em} code {font-size:.9em;overflow-wrap:anywhere}
    @media(max-width:700px) {main {margin:0;padding:28px 22px} body {font-size:16px}
    table {font-size:.78rem} th,td {padding:7px} h1 {font-size:1.8rem}}
    @media print {body {background:white;font-size:11pt} main {margin:0;max-width:none;
    padding:0;box-shadow:none;border:0} h2 {break-after:avoid} table,figure {break-inside:avoid}}
    """
    document = '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>' + html.escape(title) + '</title><style>' + style + '</style></head><body><main>' + str(soup) + '</main></body></html>'
    (OUT / (output_stem + ".html")).write_text(document)

    for name, weight, slant in [("ReviewSans", "normal", "normal"), ("ReviewSansBold", "bold", "normal"),
                                 ("ReviewSansItalic", "normal", "oblique"), ("ReviewSansBoldItalic", "bold", "oblique")]:
        pdfmetrics.registerFont(TTFont(name, findfont(FontProperties(family="DejaVu Sans", weight=weight, style=slant))))
    pdfmetrics.registerFontFamily("ReviewSans", normal="ReviewSans", bold="ReviewSansBold", italic="ReviewSansItalic", boldItalic="ReviewSansBoldItalic")
    pdfmetrics.registerFont(TTFont("ReviewMono",findfont(FontProperties(family="DejaVu Sans Mono"))))
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="ReviewBody", fontName="ReviewSans", fontSize=9.7, leading=15.2, textColor=colors.HexColor("#20313c"), spaceAfter=9))
    styles.add(ParagraphStyle(name="ReviewTitle", parent=styles["ReviewBody"], fontName="ReviewSansBold", fontSize=23, leading=28, spaceAfter=16))
    styles.add(ParagraphStyle(name="ReviewHeading", parent=styles["ReviewBody"], fontName="ReviewSansBold", fontSize=13.5, leading=18, textColor=colors.HexColor("#225b78"), spaceBefore=15, spaceAfter=9, keepWithNext=True))
    styles.add(ParagraphStyle(name="ReviewTable", parent=styles["ReviewBody"], fontSize=8, leading=11.5, spaceAfter=0))
    styles.add(ParagraphStyle(name="ReviewQuote", parent=styles["ReviewBody"], leftIndent=12, rightIndent=8, backColor=colors.HexColor("#f2f7f5"), borderPadding=8, spaceBefore=4, spaceAfter=12))
    styles.add(ParagraphStyle(name="ReviewCaption", parent=styles["ReviewBody"], fontSize=8.3, leading=12, textColor=colors.HexColor("#526975"), spaceAfter=10))
    styles.add(ParagraphStyle(name="ReviewCode", parent=styles["ReviewBody"], fontName="ReviewMono", fontSize=10.5, leading=15, leftIndent=10, spaceAfter=9))
    styles["ReviewBody"].fontSize = body_font
    if body_font > 9.7:
        styles["ReviewBody"].leading = 16.5
    source = BeautifulSoup(rendered, "html.parser")
    # A printed report cannot toggle figures: expand their contents for PDF.
    for details in source.find_all("details"):
        summary = details.find("summary", recursive=False)
        if summary:
            summary.decompose()
        details.unwrap()
    page_width, page_height = A4
    width = page_width - 88
    story = []
    def inline(element):
        clone = BeautifulSoup(str(element), "html.parser")
        for link in clone.find_all("a"):
            if not re.match(r"^https?://", link.get("href", "")):
                link.unwrap()
        for tag in clone.find_all("strong"):
            tag.name = "b"
        for tag in clone.find_all("em"):
            tag.name = "i"
        for tag in clone.find_all("code"):
            tag.unwrap()
        el = clone.find()
        return el.decode_contents() if el else html.escape(str(element))
    def paragraph(element, style_name="ReviewBody"):
        return Paragraph(inline(element), styles[style_name])
    elements = [element for element in source.children if getattr(element, "name", None)]
    consumed = set()
    for index, element in enumerate(elements):
        if index in consumed:
            continue
        if element.name == "h1":
            story.append(paragraph(element, "ReviewTitle"))
        elif element.name in ["h2", "h3"]:
            story.append(paragraph(element, "ReviewHeading"))
        elif element.name == "p":
            img = element.find("img")
            if img:
                path = source_path.parent / img["src"]
                iw, ih = ImageReader(str(path)).getSize()
                content = [Spacer(1, 7), Image(str(path), width=width, height=width*ih/iw), Spacer(1, 7)]
                if index + 1 < len(elements) and elements[index+1].name == "p" and elements[index+1].find("em"):
                    content.append(paragraph(elements[index+1], "ReviewCaption"))
                    consumed.add(index+1)
                story.append(KeepTogether(content))
            elif element.get_text(strip=True):
                story.append(paragraph(element))
        elif element.name == "blockquote":
            for p in element.find_all("p", recursive=False):
                story.append(paragraph(p, "ReviewQuote"))
        elif element.name == "pre":
            story.append(KeepTogether([Spacer(1,5),Preformatted(element.get_text().rstrip(),styles["ReviewCode"]),Spacer(1,5)]))
        elif element.name == "table":
            table_rows = []
            for row in element.find_all("tr"):
                table_rows.append([paragraph(cell, "ReviewTable") for cell in row.find_all(["td", "th"], recursive=False)])
            ncols = len(table_rows[0])
            weights = {2:[.48,.52],3:[.40,.29,.31],4:[.38,.14,.24,.24]}.get(ncols,[1/ncols]*ncols)
            table = Table(table_rows, colWidths=[width*w for w in weights], repeatRows=1, hAlign="LEFT")
            table.setStyle(TableStyle([
                ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#eaf1f4")),
                ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#f8fafb")]),
                ("LINEBELOW",(0,0),(-1,-1),.4,colors.HexColor("#d9e1e6")),
                ("VALIGN",(0,0),(-1,-1),"TOP"),("TOPPADDING",(0,0),(-1,-1),7),
                ("BOTTOMPADDING",(0,0),(-1,-1),7),("LEFTPADDING",(0,0),(-1,-1),7),
                ("RIGHTPADDING",(0,0),(-1,-1),7),
            ]))
            before_table = Spacer(1,5)
            before_table.keepWithNext = True
            story.extend([before_table,table,Spacer(1,11)])
        elif element.name in ["ul", "ol"]:
            for i, item in enumerate(element.find_all("li", recursive=False)):
                story.append(Paragraph((f"{i+1}. " if element.name=="ol" else "• ")+inline(item),styles["ReviewBody"]))
    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("ReviewSans",8)
        canvas.setFillColor(colors.HexColor("#647783"))
        canvas.drawString(44,27,"DAMP · " + report_label + " · " + review_date)
        canvas.drawRightString(page_width-44,27,str(doc.page))
        canvas.restoreState()
    doc = SimpleDocTemplate(str(OUT / (output_stem + ".pdf")),pagesize=A4,rightMargin=44,leftMargin=44,
                            topMargin=42,bottomMargin=44,title=title,author="DAMP analysis review")
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    return {"source": source_relative, "output_stem": output_stem}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--source", default="reports/FINAL_ANALYSIS.md")
    parser.add_argument("--output-stem")
    parser.add_argument("--review-date", default="10 November 2026")
    selection.add_argument("--all", action="store_true", help="Build all current reports from the catalogue")
    args = parser.parse_args()
    if args.all and args.output_stem:
        parser.error("--output-stem applies to one report, not --all")
    jobs = list(REPORTS.items()) if args.all else [(args.source, REPORTS.get(args.source, (Path(args.source).stem.upper(), "", 9.7)))]
    for source, (stem, _, _) in jobs:
        result = render_report(ROOT / source, args.output_stem or stem, args.review_date)
        print("Wrote " + result["output_stem"] + ": HTML and PDF.")


if __name__ == "__main__":
    main()
