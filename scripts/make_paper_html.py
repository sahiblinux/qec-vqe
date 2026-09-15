#!/usr/bin/env python3
"""Convert a paper markdown file to a print-ready HTML file (open in browser -> Print -> Save as PDF).

Usage:
  python scripts/make_paper_html.py                    # converts paper/research_paper.md
  python scripts/make_paper_html.py paper/paper.md     # converts any given markdown
"""
import sys
from pathlib import Path
import markdown

ROOT = Path(__file__).resolve().parent.parent
src = ROOT / sys.argv[1] if len(sys.argv) > 1 else ROOT / "paper" / "research_paper.md"
dst = src.with_suffix(".html")

text = src.read_text()
# resolve relative figure paths to absolute so images render from any location
text = text.replace("](../out/", "](file://" + str(ROOT / "out") + "/")

body = markdown.markdown(text, extensions=["tables", "fenced_code"])

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>When Does Zero-Noise Extrapolation Help? — Research Paper</title>
<style>
  @page {{ margin: 2cm; }}
  body {{ font-family: Georgia, 'Times New Roman', serif; font-size: 11pt;
         line-height: 1.55; max-width: 800px; margin: 0 auto; padding: 32px;
         color: #111; }}
  h1 {{ font-size: 19pt; line-height: 1.25; margin-bottom: 6px; }}
  h2 {{ font-size: 14pt; border-bottom: 1px solid #bbb; padding-bottom: 3px;
       margin-top: 26px; }}
  h3 {{ font-size: 12pt; margin-top: 18px; }}
  table {{ border-collapse: collapse; margin: 14px auto; font-size: 9.5pt;
          font-family: 'Helvetica Neue', Arial, sans-serif; }}
  th, td {{ border: 1px solid #999; padding: 4px 9px; text-align: center; }}
  th {{ background: #eef2f7; }}
  img {{ max-width: 100%; display: block; margin: 12px auto; }}
  code, pre {{ font-family: Menlo, Consolas, monospace; font-size: 9pt;
              background: #f5f5f5; }}
  pre {{ padding: 10px; border-radius: 4px; overflow-x: auto; }}
  blockquote {{ color: #444; border-left: 3px solid #bbb; margin-left: 0;
               padding-left: 12px; }}
  hr {{ border: none; border-top: 1px solid #ccc; margin: 22px 0; }}
</style>
</head>
<body>
{body}
</body>
</html>
"""
dst.write_text(html)
print(f"wrote {dst} ({len(html)} chars)")
