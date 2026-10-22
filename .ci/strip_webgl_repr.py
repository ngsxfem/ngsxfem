#!/usr/bin/env python3
"""Remove the text output "WebGLScene" (repr of the object returned by Draw()) from the html pages.

Only an output that consists of exactly this text is removed; other text outputs (print, warnings, ...)
of the same cell stay. Usage: strip_webgl_repr.py <html output dir>
"""
import pathlib
import re
import sys

pattern = re.compile(
    r'<div class="output text_plain[^"]*">\s*<div class="highlight">\s*<pre>\s*(?:<span></span>)?\s*WebGLScene\s*</pre>\s*</div>\s*</div>')

n = 0
for f in pathlib.Path(sys.argv[1]).rglob("*.html"):
    text = f.read_text(encoding="utf-8")
    new, c = pattern.subn("", text)
    if c:
        f.write_text(new, encoding="utf-8")
        n += c
print(f"removed {n} 'WebGLScene' outputs")
