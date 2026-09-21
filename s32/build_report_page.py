#!/usr/bin/env python
"""s32/build_report_page.py -- rebuild the published S31 report page.

Injects `s32/REPORT_S32.md` verbatim into the template `s32/report_page.html` and writes the
standalone page.  The template renders the markdown client-side, so the published artifact and the
repository's report can never drift: there is exactly one copy of the prose.

Published 2026-09-21 (S32).

    python s32/build_report_page.py [out.html]
"""
from __future__ import annotations

import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "REPORT_S32_page.html")
    md = io.open(os.path.join(HERE, "REPORT_S32.md"), encoding="utf-8").read()
    tpl = io.open(os.path.join(HERE, "report_page.html"), encoding="utf-8").read()
    if "</script" in md.lower():
        raise SystemExit("report contains a closing script tag; it cannot be embedded verbatim")
    if "<!--REPORT_MD-->" not in tpl:
        raise SystemExit("template has no <!--REPORT_MD--> placeholder")
    page = tpl.replace("<!--REPORT_MD-->", md, 1)
    io.open(out, "w", encoding="utf-8").write(page)
    n = len(page.encode("utf-8"))
    print("wrote %s  (%d bytes, %.2f MB of the 16 MB budget)" % (out, n, n / 1e6))
    print("markdown lines embedded: %d" % len(md.splitlines()))


if __name__ == "__main__":
    main()
