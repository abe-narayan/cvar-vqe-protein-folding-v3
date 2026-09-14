#!/usr/bin/env python
"""s26/w_provenance_test.py -- synthetic test of the class rule in s26/w_provenance.py."""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                   # noqa: E402

import w_provenance as P                             # noqa: E402
from core import data as D                           # noqa: E402


class Pep:
    def __init__(self, seq): self.seq = seq


def test_class_rule_on_a_fake_corpus():
    peptides = [Pep("ACDEFGHIK"), Pep("ACDEFGHIKLMNP"), Pep("WWWWACDEFGHIKWW")]
    fragments = [Pep("QQQQQQQQQQQQ")]
    orig_load, orig_frags = D.load, D.load_fragments
    D.load = lambda: peptides; D.load_fragments = lambda: fragments
    P._IDX.clear()
    try:
        idx, multi = P.substring_index(9)
        assert idx["ACDEFGHIK"] == ("whole", "peptide")          # first parent: the 9-mer itself
        assert multi >= 2                                          # the same 9-mer occurs in the two longer peptides
        assert idx["CDEFGHIKL"] == ("interior", "peptide")
        assert idx["DEFGHIKLM"] == ("interior", "peptide")
        assert idx["FGHIKLMNP"] == ("terminal", "peptide")
        assert idx["QQQQQQQQQ"] == ("fragment", "fragment")
        S = np.array([D.encode("ACDEFGHIK"), D.encode("FGHIKLMNP"), D.encode("QQQQQQQQQ"), D.encode("YYYYYYYYY")], dtype=np.int8)
        u = {"S": S, "org": np.array([True, True, False, False])}
        assert P.classify(u, "ACDEFGHIK", [0, 1, 2, 3]) == ["whole", "terminal", "fragment", "unresolved"]
    finally:
        D.load, D.load_fragments = orig_load, orig_frags
        P._IDX.clear()
    print("  class rule OK")


if __name__ == "__main__":
    test_class_rule_on_a_fake_corpus()
    print("ALL OK")
