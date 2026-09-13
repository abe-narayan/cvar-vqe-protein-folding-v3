#!/usr/bin/env python
"""s26/w_train_chain.py -- Part B of s26/PREREG_selfcopy_bound.md: the channel-B retrains, in the
order fixed in the PREREG (section 8): four carrier-out models, six control-out models, then the
reference models for folds 2, 4, 0 UNLESS lane P's `s26/models/p_ladder/pca32_fold<f>_s0.pt`
(the identical function, corpus and seed with nothing removed) already exists.  A separate
driver so that `s26/w_selfcopy.py` is not edited while its jobs run.  One .pt per model is the
checkpoint; existing models are skipped, so the job is resumable.

    python s26/jobrun.py --agent W --tag CPU --name w_train_chain --est-ram 1.5 -- python s26/w_train_chain.py
"""
from __future__ import annotations

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)

import w_selfcopy as W                               # noqa: E402


def main():
    t0 = time.time()
    for pdb, fold in W.TRAIN_ORDER:
        print("== carrier/control out: %s from fold %d (%.0f s elapsed)" % (pdb, fold, time.time() - t0), flush=True)
        W.train_without(fold, [pdb], "out_%s" % pdb)
    for fold in (2, 4, 0):
        p_ref = os.path.join(W.P_MODELS, "pca32_fold%d_s%d.pt" % (fold, 0))
        if os.path.exists(p_ref):
            print("== reference fold %d: lane P's %s exists; not retrained here" % (fold, p_ref), flush=True)
            continue
        print("== reference fold %d: training this lane's own (%.0f s elapsed)" % (fold, time.time() - t0), flush=True)
        W.train_without(fold, [], "ref")
    print("chain complete in %.0f s" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
