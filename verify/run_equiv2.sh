#!/bin/sh
# Round 2 of the two-arm equivalence, with `core.project` live.  The projection is the stage
# where the consolidated module can ship an ANALYTIC gradient in place of the reference's
# finite differences, so more than two arms are needed to separate "the consolidation is
# faithful" from "the gradient changes the science":
#
#   baseline    every backend legacy                    (s8.project, fd gradient)
#   opt_exact   consolidated, project_grad=exact        THE SHIPPED DEFAULT: the bit-exact
#                                                       builder and the reference's own finite
#                                                       difference; expected bit-identical to
#                                                       baseline
#   opt_fd      consolidated, project_grad=fd           the scan builder with the reference's
#                                                       gradient: isolates the builder
#   opt_an      consolidated, project_grad=analytic     the scan builder and the true gradient
#
# HISTORY.  As first written this script set PROJECT_GRAD in the environment, and the shipped
# default was `analytic`.  An import-time environment global cannot reach the cache key
# (verify/grad_key_collision.py: two modes, one cfg_key, the second run served the first's
# arrays), so the mode was moved into core.pipeline.Config.project_grad, hashed like every
# other science parameter, and the default became `exact`.  From then until S26 this script
# ran three arms in the SAME mode and was flagged STALE (S26 operational item 2, lane I).  The
# core.pipeline CLI has no flag for the mode, so the consolidated arms now go through
# s26/i_run_equiv2_arm.py, which builds the identical Config to `python -m core.pipeline run`
# and sets project_grad on it; it prints the same JSON, so the cfg_key grep below still works.
#
# Usage:   verify/run_equiv2.sh [--no-amber]
#          any argument is passed to every arm; --no-amber skips stage 4 in all four, which
#          leaves the projection comparison intact and keeps OpenMM out of the run
# Then:    python s26/i_equiv2_compare.py [--no-amber]     per-target comparison of the arms
cd "$(dirname "$0")/.."
EXTRA="$*"
echo "=== BASELINE (legacy) $EXTRA ==="
python -m core.pipeline run --manifest smoke8 --workers 1 --backends legacy $EXTRA \
    > verify/e2_baseline.log 2>&1
echo "baseline rc=$?"
echo "=== OPTIMISED, project_grad=exact (the shipped default) $EXTRA ==="
python s26/i_run_equiv2_arm.py --grad exact --manifest smoke8 --workers 1 $EXTRA \
    > verify/e2_opt_exact.log 2>&1
echo "opt_exact rc=$?"
echo "=== OPTIMISED, project_grad=fd $EXTRA ==="
python s26/i_run_equiv2_arm.py --grad fd --manifest smoke8 --workers 1 $EXTRA \
    > verify/e2_opt_fd.log 2>&1
echo "opt_fd rc=$?"
echo "=== OPTIMISED, project_grad=analytic $EXTRA ==="
python s26/i_run_equiv2_arm.py --grad analytic --manifest smoke8 --workers 1 $EXTRA \
    > verify/e2_opt_an.log 2>&1
echo "opt_an rc=$?"
echo "=== DONE ==="
grep -o '"cfg_key": "[a-f0-9]*"' verify/e2_baseline.log verify/e2_opt_exact.log \
    verify/e2_opt_fd.log verify/e2_opt_an.log
grep -h '"diagnostic"' verify/e2_opt_exact.log verify/e2_opt_fd.log verify/e2_opt_an.log
