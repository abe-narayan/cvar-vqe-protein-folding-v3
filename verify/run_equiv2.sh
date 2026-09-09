#!/bin/sh
# STALE as of the cache-key fix: PROJECT_GRAD is no longer read from the environment.
# `core/project.py` pins GRAD = "exact" and the mode now travels as
# `core.pipeline.Config.project_grad` so that it reaches the cache key (see that module's
# docstring, and verify/grad_key_collision.py for the collision it fixed). Running this as
# written gives three arms in the SAME mode. Kept as the record of how the three-arm
# comparison was originally driven; set Config.project_grad to reproduce it.
# Round 2, with `core.project` live. The projection is the stage where the consolidated
# module ships an ANALYTIC gradient in place of the reference's finite differences, so
# three arms are needed to separate "the consolidation is faithful" from "the analytic
# gradient changes the science":
#
#   baseline   every backend legacy                       (s8.project, fd gradient)
#   opt_fd     consolidated, PROJECT_GRAD=fd              (isolates the consolidation)
#   opt_an     consolidated, PROJECT_GRAD=analytic        (the shipped default)
cd "$(dirname "$0")/.."
echo "=== BASELINE (legacy) ==="
python -m core.pipeline run --manifest smoke8 --workers 1 --backends legacy \
    > verify/e2_baseline.log 2>&1
echo "baseline rc=$?"
echo "=== OPTIMISED, PROJECT_GRAD=fd ==="
PROJECT_GRAD=fd python -m core.pipeline run --manifest smoke8 --workers 1 \
    > verify/e2_opt_fd.log 2>&1
echo "opt_fd rc=$?"
echo "=== OPTIMISED, PROJECT_GRAD=analytic (the default) ==="
PROJECT_GRAD=analytic python -m core.pipeline run --manifest smoke8 --workers 1 \
    > verify/e2_opt_an.log 2>&1
echo "opt_an rc=$?"
echo "=== DONE ==="
grep -o '"cfg_key": "[a-f0-9]*"' verify/e2_baseline.log verify/e2_opt_fd.log verify/e2_opt_an.log
