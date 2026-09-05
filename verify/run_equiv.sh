#!/bin/sh
# Two arms of the SAME manifest, CURRENT backend set, one worker each (serialised, so the
# box never carries two heavy jobs). The baseline arm forces the legacy root modules.
cd "$(dirname "$0")/.."
echo "=== BASELINE (CORE_BACKENDS=legacy) ==="
python -m core.pipeline run --manifest smoke8 --workers 1 --backends legacy \
    > verify/equiv_baseline.log 2>&1
echo "baseline rc=$?"
echo "=== OPTIMISED (current consolidated set) ==="
python -m core.pipeline run --manifest smoke8 --workers 1 \
    > verify/equiv_optimised.log 2>&1
echo "optimised rc=$?"
echo "=== DONE ==="
