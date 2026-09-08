#!/usr/bin/env bash
# s18/run_phys_chain.sh -- the PHYSICS workstream's SEQUENTIAL queue.
#
# ONE heavy process at a time (the coordinator's constraint at 100% box CPU), BLAS pinned to a
# single thread so numpy cannot silently fan out across cores, and every stage resumable from its
# own per-target checkpoint.
#
# ORDER, revised after the coordinator closed the degree-1 branch:
#   1  phys_decorr   the ORACLE go/no-go on whether `leg_contact` corrects or amplifies the
#                    distogram's per-pair error.  Cheap, and it decides how everything else reads.
#   2  phys_lambda --bases   the contact term on the DEPLOYED objective (E_full) and on MATH's
#                    residue-additive object.  With degree-1 closed this is the LIVE form of the
#                    hypothesis and it is promoted above Phase 8 for that reason.
#   3  phys_down     PHASE 8, the causal downstream comparison.  Independent of degree-1.
#   4  phys_space    PHASE 10, secondary, pass 1 only (no AMBER).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$(dirname "$0")/.."
log() { echo "$(date +%T) $*" >> s18/results/chain.log; }
until [ "$(python -c "import json;print(json.load(open('s18/results/lam.json'))['complete'])" 2>/dev/null)" = "True" ]; do sleep 30; done
log "phase6 lam complete -> decorr"
python -m s18.phys_decorr >> s18/results/decorr.log 2>&1
log "decorr complete -> bases (the LIVE form of the hypothesis)"
python -m s18.phys_lambda --bases --out lam_bases.json >> s18/results/lam_bases.log 2>&1
log "bases complete -> down (PHASE 8)"
python -m s18.phys_down --out down.json >> s18/results/down.log 2>&1
log "down complete -> space pass 1"
python -m s18.phys_space --out space.json >> s18/results/space.log 2>&1
log "ALL PHYSICS STAGES COMPLETE"
