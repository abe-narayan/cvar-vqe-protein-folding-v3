#!/usr/bin/env bash
# s26/p_train_chain.sh -- the C2 ladder TRAINING chain, one governor job at a time, in the order
# fixed by s26/PREREG_C2.md addendum 2.  Two passes: pass 2 resumes any rung the governor killed
# (completed fold checkpoints are skipped), under a _p2 name so pass 1's jobs_done record survives.
# No evaluation happens here (phase gate; enforced in s26/p_ladder.py).
cd "$(dirname "$0")/.." || exit 1
PY=C:/Users/abena/miniforge_3/python.exe
run() { local name=$1 est=$2; shift 2; "$PY" s26/jobrun.py --agent P --tag CPU --name "$name" --est-ram "$est" -- "$PY" "$@"; }
for pass in 1 2; do
  sfx=""; [ "$pass" = 2 ] && sfx="_p2"
  run "p_train_noesm$sfx"  1.5 s26/p_ladder.py train --rung noesm  --folds 0,1,2,3,4
  run "p_train_conly$sfx"  1.5 s26/p_ladder.py train --rung conly  --folds 0,1,2,3,4
  run "p_train_pca32$sfx"  1.5 s26/p_ladder.py train --rung pca32  --folds 0,1,2,3,4
  run "p_train_wide$sfx"   1.6 s26/p_ladder.py train --rung wide   --folds 0,1,2,3,4
  run "p_train_pca32f$sfx" 1.5 s26/p_ladder.py train --rung pca32f --folds 0,1,2,3,4
  run "p_train_pca128$sfx" 1.5 s26/p_ladder.py train --rung pca128 --folds 0,1,2,3,4
  run "p_featurise_esm8m$sfx" 1.0 s26/p_ladder.py featurise-esm8m
  run "p_train_esm8m$sfx"  1.5 s26/p_ladder.py train --rung esm8m  --folds 0,1,2,3,4
  run "p_train_raw$sfx"    1.8 s26/p_ladder.py train --rung raw    --folds 0,1,2,3,4
done
echo "p_train_chain: done $(date)"
