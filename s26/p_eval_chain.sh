#!/usr/bin/env bash
# s26/p_eval_chain.sh -- lane P's post-sign-off chain: ONE governed lane-P job at a time beside the
# training chain.  Rung evals run in the registered order as their five fold checkpoints appear;
# B3, C4 and C5 fill the gaps.  Every step is resumable (evals checkpoint every 10 targets).
cd "$(dirname "$0")/.." || exit 1
PY=C:/Users/abena/miniforge_3/python.exe
run() { local name=$1 est=$2; shift 2; "$PY" s26/jobrun.py --agent P --tag CPU --name "$name" --est-ram "$est" -- "$PY" "$@"; }
wait_ckpt() { local r=$1; while :; do n=$(ls s26/models/p_ladder/${r}_fold[0-4]_s0.pt 2>/dev/null | wc -l); [ "$n" -ge 5 ] && return 0; sleep 60; done; }
run p_eval_shipped 0.8 s26/p_ladder.py eval --rung shipped
run p_eval_noesm   0.8 s26/p_ladder.py eval --rung noesm
run p_b3_run       0.6 s26/p_b3.py run
wait_ckpt conly;  run p_eval_conly 0.8 s26/p_ladder.py eval --rung conly
wait_ckpt pca32;  run p_eval_pca32 0.8 s26/p_ladder.py eval --rung pca32
run p_eval_mix     0.8 s26/p_ladder.py eval --rung mix
run p_eval_pairnet 1.0 s26/p_ladder.py eval --rung pairnet
run p_c4_run       0.8 s26/p_c4.py run
run p_c5_run       0.8 s26/p_c5.py run
for r in wide pca32f pca128 esm8m raw; do wait_ckpt "$r"; run "p_eval_$r" 1.3 s26/p_ladder.py eval --rung "$r"; done
echo "p_eval_chain: done $(date)"
