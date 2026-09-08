#!/bin/bash
# REPAIR workstream: run the passes STRICTLY ONE AT A TIME (one heavy process on the box).
# Pass A is launched separately and must already be running or finished.
cd /c/Users/abena/Protein-Folding-Algorithm
while ps -ef | grep -q "[s]16.repair --mode pass --pass A"; do sleep 30; done
python -m s16.repair --mode pass --pass B > s16/results/repair_B.log 2>&1
python -m s16.repair --mode pass --pass C > s16/results/repair_C.log 2>&1
python -m s16.repair --mode pass --pass D > s16/results/repair_D.log 2>&1
echo DONE-ALL-PASSES
