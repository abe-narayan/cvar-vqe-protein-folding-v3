#!/bin/sh
# Waits for pass 1 to exit, then runs the frame null, the flat-bottom cross-check,
# and finally the AMBER single points for the identical-candidate experiment (E4).
cd "$(dirname "$0")/.."
while ps -W 2>/dev/null | grep -q "phys_ca"; do sleep 30; done
python -m s17.phys_ca --pass 4 >> s17/results/phys_ca_p4.log 2>&1
python -m s17.phys_ident --out phys_ident_amb.json >> s17/results/ident_amb.log 2>&1
python -m s17.phys_ca --pass 3 >> s17/results/phys_ca_p3.log 2>&1
echo CHAIN2_DONE
