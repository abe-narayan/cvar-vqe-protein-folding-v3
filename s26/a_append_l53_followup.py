"""s26/a_append_l53_followup.py -- lane A: the L81 follow-up now that the pool-spread control has run."""
from __future__ import annotations
import io, re, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LED = os.path.join(ROOT, "s26", "LEDGER.md")
SENTINEL = "FOLLOW-UP TO L81"

ENTRY = """## L{a} -- FOLLOW-UP TO L81 (L53, strain_difficulty): THE POOL-SPREAD CONTROL RAN. `moved` IS THE TOP-75's OWN DISAGREEMENT IN DISGUISE (rho 0.76 WITH THE SPREAD; PARTIAL +0.08 GIVEN IT, CI SPANNING ZERO); THE CALIBRATION FLAG STANDS, THE NOVELTY SENTENCE IS VETOED AS WORDED (2026-09-14, A)

`s26/a_strain_vs_spread.py` -> `s26/results/a_strain_vs_spread.json` (job `a_strain_vs_spread`,
relaunched under the file cap per L92; ORACLE DIAGNOSTIC: `rmsd_arm` is the label only; the
signals are native-free). The obvious native-free difficulty proxy L53's prereg did not carry:
the shipped top-75's own pairwise CA-RMSD spread (`I.pairwise_rmsd` over the 75 members; no
native, and available BEFORE the relaxation runs). Spearman with `rmsd_arm`, n = 126, partial
on n and Rg as in L53, 2,000-draw iid and fold-clustered bootstraps, permutation p:

    partial | n, Rg          spread_mean  +0.452  fold CI [+0.280, +0.609]  perm p < 0.0005  5/5
                             medoid_min   +0.451  fold CI [+0.286, +0.608]  perm p < 0.0005  5/5
                             moved        +0.433  fold CI [+0.247, +0.581]  perm p < 0.0005  5/5   (L53 reproduced)
                             log_e0       +0.241  fold CI [+0.054, +0.376]  perm p 0.007     4/5   (L53 reproduced)
    partial | n, Rg, spread  moved        +0.082  iid CI [-0.103, +0.259]  fold CI [+0.011, +0.171]  perm p 0.39  4/5
                             log_e0       +0.052  iid CI [-0.122, +0.231]  perm p 0.55  3/5
    rho(moved, spread_mean) = +0.756

Reading. The relaxation's displacement carries no information about the per-target error
beyond what the pool's own disagreement already carries: given the spread, `moved` is +0.08 with
the iid CI spanning zero and a permutation p of 0.39. L53's own mechanism sentence ("a
coordinate average that the projection turned into a strained chain is one whose pool members
disagreed, and disagreement is error") is confirmed literally, and it cuts the other way for
the novelty claim: the disagreement is measurable without AMBER, at rho +0.45, and is the
better flag (5/5 folds, tighter CI, zero cost). Rulings:

- L53's calibration flag STANDS as a phenomenon (a native-free quantity at rho +0.43 to +0.45
  with the emitted chain's error, 5/5 folds); the quartile table is a fair presentable form
  of it.
- The sentence "the first native-free quantity in this programme's record with a correlation
  above 0.4 to the per-target error of the emitted structure" is VETOED as worded: the
  quantity is the pool spread, `moved` is its proxy, and "the relaxation reports when the answer
  is untrustworthy" must become "the pool's disagreement reports it, and the relaxation's move
  tracks that disagreement at rho 0.76". Owner (PH) to answer in the ledger; the slide-7 note
  and the report sentence change accordingly.
- The physics reading is unchanged in the direction the record already holds: the relaxation
  adds nothing that the pool did not already say.

Provisional label on L81 lifted; verdict STANDS WITH CAVEAT (the caveat is the veto above).

"""


def main():
    s = io.open(LED, encoding="utf-8").read()
    if SENTINEL in s:
        print("sentinel present; nothing appended"); return 0
    n = max(int(m) for m in re.findall(r"^## L(\d+)", s, flags=re.M)) + 1
    sep = "" if s.endswith("\n\n") else ("\n" if s.endswith("\n") else "\n\n")
    io.open(LED, "a", encoding="utf-8", newline="\n").write(sep + ENTRY.format(a=n) + "---\n\n")
    print("appended L%d" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
