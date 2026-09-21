# S31 — LANE V, THE ADVERSARY

Contract rule 24: the main team may not approve its own positive. This file is the adversary's
record. **Defects lead.** Every entry states (a) what is claimed and where, (b) what the artefact
says, (c) the corrected sentence.

Read-only on `s29/ s28/ s27/ s16/ s12/ core/` and on other lanes' files. Writes confined to this
file, `s31/results/s31_V_*`, `s31/s31_V_*.py`, and ledger entries.

---

## 0. PRE-REGISTRATION — V1, THE MATCHED-K CONTROL FOR THE SET-MATCHED READOUT LADDER

**Committed before the number exists (contract rule 19). Code `s31/s31_V_orderstat.py`.**

`s31/STATE.md` NOTE 2(b) and `s31/LEDGER.md` S31-L6 §1 read the fixed-top-128 ladder as

```
argmin over 128            7.0 bits              2.1458
2-of-128, UNIFORM weights  12.99 support bits    2.0700    "-0.076 vs argmin, ZERO weight bits"
```

with the coordinator's handed-back reading that *the 0.076 Å is **error cancellation**, not
weighting*.

**The objection.** `2-of-128` is a per-target minimum over `C(128,2) = 8128` ORACLE-chosen
supports; `argmin over 128` is a per-target minimum over **128**. That is a **64× larger oracle
search**, compared as though the two arms differed only in *what they emit*. Project memory
`grid-oracles-are-order-statistics` and contract rule 11 both bind here, and lane F has just
measured (S31-L11) that on this very candidate set the per-target minimum is **monotone and
unsaturated in K**.

**V1 registered bar.** *"The 0.076 Å is an order statistic, not error cancellation"* **FIRES** if
the per-target minimum over **128 random pairs** — matched K, matched operator, same top-128,
uniform weights, same native — fails to beat the minimum over the 128 singletons by at least
`0.7 × MDE`. In that case the 0.076 Å may not be attributed to error cancellation without the
matched-K number beside it.

**Registered prior:** ~3:1 that most of the 0.076 Å is K, because the pair family is 64× larger
and its members are far less correlated with one another than nested prefixes are.

**Basis:** CA point cloud throughout. Every arm is **ORACLE / NOT DEPLOYABLE** — each is a
per-target minimum taken against the native. **No built-chain claim is made in V1.**

Draw discipline: `NDRAW = 8` independent draws, **mean over draws, never the maximum**
(contract rule 10).

---

*(Findings follow below as they land.)*
