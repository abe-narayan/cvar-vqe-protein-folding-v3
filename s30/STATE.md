# S30 — RUNNING STATE

Last update: sprint open.
Endpoint: mean built-chain Cα RMSD, 126 targets. Production incumbent **3.2105 Å**.

---

## HEADLINE: THE MEAN IS A TAIL STATISTIC, AND THE TAIL IS WHERE THE PRIZE IS

Computed at the sprint open from lane O's S29 rows, before any new experiment:

```
production, n = 126        mean 3.2105    median 2.9661
the worst 18 targets       mean 6.2758
the other 108              mean 2.6997
```

**Counterfactuals on the charter's own endpoint:**

```
cap the worst 10 at 3.00 Å   ->  mean 2.9074   (-0.3031)   BEATS the primary target
cap the worst 18 at 3.00 Å   ->  mean 2.7426   (-0.4680)
cap the worst 30 at 3.00 Å   ->  mean 2.5778   (-0.6328)   approaches the ambitious target
worst 18 all the way to 2.50 ->  mean 2.1334   (-1.0772)

versus: improve EVERY ONE of the 126 by 0.20 Å  ->  mean 3.0105  (-0.2000)
```

> **Fixing ten targets is worth more than improving all 126 by 0.20 Å.**

That single comparison should govern how this sprint allocates effort. A mechanism that works on
the typical target and leaves the tail alone is close to worthless *for this endpoint*; a mechanism
that only works on hard targets can hit the primary target on its own.

**The caveat, stated up front so the number is not over-read.** Capping is an ORACLE operation — it
requires knowing which targets are bad. It bounds the prize; it does not deliver it. Realising any
of it needs either (a) a native-free regime detector, or (b) a method that is simply better on hard
pools without needing to know they are hard. **(b) does not require detection and is therefore the
stronger target.**

**And there is a live, unexploited clue.** The record says the shipped pipeline is **worse than a
sequence-blind one on its 18 hardest targets** — 5.425 blind against 6.019 shipped — while
sequence conditioning is worth 0.776 Å overall. Something the pipeline does on hard targets is
actively harmful. That is a mechanism to find, not a curiosity.

---

## WHAT I AM TREATING AS BINDING FROM S29 (until a lane breaks it)

- Every native-free operator is a displacement; its whole value is one cosine. 3.00 Å needs
  ρ = 0.358, 2.50 needs 0.628; everything measured is ≤ 0.04 **signed** — but per-target |cos| is
  0.25–0.33, so **the magnitude exists and only the sign is missing**.
- The deployed architecture's ORACLE ceiling is 2.9027 Å built chain. 2.5 Å is unreachable through
  it *with the native in hand*.
- The terminal operator consumes the set **mean**, so a perfect rank-1 is worth ≈ −0.03 Å through
  the shipped m = 75 average. **A better channel needs a different readout to be worth anything.**
- The deployed score supplies 1.442 of 7 bits against 1.405 for random; the median target is
  *worse* than chance.
- Recognition is closed across bands (theorem), within bands (F2 fails 58/70) and per target
  (incidental parameter) — **but not in the class-B2-is-empty sense**: 9 of 31 channels keep skill
  after partialling out compactness both ways.

## WHAT S29 LEFT GENUINELY OPEN

1. **LEG_torsion and 8 sibling channels** — in-band skill that survives removing monotone *and*
   V-shaped Rg dependence; LEG_torsion is a function of the structure, not the distogram, so it is
   outside the class theorem 2 bounds.
2. **A sparse weighted readout** — the only ladder class not closed by ceiling (2 members with
   oracle weights reach 1.43 Å). Needs a native-free support rule and ~18 bits, not 7.
3. **A free-energy stage** — the only native-free selector class with peptide-length precedent. It
   does not exist on disk; a rebuild, not a resume.

## LANES

| lane | remit | status |
|---|---|---|
| **R** | is nativeness recognizable from one structure at all (L11) | running |
| **L** | literature, permanent | running |
| **D** | adversary, permanent; owns the cost/RMSD meter | running |
| **T** | theory: the bit accounting (L8) and when the CVaR tail stops being a prefix | running |

## OPEN QUESTIONS I AM HOLDING

- Does the tail's difficulty have a *cause* we can name, or is it just harder sequences?
- Is there a formulation that dissolves several leads at once rather than patching each?
- The launch-gate defect (`jobrun.py` CPU_START 85 vs `governor.py` CPU_CEILING 101) is still open
  and unresolved from S29; RAM, not CPU, is currently the binding constraint so it has not bitten
  yet this sprint.
