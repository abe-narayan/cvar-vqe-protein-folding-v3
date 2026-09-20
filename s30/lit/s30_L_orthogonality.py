"""What is ORTHOGONALITY actually worth? Closed form, on this project's own numbers.

Lane P's brief is "find a source whose errors are decorrelated from the distance prior".
This script prices that request. It is arithmetic on a closed form, NOT a measurement.

SETUP (S29 section 5.1, the project's own law).  Every native-free operator is a displacement
`d`; `u` is the unit direction to the native; the whole value of the operator is
`rho = cos(d, u)`, and `RMSD = RMSD_prod * sqrt(1 - rho^2)`.

With two fields v1, v2 having cos(v1,u)=r1, cos(v2,u)=r2, cos(v1,v2)=c, the best achievable
cosine from their span is the length of u's projection onto span{v1,v2}:

    rho_max^2 = [r1 r2] G^-1 [r1 r2]'  with G = [[1,c],[c,1]]
              = (r1^2 - 2 c r1 r2 + r2^2) / (1 - c^2)

Anchors used (all from s29/REPORT_S29.md section 0 and S30's lane F entry):
    RMSD_prod (built chain)            3.2105 A
    best single field (CHAN_DISTPOT)   rho = 0.1128
    random-shape reference             rho = 0.1398
    needed for 3.00 A                  rho = 0.358
    needed for 2.50 A                  rho = 0.628
    all 21 fields, ORACLE + 1 global weighting   rho = 0.169
    all 21 fields, leave-fold-out                rho = 0.012
    Gram stable rank of the 21 fields            2.057
"""
import numpy as np

RMSD_PROD = 3.2105
R1 = 0.1128          # best single field we own
NEED_300, NEED_250 = 0.358, 0.628

def rmsd(rho):
    return RMSD_PROD * np.sqrt(max(0.0, 1.0 - rho ** 2))

def rho_max(r1, r2, c):
    return np.sqrt(max(0.0, (r1**2 - 2*c*r1*r2 + r2**2) / (1 - c**2)))

print("=" * 78)
print("1. THE HEADLINE: what does PERFECT ORTHOGONALITY discount the requirement by?")
print("=" * 78)
for need, label in ((NEED_300, "3.00 A"), (NEED_250, "2.50 A")):
    # perfectly orthogonal new channel (c = 0): rho_max^2 = r1^2 + r2^2
    r2_needed = np.sqrt(need**2 - R1**2)
    print(f"  to reach {label} (rho {need:.3f}) combined with our best field (rho {R1:.4f}):")
    print(f"     a PERFECTLY ORTHOGONAL new channel must itself have rho = {r2_needed:.4f}")
    print(f"     alone, with no help, it would need                  rho = {need:.4f}")
    print(f"     -> orthogonality is worth a discount of {100*(1-r2_needed/need):.1f}% on the requirement")
    print(f"     -> and the new channel must be {r2_needed/R1:.2f}x better than anything we own")
    print()

print("=" * 78)
print("2. WHY: the gain from a second channel is QUADRATIC in ITS OWN skill")
print("=" * 78)
print("   For c = 0 and small r2:  rho_max ~ r1 + r2^2/(2 r1).  The cross term is ABSENT.")
print(f"{'r2':>8} {'rho_max':>9} {'RMSD':>8} {'gain A':>9}")
for r2 in (0.00, 0.05, 0.10, 0.1398, 0.20, 0.30, 0.34, 0.40, 0.50):
    rm = rho_max(R1, r2, 0.0)
    print(f"{r2:8.3f} {rm:9.4f} {rmsd(rm):8.4f} {rmsd(rm)-rmsd(R1):+9.4f}")
print()
print("   A new channel as good as EVERYTHING WE OWN (r2 = 0.1128), perfectly orthogonal,")
print(f"   buys {rmsd(rho_max(R1,R1,0.0)) - rmsd(R1):+.4f} A.  A channel at the random-shape")
print(f"   reference 0.1398 buys {rmsd(rho_max(R1,0.1398,0.0)) - rmsd(R1):+.4f} A.")

print()
print("=" * 78)
print("3. THE ORACLE-COMBINATION TRAP: why 21 fields reach 0.169 with the native in hand")
print("=" * 78)
print("   With r2 = 0 EXACTLY (a channel with NO skill at all):")
for c in (0.0, 0.3, 0.5, 0.7, 0.745, 0.8, 0.9):
    rm = rho_max(R1, 0.0, c)
    print(f"     c = {c:5.3f}  ->  rho_max = {rm:.4f}   (r1 alone = {R1:.4f}, "
          f"inflation {rm/R1:5.2f}x)")
print()
print("   A ZERO-SKILL field inflates the ORACLE combination by 1/sqrt(1-c^2), purely by")
print("   cancelling the part of v1 that is orthogonal to u.  That cancellation is fitted")
print("   WITH THE NATIVE and carries NO information.")
c_implied = np.sqrt(1 - (R1 / 0.169) ** 2)
print(f"   Lane F measured 0.169 ORACLE against {R1:.4f} best-single: inflation {0.169/R1:.2f}x,")
print(f"   which is exactly what a zero-skill partner at c = {c_implied:.3f} produces.")
print(f"   And leave-fold-out the same combination gives 0.012 -- BELOW the best single field.")
print("   => the ORACLE combination gain is a cancellation artefact, not an information gain.")

print()
print("=" * 78)
print("4. THE SAME POINT IN RMSD, AND THE BAR FOR LANE P")
print("=" * 78)
print(f"   production                         {RMSD_PROD:.4f} A")
print(f"   best single field we own           {rmsd(R1):.4f} A   (rho {R1:.4f})")
print(f"   all 21 fields, leave-fold-out      {rmsd(0.012):.4f} A   (rho 0.012)")
print(f"   all 21 fields, ORACLE + 1 weight   {rmsd(0.169):.4f} A   (rho 0.169, fitted on the native)")
print(f"   target 3.00 A                      needs rho {NEED_300:.3f}")
print(f"   target 2.50 A                      needs rho {NEED_250:.3f}")
print()
print("   BAR: a new source must carry rho ~ 0.34 BY ITSELF.  Perfect orthogonality to")
print("   everything we own reduces the requirement from 0.358 to 0.340 -- a 5% discount.")
print("   DECORRELATION IS NOT THE LEVER.  SKILL IS.")

print()
print("=" * 78)
print("5. HOW MANY INDEPENDENT DIRECTIONS WOULD IT TAKE, if each is as good as our best?")
print("=" * 78)
print("   k mutually orthogonal channels each at rho = 0.1128: rho_max = sqrt(k) * 0.1128")
for k in (1, 2, 4, 8, 10, 16, 21, 32):
    rm = min(1.0, np.sqrt(k) * R1)
    print(f"     k = {k:3d}  ->  rho_max = {rm:.4f}   RMSD = {rmsd(rm):.4f} A")
need_k = (NEED_300 / R1) ** 2
print(f"   reaching 3.00 A this way needs k = {need_k:.1f} MUTUALLY ORTHOGONAL channels")
print(f"   each as good as our best.  Lane F measured the 21 we have at stable rank 2.057.")
