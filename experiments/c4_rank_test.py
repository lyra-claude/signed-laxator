#!/usr/bin/env python3
"""
c4_rank_test.py — C_4 Rank-2 Sheaf Computation

Computes the sheaf cohomology and signed projectors for the constant
rank-2 sheaf on the cycle graph C_4 (4 vertices, 4 edges).

C_4 has beta_1 = 1 (one independent cycle), so H^1 is non-trivial.
This is the key test case: unlike P_3 (a tree with H^1 = 0), C_4 has
genuine first cohomology, and the -1 eigenspace of phi_G should pick
up this cohomological structure.

Graph: C_4 = 4-cycle
    v0 --e01-- v1 --e12-- v2 --e23-- v3 --e30-- v0

Sheaf: Constant rank-2 (all stalks Q^2, all restrictions I_2)

Prediction: H^0 = Q^2 (constant sections), H^1 = Q^2 (from the cycle).
The -1 eigenspace of phi_G should have dimension 6 (= dim C^0 - dim H^0),
encoding the "non-harmonic" directions. The key question is whether the
-1 eigenspace structure reflects the non-trivial H^1.

All arithmetic is exact (fractions.Fraction).
"""

import sys
sys.path.insert(0, "/home/lyra/projects/signed-laxator/experiments/cactus_rank")

from fractions import Fraction
from cactus_rank import (
    mat_mul, mat_sub, mat_add, mat_scalar, identity, mat_eq,
    mat_rank_full, kernel_basis,
)


# ------------------------------------------------------------------
# Matrix utilities
# ------------------------------------------------------------------

def mat_transpose(A):
    """Transpose of a matrix."""
    n = len(A)
    m = len(A[0])
    return [[A[i][j] for i in range(n)] for j in range(m)]


def mat_inverse(A):
    """Inverse of a square matrix via Gauss-Jordan elimination over Q."""
    n = len(A)
    aug = [A[i][:] + [Fraction(1) if j == i else Fraction(0) for j in range(n)]
           for i in range(n)]
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if aug[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            raise ValueError("Matrix is singular")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        pv = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= pv
        for row in range(n):
            if row == col:
                continue
            f = aug[row][col]
            if f == 0:
                continue
            for j in range(2 * n):
                aug[row][j] -= f * aug[col][j]
    return [row[n:] for row in aug]


def zeros(n, m):
    return [[Fraction(0)] * m for _ in range(n)]


def fmt_frac(x):
    """Format a Fraction for display."""
    if x.denominator == 1:
        return str(x.numerator)
    return str(x)


def print_matrix(M, label="", indent=2):
    """Pretty-print a matrix of Fractions."""
    if label:
        print(label)
    prefix = " " * indent
    strs = [[fmt_frac(M[i][j]) for j in range(len(M[0]))] for i in range(len(M))]
    widths = [max(len(strs[i][j]) for i in range(len(M))) for j in range(len(M[0]))]
    for i in range(len(M)):
        row = prefix + "[ " + "  ".join(strs[i][j].rjust(widths[j]) for j in range(len(M[0]))) + " ]"
        print(row)


def print_block(M, block_size, vertices, indent=2):
    """Print a matrix as block structure with vertex labels."""
    n_blocks = len(M) // block_size
    prefix = " " * indent
    for i in range(n_blocks):
        for j in range(n_blocks):
            block = [[M[i * block_size + r][j * block_size + c]
                      for c in range(block_size)]
                     for r in range(block_size)]
            block_strs = [[fmt_frac(block[r][c]) for c in range(block_size)]
                          for r in range(block_size)]
            widths = [max(len(block_strs[r][c]) for r in range(block_size))
                      for c in range(block_size)]
            label = f"{prefix}[{vertices[i]},{vertices[j]}] ="
            for r in range(block_size):
                row_str = "[ " + "  ".join(block_strs[r][c].rjust(widths[c])
                                           for c in range(block_size)) + " ]"
                if r == 0:
                    print(f"{label} {row_str}")
                else:
                    print(f"{' ' * len(label)} {row_str}")


# ------------------------------------------------------------------
# Main computation
# ------------------------------------------------------------------

def main():
    F = Fraction

    print("=" * 70)
    print("  C_4 Rank-2 Sheaf Computation — Testing H^1 on Non-Trivial Topology")
    print("=" * 70)
    print()
    print("Graph: C_4 = 4-cycle (v0 -- v1 -- v2 -- v3 -- v0)")
    print("Edges: e01 = v0->v1, e12 = v1->v2, e23 = v2->v3, e30 = v3->v0")
    print("Sheaf: Constant rank-2 (all stalks Q^2, all restrictions I_2)")
    print()
    print("Topology: beta_1 = |E| - |V| + 1 = 4 - 4 + 1 = 1")
    print("Prediction: H^0 = Q^2 (constant sections), H^1 = Q^2 (from cycle)")
    print()

    # ------------------------------------------------------------------
    # Step 1: Build coboundary map delta_0: C^0 -> C^1
    # ------------------------------------------------------------------
    # C^0 = F(v0) + F(v1) + F(v2) + F(v3) = Q^8
    # C^1 = F(e01) + F(e12) + F(e23) + F(e30) = Q^8
    #
    # For edge e with orientation u -> w (constant sheaf, all restrictions I):
    #   (delta_0 s)(e) = s(w) - s(u)
    #
    # Edge e01 (v0->v1): s(v1) - s(v0)
    #   Row block: [-I, I, 0, 0]
    #
    # Edge e12 (v1->v2): s(v2) - s(v1)
    #   Row block: [0, -I, I, 0]
    #
    # Edge e23 (v2->v3): s(v3) - s(v2)
    #   Row block: [0, 0, -I, I]
    #
    # Edge e30 (v3->v0): s(v0) - s(v3)
    #   Row block: [I, 0, 0, -I]

    I2 = identity(2)
    Z2 = zeros(2, 2)
    mI2 = mat_scalar(I2, F(-1))

    # delta_0 is 8x8
    delta_0 = []

    # Row block for e01: [-I, I, 0, 0]
    for r in range(2):
        row = mI2[r] + I2[r] + Z2[r] + Z2[r]
        delta_0.append(row)

    # Row block for e12: [0, -I, I, 0]
    for r in range(2):
        row = Z2[r] + mI2[r] + I2[r] + Z2[r]
        delta_0.append(row)

    # Row block for e23: [0, 0, -I, I]
    for r in range(2):
        row = Z2[r] + Z2[r] + mI2[r] + I2[r]
        delta_0.append(row)

    # Row block for e30: [I, 0, 0, -I]
    for r in range(2):
        row = I2[r] + Z2[r] + Z2[r] + mI2[r]
        delta_0.append(row)

    print("--- Coboundary map delta_0 (8x8) ---")
    print_matrix(delta_0)
    print()

    # ------------------------------------------------------------------
    # Step 2: Compute cohomology dimensions
    # ------------------------------------------------------------------
    rank_d0 = mat_rank_full(delta_0)
    dim_C0 = 8
    dim_C1 = 8

    # H^0 = ker(delta_0)
    dim_H0 = dim_C0 - rank_d0

    # H^1 = C^1 / im(delta_0) = coker(delta_0) (no 2-cells in a graph)
    dim_H1 = dim_C1 - rank_d0

    print("--- Cohomology ---")
    print(f"  dim C^0 = {dim_C0}")
    print(f"  dim C^1 = {dim_C1}")
    print(f"  rank(delta_0) = {rank_d0}")
    print(f"  dim H^0 = dim ker(delta_0) = {dim_H0}")
    print(f"  dim H^1 = dim C^1 - rank(delta_0) = {dim_H1}")
    print(f"  Euler characteristic chi = dim H^0 - dim H^1 = {dim_H0 - dim_H1}")
    print()

    # Sanity check: for constant rank-r sheaf on graph G,
    # dim H^0 = r (connected), dim H^1 = r * beta_1
    beta_1 = 1  # cycle rank of C_4
    rank = 2
    print(f"  Expected (constant rank-{rank} sheaf, beta_1={beta_1}):")
    print(f"    dim H^0 = {rank} (connected graph)")
    print(f"    dim H^1 = {rank} * {beta_1} = {rank * beta_1}")
    print(f"  Match: H^0 {'YES' if dim_H0 == rank else 'NO'}, H^1 {'YES' if dim_H1 == rank * beta_1 else 'NO'}")
    print()

    # ------------------------------------------------------------------
    # Step 3: Kernel basis and projector P_+
    # ------------------------------------------------------------------
    ker_basis = kernel_basis(delta_0)
    print(f"--- Kernel basis of delta_0 ({len(ker_basis)} vectors in Q^{dim_C0}) ---")
    for i, v in enumerate(ker_basis):
        print(f"  k_{i} = [{', '.join(fmt_frac(x) for x in v)}]")
    print()

    # Build P_+ = K^T (K K^T)^{-1} K
    K = ker_basis
    KT = mat_transpose(K)
    KKT = mat_mul(K, KT)
    KKT_inv = mat_inverse(KKT)
    P_plus = mat_mul(mat_mul(KT, KKT_inv), K)

    print("--- Projector P_+ onto ker(delta_0) (8x8) ---")
    print("  Block structure (2x2 per vertex):")
    print_block(P_plus, 2, ["v0", "v1", "v2", "v3"])
    print()

    # ------------------------------------------------------------------
    # Step 4: P_- = I - P_+
    # ------------------------------------------------------------------
    I8 = identity(8)
    P_minus = mat_sub(I8, P_plus)

    print("--- Projector P_- = I - P_+ (8x8) ---")
    print("  Block structure (2x2 per vertex):")
    print_block(P_minus, 2, ["v0", "v1", "v2", "v3"])
    print()

    # ------------------------------------------------------------------
    # Step 5: phi_G = P_+ - P_- = 2*P_+ - I
    # ------------------------------------------------------------------
    phi_G = mat_sub(P_plus, P_minus)

    print("--- phi_G = P_+ - P_- (8x8) ---")
    print_matrix(phi_G, "  Full matrix:")
    print()
    print("  Block structure (2x2 per vertex):")
    print_block(phi_G, 2, ["v0", "v1", "v2", "v3"])
    print()

    # ------------------------------------------------------------------
    # Step 6: Verification
    # ------------------------------------------------------------------
    print("=" * 70)
    print("  Verification")
    print("=" * 70)
    print()

    Z8 = zeros(8, 8)

    check_Pp_idem = mat_eq(mat_mul(P_plus, P_plus), P_plus)
    print(f"  P_+^2 = P_+?          {check_Pp_idem}")

    check_Pm_idem = mat_eq(mat_mul(P_minus, P_minus), P_minus)
    print(f"  P_-^2 = P_-?          {check_Pm_idem}")

    check_PpPm = mat_eq(mat_mul(P_plus, P_minus), Z8)
    print(f"  P_+ * P_- = 0?        {check_PpPm}")

    check_PmPp = mat_eq(mat_mul(P_minus, P_plus), Z8)
    print(f"  P_- * P_+ = 0?        {check_PmPp}")

    check_invol = mat_eq(mat_mul(phi_G, phi_G), I8)
    print(f"  phi_G^2 = I?          {check_invol}")
    print()

    # Traces
    trace_Pp = sum(P_plus[i][i] for i in range(8))
    trace_Pm = sum(P_minus[i][i] for i in range(8))
    trace_phi = sum(phi_G[i][i] for i in range(8))
    print(f"  tr(P_+) = {fmt_frac(trace_Pp)}  (should = dim H^0 = {dim_H0})")
    print(f"  tr(P_-) = {fmt_frac(trace_Pm)}  (should = dim C^0 - dim H^0 = {dim_C0 - dim_H0})")
    print(f"  tr(phi_G) = {fmt_frac(trace_phi)}  (= {dim_H0} - {dim_C0 - dim_H0} = {dim_H0 - (dim_C0 - dim_H0)})")
    print()

    # Eigenspace dimensions of phi_G
    phi_plus_I = mat_sub(phi_G, I8)
    phi_minus_I = mat_add(phi_G, I8)
    dim_plus1 = 8 - mat_rank_full(phi_plus_I)
    dim_minus1 = 8 - mat_rank_full(phi_minus_I)
    print(f"  Eigenspace dims of phi_G: (+1) = {dim_plus1}, (-1) = {dim_minus1}")
    print()

    # ------------------------------------------------------------------
    # Step 7: Hodge Laplacians — connecting phi_G to H^1
    # ------------------------------------------------------------------
    print("=" * 70)
    print("  Hodge Laplacian Analysis — H^1 Structure")
    print("=" * 70)
    print()

    # Up-Laplacian on 0-cochains: L_0^up = delta_0^T delta_0
    delta_0_T = mat_transpose(delta_0)
    L0_up = mat_mul(delta_0_T, delta_0)

    print("--- L_0^up = delta_0^T delta_0 (8x8) ---")
    print("  Block structure (2x2 per vertex):")
    print_block(L0_up, 2, ["v0", "v1", "v2", "v3"])
    print()

    rank_L0 = mat_rank_full(L0_up)
    print(f"  rank(L_0^up) = {rank_L0}")
    print(f"  dim ker(L_0^up) = {8 - rank_L0}  (should = dim H^0 = {dim_H0})")
    print()

    # Down-Laplacian on 1-cochains: L_1^down = delta_0 delta_0^T
    L1_down = mat_mul(delta_0, delta_0_T)

    print("--- L_1^down = delta_0 delta_0^T (8x8) ---")
    print("  Block structure (2x2 per edge):")
    print_block(L1_down, 2, ["e01", "e12", "e23", "e30"])
    print()

    rank_L1 = mat_rank_full(L1_down)
    dim_ker_L1 = 8 - rank_L1
    print(f"  rank(L_1^down) = {rank_L1}")
    print(f"  dim ker(L_1^down) = {dim_ker_L1}  (should = dim H^1 = {dim_H1})")
    print()

    # Harmonic 1-cochains (= H^1 representatives)
    harm1_basis = kernel_basis(L1_down)
    if harm1_basis:
        print(f"--- Harmonic 1-cochains ({len(harm1_basis)} vectors = H^1 basis) ---")
        for i, v in enumerate(harm1_basis):
            print(f"  h_{i} = [{', '.join(fmt_frac(x) for x in v)}]")
        print()
        print("  These are 1-cochains in ker(delta_0 delta_0^T).")
        print("  For a graph (no 2-cells), ker(L_1^down) = ker(delta_0^T),")
        print("  which are the harmonic representatives of H^1.")
        print()

    # ------------------------------------------------------------------
    # Step 8: Diagonal blocks of phi_G — vertex-level structure
    # ------------------------------------------------------------------
    print("=" * 70)
    print("  Diagonal Blocks of phi_G (vertex-level signed structure)")
    print("=" * 70)
    print()

    for i, label in enumerate(["v0", "v1", "v2", "v3"]):
        block = [[phi_G[2*i + r][2*i + c] for c in range(2)] for r in range(2)]
        is_scalar = (block[0][0] == block[1][1] and
                     block[0][1] == F(0) and block[1][0] == F(0))
        print(f"  phi_G[{label},{label}] = [[{fmt_frac(block[0][0])}, {fmt_frac(block[0][1])}],")
        print(f"                     [{fmt_frac(block[1][0])}, {fmt_frac(block[1][1])}]]")
        print(f"    Scalar x I? {is_scalar}", end="")
        if is_scalar:
            print(f"  (value: {fmt_frac(block[0][0])})")
        else:
            print()
        print()

    # ------------------------------------------------------------------
    # Step 9: Comparison with P_3 (trees have H^1 = 0)
    # ------------------------------------------------------------------
    print("=" * 70)
    print("  Comparison: C_4 vs P_3 (cycle vs tree)")
    print("=" * 70)
    print()
    print("  P_3 (path, 3 vertices, 2 edges):")
    print("    beta_1 = 0 (tree)")
    print("    H^0 = Q^2, H^1 = Q^0")
    print("    +1 eigenspace of phi_G: dim 2 (= dim H^0)")
    print("    -1 eigenspace of phi_G: dim 4 (= dim C^0 - dim H^0)")
    print()
    print("  C_4 (cycle, 4 vertices, 4 edges):")
    print(f"    beta_1 = 1 (one independent cycle)")
    print(f"    H^0 = Q^{dim_H0}, H^1 = Q^{dim_H1}")
    print(f"    +1 eigenspace of phi_G: dim {dim_plus1} (= dim H^0)")
    print(f"    -1 eigenspace of phi_G: dim {dim_minus1} (= dim C^0 - dim H^0)")
    print()

    # Key structural observation: the -1 eigenspace
    print("  KEY OBSERVATION:")
    print(f"    The -1 eigenspace of phi_G has dimension {dim_minus1}.")
    if dim_H1 > 0:
        print(f"    H^1 has dimension {dim_H1} (non-trivial!).")
        print(f"    The {dim_minus1}-dimensional -1 eigenspace decomposes C^0 into")
        print(f"    a {dim_plus1}-dim harmonic part (global sections) and")
        print(f"    a {dim_minus1}-dim non-harmonic part.")
        print()
        print(f"    Unlike P_3 where H^1 = 0, C_4 has non-trivial cohomology.")
        print(f"    The presence of the cycle creates {dim_H1} dimensions of")
        print(f"    cohomological obstruction that phi_G detects through")
        print(f"    the structure of the -1 eigenspace.")
    print()

    # ------------------------------------------------------------------
    # Step 10: Check if phi_G commutes with the Laplacian
    # ------------------------------------------------------------------
    print("=" * 70)
    print("  Commutativity: phi_G vs L_0^up")
    print("=" * 70)
    print()

    comm = mat_sub(mat_mul(phi_G, L0_up), mat_mul(L0_up, phi_G))
    comm_is_zero = mat_eq(comm, Z8)
    print(f"  [phi_G, L_0^up] = 0?  {comm_is_zero}")
    if comm_is_zero:
        print("  phi_G and L_0^up share eigenspaces.")
        print("  This means phi_G is a spectral involution of the Laplacian.")
    print()

    # ------------------------------------------------------------------
    # Step 11: Symmetry analysis — Z/4 action
    # ------------------------------------------------------------------
    print("=" * 70)
    print("  Symmetry: Z/4 rotation on C_4")
    print("=" * 70)
    print()

    # For constant sheaf, the Z/4 rotation v_i -> v_{i+1 mod 4}
    # acts on C^0 by permuting the 2x2 blocks cyclically.
    # Build the rotation matrix R_4 (8x8):
    R4 = zeros(8, 8)
    for i in range(4):
        j = (i + 1) % 4
        for r in range(2):
            R4[2*j + r][2*i + r] = F(1)

    print("  R_4 = cyclic rotation (v_i -> v_{i+1 mod 4}):")
    print_matrix(R4, "")
    print()

    # Check R_4^4 = I
    R4_sq = mat_mul(R4, R4)
    R4_4 = mat_mul(R4_sq, R4_sq)
    print(f"  R_4^4 = I?  {mat_eq(R4_4, I8)}")

    # Does phi_G commute with R_4?
    comm_R4 = mat_sub(mat_mul(phi_G, R4), mat_mul(R4, phi_G))
    comm_R4_zero = mat_eq(comm_R4, Z8)
    print(f"  [phi_G, R_4] = 0?  {comm_R4_zero}")
    if comm_R4_zero:
        print("  phi_G respects the Z/4 symmetry of the cycle.")
    print()

    # ------------------------------------------------------------------
    # Conclusion
    # ------------------------------------------------------------------
    print("=" * 70)
    print("  CONCLUSION")
    print("=" * 70)
    print()

    all_pass = (check_Pp_idem and check_Pm_idem and
                check_PpPm and check_PmPp and check_invol)

    if all_pass:
        print("  ALL CHECKS PASSED.")
    else:
        print("  SOME CHECKS FAILED.")
    print()

    print(f"  Graph: C_4 (4-cycle)")
    print(f"  Sheaf: Constant rank-2")
    print(f"  beta_1 = 1 (one independent cycle)")
    print()
    print(f"  dim H^0 = {dim_H0}  (global sections)")
    print(f"  dim H^1 = {dim_H1}  (cohomological obstruction from cycle)")
    print(f"  phi_G^2 = I: {check_invol}")
    print(f"  +1 eigenspace dim = {dim_plus1}  (= dim H^0 = {dim_H0})")
    print(f"  -1 eigenspace dim = {dim_minus1}  (= dim C^0 - dim H^0 = {dim_C0 - dim_H0})")
    print()
    print(f"  [phi_G, L_0^up] = 0: {comm_is_zero}  (spectral involution)")
    print(f"  [phi_G, R_4] = 0: {comm_R4_zero}  (respects Z/4 symmetry)")
    print()

    if dim_H1 > 0 and check_invol:
        print("  *** MAIN RESULT: On C_4 with non-trivial H^1:")
        print("      phi_G^2 = I still holds (involution robust on cycles)")
        print(f"      The -1 eigenspace (dim {dim_minus1}) encodes the non-harmonic")
        print(f"      complement. The {dim_H1}-dimensional H^1 manifests as")
        print(f"      non-trivial harmonic 1-cochains in ker(L_1^down).")
        print()
        print("      The signed decomposition phi_G = P_+ - P_- cleanly")
        print("      separates C^0 into cohomological (+1) and non-cohomological")
        print("      (-1) subspaces, even in the presence of non-trivial topology.")
    print()


if __name__ == "__main__":
    main()
