#!/usr/bin/env python3
"""
cactus_rank.py  (pure-Python fallback; no SageMath available in container)

Computes cactus-group-generator ranks on irreducible Specht modules of S_n
for n = 4, 5.

Specht modules are constructed via Young's seminormal (orthogonal) form:
for each standard Young tableau T of shape lambda we have a basis vector v_T,
and adjacent transpositions s_i = (i, i+1) act by

    s_i v_T = (1/d) v_T + sqrt(1 - 1/d^2) v_{s_i T}      (if s_i T is standard)
    s_i v_T = (1/d) v_T                                    (if s_i T is not standard,
                                                            i.e. i, i+1 in same row/col)

where d = c_T(i+1) - c_T(i) is the axial distance from box i to box i+1
(content difference).  For the seminormal (rational) version we drop the
sqrt and use:

    s_i v_T = (1/d) v_T + v_{s_i T}                          if s_i T standard
    s_i v_{s_i T} = -(1/d) v_{s_i T} + (1 - 1/d^2) v_T

This is the classical Young's seminormal representation; the matrix of any
permutation is built by reducing it to a product of adjacent transpositions
and multiplying.  The eigenvalues of an involution are still ±1 (the
representation is similar to an orthogonal one, so eigenvalues are real and
the involutivity carries over).

For each cactus generator s_{p,q} -> w_0^{[p,q]}, we record:
  - a_lambda(p,q) = dim +1 eigenspace
  - b_lambda(p,q) = dim -1 eigenspace
For each pair we record commutator rank, joint fixed dim.
"""

from fractions import Fraction
from itertools import combinations, permutations
import csv


# ----------------------------------------------------------------------
# Partitions & SYTs
# ----------------------------------------------------------------------

def partitions(n, max_part=None):
    """All partitions of n as non-increasing tuples."""
    if n == 0:
        return [()]
    if max_part is None:
        max_part = n
    out = []
    for k in range(min(n, max_part), 0, -1):
        for rest in partitions(n - k, k):
            out.append((k,) + rest)
    return out


def shape_cells(lam):
    """Return list of (row, col) cells (0-indexed) in Young diagram."""
    cells = []
    for r, p in enumerate(lam):
        for c in range(p):
            cells.append((r, c))
    return cells


def standard_young_tableaux(lam):
    """Generate all standard Young tableaux of shape lam.
    Each SYT is returned as a dict mapping number -> (row, col)."""
    cells = shape_cells(lam)
    n = len(cells)

    # We'll build a tableau as a list of lists, filling 1..n.
    rows = [list(p for p in lam)]  # row sizes
    grid = [[0] * p for p in lam]  # 0 = unfilled
    results = []

    def addable_cells():
        """Cells where we can place the next number (corners of current shape)."""
        out = []
        for r in range(len(lam)):
            for c in range(lam[r]):
                if grid[r][c] != 0:
                    continue
                # Must be that all cells (r,c-1) and (r-1,c) are filled
                if c > 0 and grid[r][c - 1] == 0:
                    continue
                if r > 0 and grid[r - 1][c] == 0:
                    continue
                out.append((r, c))
        return out

    def recurse(k):
        if k > n:
            # Build dict
            d = {}
            for r in range(len(lam)):
                for c in range(lam[r]):
                    d[grid[r][c]] = (r, c)
            results.append(d)
            return
        for (r, c) in addable_cells():
            grid[r][c] = k
            recurse(k + 1)
            grid[r][c] = 0

    recurse(1)
    return results


def syt_to_tuple(T):
    """Hashable representation: tuple of (row,col) for each k=1..n."""
    n = len(T)
    return tuple(T[k] for k in range(1, n + 1))


# ----------------------------------------------------------------------
# Young's seminormal representation
# ----------------------------------------------------------------------

def axial_distance(T, i):
    """Content difference c(i+1) - c(i) where content(r,c) = c - r."""
    r1, c1 = T[i]
    r2, c2 = T[i + 1]
    return (c2 - r2) - (c1 - r1)


def apply_si(T, i):
    """Swap labels i and i+1 in the tableau dict T.  Returns new dict (does
    not check standardness)."""
    Tn = dict(T)
    Tn[i], Tn[i + 1] = T[i + 1], T[i]
    return Tn


def is_standard(T, lam):
    """Check that values increase along rows and columns."""
    # Build grid
    grid = [[0] * p for p in lam]
    for k, (r, c) in T.items():
        grid[r][c] = k
    for r in range(len(lam)):
        for c in range(lam[r]):
            if c + 1 < lam[r] and grid[r][c] >= grid[r][c + 1]:
                return False
            if r + 1 < len(lam) and r + 1 < len(lam) and c < lam[r + 1] and grid[r][c] >= grid[r + 1][c]:
                return False
    return True


def seminormal_si_matrix(lam):
    """Return Young's seminormal matrix of s_i for i = 1..n-1 as a list of
    matrices (each n_lam x n_lam) over Fraction.  Indexed by adjacent
    transposition index i (1-based)."""
    syts = standard_young_tableaux(lam)
    index = {syt_to_tuple(T): k for k, T in enumerate(syts)}
    d = len(syts)
    n = sum(lam)
    matrices = {}
    for i in range(1, n):
        M = [[Fraction(0)] * d for _ in range(d)]
        for k, T in enumerate(syts):
            ad = axial_distance(T, i)
            inv = Fraction(1, ad)
            # Diagonal
            M[k][k] = inv
            # Off-diagonal: T -> s_i T (if standard)
            Tn = apply_si(T, i)
            if is_standard(Tn, lam):
                kn = index[syt_to_tuple(Tn)]
                # Coefficient of v_{s_i T} in s_i v_T is 1 (seminormal, non-symmetric form)
                M[k][kn] = Fraction(1)
                # And we need to check this from the s_i T side too — the matrix
                # entry M[kn][k] gets coefficient (1 - 1/ad^2) = (ad^2 - 1)/ad^2
                # but that's filled when we iterate k over the partner.
                # Actually, axial_distance from s_i T side = -ad, so when we get
                # to k = kn we'll set M[kn][kn] = -inv (correct) and
                # M[kn][k] = 1.  But that's WRONG — the seminormal form is
                # asymmetric: s_i v_{s_iT} should produce (1-1/d^2) v_T, not 1.
                # Let me use the SYMMETRIC (orthogonal) seminormal where both
                # off-diagonal entries are sqrt(1 - 1/d^2) — but we want to
                # stay in QQ.  The classical *seminormal* (Young's first) form
                # has M[k][kn] = 1 and M[kn][k] = 1 - 1/d^2.  We should be
                # careful.
        matrices[i] = M
    # Actually, let me redo this with the classical asymmetric formula.
    matrices = {}
    for i in range(1, n):
        M = [[Fraction(0)] * d for _ in range(d)]
        for k, T in enumerate(syts):
            ad = axial_distance(T, i)
            inv = Fraction(1, ad)
            M[k][k] = inv
            Tn = apply_si(T, i)
            if is_standard(Tn, lam):
                kn = index[syt_to_tuple(Tn)]
                # Classical seminormal: s_i v_T = (1/d) v_T + v_{s_iT}  if d>0
                #                       s_i v_T = (1/d) v_T + (1 - 1/d^2) v_{s_iT}  if d<0
                # i.e. when T < s_iT lexicographically (d>0) we get coeff 1,
                # when T > s_iT we get coeff (1 - 1/d^2).
                if ad > 0:
                    M[k][kn] = Fraction(1)
                else:
                    M[k][kn] = Fraction(1) - Fraction(1, ad * ad)
        matrices[i] = M
    return matrices, syts


# ----------------------------------------------------------------------
# Linear algebra over Fraction
# ----------------------------------------------------------------------

def mat_mul(A, B):
    n = len(A)
    m = len(B[0])
    k_ = len(B)
    C = [[Fraction(0)] * m for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        Ci = C[i]
        for k in range(k_):
            a = Ai[k]
            if a == 0:
                continue
            Bk = B[k]
            for j in range(m):
                Ci[j] += a * Bk[j]
    return C


def mat_sub(A, B):
    n = len(A)
    m = len(A[0])
    return [[A[i][j] - B[i][j] for j in range(m)] for i in range(n)]


def mat_add(A, B):
    n = len(A)
    m = len(A[0])
    return [[A[i][j] + B[i][j] for j in range(m)] for i in range(n)]


def mat_scalar(A, s):
    return [[A[i][j] * s for j in range(len(A[0]))] for i in range(len(A))]


def identity(n):
    return [[Fraction(1) if i == j else Fraction(0) for j in range(n)] for i in range(n)]


def mat_eq(A, B):
    if len(A) != len(B):
        return False
    for i in range(len(A)):
        for j in range(len(A[0])):
            if A[i][j] != B[i][j]:
                return False
    return True


def mat_rank(A):
    """Rank by Gaussian elimination over Fraction."""
    M = [row[:] for row in A]
    n = len(M)
    if n == 0:
        return 0
    m = len(M[0])
    rank = 0
    pivot_col = 0
    for r in range(n):
        if pivot_col >= m:
            break
        # Find pivot
        pr = None
        for rr in range(r, n):
            if M[rr][pivot_col] != 0:
                pr = rr
                break
        if pr is None:
            pivot_col += 1
            # Decrement r so we retry this row with next column
            # Easier: implement as while loop
            return mat_rank_full(A)
        # Swap
        M[r], M[pr] = M[pr], M[r]
        # Eliminate
        pv = M[r][pivot_col]
        for rr in range(n):
            if rr == r:
                continue
            f = M[rr][pivot_col] / pv
            if f == 0:
                continue
            for cc in range(pivot_col, m):
                M[rr][cc] -= f * M[r][cc]
        rank += 1
        pivot_col += 1
    return rank


def mat_rank_full(A):
    """Robust rank via column-search Gaussian elimination."""
    M = [row[:] for row in A]
    if not M:
        return 0
    n = len(M)
    m = len(M[0])
    rank = 0
    col = 0
    row = 0
    while row < n and col < m:
        pr = None
        for rr in range(row, n):
            if M[rr][col] != 0:
                pr = rr
                break
        if pr is None:
            col += 1
            continue
        M[row], M[pr] = M[pr], M[row]
        pv = M[row][col]
        for rr in range(n):
            if rr == row:
                continue
            f = M[rr][col] / pv
            if f == 0:
                continue
            for cc in range(col, m):
                M[rr][cc] -= f * M[row][cc]
        rank += 1
        row += 1
        col += 1
    return rank


def kernel_dim(A):
    return len(A[0]) - mat_rank_full(A) if A else 0


def kernel_basis(A):
    """Return basis (list of vectors) of right null space of A using RREF."""
    if not A:
        return []
    M = [row[:] for row in A]
    n = len(M)
    m = len(M[0])
    pivots = []
    row = 0
    col = 0
    while row < n and col < m:
        pr = None
        for rr in range(row, n):
            if M[rr][col] != 0:
                pr = rr
                break
        if pr is None:
            col += 1
            continue
        M[row], M[pr] = M[pr], M[row]
        pv = M[row][col]
        # Normalize pivot row
        for cc in range(col, m):
            M[row][cc] = M[row][cc] / pv
        for rr in range(n):
            if rr == row:
                continue
            f = M[rr][col]
            if f == 0:
                continue
            for cc in range(col, m):
                M[rr][cc] -= f * M[row][cc]
        pivots.append(col)
        row += 1
        col += 1
    free = [c for c in range(m) if c not in pivots]
    basis = []
    for fc in free:
        v = [Fraction(0)] * m
        v[fc] = Fraction(1)
        for pi, pc in enumerate(pivots):
            v[pc] = -M[pi][fc]
        basis.append(v)
    return basis


def intersect_kernels(A, B):
    """Dim of common right null space of A and B."""
    # Stack A on top of B, compute kernel
    if not A:
        return kernel_dim(B)
    if not B:
        return kernel_dim(A)
    M = [row[:] for row in A] + [row[:] for row in B]
    return kernel_dim(M)


# ----------------------------------------------------------------------
# Build representation matrix for any permutation
# ----------------------------------------------------------------------

def perm_to_adj_transpositions(perm):
    """Decompose a permutation (1-indexed list) into a product of adjacent
    transpositions s_i = (i, i+1).  Returns list of i's such that
    s_{i_1} s_{i_2} ... s_{i_k} = perm (read left to right).
    Uses bubble sort: sort the inverse permutation.  Equivalent to writing
    perm as a reduced word."""
    p = list(perm)
    n = len(p)
    # We want a sequence of s_i such that applying them in order to identity
    # gives p.  Simpler: decompose by selection — find where 1 is, bubble it
    # to position 1, etc.
    swaps = []
    cur = list(p)
    for target in range(1, n):
        # Find position of `target` in cur
        pos = cur.index(target)
        while pos > target - 1:
            # Swap positions pos-1 and pos
            cur[pos - 1], cur[pos] = cur[pos], cur[pos - 1]
            swaps.append(pos)  # 1-indexed: s_{pos} swaps positions pos and pos+1
            # wait: positions are 1-indexed from outside; pos here is 0-indexed
            # so s_{pos} in 1-indexed swaps positions pos and pos+1 in 1-indexed
            # which is indices pos-1 and pos in 0-indexed.  CORRECT.
            pos -= 1
    # `swaps` brings cur to identity; so identity = swaps applied to cur = perm.
    # Thus perm = inverse-of-swaps applied to identity = reverse(swaps) since
    # adjacent transpositions are involutions.
    return list(reversed(swaps))


def perm_matrix(lam, perm, si_mats):
    """Build matrix of perm in seminormal rep of lam."""
    syts = standard_young_tableaux(lam)
    d = len(syts)
    word = perm_to_adj_transpositions(perm)
    M = identity(d)
    for i in word:
        M = mat_mul(M, si_mats[i])
    return M


def w0_interval_perm(n, p, q):
    """Permutation in S_n (1-indexed list) that reverses positions p..q."""
    one_line = list(range(1, n + 1))
    one_line[p - 1:q] = list(range(q, p - 1, -1))
    return one_line


# ----------------------------------------------------------------------
# Verification
# ----------------------------------------------------------------------

def verify_construction(n=4):
    print("=" * 70)
    print(f"VERIFICATION CROSS-CHECKS for n={n}")
    print("=" * 70)

    # Trivial: lambda = (n,)
    lam_triv = (n,)
    si_mats, syts = seminormal_si_matrix(lam_triv)
    print(f"\nTrivial rep lambda={lam_triv}, dim={len(syts)}:")
    triv_ok = True
    for p in range(1, n + 1):
        for q in range(p + 1, n + 1):
            perm = w0_interval_perm(n, p, q)
            M = perm_matrix(lam_triv, perm, si_mats)
            ok = mat_eq(M, identity(1))
            triv_ok = triv_ok and ok
            print(f"  s_({p},{q}): M={M[0][0]}  [{'OK' if ok else 'FAIL'}]")

    # Sign: lambda = (1^n)
    lam_sgn = tuple([1] * n)
    si_mats_s, syts_s = seminormal_si_matrix(lam_sgn)
    print(f"\nSign rep lambda={lam_sgn}, dim={len(syts_s)}:")
    sgn_ok = True
    for p in range(1, n + 1):
        for q in range(p + 1, n + 1):
            perm = w0_interval_perm(n, p, q)
            M = perm_matrix(lam_sgn, perm, si_mats_s)
            m = q - p + 1
            expected = (-1) ** ((m * (m - 1)) // 2)
            actual = M[0][0]
            ok = (actual == Fraction(expected))
            sgn_ok = sgn_ok and ok
            print(f"  s_({p},{q}) m={m}: expected {expected}, got {actual}  [{'OK' if ok else 'FAIL'}]")

    # Also check all si_mats have square = identity (involution)
    inv_ok = True
    for i, M in si_mats.items():
        sq = mat_mul(M, M)
        if not mat_eq(sq, identity(len(M))):
            inv_ok = False
            print(f"  !! s_{i}^2 != Id in trivial rep")
    print(f"\nTrivial rep check: {'PASSED' if triv_ok else 'FAILED'}")
    print(f"Sign rep check:    {'PASSED' if sgn_ok else 'FAILED'}")
    return triv_ok and sgn_ok


def verify_involution_general(n):
    """Spot check: for all irreps, s_i^2 = Id."""
    print(f"\nInvolution check (s_i^2 = Id) for all irreps of S_{n}:")
    all_ok = True
    for lam in partitions(n):
        si_mats, syts = seminormal_si_matrix(lam)
        for i, M in si_mats.items():
            sq = mat_mul(M, M)
            if not mat_eq(sq, identity(len(M))):
                all_ok = False
                print(f"  FAIL lambda={lam}, i={i}")
    print(f"  Involution check: {'PASSED' if all_ok else 'FAILED'}")
    return all_ok


def verify_w0_involution(n):
    """w_0^{[p,q]} should square to identity in every rep."""
    all_ok = True
    for lam in partitions(n):
        si_mats, syts = seminormal_si_matrix(lam)
        d = len(syts)
        for p in range(1, n + 1):
            for q in range(p + 1, n + 1):
                perm = w0_interval_perm(n, p, q)
                M = perm_matrix(lam, perm, si_mats)
                sq = mat_mul(M, M)
                if not mat_eq(sq, identity(d)):
                    all_ok = False
                    print(f"  FAIL lambda={lam}, (p,q)=({p},{q})")
    return all_ok


def verify_dim_sum(n):
    """Sum of squares of dims = n!"""
    total = 0
    import math
    for lam in partitions(n):
        syts = standard_young_tableaux(lam)
        total += len(syts) ** 2
    expected = math.factorial(n)
    print(f"  Sum of dim^2 over irreps of S_{n} = {total}  (expect {expected})  [{'OK' if total == expected else 'FAIL'}]")
    return total == expected


# ----------------------------------------------------------------------
# Eigenspace dimensions
# ----------------------------------------------------------------------

def eigenspace_dims(M):
    n = len(M)
    Iden = identity(n)
    Mp = mat_sub(M, Iden)  # +1 eigenspace = ker(M - I)
    Mm = mat_add(M, Iden)  # -1 eigenspace = ker(M + I)
    return kernel_dim(Mp), kernel_dim(Mm)


def joint_fixed(M1, M2):
    """Dim of joint +1 eigenspace = ker(M1 - I) intersect ker(M2 - I)."""
    n = len(M1)
    Iden = identity(n)
    return intersect_kernels(mat_sub(M1, Iden), mat_sub(M2, Iden))


def commutator_rank(M1, M2):
    A = mat_mul(M1, M2)
    B = mat_mul(M2, M1)
    return mat_rank_full(mat_sub(A, B))


def classify_arrangement(p1, q1, p2, q2):
    if q1 < p2 or q2 < p1:
        return "Separated"
    if (p1 == p2 and q1 == q2):
        return "Equal"
    if (p1 <= p2 and q2 <= q1) or (p2 <= p1 and q1 <= q2):
        return "Nested"
    endpoints1 = {p1, q1}
    endpoints2 = {p2, q2}
    if len(endpoints1 & endpoints2) == 1:
        return "Adjacent"
    return "Overlapping"


# ----------------------------------------------------------------------
# Main pipeline
# ----------------------------------------------------------------------

def run_for_n(n, csv_writer, summary):
    print("\n" + "=" * 70)
    print(f"S_{n}  --  single generator and pair scan")
    print("=" * 70)

    gens = [(p, q) for p in range(1, n + 1) for q in range(p + 1, n + 1)]
    print(f"Cactus generators ({len(gens)}): {gens}")

    for lam in partitions(n):
        si_mats, syts = seminormal_si_matrix(lam)
        d = len(syts)
        print(f"\n  lambda = {lam}, dim = {d}")
        # Build rho(s_{p,q}) for every (p,q)
        rho = {}
        for (p, q) in gens:
            perm = w0_interval_perm(n, p, q)
            M = perm_matrix(lam, perm, si_mats)
            rho[(p, q)] = M
            a, b = eigenspace_dims(M)
            print(f"    s_({p},{q}): (+1 dim, -1 dim) = ({a}, {b})")
        # Pair scan
        for (g1, g2) in combinations(gens, 2):
            p1, q1 = g1
            p2, q2 = g2
            arr = classify_arrangement(p1, q1, p2, q2)
            if arr == "Equal":
                continue
            A = rho[g1]
            B = rho[g2]
            jf = joint_fixed(A, B)
            cr = commutator_rank(A, B)
            csv_writer.writerow({
                "n": n, "lambda": str(list(lam)),
                "p1": p1, "q1": q1, "p2": p2, "q2": q2,
                "arrangement": arr, "dim_lambda": d,
                "jointfix_dim": jf, "commutator_rank": cr,
            })
            summary.setdefault((n, str(list(lam)), arr), []).append({
                "p1q1": (p1, q1), "p2q2": (p2, q2),
                "dim": d, "joint": jf, "comm_rank": cr,
            })


def summary_report(summary):
    print("\n" + "=" * 70)
    print("SUMMARY BY ARRANGEMENT (per n, per lambda)")
    print("=" * 70)
    print(f"{'n':>2} {'lambda':<14} {'arr':<13} {'count':>5} {'dim':>4} "
          f"{'mean_joint':>10} {'mean_comm':>10} {'min_comm':>9} {'max_comm':>9}")
    keys_sorted = sorted(summary.keys())
    for key in keys_sorted:
        n, lam, arr = key
        rows = summary[key]
        d = rows[0]["dim"]
        joints = [r["joint"] for r in rows]
        cranks = [r["comm_rank"] for r in rows]
        mj = sum(joints) / len(joints)
        mc = sum(cranks) / len(cranks)
        print(f"{n:>2} {lam:<14} {arr:<13} {len(rows):>5} {d:>4} "
              f"{mj:>10.3f} {mc:>10.3f} {min(cranks):>9} {max(cranks):>9}")

    # Per-irrep arrangement comparison
    print("\n" + "=" * 70)
    print("ARRANGEMENT-DEPENDENT IRREPS (dim > 1)")
    print("=" * 70)
    by_nlam = {}
    for (n, lam, arr), rows in summary.items():
        by_nlam.setdefault((n, lam), {})[arr] = rows

    striking = []
    for (n, lam), arr_dict in sorted(by_nlam.items()):
        d = list(arr_dict.values())[0][0]["dim"]
        if d == 1:
            continue
        print(f"\nn={n}, lambda={lam}, dim={d}")
        means = {}
        for arr in ("Separated", "Adjacent", "Nested", "Overlapping"):
            if arr not in arr_dict:
                continue
            rows = arr_dict[arr]
            cranks = [r["comm_rank"] for r in rows]
            joints = [r["joint"] for r in rows]
            mc = sum(cranks) / len(cranks)
            mj = sum(joints) / len(joints)
            means[arr] = (mc, mj)
            print(f"  {arr:<13} count={len(rows):>3}  "
                  f"mean(comm_rank)={mc:>6.3f}  "
                  f"mean(joint_fix)={mj:>6.3f}  "
                  f"comm_rank range=[{min(cranks)},{max(cranks)}]")
        # Striking: check if max-min of mean-commutator-rank across arrangements is large
        if means:
            mc_vals = [v[0] for v in means.values()]
            spread = max(mc_vals) - min(mc_vals)
            if spread > 0.5:
                striking.append((n, lam, d, spread, means))

    print("\n" + "=" * 70)
    print("STRIKING ORDERINGS (mean commutator rank varies across arrangement types)")
    print("=" * 70)
    for n, lam, d, spread, means in sorted(striking, key=lambda x: -x[3]):
        order = sorted(means.items(), key=lambda kv: kv[1][0])
        ordering = " < ".join(f"{a}({mc:.2f})" for a, (mc, mj) in order)
        print(f"  n={n} lambda={lam} dim={d}  spread={spread:.3f}  {ordering}")


def main():
    print("Pure-Python Specht module construction (Young's seminormal form)")
    print("Container has no SageMath; using fractions + custom linear algebra.")

    # Verifications
    ok4 = verify_construction(4)
    print()
    print("Dim sum check:")
    verify_dim_sum(4)
    verify_dim_sum(5)
    print()
    print("Involution checks:")
    verify_involution_general(4)
    print("w_0^{[p,q]} squared = Id check for n=4:",
          "PASSED" if verify_w0_involution(4) else "FAILED")
    print("w_0^{[p,q]} squared = Id check for n=5:",
          "PASSED" if verify_w0_involution(5) else "FAILED")

    csv_path = "/home/clio/projects/signed-laxator-work/cactus_rank_data.csv"
    fieldnames = ["n", "lambda", "p1", "q1", "p2", "q2",
                  "arrangement", "dim_lambda", "jointfix_dim", "commutator_rank"]
    summary = {}
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for n in (4, 5):
            run_for_n(n, writer, summary)

    summary_report(summary)
    print(f"\nCSV written to: {csv_path}")


if __name__ == "__main__":
    main()
