#!/usr/bin/env sage
"""
cactus_rank.sage

Compute cactus-group-generator ranks on irreducible Specht modules of S_n
for n = 4, 5.

For each irrep V_lambda of S_n and each cactus generator s_{p,q} (mapped to
w_0^{[p,q]} in S_n via the surjection J_n -> S_n), we record:
  - a_lambda(p,q) = dim of (+1)-eigenspace
  - b_lambda(p,q) = dim of (-1)-eigenspace

For each pair of generators we classify the relative arrangement and
record commutator rank, joint +1 eigenspace, joint fixed dim of the group
they generate.
"""

from sage.all import (
    Partitions, SymmetricGroup, SymmetricGroupRepresentation,
    Matrix, MatrixSpace, identity_matrix, QQ, ZZ, Permutation, vector
)
import csv
import itertools


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def w0_interval(n, p, q):
    """Permutation in S_n that reverses positions [p..q] and fixes the rest.
    p, q are 1-indexed.  Returns a SageMath Permutation (one-line)."""
    one_line = list(range(1, n + 1))
    rev = list(range(q, p - 1, -1))
    one_line[p - 1:q] = rev
    return Permutation(one_line)


def rho_matrix(specht, perm):
    """Apply specht representation to a Permutation, return Sage matrix over QQ."""
    M = specht(perm)
    return Matrix(QQ, M)


def eigenspace_dims(M):
    """Return (a, b) = (dim of +1 eigenspace, dim of -1 eigenspace)."""
    n = M.nrows()
    Iden = identity_matrix(QQ, n)
    a = (M - Iden).right_kernel().dimension()
    b = (M + Iden).right_kernel().dimension()
    return a, b


def joint_fixed_dim(matrices):
    """Dim of the joint +1 eigenspace (intersection of fixed spaces)."""
    if not matrices:
        return 0
    n = matrices[0].nrows()
    Iden = identity_matrix(QQ, n)
    V = (matrices[0] - Iden).right_kernel()
    for M in matrices[1:]:
        V = V.intersection((M - Iden).right_kernel())
    return V.dimension()


def joint_group_fixed_dim(M1, M2):
    """Dim of fixed space of group <M1, M2>.
    Generate the (finite) subgroup by closing under products until stable."""
    n = M1.nrows()
    Iden = identity_matrix(QQ, n)
    # Start from identity, M1, M2
    seen = {tuple(Iden.list()): Iden, tuple(M1.list()): M1, tuple(M2.list()): M2}
    queue = [Iden, M1, M2]
    while queue:
        g = queue.pop()
        for h in (M1, M2):
            r = g * h
            key = tuple(r.list())
            if key not in seen:
                seen[key] = r
                queue.append(r)
            r2 = h * g
            key2 = tuple(r2.list())
            if key2 not in seen:
                seen[key2] = r2
                queue.append(r2)
        # safety: bound the size — for our cases the group should be small
        if len(seen) > 200:
            break
    # Joint fixed dim = intersection of fixed spaces of all generators
    # (equivalently of all group elements, but generators suffice)
    V = (M1 - Iden).right_kernel().intersection((M2 - Iden).right_kernel())
    return V.dimension(), len(seen)


def classify_arrangement(p1, q1, p2, q2):
    """Classify the arrangement of two intervals."""
    # Separated: disjoint
    if q1 < p2 or q2 < p1:
        return "Separated"
    # Nested: one strictly contains the other
    if (p1 <= p2 and q2 <= q1 and (p1 < p2 or q2 < q1)):
        return "Nested"
    if (p2 <= p1 and q1 <= q2 and (p2 < p1 or q1 < q2)):
        return "Nested"
    # Equal intervals — treat as Nested degenerate (skip in pair scan since we need distinct)
    if p1 == p2 and q1 == q2:
        return "Equal"
    # Adjacent: share exactly one endpoint
    endpoints1 = {p1, q1}
    endpoints2 = {p2, q2}
    shared = endpoints1 & endpoints2
    if len(shared) == 1:
        # And not nested — must be adjacent
        return "Adjacent"
    # Otherwise partial overlap
    return "Overlapping"


# ----------------------------------------------------------------------
# Verification cross-checks (n=4)
# ----------------------------------------------------------------------

def verify_trivial_and_sign(n=4):
    print("=" * 70)
    print(f"VERIFICATION CROSS-CHECKS for n={n}")
    print("=" * 70)

    # Trivial rep: lambda = (n)
    lam_triv = [n]
    rep_triv = SymmetricGroupRepresentation(lam_triv, "specht")
    print(f"\nTrivial rep lambda = {lam_triv}:")
    triv_ok = True
    for p in range(1, n + 1):
        for q in range(p + 1, n + 1):
            M = rho_matrix(rep_triv, w0_interval(n, p, q))
            a, b = eigenspace_dims(M)
            if M != identity_matrix(QQ, M.nrows()) or b != 0 or a != 1:
                triv_ok = False
                print(f"  s_({p},{q}): FAIL  M = {M.list()}")
            else:
                print(f"  s_({p},{q}): OK (identity, +1 dim = {a})")
    print(f"  Trivial rep check: {'PASSED' if triv_ok else 'FAILED'}")

    # Sign rep: lambda = (1^n)
    lam_sgn = [1] * n
    rep_sgn = SymmetricGroupRepresentation(lam_sgn, "specht")
    print(f"\nSign rep lambda = {lam_sgn}:")
    sgn_ok = True
    for p in range(1, n + 1):
        for q in range(p + 1, n + 1):
            M = rho_matrix(rep_sgn, w0_interval(n, p, q))
            length = q - p
            # sign of w_0 on m letters is (-1)^{m(m-1)/2}; here m = q-p+1
            m = q - p + 1
            expected = (-1) ** (m * (m - 1) // 2)
            actual = M[0, 0]
            ok = (M.nrows() == 1 and actual == expected)
            sgn_ok = sgn_ok and ok
            print(f"  s_({p},{q}) m={m}: expected {expected}, got {actual}  [{'OK' if ok else 'FAIL'}]")
    print(f"  Sign rep check: {'PASSED' if sgn_ok else 'FAILED'}")

    return triv_ok and sgn_ok


# ----------------------------------------------------------------------
# Single-generator rank table
# ----------------------------------------------------------------------

def single_generator_table(n):
    print("\n" + "=" * 70)
    print(f"SINGLE-GENERATOR EIGENSPACE DIMENSIONS for n={n}")
    print("=" * 70)
    print(f"{'lambda':<14}{'dim':>5}  generators (p,q): (a_+, b_-)")
    rows = []
    for lam in Partitions(n):
        rep = SymmetricGroupRepresentation(list(lam), "specht")
        d = rep([Permutation(list(range(1, n + 1)))]).nrows() if False else None
        # easier: build identity to get dim
        Iden_perm = Permutation(list(range(1, n + 1)))
        M0 = rho_matrix(rep, Iden_perm)
        d = M0.nrows()
        line = f"{str(list(lam)):<14}{d:>5}  "
        ab_list = []
        for p in range(1, n + 1):
            for q in range(p + 1, n + 1):
                M = rho_matrix(rep, w0_interval(n, p, q))
                a, b = eigenspace_dims(M)
                ab_list.append(((p, q), a, b))
                line += f"({p},{q}):({a},{b}) "
        print(line)
        rows.append((lam, d, ab_list))
    return rows


# ----------------------------------------------------------------------
# Pair scan
# ----------------------------------------------------------------------

def pair_scan(n, csv_writer, summary):
    print("\n" + "=" * 70)
    print(f"PAIR-OF-GENERATORS SCAN for n={n}")
    print("=" * 70)

    # All cactus generators
    gens = [(p, q) for p in range(1, n + 1) for q in range(p + 1, n + 1)]

    for lam in Partitions(n):
        lam_list = list(lam)
        rep = SymmetricGroupRepresentation(lam_list, "specht")
        Iden_perm = Permutation(list(range(1, n + 1)))
        d = rho_matrix(rep, Iden_perm).nrows()

        # Cache rho(s_{p,q})
        rho_gen = {}
        for (p, q) in gens:
            rho_gen[(p, q)] = rho_matrix(rep, w0_interval(n, p, q))

        for (p1, q1), (p2, q2) in itertools.combinations(gens, 2):
            arr = classify_arrangement(p1, q1, p2, q2)
            if arr == "Equal":
                continue  # only happens if same generator — skipped by combinations
            A = rho_gen[(p1, q1)]
            B = rho_gen[(p2, q2)]
            comm = A * B - B * A
            comm_rank = comm.rank()
            jp = joint_fixed_dim([A, B])
            grp_fix, grp_size = joint_group_fixed_dim(A, B)

            csv_writer.writerow({
                "n": n,
                "lambda": str(lam_list),
                "p1": p1, "q1": q1,
                "p2": p2, "q2": q2,
                "arrangement": arr,
                "dim_lambda": d,
                "jointfix_dim": grp_fix,
                "commutator_rank": comm_rank,
            })

            key = (n, str(lam_list), arr)
            summary.setdefault(key, []).append({
                "p1q1": (p1, q1), "p2q2": (p2, q2),
                "dim": d,
                "joint": grp_fix,
                "comm_rank": comm_rank,
                "group_size": grp_size,
            })


# ----------------------------------------------------------------------
# Summary report
# ----------------------------------------------------------------------

def summary_report(summary):
    print("\n" + "=" * 70)
    print("SUMMARY BY ARRANGEMENT")
    print("=" * 70)

    # Group by (n, lambda, arrangement)
    keys_sorted = sorted(summary.keys())
    last_n = None
    last_lam = None
    print(f"{'n':>2} {'lambda':<14} {'arr':<13} {'count':>5} {'dim':>4} "
          f"{'mean_joint':>10} {'mean_comm_rk':>13} {'min_comm_rk':>11} {'max_comm_rk':>11}")
    for key in keys_sorted:
        n, lam, arr = key
        rows = summary[key]
        d = rows[0]["dim"]
        joints = [r["joint"] for r in rows]
        cranks = [r["comm_rank"] for r in rows]
        mean_joint = sum(joints) / len(joints)
        mean_crank = sum(cranks) / len(cranks)
        print(f"{n:>2} {lam:<14} {arr:<13} {len(rows):>5} {d:>4} "
              f"{float(mean_joint):>10.3f} {float(mean_crank):>13.3f} "
              f"{min(cranks):>11} {max(cranks):>11}")

    # Striking-difference detection: per (n, lambda), check if mean joint or
    # mean commutator-rank changes across arrangement types.
    print("\n" + "=" * 70)
    print("STRIKING ARRANGEMENT-DEPENDENT IRREPS")
    print("=" * 70)
    by_nlam = {}
    for (n, lam, arr), rows in summary.items():
        by_nlam.setdefault((n, lam), {})[arr] = rows

    for (n, lam), arr_dict in sorted(by_nlam.items()):
        d = list(arr_dict.values())[0][0]["dim"]
        if d == 1:
            continue  # trivial / sign — boring
        print(f"\nn={n}, lambda={lam}, dim={d}")
        for arr in ("Separated", "Adjacent", "Nested", "Overlapping"):
            if arr not in arr_dict:
                continue
            rows = arr_dict[arr]
            cranks = [r["comm_rank"] for r in rows]
            joints = [r["joint"] for r in rows]
            mean_crank = sum(cranks) / len(cranks)
            mean_joint = sum(joints) / len(joints)
            print(f"  {arr:<13} count={len(rows):>3}  "
                  f"mean(comm_rank)={float(mean_crank):>6.3f}  "
                  f"mean(joint_fix)={float(mean_joint):>6.3f}  "
                  f"comm_rank range=[{min(cranks)},{max(cranks)}]")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    ok4 = verify_trivial_and_sign(4)
    print("\nOverall verification (n=4):", "PASSED" if ok4 else "FAILED")

    csv_path = "/home/clio/projects/signed-laxator-work/cactus_rank_data.csv"
    fieldnames = ["n", "lambda", "p1", "q1", "p2", "q2",
                  "arrangement", "dim_lambda", "jointfix_dim", "commutator_rank"]
    summary = {}
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for n in (4, 5):
            single_generator_table(n)
            pair_scan(n, writer, summary)

    summary_report(summary)
    print(f"\nCSV written to: {csv_path}")


main()
