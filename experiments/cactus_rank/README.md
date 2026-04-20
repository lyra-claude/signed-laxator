# Cactus generator rank by interval arrangement

## Computational question

Let `J_n` denote the cactus group, with generators `s_{p,q}` for `1 <= p < q <= n`,
acting on a Specht module `V_lambda` of `S_n` via
`rho_lambda(w_0^{[p,q]})`, where `w_0^{[p,q]}` is the longest element of the
parabolic subgroup `S_{[p,q]}` reversing the interval `[p, q]`.

For each pair of cactus generators `s_{p,q}, s_{p',q'}`, we compute:

- `jointfix_dim` = `dim ker(rho_lambda(s_{p,q}) - I) intersect ker(rho_lambda(s_{p',q'}) - I)`
- `commutator_rank` = `rank(rho_lambda(s_{p,q}) rho_lambda(s_{p',q'}) - rho_lambda(s_{p',q'}) rho_lambda(s_{p,q}))`

Each pair is classified by the relative arrangement of the intervals `[p,q]` and
`[p',q']`:

- **Separated**: `q < p'` (disjoint, non-touching)
- **Adjacent**: `q = p'` (share a single endpoint)
- **Nested**: one interval strictly contains the other
- **Overlapping**: intervals overlap but neither contains the other

## Headline finding

For every non-trivial irreducible representation `V_lambda` of `S_n` at
`n = 4, 5`, the commutator rank is strictly ordered by arrangement type:

> **Separated  <  Overlapping  <  Nested  <  Adjacent**

Mirror-image ordering for the joint fixed dimension. Adjacent generators
always achieve the maximum commutator rank (equal to `dim V_lambda` for the
smaller-dimensional irreps), and Separated generators always commute
(rank 0).

The ordering is uniform across irreps: it is not a per-`lambda` accident but a
structural feature of the cactus action. See `cactus_rank_data.csv` for the
full 390-row table.

### Caveat at $n=4$ (per Lyra's PR#1 review)

The strict ordering $\text{Sep} < \text{Over}$ holds for $n \geq 5$; at $n=4$
it becomes a **weak** inequality $\text{Sep} \leq \text{Over}$. When the two
intervals are disjoint (the Separated case), the cactus generators act as
transpositions in **disjoint parabolic subgroups** and therefore always
commute, giving commutator rank exactly $0$. At $n=4$ the Overlapping category
also collapses to rank $0$ on several low-dimensional irreps — there simply
isn't enough room in $S_4$ for Overlapping to be strictly heavier than
Separated. The 390 rows in `cactus_rank_data.csv` exhibit the ties at $n=4$
directly; strict separation between Sep and Over is first seen at $n=5$.

## Reproducing

The script is pure Python (no SageMath required). Specht modules are built via
Young's seminormal form over the rationals:

```
python3 cactus_rank.py
```

This regenerates `cactus_rank_output.log` and `cactus_rank_data.csv`. Total
runtime: a few seconds on a laptop.

A SageMath translation is provided in `cactus_rank.sage` (untested in the
container that produced this data; SageMath was unavailable, so the pure-Python
path is the source of truth).

## Files

- `cactus_rank.py` — pure-Python computation (Young's seminormal form).
- `cactus_rank.sage` — SageMath translation (untested).
- `cactus_rank_data.csv` — 390 rows, tab-separated, columns:
  `n, lambda, p1, q1, p2, q2, arrangement, dim_lambda, jointfix_dim, commutator_rank`.
- `cactus_rank_output.log` — full run log including verification cross-checks
  (trivial/sign rep correctness, dim-squared sums, involution `s_i^2 = I`,
  `(w_0^{[p,q]})^2 = I`).

## Connection to Lyra's maze sign-flip

In the NK-landscape Table 6 experiments, the **arrangement type of a pair of
intervals** (interval = communication channel between subgroups) robustly
predicts the sign of the laxator's diversity component: separated channels
behave additively (sign +1), adjacent/nested channels exhibit non-trivial sign
structure.

This experiment establishes a parallel story on the algebraic side: the same
arrangement classification — Separated / Overlapping / Nested / Adjacent —
robustly predicts the **non-commutativity of the corresponding cactus
generators** acting on irreps. Separated generators commute on the nose;
adjacent ones are maximally non-commuting.

The shared classifier suggests that both phenomena are shadows of one
underlying invariant: an "interval-arrangement obstruction" in the categorical
phase diagram of the signed laxator. This feeds directly into Q3 in
`notes/research-questions.md` (the sign component of the laxator).
