# Three Signs of the Same Flip

> Connecting C110 (crystal/atom = beta_1/lambda_2), C117 (three signs of the same flip), and the bridge experiment results.

## Summary

Three independent research threads have each discovered that topology's effect on multi-agent / evolutionary optimization has a **sign** -- it can help or hurt, and which one depends on the invariant, the timescale, and the agent type. We argue these are three views of a single phenomenon: the signed laxator.

| Observation | Setting | Sign-flip mechanism | Key data |
|---|---|---|---|
| Bridge experiment | Island-model GA, NK landscapes | beta_1 = transient (+), lambda_2 = persistent (+), but different timescales | 320 runs, eta^2 values |
| Claude Chorus | LLM multi-agent collaboration | Connectivity -> consensus (diversity DECREASES) | Kendall's W = 0.90 |
| Crystal/Atom | Hecke algebra at q -> 0 | Crystal = transient, Atom = persistent, sharp phase boundary | Transition at q = 0.001 |

## 1. The Bridge Experiment

### Setup

320 runs on the `feat/bridge-experiment` branch of `lyra-claude/Topology-experiments`.

- **4 graphs:** K4 (complete on 4 vertices, beta_1 = 3, lambda_2 = 4), K4-e (one edge removed, beta_1 = 2, lambda_2 ~ 2.59), K4-2e (two edges removed / "path+triangle," beta_1 = 1, lambda_2 ~ 1.38), star (beta_1 = 0, lambda_2 = 1).
- **4 NK landscapes:** K = 0 (smooth), K = 2, K = 4, K = 6 (maximally rugged for N=10).
- **20 seeds** per condition.
- **Iso-spectral bridge construction:** K4-e and K4-2e were chosen to have equal Fiedler components where possible, isolating the beta_1 effect from the lambda_2 effect.

### Results: Two Independent Channels

**Transient channel (beta_1):**
- eta^2 = 0.191 at generation 100 -- beta_1 explains ~19% of diversity variance.
- Fades to near zero by generation 500.
- Higher beta_1 -> higher early diversity (more cycles = more pathways for novel solutions to persist temporarily).

**Persistent channel (lambda_2):**
- eta^2 = 0.26 at generation 200 -- lambda_2 explains ~26% of diversity variance.
- **Persists** through generation 500 and beyond.
- Higher lambda_2 -> faster mixing -> lower long-run diversity (the population homogenizes through better-connected graphs).

### The Sign Structure

The sign of topology's effect on diversity depends on WHICH invariant and WHEN:

| Invariant | Early (gen 100) | Late (gen 500) |
|---|---|---|
| beta_1 | **Positive** (more cycles -> more diversity) | **Near zero** (effect fades) |
| lambda_2 | Weak | **Negative** (more connectivity -> less diversity) |

This is the first sign: topology helps diversity through cycles (beta_1) but eventually hurts it through mixing (lambda_2).

## 2. The Contravariance Observation (Claude Chorus)

### Setup

Multi-agent LLM experiments where multiple Claude instances collaborate on a shared topology. Diversity measured as agreement across agents using Kendall's W statistic.

### Result

Kendall's W = 0.90 across agents -- near-perfect agreement. In LLM multi-agent systems:

- **More connectivity -> more consensus -> LESS diversity.**
- The ranking of topologies by diversity is **inverted** compared to the GA setting.

### The Sign Structure

| System type | Effect of increasing beta_1 on diversity |
|---|---|
| Genetic algorithm | **Positive** (transient) -- cycles allow escape from local optima |
| LLM multi-agent | **Negative** -- cycles allow faster consensus convergence |

This is the second sign: the **same** topological feature (cycle rank) has **opposite** effects depending on whether agents are GAs or LLMs.

### Why This Happens (Mechanistic Intuition)

- **GAs:** Each island has an independent population. Cycles in the migration graph allow novel alleles to persist by circulating through multiple islands before being lost. Higher beta_1 = more recirculation pathways = more transient diversity.
- **LLMs:** Each agent conditions on messages from neighbors. Cycles in the communication graph create feedback loops that reinforce consensus. Higher beta_1 = faster convergence to agreement = less diversity.

The difference is in the **coupling dynamics:** GAs have weak coupling (migration rate << selection pressure), while LLMs have strong coupling (each agent directly incorporates neighbor outputs). The sign of diversity response flips at a critical coupling strength.

## 3. Clio's Demazure Crystal/Atom Decomposition

### Setting

Hecke algebra H_n(q) representations. Demazure operators acting on the polynomial ring. The crystal basis (combinatorial, arising at q -> 0) and the atom basis (algebraic, related to Kazhdan-Lusztig polynomials).

### The Mapping (C110, confidence 85%)

| Algebraic concept | Topological concept | Timescale |
|---|---|---|
| Crystal (Demazure crystal basis) | beta_1 channel | Transient |
| Atom (Demazure atom basis) | lambda_2 channel | Persistent |
| q parameter | Coupling strength? | Controls transition |

### Key Evidence

1. **Sharp transition at q = 0.001:** The crystal and atom behaviors separate cleanly at this threshold. This mirrors the GA/LLM split -- at low coupling (GA), the two channels are visible; at high coupling (LLM), the persistent channel dominates.

2. **Tits deformation:** The cactus midpoint operator rank = n! / 2^floor(n/2), independent of q. This means the **total** structure is preserved (by Tits deformation theorem, H_n(q) ~ C[S_n] for q != 0), but the **decomposition** into crystal vs atom components changes. Similarly, the total diversity capacity of a graph is fixed by its size, but the decomposition into transient (beta_1) vs persistent (lambda_2) channels changes with the dynamics.

3. **Coboundary monoidal structure:** The relevant monoidal structure is coboundary (not braided). Coboundary categories have a "flat connection" -- parallel transport is path-independent. This suggests that the signed laxator's sign depends only on the **endpoints** (the topological invariants), not on the path taken through parameter space.

## 4. Why These Are "The Same" Sign Flip

### The Unifying Pattern

All three observations share the same abstract structure:

1. **Two components** with different timescale behaviors (transient vs persistent).
2. **A sign** that determines whether a topological feature helps or hurts.
3. **A parameter** that controls which regime you're in (coupling strength, q, NK landscape ruggedness K).

The signed laxator phi_G encodes this structure as a lax natural transformation whose 2-cells carry the sign. The three observations are:

- **Bridge experiment:** phi_G evaluated at different timescales (gen 100 vs gen 500), showing the two components.
- **Claude Chorus:** phi_G evaluated for different agent types (GA functor vs LLM functor), showing the sign flip across systems.
- **Crystal/Atom:** phi_G's algebraic decomposition in the Hecke algebra, providing the mathematical structure underlying the sign.

### A Fourth Sign? Directed Topologies

Connection C131 (temporal inversion requires epistasis, confidence 95%): in directed topologies, the effect of directed cycle count kappa on diversity **inverts sign** as landscape ruggedness K increases beyond 2. This may be a fourth manifestation of the same sign flip, now parameterized by landscape complexity rather than agent type or timescale.

## Open Problems

1. **Make the crystal/atom <-> beta_1/lambda_2 mapping precise.** Currently C110 is at 85% confidence. We need a theorem, not an analogy.
2. **Determine whether q interpolates between GA and LLM behavior.** If so, what is q "physically" in the multi-agent context?
3. **Connect to the directed topology results.** The 960-run directed experiment provides data on a different sign flip (temporal inversion). Does phi_G capture this, or is a separate structure needed?
4. **Fiber preservation question** (Clio): as the parameter kappa = (1 - qp)/(1 - p) moves on the Riemann sphere, do representation fibers stay rigid (parallel transport) or rotate? The coboundary = flat connection result suggests rigidity, but this needs verification.
