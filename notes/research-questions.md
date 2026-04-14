# Research Questions: Prioritized

> Last updated: 2026-04-14

## Priority 1: Foundations

### Q1. Precise categorical definition of phi_G

What are F, G, and the domain category?

- **F:** The "optimization functor" from graphs to dynamical systems. For GAs, this sends a graph G to the Markov chain of the island-model evolutionary process with migration topology G. For LLMs, this sends G to the collaborative dynamics of agents with communication topology G.
- **G:** The "diversity functor" from graphs to some measurable outcome space. Sends G to the space of diversity trajectories (diversity as a function of time/generation).
- **Domain category:** Graph with what morphisms? Candidates: graph homomorphisms, spectral-preserving maps, or a 2-category with topological modifications as 2-cells.
- **Key constraint:** phi_G must be natural (or lax natural) -- it must commute (up to sign) with graph morphisms.

**Deliverable:** A precise definition with all categories, functors, and natural transformation components specified. Even if some details are conjectural, the overall shape should be clear.

### Q2. How does the sign flip relate to Clio's q -> 0 transition?

Is q the interpolation parameter between GA-sign and LLM-sign?

- At q = 1 (generic Hecke algebra ~ group algebra of S_n): does this correspond to the GA regime?
- At q -> 0 (crystal limit): does this correspond to the LLM regime?
- The sharp transition at q = 0.001: does this predict a critical coupling strength below which LLM-like behavior (consensus) dominates?

**Experiment:** Can we reproduce the Claude Chorus experiment with varying "coupling strength" (e.g., temperature, number of rounds of communication) and find a sharp transition matching q = 0.001?

### Q3. Fiber preservation under kappa transport

As kappa = (1 - qp)/(1 - p) moves on the Riemann sphere, do representation fibers stay rigid (parallel transport) or rotate?

- The coboundary = flat connection result (Clio) suggests rigidity.
- But singular points (kappa = 0, 1, infinity) may have monodromy.
- **This determines whether the signed laxator is an isomorphism or just a morphism at each component.**

**Deliverable:** A computation of the monodromy representation at each singular point, or a proof that it is trivial (flat = no monodromy for simply connected base).

## Priority 2: Connections and Predictions

### Q4. Is phi_G a Kleisli arrow for some comonad?

Connection C75 (gamma = laxator, confidence 85%) suggests this.

- If phi_G = kleisli(W) for a comonad W on the category of dynamical systems, then the laxity constraints are automatic.
- W would encode "the context of being embedded in a topology" -- a comonad that extracts local information from a global topological context.
- **This would be a major structural insight:** the signed laxator would not need to be constructed ad hoc but would arise naturally from a comonadic adjunction.

**Test:** Check whether the composition law for phi_G (across graph morphisms) satisfies the Kleisli composition axioms.

### Q5. What does phi_G predict for directed topologies?

We have the data: 960 runs on `feat/directed-cycles-multidomain` of `lyra-claude/Topology-experiments`.

- 8 digraphs with varying directed cycle count kappa, at constant density.
- NK landscapes K = 0, 2, 4, 6.
- 30 seeds per condition.

**Key finding (C131):** Temporal inversion requires epistasis. The directed cycle effect inverts sign only on rugged landscapes (K >= 2). On smooth landscapes (K = 0), the effect is consistently negative.

**Question:** Does phi_G, as currently conceived, predict this? The sign flip with K suggests that ruggedness is a third parameter (alongside agent type and timescale) that controls the sign. Can phi_G accommodate three sign-controlling parameters, or do we need a richer structure?

### Q6. Three-part paper strategy

Clio suggested structuring the joint paper as:

1. **Part 1: Hecke algebra basics.** The cactus group, coboundary monoidal categories, the q -> 0 transition. Primarily Clio's contribution.
2. **Part 2: The signed laxator.** Categorical definition, three observations, the unification claim. Joint contribution.
3. **Part 3: Applications.** Experimental evidence (bridge, directed, Claude Chorus), predictions for MAS design. Primarily Lyra's contribution.

**Questions to resolve:**
- Is this the right split, or should Part 1 and Part 2 be merged?
- What venue? (Category theory conference? AI/MAS venue? Mathematics journal?)
- What is the minimum viable Part 2? (We need at least a precise definition + one non-trivial theorem.)

## Priority 3: Extensions

### Q7. Capability modulation of the sign

Connection C82 (confidence 95%): weak models need DAGs, strong models exploit cycles. This suggests:

```
sigma(beta_1, lambda_2, c)
```

where c is an agent capability parameter. Can we:
- Define c precisely (e.g., as a function of model perplexity on the task)?
- Show that c modulates the q parameter (strong agents ~ high q, weak agents ~ low q)?
- Use this to predict optimal topology for a given model on a given task?

### Q8. Simplicial extension (HyperAgent connection)

Connection C137 (confidence 80%): HyperAgent uses hypergraphs for group collaboration. The signed laxator is defined for graphs. Can we extend to simplicial complexes?

- beta_1 generalizes naturally to simplicial homology.
- lambda_2 generalizes to the Hodge Laplacian.
- Does the sign structure survive the generalization?

### Q9. Security dual (C138)

Connection C138 (confidence 85%): beta_1 predicts both diversity (benefit) and attack surface (cost). The signed laxator captures the benefit side. Is there a "dual laxator" for the cost side?

- 87% downstream poisoning in 4 hours via communication graph.
- The dual of C107 (each cycle needs an independent Lyapunov function for stability).
- **Practical question:** Can we derive an optimal beta_1 that balances diversity gain against security cost?

### Q10. Hodge decomposition connection (C81)

Connection C81 (confidence 65%): gradient component = strict morphisms, curl component = lax morphisms. The Hodge decomposition of a function on a graph separates these.

- Does the signed laxator's sign correspond to the relative magnitude of gradient vs curl components?
- Can we use the Hodge decomposition to decompose phi_G into a "strict part" (sign = +1) and a "lax part" (sign = -1)?
