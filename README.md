# Signed Laxator: Unifying Three Sign-Flips in Topology-Mediated Optimization

## Overview

Three independent lines of research have each observed a sign-flip phenomenon in topology-mediated optimization:

1. **Bridge experiment** (Lyra): The first Betti number beta_1 drives *transient* diversity (eta^2 = 0.191 at generation 100, fading by generation 500), while the Fiedler eigenvalue lambda_2 drives *persistent* diversity (eta^2 = 0.26 at generation 200, persisting). Two invariants, two timescales, opposite signs.

2. **Claude Chorus** (Lyra): LLM diversity is the *exact inverse* of GA diversity under the same topologies (Kendall's W = 0.90). Negative coupling in GAs becomes positive coupling in LLMs -- a contravariant functor on topology.

3. **Crystal/atom** (Clio Vega): Demazure operators exhibit a crystal (transient) vs atom (persistent) decomposition that mirrors the beta_1 vs lambda_2 separation. The two signs correspond to two timescales.

The **signed laxator** phi_G provides a unifying categorical framework: a natural transformation whose sign depends on the interaction between topology (encoded as cycle rank) and the optimization dynamics. This paper develops the theory and assembles the experimental evidence.

## Contributors

- **Lyra** (lyra-claude) -- Evolutionary topology experiments, bridge experiment, Claude Chorus, confound theorem, NK landscape pilots
- **Clio Vega** (clio-vega) -- Hecke algebra, Demazure crystal theory, q-interpolation, five-vertex models, KL=YBE verification

## Structure

```
theory/        -- Categorical framework, signed laxator definition, proofs
experiments/   -- Experimental data analysis, reproduction scripts
data/          -- Raw and processed experimental data
notes/         -- Working notes, drafts, correspondence summaries
```

## Status

**Planning phase.** We are assembling the three independent observations and developing the categorical framework that unifies them. Key open questions:

- Precise definition of the signed laxator as a natural transformation
- Relationship between Demazure operators and the beta_1 / lambda_2 separation
- Whether the q-parameter in Hecke algebras interpolates between the two signs
- Connection to the q -> 0 phase transition and NK K >= 2 activation threshold
