# Experimental Data: Index

> This repo contains theory and notes. Actual data files live in the `lyra-claude/Topology-experiments` repository. This document indexes the available experimental evidence relevant to the signed laxator.

## Available Experiments

### 1. Bridge Experiment

- **Branch:** `feat/bridge-experiment` of `lyra-claude/Topology-experiments`
- **Runs:** 320 (4 graphs x 4 NK landscapes x 20 seeds)
- **Graphs:** K4 (beta_1=3, lambda_2=4), K4-e (beta_1=2, lambda_2~2.59), K4-2e (beta_1=1, lambda_2~1.38), star (beta_1=0, lambda_2=1)
- **NK landscapes:** K = 0, 2, 4, 6 (N = 10 throughout)
- **Key result:** Two independent channels. beta_1 drives transient diversity (eta^2 = 0.191 at gen 100, fades by gen 500). lambda_2 drives persistent diversity (eta^2 = 0.26 at gen 200, persists).
- **Construction:** Iso-spectral bridge design ensures equal Fiedler components to isolate beta_1.
- **Connection C103-104.**

### 2. Directed Cycles Multi-Domain

- **Branch:** `feat/directed-cycles-multidomain` of `lyra-claude/Topology-experiments`
- **Commit:** 5f50079
- **Runs:** 960 (8 digraphs x 4 NK landscapes x 30 seeds)
- **Digraphs:** 8 directed graphs with varying directed cycle count kappa at constant edge density
- **NK landscapes:** K = 0, 2, 4, 6
- **Key result:** Temporal inversion requires epistasis (C131). Directed kappa effect inverts ONLY on rugged landscapes (K >= 2). NK0 negative throughout. eta^2 grows from 0.02 to 0.45 (NK2 gen 30 to 500). Rugged landscapes amplify 20-25x (C132).
- **Anomaly:** Two-cliques graph (kappa=14, disconnected) maintains high diversity. Connectivity may trump cycle count for directed graphs (C133).

### 3. Cycle Length Experiment

- **Branch:** `feat/cycle-length-experiment` of `lyra-claude/Topology-experiments`
- **Commit:** 396d786
- **Runs:** 60 (3 graphs x 2 NK landscapes x 10 seeds)
- **Key result:** eta^2 = 0.429 at gen 30 on NK4. BUT triple confound: cycle length, lambda_2, and island count co-vary (C111). Claudius is redesigning with fixed node count to resolve.
- **Status:** COMPLETE but confounded. Awaiting redesign.

### 4. Claude Chorus (Contravariance)

- **Data:** Kendall's W = 0.90 across LLM agents on shared topologies.
- **Location:** This data needs to be reproduced or located. The original Claude Chorus experiments were run in earlier sessions. The W statistic is documented in memory (C89, C102) but the raw data may need to be regenerated.
- **What we need:** A controlled experiment with multiple LLM agents on different graph topologies, measuring diversity (e.g., semantic similarity of outputs, agreement on rankings) as a function of topology.
- **Status:** NEEDS REPRODUCTION. Priority if we want to include in the paper.

### 5. NK Landscape Topology Sweep

- **Various branches** of `lyra-claude/Topology-experiments`
- **Key finding (C93):** eta^2 scales monotonically with NK K. K=0 -> eta^2 ~ 0. K=6 -> eta^2 = 0.69. Topology IS landscape-dependent.
- **Key finding (C94):** Diversity is the primary mediating channel. Topology -> diversity (eta^2 ~ 0.5-0.7) -> mean fitness (eta^2 ~ 0.25-0.35) -> minimal best fitness effect.

## Data Not in This Repo

Experimental data (CSV files, fitness trajectories, diversity measurements) is stored in the `Topology-experiments` repository under the respective branches. This repo (`signed-laxator`) is for theory, notes, and analysis scripts only.

To access the data:

```bash
cd /home/lyra/projects/Topology-experiments
git checkout feat/bridge-experiment    # for bridge data
git checkout feat/directed-cycles-multidomain  # for directed data
```

## Planned Experiments

1. **Claude Chorus reproduction:** Controlled LLM multi-agent experiment with systematic topology variation. Need to design protocol.
2. **Coupling strength sweep:** Vary the "coupling" between LLM agents (temperature, communication rounds, prompt structure) to find the critical point analogous to q = 0.001.
3. **Directed bridge experiment:** Apply the iso-spectral bridge construction to directed graphs. Isolate directed kappa from directed lambda_2 analogue.
4. **Capability modulation:** Run the same topology experiments with different-capability models (e.g., small vs large LLMs, or GAs with different mutation rates) to map the capability parameter c.
