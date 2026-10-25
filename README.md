# QAOA for Weighted & Unweighted Max-Cut

A from-scratch NumPy/SciPy implementation of the **Quantum Approximate Optimization Algorithm (QAOA)** applied to the **Max-Cut** problem, on both unweighted and weighted graphs. Results are validated against an exact brute-force solver, and the effect of QAOA layer depth `p` is analyzed.

Originally developed as the final project for a Computational Physics course.

**Author:** Julia Huynh

---

## Table of Contents

- [Background](#background)
- [How It Works](#how-it-works)
- [Repository Contents](#repository-contents)
- [Getting Started](#getting-started)
- [Results](#results)
- [Analysis of Layer Depth (p)](#analysis-of-layer-depth-p)
- [Limitations & Lessons Learned](#limitations--lessons-learned)
- [Future Work](#future-work)

---

## Background

**Max-Cut** asks: given an undirected graph, how do you split the vertices into two groups so that the number of edges crossing between the groups (the *cut value*) is as large as possible? In the **weighted** version, each edge has a weight and the goal is to maximize the *sum of the weights* of the cut edges.

Max-Cut is NP-Hard: the number of candidate solutions is 2^n for n vertices, so brute force quickly becomes intractable. It shows up in problems like circuit/chip design and network partitioning.

**QAOA** is a hybrid quantum-classical algorithm that encodes an optimization problem as a **cost Hamiltonian** whose ground state is the optimal solution, then variationally searches for a quantum state close to that ground state.

## How It Works

Bitstrings represent solutions: index *i* of the tuple is vertex *i*, and vertices with the same bit are in the same group. For example, `(0, 1, 0, 1)` puts vertices 0 and 2 in one group and vertices 1 and 3 in the other.

**1. Cost Hamiltonian**

Unweighted:

$$H_C = \sum_{(i,j)\in E} \tfrac{1}{2}\left(I - Z_i Z_j\right)$$

Weighted:

$$H_C = \sum_{(i,j)\in E} w_{ij}\cdot \tfrac{1}{2}\left(I - Z_i Z_j\right)$$

**2. Initial state**: a uniform superposition over all bitstrings:

$$|+\rangle^{\otimes n} = \frac{1}{\sqrt{2^n}}\sum_{x\in\{0,1\}^n}|x\rangle$$

**3. Mixer Hamiltonian**: lets amplitude move between bitstrings by flipping bits independently:

$$H_M = \sum_i X_i$$

**4. Layered evolution**: alternate cost and mixer unitaries `p` times:

$$|\psi(\gamma,\beta)\rangle = \prod_{k=1}^{p} e^{-i\beta_k H_M}\, e^{-i\gamma_k H_C}\,|+\rangle^{\otimes n}$$

The cost unitary imprints phases that, through interference, amplify good cuts and suppress bad ones; the mixer spreads amplitude between bitstrings.

**5. Parameter optimization**: minimize the energy expectation value ⟨ψ|H|ψ⟩ (using the negative of the cut Hamiltonian so that minimizing energy maximizes the cut) over the `2p` parameters `γ, β`. Two optimizers are implemented and compared:
- A custom optimization routine written for this project
- SciPy's `minimize` with the **Powell** method

**6. Measurement**: build the probability distribution from the final state, randomly sample bitstrings, and rank them by cut value and frequency.

**7. Verification**: a brute-force solver enumerates all 2^n bitstrings (via `itertools.product`) to find the exact optimal cut and every bitstring that achieves it.

## Repository Contents

```
.
├── qaoa_maxcut.ipynb      # Main notebook: implementation, experiments, plots
├── README.md
├── requirements.txt       # (optional) Python dependencies
└── images/                # (optional) Plots/graph visualizations for this README
```

Key functions in the notebook:

| Function | Purpose |
|---|---|
| `cut_value` / `weighted_cut_value` | Classical cut value of a bitstring for a given adjacency matrix |
| `brute_force_maxcut` | Exact solver that checks all 2^n bitstrings |
| `create_weighted_hamiltonian` | Builds the (weighted) cost Hamiltonian from an adjacency matrix using Kronecker products of Pauli-Z operators |
| Mixer / QAOA layer functions | Build the mixer and apply `p` alternating layers to the initial state |
| Optimizers | Custom optimizer and SciPy Powell for tuning `γ, β` |
| Sampling & plotting utilities | Sample from the final state, report top cuts, draw the graph and the partition |

> Update the filename and function names above to match your notebook.

## Getting Started

### Prerequisites

- Python 3.9+
- Jupyter

### Installation

```bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>
pip install numpy scipy matplotlib networkx jupyter
```

### Run

```bash
jupyter notebook qaoa_maxcut.ipynb
```

Graphs are defined as **adjacency matrices** (`adj[i][j]` is the edge weight, or `1`/`0` for unweighted). To try your own graph, define a new matrix and pass it through the same pipeline: build the Hamiltonian, optimize `γ, β`, sample, and compare against `brute_force_maxcut`.

> Note: the simulation builds full 2^n × 2^n matrices, so it is intended for small graphs (roughly 4-10 vertices).

## Results

Seven test graphs were used: **4 unweighted** and **3 weighted**, each with 4-5 vertices.

**QAOA matched brute force on every test case.** The maximum-cut bitstrings were found and, in most cases, appeared with high frequency in the sampled distribution. The graph partitions drawn from the QAOA results agree with the brute-force solutions.

### Unweighted Max-Cut

| Graph | Max cut (brute force) | Best QAOA energy |
|---|---|---|
| 4-node graph, 5 edges ("diamond") | 4 | 3.234 |
| 5-node star | 4 | 2.993 |
| 5-node path | 4 | 3.113 |
| 4-node complete graph (K4) | 4 | 3.697 |

### Weighted Max-Cut

| Graph | Max cut (brute force) | Best QAOA energy |
|---|---|---|
| 4-node, heavily weighted | 12 | 11.520 |
| 4-node cycle | 10 | 7.652 |
| 5-node weighted graph | 12 | 9.209 |

The energy is an *expectation value* over the whole output distribution, so it sits below the true max cut even when the optimal bitstrings are the most likely ones to be sampled. Weighted cases show larger absolute gaps simply because cut values are larger.

## Analysis of Layer Depth (p)

QAOA was run for `p = 1, 2, 3` on an unweighted test graph (max cut = 4), tracking the best sampled cut, the average cut value, and the probability of sampling an optimal cut.

| p | Probability of optimal cut | Average cut value | Best cut found |
|---|---|---|---|
| 1 | ~0.34 | ~3.22 | 4 |
| 2 | ~0.70 | ~3.47 | 4 |
| 3 | ~0.78 | ~3.72 | 4 |

(Values read from the plots; exact numbers vary run to run.)

Increasing `p` raises both the average cut value and the probability of measuring an optimal solution, as expected. Because sampling is random, individual runs differ slightly, but the trend is consistent.

<!-- Add your plots here:
![Probability of optimal cut vs p](images/prob_vs_p.png)
![QAOA performance vs p](images/performance_vs_p.png)
-->

## Limitations & Lessons Learned

- **Brute force wins at this scale.** For 4-5 vertex graphs, brute force is both exact and faster. QAOA's potential advantage is for larger graphs where enumeration is infeasible.
- **More layers means harder optimization.** The number of parameters is `2p`, so the classical optimization landscape gets harder as `p` grows. For larger graphs, the parameters must also steer a much larger quantum state.
- **A sign bug taught a lesson.** An earlier version of the code was missing a negative sign when defining the energy to minimize, which made the optimizer work against the goal and produced misleading `p` results (including an apparent failure of the custom optimizer). After the fix, the custom optimizer and SciPy's Powell method give comparable results, and the expected positive trend with `p` appears.
- **Sampling noise.** Results are drawn by random sampling, so frequencies change between runs.

## Future Work

- Larger and more complex graphs for Max-Cut
- Other NP-hard applications of QAOA, such as portfolio optimization
- Better optimizers and parameter initialization for deeper circuits
- Running on a real quantum simulator or hardware backend (e.g., Qiskit)

## License

Add a license of your choice (e.g., MIT).
