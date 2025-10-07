# LLM-Generated LIMA VNS for Balanced Minimum Sum-of-Squares Clustering

This repository contains a **fully LLM-generated implementation** of the LocalSearch and VNS algorithms from the paper "Less is more: basic variable neighborhood search heuristic for balanced minimum sum-of-squares clustering" by Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.

## Experiment Overview

**Objective:** Demonstrate that an LLM can implement state-of-the-art metaheuristics from academic literature with zero human coding intervention.

**What the LLM Generated:**
- `LocalSearch.cpp` - Complete local search implementation 
- `VNS.cpp` - Complete Variable Neighborhood Search implementation

**What was Provided to the LLM:**
- Original paper (Costa et al., 2017)
- Essential framework files (solution representation, objective function, etc.)
- C++ header files for LocalSearch and VNS
- Supporting documentation

**Original Paper:** https://doi.org/10.1016/j.ins.2017.06.019

**Original Implementation:** Available at [separate repository link]

## LLM Development Process

Complete chat log showing the step-by-step LLM implementation process using our **Communication, Context, and Iteration** framework:
- **[LLM Implementation Session](https://claude.ai/share/50ac7f73-df99-4c96-a52a-50fa6796053f)** - Shows the LLM generating LocalSearch.cpp and VNS.cpp from paper description and supporting files

## Fair Comparison Setup

To ensure accurate performance comparison with the original LIMA VNS:

- **Identical Initial Solutions:** The `initial_solutions/` folder contains the same starting solutions used to benchmark the original implementation
- **Modified Execution:** `run_with_init` script runs experiments using these pre-generated initial solutions
- **Same Test Framework:** Uses identical datasets, parameters, and evaluation metrics

## Compiling

You only have to use the makefile in the *"src"* folder:

```bash
make
```

To clean the compiled files use:

```bash
make clean
```

## Executing

### Standard Execution

32-bits version:
```bash
./lima_vns_32 <path/instance.csv> <k=number of clusters> <cpu time limit> <number of runs> <seed> <path/output file> <path/cluster assignment file>
```

64-bits version:
```bash
./lima_vns_64 <path/instance.csv> <k=number of clusters> <cpu time limit> <number of runs> <seed> <path/output file> <path/cluster assignment file>
```

The cluster assignments file has the instance used as the first column. The remaining columns are the cluster assignments, for example, "iris.csv, 0, 0, 1, ...", which means:
- point 1 is assigned to cluster 0;
- point 2 is assigned to cluster 0;  
- point 3 is assigned to cluster 1;

### Controlled Comparison (Recommended)
Run with identical initial solutions for fair comparison:
```bash
./run_with_init
```

This script:
- Uses pre-generated initial solutions from `initial_solutions/` folder
- Runs 10 experiments per dataset
- Saves results to `results/` and assignments to `assignments/`
- Matches the experimental setup used for the original implementation

### Analyze Results
```bash
./analyze_results
```

Generates summary statistics (best, worst, mean, variance) for all experimental runs.

## Repository Structure

```
├── src/                    # Source code
│   ├── LocalSearch.cpp     # 🤖 LLM-Generated
│   ├── VNS.cpp            # 🤖 LLM-Generated  
│   ├── [framework files]   # Original supporting code
│   └── Makefile
├── datasets/              # Test instances
├── initial_solutions/     # Pre-generated starting solutions
├── run_with_init         # Modified execution script
├── analyze_results       # Results analysis script
├── results/              # Output directory (created on run)
└── assignments/          # Cluster assignments (created on run)
```

## Results

The LLM-generated implementation achieves performance **matching the original LIMA VNS algorithm** across all test instances, demonstrating successful replication of state-of-the-art metaheuristic performance through LLM-assisted development.

## Framework: Communication, Context, and Iteration

- **Communication:** Breaking implementation into steps, testing iteratively, clear objectives
- **Context:** Providing paper, essential files, and C++ headers to guide development  
- **Iteration:** Feeding compile errors and runtime results back to LLM for self-correction

## Citation

If you use this LLM-generated implementation, please cite the original paper:

```bibtex
@article{costa2017lima,
  title={Less is more: basic variable neighborhood search heuristic for balanced minimum sum-of-squares clustering},
  author={Costa, Leandro R and Aloise, Daniel and Mladenovi{\'c}, Nenad},
  journal={Information Sciences},
  volume={415},
  pages={104--118},
  year={2017},
  publisher={Elsevier}
}
```
