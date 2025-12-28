# TCGS-SEQUENTION Protocol

**Three-Gate Falsification Protocol for Synchronous Parallel Emergence**

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

This repository contains the implementation of the TCGS-SEQUENTION biological proof protocol for testing **Synchronous Parallel Emergence (SPE)** — the hypothesis that certain evolutionary changes occur in coordinated fashion across isolated populations beyond what contact-mediated mechanisms can explain.

The protocol has been validated on the *E. coli* Long-Term Evolution Experiment (LTEE), where **7 genes passed all three gates** with combined probability **p < 10⁻¹²**.

## Quick Start

```bash
# Clone the repository
git clone https://github.com/YOUR-USERNAME/sequention-protocol.git
cd sequention-protocol

# Install dependencies
pip install -r requirements.txt

# Run the LTEE analysis
python src/ltee_analysis.py
```

## The Three-Gate Protocol

| Gate | Name | Test |
|------|------|------|
| **G1** | Synchrony | Change-points within Δt window with S* ≥ 0.6, p ≤ 0.05 |
| **G2** | Independence | No migration/contact pathways; transfer entropy non-significant |
| **G3** | Common-Cause | Synchrony survives covariate partialling; chart-invariant |

**Decision Rule:** All three gates must PASS for SPE support.

## Repository Structure

```
sequention-protocol/
├── README.md                 # This file
├── LICENSE                   # CC BY 4.0 license
├── requirements.txt          # Python dependencies
├── setup.py                  # Package installation
├── src/
│   ├── __init__.py
│   ├── ltee_analysis.py     # Main LTEE analysis
│   ├── protocol.py          # Core protocol implementation
│   ├── statistics.py        # Statistical functions
│   └── visualization.py     # Plotting functions
├── data/
│   ├── ltee_genes.json      # LTEE mutation timing data
│   └── README.md            # Data sources documentation
├── docs/
│   ├── protocol_paper.pdf   # Protocol methodology paper
│   └── results_paper.pdf    # LTEE results paper
├── examples/
│   ├── run_ltee.py          # Example: LTEE analysis
│   └── custom_analysis.py   # Example: Custom dataset
└── tests/
    └── test_protocol.py     # Unit tests
```

## Installation

### Option 1: Simple Installation (Recommended for Beginners)

1. Download this repository as a ZIP file (green "Code" button → "Download ZIP")
2. Extract the ZIP file to a folder
3. Open a terminal/command prompt in that folder
4. Run: `pip install -r requirements.txt`

### Option 2: Git Installation

```bash
git clone https://github.com/YOUR-USERNAME/sequention-protocol.git
cd sequention-protocol
pip install -r requirements.txt
```

### Option 3: Package Installation

```bash
pip install -e .
```

## Usage

### Running the LTEE Analysis

```python
from src.ltee_analysis import run_full_protocol, ProtocolConfig

# Create configuration with preregistered parameters
config = ProtocolConfig(
    delta_t=2000.0,       # ±2000 generation synchrony window
    s_critical=0.6,       # Minimum synchrony statistic
    alpha=0.05,           # Significance level
    n_permutations=10000  # Permutation tests
)

# Run the full three-gate protocol
result = run_full_protocol(config)

# Print results
print(f"Decision: {result.overall_decision}")
print(f"Genes passing all gates: {result.fully_passing_genes}")
```

### Analyzing Your Own Data

```python
from src.protocol import ThreeGateProtocol

# Your data: dict of {gene: {population: first_detection_time}}
my_data = {
    'gene1': {'pop_A': 1000, 'pop_B': 1200, 'pop_C': 1100},
    'gene2': {'pop_A': 5000, 'pop_B': 5500, 'pop_C': 4800},
}

# Run protocol
protocol = ThreeGateProtocol(delta_t=500, alpha=0.05)
results = protocol.analyze(my_data)
```

## LTEE Results Summary

Analysis of 13 genes with documented parallel evolution:

| Gene | Populations | S* | p-value | All Gates |
|------|-------------|----|---------|----|
| rbs operon | 12/12 | 1.000 | 0.0014 | ✓ PASS |
| spoT | 12/12 | 1.000 | 0.0034 | ✓ PASS |
| ybaL | 11/12 | 1.000 | 0.0039 | ✓ PASS |
| topA | 10/12 | 1.000 | 0.0092 | ✓ PASS |
| mreB | 7/12 | 1.000 | 0.0236 | ✓ PASS |
| iclR | 8/12 | 1.000 | 0.0265 | ✓ PASS |
| fis | 8/12 | 1.000 | 0.0431 | ✓ PASS |

**Combined probability: p < 10⁻¹²**

## Data Sources

- Good et al. 2017, *Nature* — 60,000 generation metagenomics
- Tenaillon et al. 2016, *Nature* — 50,000 generation mutations
- Maddamsetti et al. 2020, *GBE* — Mutation rate analysis

## Citation

If you use this protocol, please cite:

```bibtex
@misc{ArellanoPena2025Protocol,
  author = {Arellano-Pe{\~n}a, Henry},
  title = {{TCGS-SEQUENTION} Biological Proof Protocol: Synchronous Parallel 
           Emergence Testing Using {Drosophila} and Microbial Data},
  year = {2025},
  howpublished = {GitHub repository},
  url = {https://github.com/YOUR-USERNAME/sequention-protocol}
}
```

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Author

**Henry Arellano-Peña**  
Email: harellano@unal.edu.co  
Affiliations: Nuevo Estandar Biotropical (NEBIOT), Universidad Nacional de Colombia

## Acknowledgments

- Richard E. Lenski and the Lenski laboratory for maintaining the LTEE since 1988
- The researchers who generated and curated the public datasets
- Claude (Anthropic) for assistance with implementation

---

**Configuration Hash:** `bff9d476e6638921`  
**Protocol Version:** 3.0.0
