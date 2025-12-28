# LTEE Analysis Results

## Protocol Execution Summary

**Date:** December 2025  
**Configuration Hash:** `bff9d476e6638921`  
**Framework:** TCGS-SEQUENTION  
**Data Source:** E. coli LTEE (Good et al. 2017, Tenaillon et al. 2016)

---

## Overall Decision: **SUPPORT**

Significant synchrony detected in **7 genes** across **strictly independent** populations, not explained by mutator status or environmental covariates.

---

## Results Summary

### Genes Passing ALL Three Gates

| Gene | Populations | S* | p-value | Spread (gen) |
|------|-------------|-----|---------|--------------|
| **rbs_operon** | 12/12 | 1.000 | 0.0014 | 3,500 |
| **spoT** | 12/12 | 1.000 | 0.0034 | 3,000 |
| **ybaL** | 11/12 | 1.000 | 0.0039 | 3,500 |
| **topA** | 10/12 | 1.000 | 0.0092 | 3,000 |
| **mreB** | 7/12 | 1.000 | 0.0236 | 4,000 |
| **iclR** | 8/12 | 1.000 | 0.0265 | 3,000 |
| **fis** | 8/12 | 1.000 | 0.0431 | 2,500 |

**Combined probability:** p < 10⁻¹²

### Genes NOT Passing (for comparison)

| Gene | Populations | S* | p-value | Reason |
|------|-------------|-----|---------|--------|
| pykF | 12/12 | 0.583 | 0.0984 | S* < 0.6 |
| nadR | 9/12 | 0.556 | 0.3477 | S* < 0.6, p > 0.05 |
| hslU | 7/12 | 0.857 | 0.1528 | p > 0.05 |
| mrdA | 6/12 | 0.833 | 0.2688 | p > 0.05 |
| infB | 5/12 | 0.600 | 0.6033 | p > 0.05 |
| arcA | 6/12 | 0.667 | 0.4515 | p > 0.05 |

---

## Three-Gate Results

### Gate G1: Synchrony ✅
- **7 of 13 genes pass** (S* ≥ 0.6 AND p ≤ 0.05)
- Synchrony window: ±2,000 generations
- All passing genes show S* = 1.0 (perfect synchrony)

### Gate G2: Independence ✅
- **AUTOMATICALLY PASSES** by LTEE experimental design
- 12 populations in separate flasks since 1988
- Zero physical contact between populations
- Max IBD: 0.0 | Min FST: 0.95

### Gate G3: Common-Cause Exclusion ✅
- **All 7 G1-passing genes also pass G3**
- Synchrony retained: 100% for all genes
- ΔAIC: -20 to -24 (strongly favors Sequention model)
- Chart drift: < 0.6 for all genes

---

## Files in This Folder

| File | Description |
|------|-------------|
| `ltee_protocol_result.json` | Complete results in JSON format |
| `ltee_protocol_results.png` | Visualization of results |
| `README.md` | This file |

---

## Interpretation

The detection of perfect synchrony (S* = 1.0) across 7 genes in populations that have been **physically isolated for over 35 years** suggests coordination beyond what standard evolutionary mechanisms can explain. Gate G2 (independence) is guaranteed by the LTEE experimental design, eliminating the possibility of gene flow or contamination.

This pattern is consistent with the TCGS-SEQUENTION prediction of Synchronous Parallel Emergence (SPE) — coordinated evolutionary changes mediated by geometric constraints in counterspace.

---

## Citation

```bibtex
@misc{ArellanoPena2025LTEEResults,
  author = {Arellano-Peña, Henry},
  title = {{TCGS-SEQUENTION}: Complete Three-Gate Analysis of the 
           {E. coli} Long-Term Evolution Experiment},
  year = {2025},
  howpublished = {GitHub repository},
  url = {https://github.com/hyrucanji77/TCGS-SEQUENTION-protocol}
}
```
