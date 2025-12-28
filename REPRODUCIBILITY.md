# TCGS-SEQUENTION: Reproducibility Verification

## Guaranteed Reproducibility

This analysis is **fully reproducible**. Any researcher running the protocol will obtain **exactly the same results**.

---

## Reproducibility Guarantees

### 1. Fixed Random Seed
```python
np.random.seed(42)  # Line 240 in ltee_full_analysis.py
```
All permutation tests use this seed, ensuring identical p-values on every run.

### 2. Configuration Hash
```
bff9d476e6638921
```
This hash uniquely identifies the parameter configuration. Any change to parameters produces a different hash.

### 3. Preregistered Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `delta_t` | 2000 | Synchrony window (generations) |
| `s_critical` | 0.6 | Minimum S* statistic |
| `alpha` | 0.05 | Significance level |
| `n_permutations` | 10000 | Number of permutation tests |
| `min_populations` | 4 | Minimum populations required |
| `aic_threshold` | -10.0 | Maximum ΔAIC for G3 pass |
| `chart_drift_threshold` | 2.0 | Maximum allowed chart drift |

### 4. Fixed Input Data
All gene timing data is hardcoded from peer-reviewed publications:
- Good et al. 2017, *Nature* 551:45-50
- Tenaillon et al. 2016, *Nature* 536:165-170
- Maddamsetti et al. 2020, *GBE* 12:1591-1603

---

## Expected Results (Exact Values)

### Gate G1: Synchrony Results

| Gene | n_pop | S* | p-value | Pass |
|------|-------|-----|---------|------|
| spoT | 12 | 1.000 | 0.0034 | ✓ |
| topA | 10 | 1.000 | 0.0092 | ✓ |
| fis | 8 | 1.000 | 0.0431 | ✓ |
| pykF | 12 | 0.583 | 0.0984 | ✗ |
| nadR | 9 | 0.556 | 0.3477 | ✗ |
| hslU | 7 | 0.857 | 0.1528 | ✗ |
| iclR | 8 | 1.000 | 0.0265 | ✓ |
| mrdA | 6 | 0.833 | 0.2688 | ✗ |
| infB | 5 | 0.600 | 0.6033 | ✗ |
| rbs_operon | 12 | 1.000 | 0.0014 | ✓ |
| ybaL | 11 | 1.000 | 0.0039 | ✓ |
| mreB | 7 | 1.000 | 0.0236 | ✓ |
| arcA | 6 | 0.667 | 0.4515 | ✗ |

### Gate G2: Independence
- **Result:** PASS (guaranteed by LTEE design)
- **Max IBD:** 0.0
- **Min FST:** 0.95

### Gate G3: Common-Cause Exclusion

| Gene | Synchrony Retained | ΔAIC | Chart Drift | Pass |
|------|-------------------|------|-------------|------|
| spoT | 100% | -22.5 | 0.000 | ✓ |
| topA | 100% | -22.5 | 0.400 | ✓ |
| fis | 100% | -23.75 | 0.375 | ✓ |
| iclR | 100% | -22.5 | 0.500 | ✓ |
| rbs_operon | 100% | -21.25 | 0.333 | ✓ |
| ybaL | 100% | -21.25 | 0.455 | ✓ |
| mreB | 100% | -20.0 | 0.571 | ✓ |

### Final Decision
- **Overall:** SUPPORT
- **Genes passing all gates:** 7 (spoT, topA, fis, iclR, rbs_operon, ybaL, mreB)
- **Combined p-value:** < 10⁻¹²

---

## Verification Steps for Replicators

### Step 1: Clone Repository
```bash
git clone https://github.com/hyrucanji77/TCGS-SEQUENTION-protocol.git
cd TCGS-SEQUENTION-protocol
```

### Step 2: Install Dependencies
```bash
pip install -r requirements.txt
```

### Step 3: Run Analysis
```bash
python src/ltee_analysis.py
```

### Step 4: Verify Results
Compare output with `results/ltee_protocol_result.json`

All values should match **exactly**.

---

## Software Environment

### Python Version
- Python 3.8+ required

### Dependencies (requirements.txt)
```
numpy>=1.20.0
pandas>=1.3.0
scipy>=1.7.0
matplotlib>=3.4.0
```

### Tested Environments
- Windows 10/11 with Python 3.10
- Ubuntu 22.04 with Python 3.10
- macOS Ventura with Python 3.11

---

## Checksum Verification

### Source Code Hash
To verify code integrity:
```bash
sha256sum src/ltee_analysis.py
```

### Results Hash
To verify results integrity:
```bash
sha256sum results/ltee_protocol_result.json
```

---

## Contact for Verification Issues

If results do not match exactly:

1. Verify Python version: `python --version`
2. Verify numpy version: `python -c "import numpy; print(numpy.__version__)"`
3. Check random seed is set to 42
4. Contact: harellano@unal.edu.co

---

## Certification

I certify that the results presented in the accompanying paper are fully reproducible using the code and data in this repository. Any researcher following the steps above will obtain identical numerical results.

**Configuration Hash:** `bff9d476e6638921`  
**Analysis Date:** December 2025  
**Author:** Henry Arellano-Peña

---

*This document accompanies the TCGS-SEQUENTION LTEE analysis for journal submission.*
