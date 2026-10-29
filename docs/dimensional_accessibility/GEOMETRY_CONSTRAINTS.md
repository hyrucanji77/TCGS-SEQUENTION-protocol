# Hyperbolic geometry and power-law constraint supplement

Henry Arellano-Peña. Companion to Sections 6.3 and 14.4 and Appendices E.4 and F.3 of *Dimensional Accessibility and Gravitational Localization*, 29 October 2026.

## Run

From the repository root:

```sh
python -m pip install -r docs/dimensional_accessibility/requirements.txt
python code/verify_geometry_constraints.py --output results/dimensional_accessibility
python -m unittest discover -s tests -p 'test_geometry_constraints.py' -v
```

The article package uses `--output results`. This script retains the 118 formula/unit-conversion comparisons. The two earlier suites supply 186 and 225 comparisons, for **529 numerical checks**. There are **nine separate software tests** across the two paper-specific test files. Numerical checks do not certify general proofs or empirical interpretations.

## Geometry

For A_gamma(z)=(1+z/ell)^(-gamma), proper-normal differentiation gives tangential and normal sectional curvatures -gamma^2/ell^2*(1+z/ell)^(2*gamma-2) and -gamma/ell^2*(1+z/ell)^(2*gamma-2). With n>=2, equality forces gamma=0 or 1. Only gamma=1 localizes a constant scalar mode, and only for n>2. The selected domain is a horoball in H^(n+1); the boundary-value problem is not the entire maximally symmetric space.

## Experimental-template conversion and history

The exact potential is -G*m1*m2/r*[1+2*ell/(pi*r)]. In the published k=2 convention, ell=pi*beta2*r_ref/2. The fractional extra force is twice the fractional extra potential. For r_ref=1 mm:

| Experimental attribution | beta2 magnitude limit | Conditional ell limit | Coefficient provenance |
| --- | ---: | ---: | --- |
| Hoskins et al. (1985) | 1.3e-3 | 2.04204 micrometres | Adelberger 2007 Table I, previous k=2 entry, reference 15; also Hoyle 2004 Table XIV, reference 9 |
| Hoyle et al. (2004) | 3.6e-3 | 5.65487 micrometres | Hoyle 2004 Table XIV, this work |
| Adelberger et al. (2007) | 4.5e-4 | 0.706858 micrometres | Adelberger 2007 Table I, this work |

All three entries use the **68% convention printed in the cited comparison tables**. No original Hoskins likelihood was reconstructed. Hoyle's value was not the tightest pre-2007 k=2 bound. The unchanged 2007 benchmark is not asserted to be the strongest current result. Bsaibes and Decca (2024), DOI 10.3390/metrology4020014, report a later power-law reanalysis without improvement over existing constraints. Murata et al., arXiv:2605.18212v2, Section 5, distinguish reconstructed estimates from direct experimental confidence limits. Neither supplies a new numerical input here.

The transfer requires the same critical, sharp-readout, universal, unscreened, linear pairwise template. A single-Yukawa range limit is not a power-law coefficient limit. Published inputs and table/page locators are in `published_powerlaw_limits.json` in the selected output directory.

## Conditional cross-sector hierarchy

If the identical horosphere sets a_star=xi*c^2*|H| with |H|=eta_H/ell, then xi*eta_H=a_star*ell/c^2. For ell at the 2007 upper benchmark and the illustrative central calibration a_star=1.2e-10 m/s^2 (McGaugh et al. 2016, DOI 10.1103/PhysRevLett.117.201101), the product is at most about 9.44e-34, while c^2/ell is at least about 1.27e23 m/s^2. Here eta_H=1 for averaged principal curvature or 3 for the trace. This is a conditional dimensional comparison, not a joint-confidence interval, a fit, or a derivation of the transition rule. The machine-readable calculation is `cross_sector_hierarchy.json`.

## Release

This suite and its provenance outputs are included in the complete public commit pinned in the manuscript. The earlier `e6257110b8fd2c913f4ce788ced85690204fa158` release continues to identify the unchanged first two suites, not this geometry suite.
