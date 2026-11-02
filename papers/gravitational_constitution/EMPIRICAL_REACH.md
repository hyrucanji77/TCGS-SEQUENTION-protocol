# Empirical-reach calculations: scope and provenance

The five suites in computational release 4 accompany manuscript version 5. The first four programs retain their release-3 bytes. The fifth preserves all 80 empirical-reach checks from manuscript version 4 and adds 20 threshold tests, giving 276 checks in total.

## Tolerance thresholds

The diagnostic interval is r_E (C_p/tau_n)^(1/(2p)) <= r <= r_J tau_t^(1/3). Its endpoints touch at tau_t,crit = (r_E/r_J)^3 (C_p/tau_n)^(3/(2p)), or equivalently tau_n,crit = C_p (r_E/r_J)^(2p) tau_t^(-2p/3). A finite interval requires a strict inequality on the varied tolerance. `tolerance_thresholds.json` records both the p=1 and p=2 scenarios. The tests verify the equal-endpoint identities and open/closed cases on both sides.

For p=2 and C_2=1, the endpoint ratio under both 0.001 allocations is 1.27101210. With tau_n held at 0.001 the opening threshold is tau_t=0.00205328414; with tau_t held at 0.001 it is tau_n=0.00260974898. The quadratic window is therefore tolerance-sensitive. The local parity argument motivates examining p=2 but does not prove its coefficient or a uniform error bound for stellar interiors. Relaxing separate budgets must still be consistent with the covariance and bias requirements of the matched-gradient residual.

## Inputs and limitations

The nominal solar mass parameter and astronomical lengths follow Prša et al., AJ 152, 41 (2016), DOI 10.3847/0004-6256/152/2/41, arXiv:1605.09788v1. Adopted values are GM_sun=1.3271244e20 m^3 s^-2, au=149597870700 m, pc=648000 au/pi, a_star=1.2e-10 m s^-2 and y_b=1.8. The latter two are sensitivity inputs, not fits.

The 1.7 pc tidal benchmark follows Jiang and Tremaine, MNRAS 401, 977–994 (2010), DOI 10.1111/j.1365-2966.2009.15744.x, arXiv:0907.2952. It is a Solar-neighborhood Newtonian Jacobi-radius proxy for two solar masses, not a universal MOND orbit radius or irreducible noise floor.

Banik and Zhao, arXiv:1805.12273v9, section 7.3, discuss cancellation of certain first-order mass-ratio effects. Neither those cancellations nor the local analytic-branch parity argument provides a uniform two-body error estimate for compact stars. No feasible binary observing window is certified by this supplement.

The NGC 2419 environmental estimate of order 0.1 a_0 follows Ibata et al., arXiv:1106.4909v1. No mass slope at the selected matching radius is inferred from this qualitative estimate. Environmental and mass-profile uncertainties must be propagated to the actual residual.

## Source checks

Partanen–Tulkki, Reports on Progress in Physics 88, 057802 (2025), printed p. 37, equation (184), prints the ultraviolet coefficient -kappa c p^2/(24 pi^2 hbar epsilon_UV). This is a source transcription check, not an independent loop calculation. The May corrigendum changes interaction factors; the operator obstruction requires only a nonzero extra-momentum coefficient.

The rendered page 4 of arXiv:2505.14446v6 confirms equations (26)–(28): second-order coefficients 23 pi/8, 11 pi/4 and 15 pi/4. The pinned version is dated 16 September 2025.

Lee's dissertation Table 11.4, printed p. 134, gives absolute-value 95% limits 2.08 and 1.45 at 32 and 35 micrometres, and positive delta-chi-squared=6.17 boundaries 2.56 and 1.77. The university repository's extracted text reproduces these rows and the distinct confidence conventions. A direct PDF screenshot was not retrievable during this release check; the values were checked in that repository text, not in an independently refitted likelihood. Source: https://digital.lib.washington.edu/researchworks/bitstreams/6013936a-da00-42b7-ba0a-4dff0cb05bf8/download

The width is a historical single-Yukawa benchmark; none of the suites performs a full-tower apparatus likelihood, an empirical data fit, or an ontological test.
