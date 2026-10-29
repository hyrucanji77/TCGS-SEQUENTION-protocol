#!/usr/bin/env python3
"""Check the hyperbolic candidate and published power-law-limit conversions.

No experimental data are fitted. Published beta_2 limits retain their stated
68% confidence convention. Checks concern formulas and unit conversions, not
statistical reanalysis, completeness of the literature, or source ontology.
Copyright 2026 Henry Arellano-Pena. CC BY 4.0; see accompanying license.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import mpmath as mp
from scipy.integrate import quad


def convert_beta2(beta2: float, reference_m: float = 1e-3) -> float:
    """Map a nonnegative potential coefficient to ell in metres."""
    if not math.isfinite(beta2) or beta2 < 0:
        raise ValueError('beta2 must be finite and nonnegative')
    if not math.isfinite(reference_m) or reference_m <= 0:
        raise ValueError('reference length must be finite and positive')
    return math.pi * beta2 * reference_m / 2


def sectional_curvatures(gamma: float, ell: float, z: float) -> tuple[float, float]:
    """Tangential and normal sectional curvature from the closed formulas."""
    if not all(math.isfinite(v) for v in (gamma, ell, z)) or ell <= 0 or z < 0:
        raise ValueError('finite parameters, ell > 0 and z >= 0 required')
    scale = (1 + z/ell)**(2*gamma-2) / ell**2
    return -gamma**2*scale, -gamma*scale


class Checks:
    def __init__(self) -> None:
        self.rows: list[dict] = []

    def close(self, name: str, actual: float, expected: float,
              rtol: float = 2e-10, atol: float = 1e-12) -> None:
        a, e = float(actual), float(expected)
        self.rows.append(dict(name=name, actual=a, expected=e, rtol=rtol,
                              atol=atol, passed=math.isfinite(a) and
                              math.isfinite(e) and math.isclose(a, e,
                                  rel_tol=rtol, abs_tol=atol)))


def run_checks() -> list[dict]:
    c = Checks()
    with mp.workdps(50):
        # Numerical differentiation of the warp, independent of the closed
        # gamma-dependent sectional-curvature formulas being tested.
        for gamma in (0., .5, 1., 1.5, 2.):
            for ell in (.3, 2.):
                for z_over_ell in (0., .7, 4.):
                    z = mp.mpf(str(z_over_ell)) * mp.mpf(str(ell))
                    length, power = mp.mpf(str(ell)), mp.mpf(str(gamma))
                    warp = lambda v: (1 + v/length)**(-power)
                    a, da, dda = warp(z), mp.diff(warp, z), mp.diff(warp, z, 2)
                    tangential = -da**2/a**4
                    normal = -(dda*a-da**2)/a**4
                    kt, kn = sectional_curvatures(gamma, ell, float(z))
                    label = f'g{gamma}_ell{ell}_z{z_over_ell}'
                    c.close('section_tangent_'+label, kt, tangential)
                    c.close('section_normal_'+label, kn, normal)
        for n in (3, 4, 6):
            for ell in (.3, 2.):
                # Ricci components and scalar trace from plane counts.
                kt, kn = sectional_curvatures(1., ell, .7*ell)
                c.close(f'ricci_tangent_n{n}_ell{ell}',
                        (n-1)*kt+kn, -n/ell**2)
                c.close(f'ricci_normal_n{n}_ell{ell}', n*kn, -n/ell**2)
                c.close(f'scalar_curvature_n{n}_ell{ell}',
                        n*(n-1)*kt+2*n*kn, -n*(n+1)/ell**2)
                norm = quad(lambda t: (1+t)**(-(n-1)), 0, math.inf,
                            epsabs=1e-12)[0] * ell
                c.close(f'normal_norm_n{n}_ell{ell}', norm, ell/(n-2))
        for ell in (.3, 2.):
            length = mp.mpf(str(ell))
            # (1/2) Lie_N h at w=ell, with h_ii=ell^2/w^2.
            k = length**2/2 * (length/length) * mp.diff(lambda w: w**-2, length)
            c.close(f'horosphere_second_form_ell{ell}', k, -1/ell)
            c.close(f'horosphere_gauss_flat_ell{ell}', -1/ell**2+float(k)**2, 0.)
            c.close(f'horosphere_induced_metric_ell{ell}', length**2/length**2, 1.)
        # Exact derivative, not a large-r approximation. G*m1*m2 is set to 1.
        for ell in (.1, 1., 3.):
            for r in (.02, 2., 80.):
                length, radius = mp.mpf(str(ell)), mp.mpf(str(r))
                potential = lambda v: -(1/v+2*length/(mp.pi*v**2))
                excess = mp.diff(potential, radius)*radius**2 - 1
                c.close(f'force_derivative_ell{ell}_r{r}', excess,
                        4*length/(mp.pi*radius))
                beta = 2*ell/(math.pi*1e-3)
                c.close(f'potential_template_ell{ell}_r{r}',
                        -(1/r)*(1+beta*1e-3/r), potential(radius))
        for n in (1, 2):
            for cutoff in (2., 20.):
                integral = quad(lambda t: (1+t)**(-(n-1)), 0, cutoff)[0]
                expected = cutoff if n == 1 else math.log1p(cutoff)
                c.close(f'nonlocalizing_truncated_norm_n{n}_T{cutoff}', integral, expected)
    # The reference length is an arbitrary convention: beta rescales inversely.
    for beta in (3.6e-3, 4.5e-4):
        length = convert_beta2(beta)
        c.close(f'reference_length_rescaling_{beta}',
                convert_beta2(beta/10, 1e-2), length)
    c.close('hoyle2004_ell_micrometres', convert_beta2(3.6e-3)*1e6,
            5.654866776461628)
    c.close('adelberger2007_ell_micrometres', convert_beta2(4.5e-4)*1e6,
            .7068583470577035)
    length = convert_beta2(4.5e-4)
    c.close('curvature_bound_ratio', (12/length**2)/(1/length**2), 12.)
    c.close('force_potential_equal_contribution_radius_ratio',
            (4*length/math.pi)/(2*length/math.pi), 2.)
    return c.rows


def published_limits() -> dict:
    """Published inputs with explicit table provenance and historical ranking.

    Hoskins' coefficient and confidence convention are transcribed from the
    later comparison tables, not claimed to have been independently extracted
    from the 1985 likelihood. No strongest-current-bound claim is made.
    """
    records = (
        ('hoyle2004', 'Hoyle et al. (2004)', 3.6e-3,
         '10.1103/PhysRevD.70.042004', 'XIV',
         'https://arxiv.org/abs/hep-ph/0405262',
         'Weaker k=2 comparison; not the strongest pre-2007 k=2 bound.'),
        ('adelberger2007', 'Adelberger et al. (2007)', 4.5e-4,
         '10.1103/PhysRevLett.98.131104', 'I',
         'https://arxiv.org/abs/hep-ph/0611223v3',
         'Main matched-template bound used in the manuscript.'),
        ('hoskins1985', 'Hoskins et al. (1985)', 1.3e-3,
         '10.1103/PhysRevD.32.3084', 'I (Adelberger 2007), XIV (Hoyle 2004)',
         'https://dcc.ligo.org/LIGO-T850003/public',
         'Tighter pre-2007 k=2 comparator identified in both later tables.'))
    rows = []
    for ident, name, beta, doi, table, url, role in records:
        ell = convert_beta2(beta)
        provenance = ({'doi': '10.1103/PhysRevLett.98.131104',
                       'table': 'I, k=2, previous constraints, reference 15',
                       'pdf_page_1based': 3,
                       'pdf': 'https://arxiv.org/pdf/hep-ph/0611223v3',
                       'corroboration': 'Hoyle 2004 Table XIV, k=2, previous work, reference 9'}
                      if ident == 'hoskins1985' else
                      {'doi': doi, 'table': table + ', k=2, this work',
                       'pdf_page_1based': 29 if ident == 'hoyle2004' else 3,
                       'pdf': ('https://arxiv.org/pdf/hep-ph/0405262' if ident == 'hoyle2004'
                               else 'https://arxiv.org/pdf/hep-ph/0611223v3')})
        rows.append(dict(id=ident, source=name, doi=doi, table=table,
                         source_url=url, historical_role=role,
                         coefficient_and_confidence_provenance=provenance,
                         published_confidence_percent=68, beta2_abs_limit=beta,
                         reference_length_m=1e-3, ell_limit_m=ell,
                         ell_limit_micrometres=ell*1e6,
                         comparison_sectional_curvature_magnitude_min_m2=1/ell**2,
                         comparison_H4_scalar_curvature_magnitude_min_m2=12/ell**2))
    return dict(template='U=-G*m1*m2/r*(1+beta2*r_ref/r)',
                mapping='ell=pi*beta2*r_ref/2; E=2*beta2*r_ref/r',
                assumptions=['critical nu=3/2', 'sharp boundary readout',
                             'universal unscreened pairwise mass coupling',
                             'no independent continuum amplitude'],
                confidence_note='Inherited published magnitude constraints; no new likelihood, no 95% claim.',
                scope='Historical comparisons, not strongest-current-bound claims.',
                sources=rows,
                later_literature=[
                    dict(source='Bsaibes and Decca (2024)', doi='10.3390/metrology4020014',
                         role='Power-law reanalysis; authors report no improvement over existing limits.',
                         quantitative_beta2_input_used=False),
                    dict(source='Murata, Fujiie, and Suzuki (2026)', arxiv='2605.18212v2',
                         role='Section 5 distinguishes reconstructed estimates from directly fitted confidence limits.',
                         quantitative_beta2_input_used=False)],
                last_literature_check='2026-10-29')


def cross_sector_hierarchy(ell_m: float, a_star: float) -> dict:
    """Conditional scale comparison, not a fitted or predicted AQUAL scale.

    For |H|=eta_H/ell and a_star=xi*c^2*|H|, return xi*eta_H.
    With an upper bound on ell, this product is an upper value, and c^2/ell
    is a lower value. The caller supplies the illustrative acceleration.
    """
    if not all(math.isfinite(v) and v > 0 for v in (ell_m, a_star)):
        raise ValueError('ell and a_star must be finite and positive')
    c = 299792458.0
    return dict(ell_m=ell_m, a_star_m_s2=a_star, c_m_s=c,
                comparison_acceleration_m_s2=c*c/ell_m,
                xi_times_eta_H=a_star*ell_m/(c*c),
                unsuppressed_curvature_radius_m=c*c/a_star,
                scale_ratio=c*c/(a_star*ell_m))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'results')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    checks = run_checks()
    failed = sum(not row['passed'] for row in checks)
    report = dict(checks=len(checks), passed=len(checks)-failed, failed=failed,
                  python=platform.python_version(), mpmath=mp.__version__,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  note='Finite numerical consistency checks; not proof or empirical certification.',
                  results=checks)
    (args.output / 'geometry_constraint_verification.json').write_text(
        json.dumps(report, indent=2)+'\n', encoding='utf-8')
    (args.output / 'published_powerlaw_limits.json').write_text(
        json.dumps(published_limits(), indent=2)+'\n', encoding='utf-8')
    hierarchy = cross_sector_hierarchy(convert_beta2(4.5e-4), 1.2e-10)
    hierarchy.update(
        status='Conditional dimensional comparison, not an empirical likelihood or derived ECL rule.',
        acceleration_source='McGaugh, Lelli, and Schombert (2016), DOI 10.1103/PhysRevLett.117.201101',
        acceleration_role='Illustrative central calibration only; uncertainties are not combined with the laboratory bound.',
        geometry_convention='|H|=eta_H/ell; eta_H=1 for averaged principal curvature, eta_H=3 for its trace.',
        ell_source='Adelberger et al. (2007), Table I, k=2; inherited 68% confidence convention.')
    (args.output / 'cross_sector_hierarchy.json').write_text(
        json.dumps(hierarchy, indent=2)+'\n', encoding='utf-8')
    print(f'Geometry and constraint checks: {len(checks)-failed}/{len(checks)} passed; {failed} failed')
    return int(failed != 0)


if __name__ == '__main__':
    raise SystemExit(main())
