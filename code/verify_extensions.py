#!/usr/bin/env python3
"""Numerical checks for the additional conditional results in the manuscript.

Run alongside the unchanged baseline: python code/reproduce.py --output results
Then: python code/verify_extensions.py --output results
No data download or fitting is performed. Numerical checks do not certify proofs.
Copyright 2026 Henry Arellano-Pena. CC BY 4.0 (see the accompanying license).
"""
from __future__ import annotations

import argparse
import cmath
import json
import math
import platform
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy import integrate, special

import reproduce as base


class Checks:
    """Record each finite numerical comparison, its tolerance and its outcome."""
    def __init__(self) -> None:
        self.rows: list[dict] = []

    def close(self, name: str, actual: float, expected: float,
              rtol: float = 2e-8, atol: float = 1e-11) -> None:
        actual, expected = float(actual), float(expected)
        ok = math.isfinite(actual) and math.isfinite(expected) and math.isclose(
            actual, expected, rel_tol=rtol, abs_tol=atol)
        self.rows.append(dict(name=name, actual=actual, expected=expected,
                              rtol=rtol, atol=atol, passed=ok))


def aqual_flux(p: np.ndarray, a_star: float = 1.0) -> np.ndarray:
    """Flux for the illustrative mu(y)=y/sqrt(1+y^2), not a derived ECL."""
    if a_star <= 0:
        raise ValueError('a_star must be positive')
    y = np.linalg.norm(p) / a_star
    return (y / math.hypot(1.0, y)) * p


def aqual_tensor(p: np.ndarray, a_star: float = 1.0) -> np.ndarray:
    g = float(np.linalg.norm(p))
    if g <= 0 or a_star <= 0:
        raise ValueError('nonzero background and positive a_star required')
    y, direction = g / a_star, p / g
    mu, derivative = y / math.hypot(1.0, y), (1.0 + y*y)**(-1.5)
    return mu*np.eye(3) + y*derivative*np.outer(direction, direction)


def run_checks() -> list[dict]:
    c = Checks()
    # Weighted scalar action -> canonical inverse-square form, ell not fixed to 1.
    for n in (3, 4, 6):
        for gamma in (0.75, 1.0, 1.25):
            s = gamma*(n-1)/2
            nu = s + 0.5
            for ell in (0.4, 2.0):
                key = f'n{n}_gamma{gamma}_ell{ell}'
                b = 0.8/ell
                f2 = lambda z: math.exp(-2*b*z)
                fac = integrate.quad(lambda z: f2(z)*(-b+s/(z+ell))**2,
                                     0, np.inf, epsabs=1e-11)[0]
                direct = integrate.quad(
                    lambda z: f2(z)*(b*b+s*(s+1)/(z+ell)**2),
                    0, np.inf, epsabs=1e-11)[0] - s/ell
                c.close('warped_form_'+key, direct, fac)
                c.close('inverse_square_coefficient_'+key, s*(s+1), nu*nu-.25)
                z0 = (2*s-1)/ell
                norm = integrate.quad(lambda z: z0*(1+z/ell)**(-2*s),
                                      0, np.inf, epsabs=1e-11)[0]
                c.close('warped_zero_norm_'+key, norm, 1.0)
                z = 0.37*ell
                weighted = (1+z/ell)**(-2*s) * (
                    (1+z/ell)**s*math.exp(-b*z)*(-b+s/(z+ell)))**2
                c.close('weighted_canonical_integrand_'+key,
                        weighted, f2(z)*(-b+s/(z+ell))**2)
    for ell in (0.3, 1.0, 4.0):
        for q in (0.05, 0.7, 6.0):
            c.close(f'candidate_kernel_{ell}_{q}', ell*ell*base.kernel(1.5,q*ell),
                    1/q**2+ell/q, rtol=2e-12)
        for r in (0.2, 3.0, 80.0):
            x = r/ell
            e = base.response(1.5,x,True)[0]
            c.close(f'candidate_force_{ell}_{r}', e, 4*ell/(math.pi*r))
            c.close(f'candidate_ratio_{ell}_{r}',
                    base.response(1.5,3*x,True)[0]/e, 1/3)
        r0, epsilon = 2.3*ell, 0.004
        ell_bound = math.pi*epsilon*r0/4
        c.close(f'candidate_bound_{ell}', 4*ell_bound/(math.pi*r0), epsilon)
    # Finite differences computed at 60 digits; alternating sums can cancel.
    with mp.workdps(60):
        atoms = [(mp.mpf('0.4'),mp.mpf('0.7')),
                 (mp.mpf('2.1'),mp.mpf('0.3'))]
        density = 2/mp.pi
        def d(r):
            return density/r + sum(w*mp.exp(-m*r) for m,w in atoms)
        for r,h in ((mp.mpf('.3'),mp.mpf('.1')),
                    (mp.mpf('3'),mp.mpf('.7'))):
            for order in range(9):
                lhs = sum((-1)**j*mp.binomial(order,j)*d(r+j*h)
                          for j in range(order+1))
                rhs = mp.quad(lambda m: density*mp.exp(-m*r)*
                              (1-mp.exp(-m*h))**order, [0,1,mp.inf])
                rhs += sum(w*mp.exp(-m*r)*(1-mp.exp(-m*h))**order
                           for m,w in atoms)
                c.close(f'finite_difference_{r}_{h}_{order}',lhs,rhs,
                        rtol=1e-12,atol=1e-14)
        for r in (0.2,1.0,5.0):
            # A discrete measure has a finite resolvent without having a density.
            c.close(f'atomic_integrability_{r}',
                    integrate.quad(lambda t:t*math.exp(-0.7*t),0,r)[0],
                    (1-(1+0.7*r)*math.exp(-0.7*r))/0.7**2)
            # D(r)=r^-2 is CM but rho(mu)=mu fails unsubtracted resolvent integrability.
            c.close(f'nonintegrable_density_laplace_{r}',
                    integrate.quad(lambda m:m*math.exp(-m*r),0,np.inf)[0],r**-2)
        for cutoff in (1.,10.,1000.):
            c.close(f'log_resolvent_cutoff_{cutoff}',
                    integrate.quad(lambda m:m/(1+m*m),0,cutoff)[0],
                    0.5*math.log1p(cutoff*cutoff))
        for rstar in (0.3,2.0):
            # Tonelli criterion for constant density: integral r*D(r) = c*rstar.
            # Lower incomplete gamma avoids cancellation of 1-(1+x)*exp(-x).
            k = lambda m: mp.gammainc(2, 0, m*rstar)/m**2 if m else rstar*rstar/2
            spectral = mp.quad(lambda m:density*k(m),[0,1,mp.inf])
            c.close(f'tonelli_integrability_{rstar}',spectral,density*rstar,rtol=1e-12)
    # Exact continued resolvent in the solvable case, plus positive detuning threshold.
    for ell in (0.5,2.0):
        for delta in (1e-4,0.05,0.5):
            for q in (0.02,0.7,4.0):
                g = 1/q+1/(ell*q*q)
                exact = (q+1/ell)/(q*q+delta*q+delta/ell)
                c.close(f'detuned_exact_{ell}_{delta}_{q}',g/(1+delta*g),exact,rtol=1e-12)
            pole = cmath.sqrt(delta/ell-delta*delta/4)-1j*delta/2
            qpole = -1j*pole
            residual = qpole*qpole+delta*qpole+delta/ell
            c.close(f'resonance_pole_{ell}_{delta}',abs(residual),0,atol=1e-14)
            ratio = abs(pole.imag)/pole.real
            c.close(f'resonance_width_{ell}_{delta}',ratio,
                    0.5*math.sqrt(delta*ell)/math.sqrt(1-delta*ell/4),rtol=1e-12)
    for nu in (1.5,2.0,2.5,3.0):
        mu, delta, z0 = 1e-6, 0.03, 2*(nu-1)
        # K_nu(-i mu)/((-i mu)K_{nu-1}(-i mu)) = -H_nu^(1)/(mu H_{nu-1}^(1)).
        gplus = -special.hankel1(nu,mu)/(mu*special.hankel1(nu-1,mu))
        rho = base.weight(nu,mu)
        measured = rho/abs(1+delta*gplus)**2
        expected = rho*mu**4/(delta*delta*z0*z0)
        c.close(f'detuned_threshold_{nu}',measured,expected,rtol=2e-8,atol=0)
    # Analytic flux Jacobian checked by central differences on nonzero backgrounds.
    direction = np.array([1.,2.,3.])/math.sqrt(14)
    for y in (0.01,0.1,1.,10.,100.):
        p = y*direction
        tensor = aqual_tensor(p)
        h = max(1e-8,y*1e-5)
        numerical = np.column_stack([
            (aqual_flux(p+h*np.eye(3)[j])-aqual_flux(p-h*np.eye(3)[j]))/(2*h)
            for j in range(3)])
        for i in range(3):
            for j in range(3):
                c.close(f'aqual_jacobian_{y}_{i}_{j}',tensor[i,j],numerical[i,j],
                        rtol=4e-6,atol=2e-10)
        expected = [y/math.hypot(1,y)]*2 + [y/math.hypot(1,y)+y*(1+y*y)**(-1.5)]
        for j, ev in enumerate(np.linalg.eigvalsh(tensor)):
            c.close(f'aqual_eigenvalue_{y}_{j}',ev,expected[j],rtol=1e-12)
    return c.rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,
                        default=Path(__file__).resolve().parents[1]/'results')
    args = parser.parse_args()
    rows = run_checks()
    failed = sum(not row['passed'] for row in rows)
    report = dict(passed=len(rows)-failed, failed=failed,
                  python=platform.python_version(), numpy=np.__version__,
                  scipy=scipy.__version__, mpmath=mp.__version__, checks=rows,
                  note='Additional numerical consistency tests; not experimental data or proof certification.')
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'extension_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('passed','failed')},indent=2))
    if failed:
        raise SystemExit('Failed checks: '+', '.join(r['name'] for r in rows if not r['passed']))


if __name__ == '__main__':
    main()
