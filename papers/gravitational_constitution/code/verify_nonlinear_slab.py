#!/usr/bin/env python3
"""Variational check of the smooth nonlinear slab in Appendix D.

This is a finite-dimensional numerical check, not an existence proof or a fit.
The trace is u(x,L)=x+sin(x)/4, a_*=1, and gamma*L=1. Energies
are averaged over x and normalized by gamma*L. The two suppressed tangential
coordinates have unit volume. Even Legendre polynomials impose both face
conditions exactly. Run from any directory with --output PATH.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from typing import Any
import numpy as np
import scipy
from numpy.polynomial.legendre import Legendre, leggauss
from scipy.optimize import minimize


def solve(length: float, nx: int, n_basis: int, n_quad: int) -> dict[str, Any]:
    """Minimize the convex energy on a fixed spectral trial space."""
    if length <= 0 or nx % 2 or min(nx, n_basis, n_quad) < 2:
        raise ValueError('Positive length, even nx, and adequate resolutions required.')
    x = 2*np.pi*np.arange(nx)/nx
    freq = np.fft.fftfreq(nx, d=1/nx)
    # The even-grid Nyquist mode has no consistent real collocation derivative.
    freq[nx//2] = 0.0
    def dx(values: np.ndarray) -> np.ndarray:
        return np.fft.ifft(1j*freq[:, None]*np.fft.fft(values, axis=0), axis=0).real
    nodes, weights = leggauss(n_quad)
    t, weights = (nodes+1)/2, weights/2
    polys = [Legendre.basis(2*j) for j in range(1, n_basis+1)]
    B = np.array([poly(t)-1 for poly in polys])
    Bp = np.array([poly.deriv()(t) for poly in polys])
    p0 = 1+0.25*np.cos(x)
    p0x = -0.25*np.sin(x)
    mu0 = p0/np.sqrt(1+p0*p0)
    ratio0 = 1+1/(1+p0*p0)
    Bphi = mu0*ratio0*p0x
    Q = Bphi/mu0
    def density(s: np.ndarray) -> np.ndarray:
        return 0.5*(np.sqrt(s*(1+s))-np.arcsinh(np.sqrt(s)))
    E0 = float(np.mean(density(p0*p0)))
    coefficient = -float(np.mean(Bphi*Bphi/mu0))/6
    def energy_gradient(flat: np.ndarray) -> tuple[float, np.ndarray]:
        coeff = flat.reshape(nx, n_basis)
        px = p0[:, None]+length**2*(dx(coeff)@B)
        normal_scaled = coeff@Bp
        pw = length*normal_scaled
        ss = px*px+pw*pw
        mu = np.sqrt(ss/(1+ss))
        energy = float(np.mean(density(ss)@weights))
        # Gradient of (E-E0)/L^2; tangential integration by parts is spectral.
        grad = ((-dx(mu*px)*weights)@B.T
                +(mu*normal_scaled*weights)@Bp.T)/nx
        return (energy-E0)/length**2, grad.ravel()
    c0 = np.zeros((nx, n_basis))
    c0[:, 0] = -Q/3  # (1-t^2)Q/2 = -(P_2(t)-1)Q/3.
    result = minimize(energy_gradient, c0.ravel(), jac=True, method='L-BFGS-B',
                      options={'maxiter': 3000, 'ftol': 1e-15, 'gtol': 2e-11,
                               'maxls': 50, 'maxcor': 25})
    value, grad = energy_gradient(result.x)
    coeff = result.x.reshape(nx, n_basis)
    px = p0[:, None]+length**2*(dx(coeff)@B)
    pw = length*(coeff@Bp)
    energy = E0+length**2*value
    approx = E0+length**2*coefficient
    gradient_norm = float(np.max(np.abs(grad)))
    if not np.isfinite(energy) or gradient_norm > 2e-7:
        raise RuntimeError(f'Insufficient minimization: L={length}, grad={gradient_norm}, {result.message}')
    return {'L': length, 'nx': nx, 'normal_basis': n_basis, 'normal_quadrature': n_quad,
            'energy': energy, 'leading_energy': E0, 'L2_coefficient': coefficient,
            'L2_approximation': approx, 'absolute_remainder': abs(energy-approx),
            'remainder_over_L4': (energy-approx)/length**4,
            'max_scaled_gradient': gradient_norm,
            'minimum_gradient_magnitude': float(np.min(np.sqrt(px*px+pw*pw))),
            'optimizer_success': bool(result.success), 'optimizer_message': str(result.message),
            'iterations': int(result.nit)}


def main(output: Path) -> None:
    lengths = [0.4, 0.2, 0.1, 0.05]
    low = [solve(L, 48, 4, 20) for L in lengths]
    high = [solve(L, 96, 6, 32) for L in lengths]
    checks: dict[str, bool] = {}
    for a, b in zip(low, high):
        L = a['L']
        checks[f'energy_resolution_L={L}'] = abs(a['energy']-b['energy']) < 2e-10
        checks[f'regular_branch_L={L}'] = b['minimum_gradient_magnitude'] > 0.7
        checks[f'energy_below_constant_extension_L={L}'] = b['energy'] < b['leading_energy']
        checks[f'energy_above_L2_truncation_L={L}'] = b['energy'] > b['L2_approximation']
        checks[f'optimizer_stationarity_L={L}'] = b['max_scaled_gradient'] < 2e-7
    orders=[]
    for a,b in zip(high[:-1],high[1:]):
        order=float(np.log(a['absolute_remainder']/b['absolute_remainder'])/np.log(2))
        orders.append({'L_from':a['L'],'L_to':b['L'],'observed_order':order})
        checks[f'fourth_order_remainder_L={a["L"]}->{b["L"]}'] = 3.7 < order < 4.1
    report={'description':__doc__, 'normalization':'mean_x integral_0^1 F/2 dt; gamma L = a_* = 1',
            'environment':{'numpy':np.__version__,'scipy':scipy.__version__},
            'low_resolution':low,'high_resolution':high,'convergence_orders':orders,
            'checks':checks,'number_of_checks':len(checks),'all_checks_passed':all(checks.values())}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'number_of_checks':len(checks),'all_checks_passed':all(checks.values()),
                      'convergence_orders':orders,'high_resolution':high},indent=2))
    if not all(checks.values()):
        raise AssertionError([name for name, ok in checks.items() if not ok])

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('nonlinear_slab_results.json'))
    main(parser.parse_args().output)
