#!/usr/bin/env python3
"""Verify the specified AQUAL--QUMOND matched-gradient discriminator.

The real-space squared axis-force contrast is not a Fourier-axis susceptibility.
All values are model consequences; this program fits no astronomical catalogue.
The source equations are pinned to Milgrom (2010), arXiv:0911.5464v2, and
Banik & Zhao, arXiv:1509.08457v3, equations (18), (24), and (36)--(38).
Run from any directory; --output selects the JSON destination.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from typing import Any
import mpmath as mp
import sympy as sp


def main(output: Path) -> None:
    mp.mp.dps = 60
    checks: dict[str, bool] = {}

    def check(name: str, condition: bool) -> None:
        checks[name] = bool(condition)
        if not condition:
            raise AssertionError(name)

    def symbolic(name: str, expr: sp.Expr) -> None:
        check(name, sp.simplify(expr) == 0)

    def close(name: str, a: Any, b: Any, tol: str = '1e-48') -> None:
        aa, bb = mp.mpf(a), mp.mpf(b)
        check(name, abs(aa-bb) <= mp.mpf(tol)*max(mp.mpf(1), abs(bb)))

    ell, eta, t = sp.symbols('ell eta t', real=True)
    K = -ell/(1+ell)
    RA = 1+ell
    RQ = (1+K/2)**-2
    slope = (2-eta)/(1+ell)
    symbolic('QUMOND axis response', RQ-4*(1+ell)**2/(2+ell)**2)
    symbolic('AQUAL residual', slope*RA-(2-eta))
    symbolic('QUMOND residual', slope*RQ-(2-eta)+(2-eta)*ell**2/(2+ell)**2)
    symbolic('relative response deficit', RQ/RA-(1-ell**2/(2+ell)**2))
    symbolic('Newton null with general mass profile', (slope*RQ-(2-eta)).subs(ell, 0))
    symbolic('zero radial slope blind case', (slope*RQ-(2-eta)).subs(eta, 2))
    symbolic('deep limit axis response', RQ.subs(ell, 1)-sp.Rational(16,9))
    symbolic('deep limit residual', (slope*RQ-(2-eta)).subs({ell:1, eta:0})+sp.Rational(2,9))
    symbolic('Fourier transfer ratio', (1+K*t)*(1+ell*t)-1-ell**2*t*(1-t)/(1+ell))
    symbolic('parallel Fourier equality', (1+K)*(1+ell)-1)
    symbolic('transverse Fourier equality', ((1+K*t)*(1+ell*t)-1).subs(t, 0))
    symbolic('oblique Fourier difference', ((1+K*t)*(1+ell*t)-1).subs(t,sp.Rational(1,2))-ell**2/(4*(1+ell)))

    y = sp.symbols('y', positive=True)
    mu = y/sp.sqrt(1+y*y)
    log_slope = sp.simplify(y*sp.diff(mu,y)/mu)
    xi = y*mu
    # Chain rule d ln N / d ln xi, with N(xi(y))=1/mu(y).
    K_mu = sp.diff(sp.log(1/mu), y)/sp.diff(sp.log(xi), y)
    symbolic('paired interpolation derivative', K_mu+log_slope/(1+log_slope))
    symbolic('mu0 logarithmic slope', log_slope-1/(1+y*y))
    deriv_R = sp.simplify(y*sp.diff(1+log_slope,y)/(1+log_slope))
    symbolic('field matching derivative', deriv_R+2*y*y/((1+y*y)*(2+y*y)))

    r, theta, kk = sp.symbols('r theta K', positive=True)
    # Set G*m*N=1. The radial and angular operators verify the PDE away from 0.
    phi = -(1+kk*sp.sin(theta)**2/2)/r
    lap_phi = (sp.diff(r*r*sp.diff(phi,r),r)/r**2
               +sp.diff(sp.sin(theta)*sp.diff(phi,theta),theta)/(r*r*sp.sin(theta)))
    symbolic('QUMOND exterior PDE', sp.trigsimp(lap_phi-kk*(1-3*sp.cos(theta)**2)/r**3))
    symbolic('QUMOND radial force', sp.diff(phi,r)-(1+kk*sp.sin(theta)**2/2)/r**2)
    symbolic('QUMOND angular force', -sp.diff(phi,theta)/r-kk*sp.sin(theta)*sp.cos(theta)/r**2)
    flux = sp.integrate((1+kk*sp.sin(theta)**2/2)*sp.sin(theta),(theta,0,sp.pi))*2*sp.pi
    symbolic('distributional source flux', flux-4*sp.pi*(1+kk/3))

    u = sp.symbols('u', real=True)
    fA = (1+ell*u)**(-sp.Rational(1,2))
    fQ = 1-ell*u/(2*(1+ell))
    symbolic('angular difference starts at second order', sp.series(fQ-fA,ell,0,3).removeO()-ell**2*(u/2-3*u*u/8))
    symbolic('residual starts at second order', sp.series(slope*RQ-(2-eta),ell,0,3).removeO()+(2-eta)*ell**2/4)

    rows = []
    for text in ['0.1','1','1.8','10']:
        yy = mp.mpf(text)
        mm = yy/mp.sqrt(1+yy*yy)
        ll = 1/(1+yy*yy)
        nn = 2/(1+ll)
        kg = -ll/(1+ll)
        ra, rq = 1+ll, (1+kg/2)**-2
        aa, aq = mp.sqrt(ra), mp.sqrt(rq)
        tq = nn*rq-2
        close(f'AQUAL product y={text}', nn*ra, 2)
        close(f'QUMOND residual y={text}', tq, -2*ll*ll/(2+ll)**2)
        close(f'axis contrast from potential y={text}', aq, 1/(1+kg/2))
        close(f'angular monopole integral y={text}', mp.quad(lambda z: (1+kg*mp.sin(z)**2/2)*mp.sin(z),[0,mp.pi]), 2*(1+kg/3))
        close(f'Fourier axial transfer y={text}', (1+kg)*(1+ll), 1)
        check(f'Fourier oblique separation y={text}', (1+kg/2)*(1+ll/2)>1)
        check(f'QUMOND axis contrast below AQUAL y={text}', aq<aa)
        # Exact positive-root inversion for this mu0, evaluated independently.
        xx = yy*mm
        yy_back = mp.sqrt((xx*xx+xx*mp.sqrt(xx*xx+4))/2)
        close(f'spherical inverse relation y={text}', yy_back, yy)
        rows.append({k:mp.nstr(v,30) for k,v in {
            'y':yy,'mu':mm,'ell':ll,'nu_exterior':nn,'R_A':ra,'R_Q':rq,
            'axis_A':aa,'axis_Q':aq,'T_Q_exterior':tq}.items()})

    row = rows[2]
    nn, rq, ra = (mp.mpf(row[k]) for k in ['nu_exterior','R_Q','R_A'])
    # Re-evaluate rather than truncate to the serialized decimal precision.
    yy = mp.mpf('1.8'); ll = 1/(1+yy*yy)
    nn = 2/(1+ll); ra=1+ll; rq=4*(1+ll)**2/(2+ll)**2
    tq = nn*rq-2
    close('local exact ell', ll, mp.mpf(25)/106)
    close('local exact spherical slope', nn, mp.mpf(212)/131)
    close('local exact QUMOND residual', tq, -mp.mpf(1250)/56169)
    close('printed local residual', tq, '-0.02225427','5e-9')
    close('printed local AQUAL response', ra, '1.23584906','5e-9')
    close('printed local QUMOND response', rq, '1.22209760','5e-9')
    signal = abs(tq)
    sigT = signal/3
    sigRfrac = sigT/(nn*rq)
    sigAxisfrac = sigRfrac/2
    sigBothfrac = sigRfrac/mp.sqrt(2)
    dlnR = -2*yy*yy/((1+yy*yy)*(2+yy*yy))
    mismatch_exact = 2*((1+1/(1+(yy*mp.exp(mp.mpf('.01')))**2)))/ra-2
    mismatch_linear = 2*dlnR*mp.mpf('.01')
    check('one percent mismatch linear approximation', abs(mismatch_exact-mismatch_linear)<mp.mpf('1e-4'))
    close('quoted three sigma residual precision', sigT,'0.00742','3e-6')
    close('quoted fractional R precision',sigRfrac,'0.00375','1e-6')
    close('quoted fractional axis precision',sigAxisfrac,'0.00188','5e-6')
    close('quoted equal fractional precisions',sigBothfrac,'0.00265','3e-6')

    # Check the direction of the published absolute-value and positive limits.
    for convention, lower, upper in [('absolute_95','2.08','1.45'),('positive_delta_chi2_6.17','2.56','1.77')]:
        check(f'Lee alpha=2 not excluded at 32 um: {convention}',mp.mpf(2)<mp.mpf(lower))
        check(f'Lee alpha=2 excluded at 35 um: {convention}',mp.mpf(2)>mp.mpf(upper))
    close('Lee lower edge rounded micrometres',32*mp.pi,'100.5','5e-4')
    close('Lee upper edge rounded micrometres',35*mp.pi,'110.0','5e-4')

    sensitivity = {k:mp.nstr(v,30) for k,v in {
        'relative_R_separation':1-rq/ra,
        'relative_axis_separation':1-mp.sqrt(rq/ra),
        'sigma_T_for_nominal_3_sigma':sigT,
        'fractional_sigma_R_if_slope_exact':sigRfrac,
        'fractional_sigma_axis_if_slope_exact':sigAxisfrac,
        'equal_independent_fractional_sigma_nu_R':sigBothfrac,
        'dlogR_dlogg':dlnR,'one_percent_mismatch_exact_T':mismatch_exact,
        'one_percent_mismatch_linear_T':mismatch_linear}.items()}
    payload = {
        'scope':'Conditional field-theory calculations; no observational fit or quantum-loop recomputation.',
        'statistic':'R_WB = (absolute internal radial force parallel / transverse)^2, same separation and mass.',
        'background_assumption':'Paired collinear constant fields g_N,b = mu(g_b/a_*) g_b; not assumed for every galaxy.',
        'references':['arXiv:0911.5464v2','arXiv:1509.08457v3','Lee 2020 dissertation Table 11.4, printed page 134'],
        'n_checks':len(checks),'all_passed':all(checks.values()),'checks':checks,
        'table_6':rows,'sensitivity_y_1_8':sensitivity,
        'lee_conventions':{
            'absolute_95':{'lambda_um':[32,35],'alpha_limit':['2.08','1.45'],
                'definition':'Integral of fitted normal density from -a to +a is 0.95.'},
            'positive_delta_chi2_6.17':{'lambda_um':[32,35],'alpha_limit':['2.56','1.77']}}
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(f'{len(checks)}/{len(checks)} checks passed; {output}')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('model_discrimination_results.json'))
    main(parser.parse_args().output)
