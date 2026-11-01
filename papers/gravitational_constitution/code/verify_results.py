#!/usr/bin/env python3
"""Reproduce the analytic checks and dimensionless tables accompanying the article.

Requires Python >= 3.10, mpmath, and sympy. No experimental data are used.
This program does not recompute the loop integrals in Partanen and Tulkki.
Run from any directory: python verify_results.py --output results.json
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import mpmath as mp
import sympy as sp

mp.mp.dps = 60
checks: dict[str, bool] = {}

def symbolic(name: str, expression: sp.Expr) -> None:
    ok = sp.simplify(expression) == 0
    checks[name] = bool(ok)
    if not ok:
        raise AssertionError(f"{name}: {sp.simplify(expression)}")

def close(name: str, observed: mp.mpf, expected: mp.mpf, tol: str = '1e-45') -> None:
    error = abs(observed - expected)
    ok = error <= mp.mpf(tol) * max(mp.mpf(1), abs(expected))
    checks[name] = bool(ok)
    if not ok:
        raise AssertionError(f"{name}: error = {error}")

def enhancement(z: mp.mpf) -> mp.mpf:
    return mp.coth(z) + z / mp.sinh(z)**2

def bisect_root(fun, left: mp.mpf, right: mp.mpf) -> mp.mpf:
    if fun(left) * fun(right) >= 0:
        raise ValueError('Root must be bracketed')
    for _ in range(220):
        mid = (left + right) / 2
        if fun(left) * fun(mid) > 0:
            left = mid
        else:
            right = mid
    return (left + right) / 2

def main(output: Path) -> None:
    s, y, x, k, L, w, r = sp.symbols('s y x k L w r', positive=True)
    mu = y/sp.sqrt(1+y*y)
    f = sp.sqrt(s*(1+s))-sp.asinh(sp.sqrt(s))
    symbolic('AQUAL antiderivative', sp.diff(f,s)-sp.sqrt(s)/sp.sqrt(1+s))
    ell = sp.simplify(y*sp.diff(mu,y)/mu)
    symbolic('logarithmic susceptibility', ell-1/(1+y*y))
    q = (x*x+x*sp.sqrt(x*x+4))/2
    symbolic('spherical root squared', q*q-x*x*(1+q))
    symbolic('response-slope product', (2/(1+ell))*(1+ell)-2)
    u = sp.cosh(k*w)/sp.cosh(k*L)
    symbolic('slab harmonic extension', sp.diff(u,w,2)-k*k*u)
    symbolic('reflecting face', sp.diff(u,w).subs(w,0))
    symbolic('Dirichlet face', u.subs(w,L)-1)
    symbolic('Neumann response', sp.diff(u,w).subs(w,L)-k*sp.tanh(k*L))
    ud = sp.sinh(k*w)/sp.sinh(k*L)
    symbolic('alternative boundary response', sp.diff(ud,w).subs(w,L)-k/sp.tanh(k*L))
    potential = -sp.coth(sp.pi*r/(2*L))/r
    symbolic('point-source force', sp.diff(potential,r) - (sp.coth(sp.pi*r/(2*L))+(sp.pi*r/(2*L))/sp.sinh(sp.pi*r/(2*L))**2)/r**2)
    z = sp.symbols('z',positive=True)
    R = sp.coth(z)+z/sp.sinh(z)**2
    symbolic('force monotonicity', sp.diff(R,z)+2*z*sp.coth(z)/sp.sinh(z)**2)
    for kval in [mp.mpf('0.1'),mp.mpf('1'),mp.mpf('4')]:
        length=mp.mpf('1.3')
        kernel=kval*mp.tanh(kval*length)
        integrand=lambda t: (kval*mp.sinh(kval*t)/mp.cosh(kval*length))**2 + kval**2*(mp.cosh(kval*t)/mp.cosh(kval*length))**2
        close(f'slab energy identity k={kval}', mp.quad(integrand,[0,length]), kernel)
        spectral=1/kval**2 + 2*mp.nsum(lambda n:1/(kval**2+(n*mp.pi/length)**2),[1,mp.inf])
        close(f'spectral resolvent k={kval}', spectral/length, 1/kernel)
    for rv in [mp.mpf('.25'),mp.mpf('1'),mp.mpf('3')]:
        zv=mp.pi*rv/2
        image=1+2*mp.nsum(lambda n:mp.exp(-n*mp.pi*rv),[1,mp.inf])
        close(f'Yukawa sum r/L={rv}', image, mp.coth(zv))
    # Fundamental solution of an anisotropic, constant-coefficient operator.
    xx, yy, zz = sp.symbols('xx yy zz',real=True)
    a, b = sp.symbols('a b',positive=True)
    green = 1/(a*sp.sqrt(b)*sp.sqrt(xx*xx/a+yy*yy/a+zz*zz/b))
    residual = a*(sp.diff(green,xx,2)+sp.diff(green,yy,2))+b*sp.diff(green,zz,2)
    symbolic('anisotropic Green equation away from source', residual)
    # Positive symmetric matrices need not be Hessians.
    f1,f2=sp.symbols('f1 f2')
    integrability=sp.diff(1+f2,f2)-sp.diff(sp.Integer(0),f1)
    symbolic('nonintegrable susceptibility obstruction equals one',integrability-1)
    # Mostly-plus convention: one index raised in the projector identity.
    eta=sp.diag(-1,1,1,1)
    p=sp.Matrix([2,3,5,7]); pl=eta*p; p2=(pl.T*p)[0]
    P=p2*eta-p*p.T
    symbolic('transverse tensor',sum((pl.T*P)[j]**2 for j in range(4)))
    symbolic('projector composition',sum(v*v for v in (P*eta*P-p2*P)))
    cgr=sp.Rational(15,4); cs=sp.Rational(23,8); cp=sp.Rational(11,4)
    symbolic('out-of-plane 2PN deficit', (cgr-cs)/cgr-sp.Rational(7,30))
    symbolic('in-plane 2PN deficit', (cgr-cp)/cgr-sp.Rational(4,15))
    symbolic('2PN polarization splitting',cs-cp-sp.Rational(1,8))
    slab=[]
    for rv in ['0.25','0.5','1','2','4']:
        zzv=mp.pi*mp.mpf(rv)/2
        slab.append({'r_over_L':rv,'potential_ratio':float(mp.coth(zzv)), 'force_ratio':float(enhancement(zzv))})
    bounds=[]
    for tolerance in ['0.1','0.01','0.001']:
        root=bisect_root(lambda z0:enhancement(z0)-1-mp.mpf(tolerance),mp.mpf('.1'),mp.mpf(20))
        close(f'illustrative force bound eps={tolerance}',enhancement(root)-1,mp.mpf(tolerance))
        bounds.append({'epsilon':tolerance,'z_min':float(root),'L_over_r_max':float(mp.pi/(2*root))})
    aqual=[]
    for yv in ['0.1','1','10']:
        v=mp.mpf(yv); muv=v/mp.sqrt(1+v*v); lv=1/(1+v*v)
        aqual.append({'y':yv,'mu':float(muv),'longitudinal_over_transverse':float(1+lv),'radial_slope':float(2/(1+lv))})
    # Explicit convergence estimate for the first derivative correction.
    series=[]
    for xv in ['0.1','0.25','0.5']:
        v=mp.mpf(xv); exact=mp.tanh(v)/v; approx=1-v*v/3
        series.append({'kL':xv,'exact_normalized_kernel':float(exact),'quadratic_relative_error':float(abs(approx/exact-1))})
    payload={'description':'Analytic identities and dimensionless model illustrations; no empirical fit.',
             'all_checks_passed':all(checks.values()),'number_of_checks':len(checks),'checks':checks,
             'slab_response':slab,'illustrative_bounds_not_measurements':bounds,'aqual_response':aqual,'slab_derivative_expansion':series}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(payload,indent=2))

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('results.json'))
    args=parser.parse_args()
    main(args.output)
