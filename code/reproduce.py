#!/usr/bin/env python3
"""Reproduce dimensionless spectral-response tables and analytical consistency tests.

All examples are mathematical, not fits to experimental data.  The inverse-square
half-line family is a scalar constitutive comparison model, not a tensor theory.
Run: python code/reproduce.py --output results
Dependencies: numpy, scipy, mpmath (Python 3.10+).
"""
from __future__ import annotations
import argparse, csv, json, math, platform
from pathlib import Path
import numpy as np
import scipy
from scipy import integrate, special, optimize
import mpmath as mp


def _validate(nu: float, x: float) -> None:
    if not math.isfinite(nu) or nu <= 1 or not math.isfinite(x) or x <= 0:
        raise ValueError('nu must exceed 1 and the dimensionless argument must be positive')


def weight(nu: float, u: float) -> float:
    """W_nu(u) = rho_nu(mu)/ell, u = mu*ell."""
    _validate(nu, u)
    a = nu - 1.0
    if u < 1e-12:
        return math.exp((1-2*a)*math.log(2)-math.log(a)-2*special.gammaln(a)+(2*a-1)*math.log(u))
    j, y = special.jv(a,u), special.yv(a,u)
    return 2.0/(math.pi**2*a*u*(j*j+y*y))


def coefficients(nu: float) -> tuple[float,float]:
    _validate(nu,1.0)
    p=2*nu-2
    c=math.exp((3-2*nu)*math.log(2)+special.gammaln(p)-math.log(nu-1)-2*special.gammaln(nu-1))
    return p,c


def response(nu: float, x: float, force: bool=False) -> tuple[float,float]:
    """Return potential correction D or force excess E at x=r/ell."""
    _validate(nu,x)
    def f(t: float) -> float:
        if t == 0: return 0.0  # isolated endpoint value does not affect integration
        return weight(nu,t/x)*math.exp(-t)*(1+t if force else 1)/x
    val,err=integrate.quad(f,0,np.inf,epsabs=2e-12,epsrel=3e-11,limit=400)
    return float(val),float(err)


def kernel(nu: float, s: float) -> float:
    _validate(nu,s)
    # Exponentially scaled K functions avoid underflow at large s.
    return special.kve(nu,s)/(2*(nu-1)*s*special.kve(nu-1,s))


def continuum_kernel(nu: float, s: float) -> float:
    _validate(nu,s)
    # Recurrence isolates the pole without subtracting two nearly equal terms.
    return special.kve(nu-2,s)/(2*(nu-1)*s*special.kve(nu-1,s))


def spectral_kernel(nu: float, s: float) -> float:
    _validate(nu,s)
    return integrate.quad(lambda u: weight(nu,u)/(s*s+u*u),0,np.inf,
                          epsabs=1e-11,epsrel=3e-10,limit=400)[0]


def detuning_root(nu: float, d: float) -> float:
    """Dimensionless negative eigenvalue is -s**2 for d=delta*ell<0."""
    _validate(nu,1.0)
    if d >= 0: raise ValueError('Negative detuning required')
    def f(s: float) -> float:
        return s*special.kve(nu-1,s)/special.kve(nu,s)+d
    hi=max(1,abs(d)*2+nu)
    while f(hi)<0: hi*=2
    return float(optimize.brentq(f,1e-14,hi,xtol=2e-14))


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'results')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    checks=[]
    def check(name: str, actual: float, expected: float, rtol: float=1e-8, atol: float=1e-11):
        ok=math.isclose(actual,expected,rel_tol=rtol,abs_tol=atol)
        checks.append(dict(name=name,actual=float(actual),expected=float(expected),rtol=rtol,atol=atol,passed=ok))
        if not ok: raise AssertionError(f'{name}: {actual} versus {expected}')
    nus=[1.25,1.5,2.,2.5,3.]
    for nu in nus:
        a=nu-.5
        z0=2*(nu-1)
        norm=integrate.quad(lambda z:z0*(1+z)**(1-2*nu),0,np.inf)[0]
        check(f'zero_mode_norm_{nu}',norm,1)
        # A test function exp(-z); beta_*=-a for ell=1.
        direct=integrate.quad(lambda z:math.exp(-2*z)*(1+(nu*nu-.25)/(1+z)**2),0,np.inf)[0]-a
        fact=integrate.quad(lambda z:math.exp(-2*z)*(-1+a/(1+z))**2,0,np.inf)[0]
        check(f'form_factorization_{nu}',direct,fact)
        for s in [.03,.3,1.,3.,30.]:
            check(f'spectral_resolvent_{nu}_{s}',spectral_kernel(nu,s),continuum_kernel(nu,s),rtol=2e-8)
            check(f'bessel_recurrence_{nu}_{s}',kernel(nu,s),1/s**2+continuum_kernel(nu,s),rtol=1e-12)
        p,c=coefficients(nu)
        for x in [.3,3,30]:
            D=response(nu,x)[0]; E=response(nu,x,True)[0]
            h=x*1e-4
            derivative=(response(nu,x+h)[0]-response(nu,x-h)[0])/(2*h)
            check(f'force_derivative_{nu}_{x}',E,D-x*derivative,rtol=2e-7)
        previous=float('inf')
        for x in np.logspace(-2,3,18):
            E=response(nu,float(x),True)[0]
            if E <= 0 or E >= previous: raise AssertionError('Force monotonicity failed')
            checks.append(dict(name=f'force_monotonicity_{nu}_{x:.6g}',passed=True))
            previous=E
        s=detuning_root(nu,-1e-7)
        check(f'small_detuning_{nu}',s*s,2*(nu-1)*1e-7,rtol=.03,atol=0)
    for x in [.01,.1,1,10,100]:
        check(f'exact_nu_1p5_potential_{x}',response(1.5,x)[0],2/(math.pi*x))
        check(f'exact_nu_1p5_force_{x}',response(1.5,x,True)[0],4/(math.pi*x))
    # Independent 40-digit quadrature on the original u integration variable.
    mp.mp.dps=40
    independent=[]
    for nu,x in [(1.5,1),(2,1),(2,10),(2.5,10),(3,1),(3,100)]:
        a=mp.mpf(str(nu-1));xx=mp.mpf(str(x))
        def mf(u):
            return 2*mp.exp(-xx*u)/(mp.pi**2*a*u*(mp.besselj(a,u)**2+mp.bessely(a,u)**2))
        result=mp.quad(mf,[0,1/xx,10/xx,mp.inf])
        actual=response(nu,x)[0]
        check(f'mpmath_crosscheck_{nu}_{x}',actual,float(result),rtol=2e-8,atol=1e-12)
        independent.append(dict(nu=nu,x=x,mpmath=str(result),scipy=actual))
    rows=[]
    for nu in [1.5,2.,2.5,3.]:
        p,c=coefficients(nu)
        for x in [1.,10.,100.]:
            D,err=response(nu,x);E,eerr=response(nu,x,True)
            rows.append(dict(nu=nu,x=x,D=D,force_excess=E,asymptotic_D=c*x**(-p),
                             asymptotic_ratio=D/(c*x**(-p)),quad_error=err,force_quad_error=eerr))
    with (args.output/'response_table.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    tex=[]
    for row in rows:
        tex.append(f"{row['nu']:g} & {row['x']:g} & {row['D']:.8g} & {row['force_excess']:.8g} & {row['asymptotic_ratio']:.6f} \\\\")
    (args.output/'response_rows.tex').write_text('\n'.join(tex)+'\n')
    report=dict(passed=len(checks),failed=0,python=platform.python_version(),numpy=np.__version__,
                scipy=scipy.__version__,mpmath=mp.__version__,checks=checks,independent_quadrature=independent,
                note='Mathematical consistency tests; not experimental observations or proof certification.')
    (args.output/'verification.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'passed':len(checks),'failed':0,'rows':len(rows)},indent=2))
    print((args.output/'response_rows.tex').read_text())

if __name__=='__main__': main()
