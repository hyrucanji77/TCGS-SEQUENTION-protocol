#!/usr/bin/env python3
"""Reproduce the additional algebraic and spectral checks in the article.

No experiment is fitted and no quantum loop integral is recomputed. The Lee
2020 entries are transcribed historical two-sided 95% confidence limits, not a
likelihood or an exclusion bound on the complete slab tower.
"""
from __future__ import annotations
import argparse
import json
from itertools import combinations
from pathlib import Path
import mpmath as mp
import sympy as sp

mp.mp.dps=50
checks: dict[str, bool]={}

def check(name: str, ok: bool) -> None:
    checks[name]=bool(ok)
    if not ok:
        raise AssertionError(name)

def symbolic(name: str, expr: sp.Expr) -> None:
    check(name, sp.simplify(expr)==0)

def close(name: str, observed: mp.mpf, expected: mp.mpf, tol: str='1e-38') -> None:
    check(name, abs(observed-expected) <= mp.mpf(tol)*max(mp.mpf(1),abs(expected)))

def constitutive_rank(symmetric: bool) -> int:
    """Rank of chi -> B(k)^T chi B(k), for all quadratic monomials in k.

    F_mu,nu = k_mu A_nu - k_nu A_mu. Evaluating at e_i and e_i+e_j
    spans all quadratic monomials, so this is an exact integer rank test.
    A symmetric diagonal internal block has 21 coefficients; an off-diagonal
    internal block has 36, its transpose residing in the opposite block.
    """
    pairs=list(combinations(range(4),2))
    kvecs=[sp.eye(4)[:,i] for i in range(4)]
    kvecs += [sp.eye(4)[:,i]+sp.eye(4)[:,j] for i,j in pairs]
    bases=[]
    for i in range(6):
        for j in range(i if symmetric else 0,6):
            E=sp.zeros(6);E[i,j]=1
            if symmetric and i!=j:E[j,i]=1
            bases.append(E)
    columns=[[] for _ in bases]
    for k in kvecs:
        B=sp.zeros(6,4)
        for a,(mu,nu) in enumerate(pairs):
            B[a,nu]=k[mu];B[a,mu]=-k[nu]
        for col,E in zip(columns,bases):
            col.extend(list(B.T*E*B))
    return sp.Matrix.hstack(*(sp.Matrix(col) for col in columns)).rank()


def main(output: Path) -> None:
    rank_sym=constitutive_rank(True);rank_cross=constitutive_rank(False)
    check('single internal diagonal block rank = 20',rank_sym==20)
    check('single internal off-diagonal block rank = 35',rank_cross==35)
    check('four-field dynamic rank = 290',4*rank_sym+6*rank_cross==290)
    check('four-field constant null sector = 10',300-(4*rank_sym+6*rank_cross)==10)
    for rv in ['0.25','0.5','1','2','4']:
        r=mp.mpf(rv);L=mp.mpf(1);z=mp.pi*r/2
        images=1/r**2+2*mp.nsum(lambda j:1/(r*r+(2*L*j)**2),[1,mp.inf])
        close(f'compact image sum r/L={rv}',images,mp.pi*mp.coth(z)/(2*L*r))
        half=2*mp.nsum(lambda n:mp.exp(-(n+mp.mpf('.5'))*mp.pi*r/L),[0,mp.inf])
        close(f'Dirichlet half-integer potential r/L={rv}',half,mp.csch(z))
        tail=2*mp.nsum(lambda n:mp.exp(-n*mp.pi*r/L),[2,mp.inf])
        first=2*mp.exp(-mp.pi*r/L)
        close(f'first-Yukawa remainder ratio r/L={rv}',tail/first,mp.exp(-mp.pi*r/L)/(1-mp.exp(-mp.pi*r/L)))
    for kv in ['0.1','1','4']:
        k=mp.mpf(kv);L=mp.mpf('1.3')
        spec=2/L*mp.nsum(lambda n:1/(k*k+((n+mp.mpf('.5'))*mp.pi/L)**2),[0,mp.inf])
        close(f'Dirichlet resolvent k={kv}',spec,mp.tanh(k*L)/k)
    r,L,z,w=sp.symbols('r L z w',positive=True)
    PhiD=-sp.csch(sp.pi*r/(2*L))/r
    Z=sp.pi*r/(2*L)
    symbolic('Dirichlet force derivative',sp.diff(PhiD,r)-(sp.csch(Z)+Z*sp.csch(Z)*sp.coth(Z))/r**2)
    symbolic('common ultraviolet pole coefficient',sp.limit(z*sp.csch(z),z,0)-sp.limit(z*sp.coth(z),z,0))
    symbolic('Dirichlet exponential infrared coefficient',sp.limit(sp.exp(z)*sp.csch(z),z,sp.oo)-2)
    q=sp.symbols('q',positive=True)
    tangent=sp.cosh(q*w)/sp.cosh(q*L)
    symbolic('anisotropic extension equation',sp.diff(tangent,w,2)-q*q*tangent)
    symbolic('anisotropic reflecting face',sp.diff(tangent,w).subs(w,0))
    symbolic('anisotropic upper trace',tangent.subs(w,L)-1)
    symbolic('anisotropic upper flux',sp.diff(tangent,w).subs(w,L)-q*sp.tanh(q*L))
    symbolic('normalized tangent kernel L2 coefficient',sp.limit((sp.tanh(z)/z-1)/z**2,z,0)+sp.Rational(1,3))
    symbolic('normalized tangent kernel L4 coefficient',sp.limit((sp.tanh(z)/z-1+z**2/3)/z**4,z,0)-sp.Rational(2,15))
    symbolic('thin-slab tangential coefficient',sp.integrate((L*L-w*w)/2,(w,0,L))-L**3/3)
    symbolic('thin-slab normal coefficient',sp.integrate(w*w/2,(w,0,L))-L**3/6)
    symbolic('thin-slab net energy coefficient',-sp.Rational(1,3)+sp.Rational(1,6)+sp.Rational(1,6))
    for yv in ['0.1','1','10']:
        y=mp.mpf(yv);ratio=1+1/(1+y*y)
        nu=2/ratio
        close(f'matched isolated-gradient residual y={yv}',nu*ratio,mp.mpf(2))
        mass_slope=mp.mpf('.3')
        extended_nu=(2-mass_slope)/ratio
        close(f'matched extended-source residual y={yv}',extended_nu*ratio,2-mass_slope)
        mu=y/mp.sqrt(1+y*y)
        gpar=1/mu;gperp=1/(mu*mp.sqrt(ratio))
        close(f'axis-force ratio squared y={yv}',(gpar/gperp)**2,ratio)
        kpar,kperp,length=mp.mpf('.4'),mp.mpf('.7'),mp.mpf('.3')
        wave=mp.sqrt(kperp*kperp+ratio*kpar*kpar)
        # Set 4*pi*G_N=1, gamma*L=1 for this dimensionless check.
        kernel=mu/length*wave*mp.tanh(wave*length)
        close(f'slab-sector transfer y={yv}',kernel/(mu*wave*wave),mp.tanh(wave*length)/(wave*length))
    cgr,cs,cp=sp.Rational(15,4),sp.Rational(23,8),sp.Rational(11,4)
    symbolic('pinned-v6 s deficit',(cgr-cs)/cgr-sp.Rational(7,30))
    symbolic('pinned-v6 p deficit',(cgr-cp)/cgr-sp.Rational(4,15))
    symbolic('pinned-v6 polarization difference',cs-cp-sp.Rational(1,8))
    symbolic('pinned-v6 total-angle s prefactor',(cgr-cs)/4-sp.Rational(7,32))
    historical={'source':'John G. Lee, A Fourier-Bessel Test of the Gravitational Inverse-Square Law (2020), Table 11.4, printed p. 134',
                'url':'https://digital.lib.washington.edu/researchworks/bitstreams/6013936a-da00-42b7-ba0a-4dff0cb05bf8/download',
                'convention':'two-sided 95% |alpha| limit for one Yukawa, not the full tower',
                'rows':[{'lambda_um':32,'abs_alpha_95':2.08},{'lambda_um':35,'abs_alpha_95':1.45}],
                'alpha_compared':2,'L_crossing_bracket_um':[float(32*mp.pi),float(35*mp.pi)],
                'claimed_tower_exclusion':False}
    check('alpha-two crossing bracket',historical['rows'][0]['abs_alpha_95']>2>historical['rows'][1]['abs_alpha_95'])
    check('quoted L lower bracket rounding',round(float(32*mp.pi),1)==100.5)
    check('quoted L upper bracket rounding',round(float(35*mp.pi),1)==110.0)
    payload={'description':__doc__,'number_of_checks':len(checks),'all_checks_passed':all(checks.values()),
             'checks':checks,'constitutive_block_ranks':{'symmetric':rank_sym,'cross':rank_cross},
             'historical_yukawa_input':historical,
             'light_coefficients':{'version':'2505.14446v6','date':'2025-09-16','printed_page':4,
                                   'source_equations':[26,27,28],'pi_coefficients':['23/8','11/4','15/4']}}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'number_of_checks':len(checks),'all_checks_passed':all(checks.values())},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('extension_results.json'))
    main(parser.parse_args().output)
