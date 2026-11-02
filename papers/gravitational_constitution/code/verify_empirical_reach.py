#!/usr/bin/env python3
"""Check scale separation and conditional observing-window calculations.

No likelihood, astronomical catalogue, or nonlinear stellar two-body solution is
fitted. The 110 micrometre width is a historical single-Yukawa benchmark, not a
new full-tower exclusion. C_p and the assigned error budgets are assumptions.
The four release-3 programs are unchanged; this fifth suite is included in
computational release 4 with the tolerance-opening checks. Copyright 2026 Henry Arellano-Pena. CC BY 4.0.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import mpmath as mp
import sympy as sp


def run(output: Path) -> None:
    mp.mp.dps = 90
    rows: list[dict[str, Any]] = []

    def check(name: str, passed: bool, **detail: Any) -> None:
        row = {"name": name, "passed": bool(passed), **detail}
        rows.append(row)
        if not passed:
            raise AssertionError(name)

    def close(name: str, actual: Any, expected: Any,
              rtol: str = "1e-65", atol: str = "0") -> None:
        a, e = mp.mpf(actual), mp.mpf(expected)
        check(name, mp.isfinite(a) and abs(a-e) <= mp.mpf(atol)+mp.mpf(rtol)*abs(e),
              actual=mp.nstr(a, 25), expected=mp.nstr(e, 25), rtol=rtol, atol=atol)

    def symbolic(name: str, expression: sp.Expr) -> None:
        check(name, sp.simplify(expression.doit()) == 0)

    c, astar = mp.mpf(299792458), mp.mpf("1.2e-10")
    au = mp.mpf(149597870700)
    pc = 648000*au/mp.pi
    gm = 2*mp.mpf("1.3271244e20")  # IAU nominal solar mass parameter
    width, yb, rj = mp.mpf("110e-6"), mp.mpf("1.8"), mp.mpf("1.7")*pc
    gb = yb*astar
    mub, ellb = yb/mp.sqrt(1+yb**2), 1/(1+yb**2)
    ra, rq = 1+ellb, 4*(1+ellb)**2/(2+ellb)**2
    re = mp.sqrt(gm/(mub*gb))

    close("parsec conversion", pc, "3.0856775814913673e16", rtol="3e-17")
    close("Newtonian force at 10 kAU", gm/(10000*au)**2,
          "1.1860167036305601e-10", rtol="1e-15")
    close("external transition radius pc", re/pc, "0.03842364893548668", rtol="3e-16")
    close("single-scale cosmological width m", c**2/astar,
          "7.48962648947348e26", rtol="2e-16")
    close("dimensionless scale hierarchy", astar*width/c**2,
          "1.4686980739907764e-31", rtol="2e-16")
    check("cosmological width exceeds lab benchmark by 30 orders",
          (c**2/astar)/width > mp.mpf("1e30"))

    # Stable high-background calculation (avoid cancellation at machine precision).
    ylab = mp.mpf("9.81")/astar
    mulab = ylab/mp.sqrt(1+ylab**2)
    elllab = 1/(1+ylab**2)
    close("laboratory 1-mu", 1-mulab, mp.mpf("0.5")/ylab**2, rtol="2e-22")
    close("laboratory log slope", elllab, 1/ylab**2, rtol="2e-22")
    check("laboratory anisotropy below 2e-22", elllab < mp.mpf("2e-22"))
    z = sp.symbols("z", positive=True)
    symbolic("high-field mu expansion", sp.series((1+z*z)**(-sp.Rational(1,2)), z, 0, 6).removeO()
             -(1-z*z/2+3*z**4/8))
    symbolic("high-field slope expansion", sp.series(z*z/(1+z*z), z, 0, 6).removeO()
             -(z*z-z**4))
    symbolic("finite-width expansion", sp.series(1-sp.tanh(z)/z, z, 0, 6).removeO()
             -(z*z/3-2*z**4/15))

    width_rows = []
    for label, distance, expected in [
        ("10 kAU", 10000*au, "1.802242768139618e-39"),
        ("20 pc", 20*pc, "1.0590172489582289e-44"),
        ("1 kpc", 1000*pc, "4.236068995832915e-48"),
    ]:
        ql = width/distance
        exact = 1-mp.tanh(ql)/ql
        leading = ql*ql/3
        # At 1 kpc, subtraction consumes 48 decimal places; 90 digits still
        # leave more than 40 reliable relative digits.
        close(f"exact/leading width correction {label}", exact, leading, rtol="1e-35")
        close(f"printed width correction {label}", leading, expected, rtol="4e-16")
        check(f"width is negligible versus 1e-3 at {label}", exact < mp.mpf("1e-30"))
        width_rows.append(dict(scale=label, qL=mp.nstr(ql,25), qL_squared=mp.nstr(ql**2,25),
                               relative_kernel_deficit=mp.nstr(exact,25)))

    # Independent finite-separation diagnostics; not nonlinear binary solutions.
    binary_rows = []
    expected_rows = [
        ("10 kAU",10000*au,"0.5490818072363705","0.6281270227167100","2.3194116613155793e-5"),
        ("0.1 pc",mp.mpf(".1")*pc,"0.1290585519846960","0.1476376797517527","0.0002035416242621617"),
        ("0.2 pc",mp.mpf(".2")*pc,"0.0322646379961740","0.03690941993793816","0.001628332994097293"),
        ("0.5 pc",mp.mpf(".5")*pc,"0.005162342079387841","0.005905507190070107","0.02544270303277020"),
        ("1 pc",pc,"0.0012905855198469603","0.001476376797517527","0.2035416242621616"),
    ]
    for label,r,en0,ea0,t0 in expected_rows:
        en, ea, tide = gm/(r*r*gb), gm/(r*r*mub*gb), (r/rj)**3
        close(f"Newtonian internal/background {label}",en,en0,rtol="8e-16")
        close(f"linear maximum-axis/background {label}",ea,ea0,rtol="8e-16")
        close(f"tidal proxy {label}",tide,t0,rtol="8e-16")
        close(f"transition-radius scaling {label}",ea,(re/r)**2)
        binary_rows.append(dict(separation=label, separation_pc=mp.nstr(r/pc,25),
                                epsilon_N=mp.nstr(en,25),epsilon_A=mp.nstr(ea,25),
                                Newtonian_tidal_proxy=mp.nstr(tide,25)))
    tau = mp.mpf(".001")
    upper = rj*mp.root(tau,3)
    close("tidal upper separation pc",upper/pc,".17")
    windows = []
    for order, expected in [(1,"1.2150624665084204"),(2,"0.21607205661030326")]:
        lower = re*tau**(-mp.mpf(1)/(2*order))
        close(f"conditional lower separation p={order}",lower/pc,expected,rtol="5e-16")
        close(f"nonlinear budget at lower p={order}",(re/lower)**(2*order),tau)
        check(f"unit-coefficient uncorrected window empty p={order}",lower>upper)
        # Maximum admissible coefficient for this chosen pair of budgets.
        cmax = tau*(upper/re)**(2*order)
        close(f"window threshold coefficient p={order}",re*(cmax/tau)**(mp.mpf(1)/(2*order)),upper)
        windows.append(dict(order=order, C_p="1", nonlinear_budget=".001",tidal_budget=".001",
                            lower_pc=mp.nstr(lower/pc,25),upper_pc=mp.nstr(upper/pc,25),
                            maximum_C_p_for_nonempty_window=mp.nstr(cmax,25)))

    # Invert the endpoint condition in each tolerance. Equality is a touching
    # interval, not a finite observing window; test both sides independently.
    threshold_rows = []
    for order in (1, 2):
        exponent = mp.mpf(1)/(2*order)
        lower_at_fixed_n = re*tau**(-exponent)
        critical_t = (re/rj)**3*tau**(-3*exponent)
        critical_n = (re/rj)**(2*order)*tau**(-mp.mpf(2*order)/3)
        close(f"tidal-budget touching endpoints p={order}",
              rj*mp.root(critical_t,3), lower_at_fixed_n)
        close(f"nonlinear-budget touching endpoints p={order}",
              re*critical_n**(-exponent), upper)
        check(f"below tidal threshold closed p={order}",
              rj*mp.root(critical_t*mp.mpf('.99'),3) < lower_at_fixed_n)
        check(f"above tidal threshold open p={order}",
              rj*mp.root(critical_t*mp.mpf('1.01'),3) > lower_at_fixed_n)
        check(f"below nonlinear threshold closed p={order}",
              re*(critical_n*mp.mpf('.99'))**(-exponent) > upper)
        check(f"above nonlinear threshold open p={order}",
              re*(critical_n*mp.mpf('1.01'))**(-exponent) < upper)
        threshold_rows.append(dict(order=order, C_p="1", held_tolerance=".001",
                                   tidal_threshold=mp.nstr(critical_t,25),
                                   nonlinear_threshold=mp.nstr(critical_n,25),
                                   lower_to_upper_ratio=mp.nstr(lower_at_fixed_n/upper,25)))
    close("p2 endpoint ratio", threshold_rows[1]["lower_to_upper_ratio"],
          "1.2710120977076662", rtol="4e-16")
    close("p2 critical tidal tolerance", threshold_rows[1]["tidal_threshold"],
          "0.002053284140956962", rtol="5e-16")
    close("p2 critical nonlinear tolerance", threshold_rows[1]["nonlinear_threshold"],
          "0.002609748983187591", rtol="7e-16")
    x, cp, tn, pp = sp.symbols("x cp tn pp", positive=True)
    tt = x**3*(cp/tn)**(sp.Rational(3,2)/pp)
    symbolic("general threshold reciprocal identity",
             sp.powdenest(cp*x**(2*pp)*tt**(-2*pp/3), force=True)-tn)
    p2lower = re*tau**(-mp.mpf(1)/4)
    check("p2 tidal tolerance 0.0021 opens interval",
          rj*mp.root(mp.mpf('.0021'),3) > p2lower)
    check("p2 nonlinear tolerance 0.0027 opens interval",
          re*mp.mpf('.0027')**(-mp.mpf(1)/4) < upper)
    check("p2 tidal tolerance 0.002 remains closed",
          rj*mp.root(mp.mpf('.002'),3) < p2lower)
    check("p2 nonlinear tolerance 0.0026 remains closed",
          re*mp.mpf('.0026')**(-mp.mpf(1)/4) > upper)

    # General isotropic-flux Hessian tested at a positive z-directed background.
    # Its rotational covariance makes this an arbitrary background direction.
    px,py,pz,b,a = sp.symbols("px py pz b a", real=True, positive=True)
    hx,hy,hz = sp.symbols("hx hy hz",real=True)
    length = sp.sqrt(px*px+py*py+pz*pz)
    mufun = sp.Function("m")
    coords, perturb = [px,py,pz], sp.Matrix([hx,hy,hz])
    flux = sp.Matrix([mufun(length)*v for v in coords])
    n = sp.Matrix([0,0,1])
    m1,m2 = sp.diff(mufun(b),b),sp.diff(mufun(b),b,2)
    Q = m1*hz*perturb + m1*(hx*hx+hy*hy+hz*hz)*n/2 + (b*m2-m1)*hz*hz*n/2
    for i in range(3):
        hessian = sp.hessian(flux[i],coords).subs({px:0,py:0,pz:b})
        symbolic(f"quadratic constitutive Taylor component {i}",
                 (perturb.T*hessian*perturb)[0]/2-Q[i])
        symbolic(f"quadratic flux even under inversion component {i}",
                 Q[i].subs({hx:-hx,hy:-hy,hz:-hz}, simultaneous=True)-Q[i])

    y=sp.symbols("y",positive=True)
    mu=y/sp.sqrt(1+y*y)
    ell=sp.simplify(y*sp.diff(mu,y)/mu)
    cpar=sp.simplify((ell+ell*ell+y*sp.diff(ell,y))/(2*(1+ell)))
    closed=(2-y*y)/(2*(1+y*y)*(2+y*y))
    symbolic("parallel coefficient closed form",cpar-closed)
    f=y*mu
    symbolic("parallel coefficient from direct second derivative",cpar-y*sp.diff(f,y,2)/(2*sp.diff(f,y)))
    cpar_value=mp.mpf(str(sp.N(cpar.subs(y,sp.Rational(9,5)),85)))
    close("parallel local coefficient rounded",cpar_value,"-0.02790580440731672",rtol="1e-15")
    check("local coefficient not unity",abs(cpar_value)<mp.mpf(".03"))
    check("parallel coefficient crosses sign",closed.subs(y,1)>0 and closed.subs(y,2)<0)

    eta,rr,ss=sp.symbols("eta R s",real=True)
    residual=ss*rr-(2-eta)
    symbolic("observed residual mass-slope derivative",sp.diff(residual,eta)-1)
    symbolic("spherical mass-profile slope",sp.diff(sp.log(sp.Function('M')(y)),y)*y
             -y*sp.diff(sp.Function('M')(y),y)/sp.Function('M')(y))
    residual_rows=[]
    for ee in [0,1,2]:
        slope=(2-ee)/ra
        ta=slope*ra-(2-ee)
        tq=slope*rq-(2-ee)
        close(f"local flux residual eta={ee}",ta,0,atol="1e-80")
        close(f"QUMOND mass-profile residual eta={ee}",tq,-(2-ee)*ellb**2/(2+ellb)**2,
              atol="1e-80")
        residual_rows.append(dict(eta_M=ee, slope=mp.nstr(slope,25),
                                  T_Q=mp.nstr(tq,25),three_sigma_allowance=mp.nstr(abs(tq)/3,25)))
    close("eta one QUMOND residual",mp.mpf(residual_rows[1]['T_Q']),"-.01112713418433655575",rtol="2e-15")
    close("eta one 3-sigma budget",abs(mp.mpf(residual_rows[1]['T_Q']))/3,
          ".00370904472811218525",rtol="2e-15")
    close("environment/selected acceleration",mp.mpf('.1')/yb,mp.mpf(1)/18)
    close("squared environmental ratio",(mp.mpf('.1')/yb)**2,mp.mpf(1)/324)

    # Full covariance expression, including the mass-profile derivative.
    v1,v2,v3,c12,c13,c23=sp.symbols("v1 v2 v3 c12 c13 c23", real=True)
    covariance=sp.Matrix([[v1,c12,c13],[c12,v2,c23],[c13,c23,v3]])
    jac=sp.Matrix([[rr,ss,1]])
    expanded=rr**2*v1+ss**2*v2+v3+2*rr*ss*c12+2*rr*c13+2*ss*c23
    symbolic("full covariance propagation",(jac*covariance*jac.T)[0]-expanded)
    # Positive-definite examples with nonzero nuisance correlations.
    for k in [1,2,3]:
        lower=sp.Matrix([[sp.Rational(k,1000),0,0],
                         [sp.Rational(1,3000),sp.Rational(k,2000),0],
                         [-sp.Rational(1,4000),sp.Rational(1,5000),sp.Rational(k,3000)]])
        cov=lower*lower.T
        jnum=sp.Matrix([[sp.Rational(131,106),sp.Rational(106,131),1]])
        check(f"nonnegative correlated variance example {k}",(jnum*cov*jnum.T)[0]>0)
        symbolic(f"covariance as squared projection example {k}",
                 (jnum*cov*jnum.T)[0]-sum(x*x for x in (jnum*lower)))

    report={
        "scope":"Finite algebraic and numerical checks; no catalogue fit or observing-window certification.",
        "inputs": {"GM_total_m3_s2":str(gm),"a_star_m_s2":str(astar),"y_background":str(yb),
                   "width_benchmark_m":str(width),"tidal_benchmark_pc":"1.7",
                   "au_m":str(au),"pc_m":mp.nstr(pc,25),"c_m_s":str(c)},
        "checks":rows,"total_checks":len(rows),"passed":sum(row['passed'] for row in rows),
        "width_suppression":width_rows,"binary_regime_diagnostics":binary_rows,
        "conditional_windows":windows,"tolerance_opening_thresholds":threshold_rows,
        "mass_profile_sensitivity":residual_rows,
        "parallel_local_flux_coefficient":mp.nstr(cpar_value,25),
        "limitations":["The width is a single-Yukawa historical benchmark, not a full-tower likelihood.",
                       "The local flux Taylor coefficient does not bound a nonlinear two-star solution.",
                       "The tidal scale and C_p are diagnostic assumptions, not irreducible noise models.",
                       "The environmental field ratios are not calculated biases of the matched residual."]}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(f"{len(rows)}/{len(rows)} checks passed; report: {output}")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name("empirical_reach_results.json"))
    args=parser.parse_args()
    run(args.output)
