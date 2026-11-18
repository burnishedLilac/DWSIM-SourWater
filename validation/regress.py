"""Regress SWEQ Pitzer parameters against Rumpf et al. (1999) NH3-H2S-H2O data.

Fit set:   Rumpf 1999 isotherms 313 K and 393 K
Hold-out:  Rumpf 1999 isotherm 353 K (temperature interpolation), plus the independent
           data sets of Oratovskii 1964 and Leyko 1959 (different laboratories)
Objective: least squares on ln(p_calc/p_exp) for p(NH3) and p(H2S); points with p <= 0.01 bar
           (below detection or nearly so) are skipped.

Run: python validation/regress.py [--write]   (--write stores the result in validation/fitted_params.json)
"""
import os, sys, math, json, csv
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'dwsim'))
import SWEQ_Solver as S

PSI_TO_BAR = 0.0689476
MMHG = 1.01325 / 760.0
MW = S.CONST['MW']

def read(name):
    with open(os.path.join(HERE, 'data', name)) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(l for l in f if not l.startswith('#'))]

def rumpf(Tsel):
    return [(r['T_K'], r['m_NH3'], r['m_H2S'], r['p_NH3_bar'], r['p_H2S_bar'])
            for r in read('rumpf_1999.csv') if round(r['T_K'] / 10) * 10 in Tsel and r['m_H2S'] > 0]

def oratovskii():
    out = []
    for r in read('oratovskii_1964.csv'):
        kw = (100 - r['w_H2S_pct'] - r['w_NH3_pct']) / 1000
        out.append((r['T_K'], r['w_NH3_pct'] / MW['NH3'] / kw, r['w_H2S_pct'] / MW['H2S'] / kw, r['p_NH3_bar'], r['p_H2S_bar']))
    return out

def leyko():
    out = []
    for r in read('leyko_1959.csv'):
        kw = 1.0 - (r['c_H2S_molL'] * MW['H2S'] + r['c_NH3_molL'] * MW['NH3']) / 1000  # density assumed 1.0
        out.append((r['T_K'], r['c_NH3_molL'] / kw, r['c_H2S_molL'] / kw, r['p_NH3_mmHg'] * MMHG, r['p_H2S_mmHg'] * MMHG))
    return out

SETS = {'fit: Rumpf 313+393 K': rumpf({310, 390}), 'hold-out: Rumpf 353 K': rumpf({350}),
        'external: Oratovskii 1964': oratovskii(), 'external: Leyko 1959': leyko()}

NAMES = [('beta0_NH4_HS', 0), ('beta0_NH4_HS', 1), ('beta1_NH4_HS', 0),
         ('lam_H2S_salt', 0), ('lam_H2S_salt', 1), ('lam_NH3_salt', 0)]
SCALE = np.array([1, 300, 1, 1, 300, 1], float)

def apply(x):
    for (k, i), v in zip(NAMES, x * SCALE):
        S.PITZER_PARAMS[k][i] = float(v)

def residuals(x, pts):
    apply(x)
    out = []
    for T, mN, mS, eN, eS in pts:
        pp = S.calculate_equilibrium(T, {'NH3': mN, 'H2S': mS, 'CO2': 0.0})['pp']
        for calc, exp in ((pp['NH3'] * PSI_TO_BAR, eN), (pp['H2S'] * PSI_TO_BAR, eS)):
            if exp > 0.01:
                out.append(math.log(max(calc, 1e-12) / exp))
    return np.array(out)

def aad(x, pts):
    apply(x)
    dN, dS = [], []
    for T, mN, mS, eN, eS in pts:
        pp = S.calculate_equilibrium(T, {'NH3': mN, 'H2S': mS, 'CO2': 0.0})['pp']
        if eN > 0.01: dN.append(abs(pp['NH3'] * PSI_TO_BAR / eN - 1) * 100)
        if eS > 0.01: dS.append(abs(pp['H2S'] * PSI_TO_BAR / eS - 1) * 100)
    return np.mean(dN), np.mean(dS)

def table(x, title):
    print(title)
    for name, pts in SETS.items():
        print("  %-28s p(NH3) %6.1f%%   p(H2S) %6.1f%%" % ((name,) + aad(x, pts)))

if __name__ == '__main__':
    x0 = np.array([S.PITZER_PARAMS[k][i] for k, i in NAMES]) / SCALE
    table(x0, "Before (current defaults):")
    sol = least_squares(residuals, x0, args=(SETS['fit: Rumpf 313+393 K'],), x_scale='jac', max_nfev=500)
    table(sol.x, "After regression:")
    J = sol.jac; s2 = (sol.fun ** 2).sum() / max(len(sol.fun) - len(sol.x), 1)
    se = np.sqrt(np.diag(np.linalg.pinv(J.T @ J) * s2))
    fitted = {}
    print("Parameters (value ± 1 s.e.):")
    for (k, i), v, e in zip(NAMES, sol.x * SCALE, se * SCALE):
        print("  %-14s[%d] = %+.6g ± %.2g" % (k, i, v, e))
        fitted.setdefault(k, list(S.PITZER_PARAMS[k]))[i] = float(v)
    if "--write" in sys.argv:  # optional: store the regressed values for inspection
        json.dump(fitted, open(os.path.join(HERE, 'fitted_params.json'), 'w'), indent=2)
        print("Wrote validation/fitted_params.json")
