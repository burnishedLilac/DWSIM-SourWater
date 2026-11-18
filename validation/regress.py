"""Regress Pitzer parameters of SWEQ against experimental NH3-H2S-H2O data.

Fit set:      Oratovskii et al. 1964 (exact wt% -> molality conversion)
Hold-out set: Leyko 1959 (not used in the fit; density assumed 1.0 kg/L)
Objective:    least squares on ln(p_calc / p_exp) for p(H2S) and p(NH3)

Run: python validation/regress.py   -> prints fitted values to paste into PITZER_PARAMS
"""
import os, sys, math
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'dwsim'))
import SWEQ_Solver as S
import validate as V   # reuses data readers / conversions (also regenerates RESULTS.md on import)

NAMES = [('beta0_NH4_HS', 0), ('beta0_NH4_HS', 1), ('beta1_NH4_HS', 0),
         ('lam_NH3_salt', 0), ('lam_NH3_salt', 1), ('lam_H2S_salt', 0), ('lam_H2S_salt', 1)]
SCALE = np.array([1, 300, 1, 1, 300, 1, 300], float)   # 1/T terms scaled for conditioning

def load(name, conv):
    pts = []
    for r in V.read(name):
        T = float(r['T_K'])
        if conv == 'wt':
            mN, mS = V.from_wt(float(r['w_H2S_pct']), float(r['w_NH3_pct']))
            eS, eN = float(r['p_H2S_bar']), float(r['p_NH3_bar'])
        else:
            mN, mS = V.from_molar(float(r['c_H2S_molL']), float(r['c_NH3_molL']))
            eS, eN = float(r['p_H2S_mmHg']) * V.MMHG_TO_BAR, float(r['p_NH3_mmHg']) * V.MMHG_TO_BAR
        pts.append((T, mN, mS, eS, eN))
    return pts

FIT = load('oratovskii_1964.csv', 'wt')
HOLD = load('leyko_1959.csv', 'molar')

def apply(x):
    for (k, i), v in zip(NAMES, x * SCALE):
        S.PITZER_PARAMS[k][i] = float(v)

def residuals(x, pts):
    apply(x)
    out = []
    for T, mN, mS, eS, eN in pts:
        pp, _ = V.predict(T, mN, mS)
        out += [math.log(pp['H2S'] / eS), math.log(pp['NH3'] / eN)]
    return np.array(out)

def aad(x, pts):
    r = residuals(x, pts)
    d = np.abs(np.exp(r) - 1) * 100
    return d[0::2].mean(), d[1::2].mean()

x0 = np.array([S.PITZER_PARAMS[k][i] for k, i in NAMES]) / SCALE
print("Before  fit-set AAD H2S %.1f%% NH3 %.1f%% | hold-out H2S %.1f%% NH3 %.1f%%" % (*aad(x0, FIT), *aad(x0, HOLD)))
sol = least_squares(residuals, x0, args=(FIT,), method='trf', x_scale='jac', max_nfev=400)
x = sol.x
print("After   fit-set AAD H2S %.1f%% NH3 %.1f%% | hold-out H2S %.1f%% NH3 %.1f%%" % (*aad(x, FIT), *aad(x, HOLD)))

J = sol.jac; dof = max(len(sol.fun) - len(x), 1)
s2 = (sol.fun ** 2).sum() / dof
try:
    se = np.sqrt(np.diag(np.linalg.inv(J.T @ J) * s2))
except np.linalg.LinAlgError:
    se = np.full_like(x, np.nan)
print("\nFitted parameters (value ± 1 std. error):")
for (k, i), v, e in zip(NAMES, x * SCALE, se * SCALE):
    print("  %-14s[%d] = %+.6g  ± %.2g" % (k, i, v, e))
