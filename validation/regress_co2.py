"""Regress NH3-CO2 molecule-anion Pitzer parameters to Kurz et al. (1995).
Fit: 313 K + 353 K isotherms. Hold-out: 333 K isotherm. VLE points only; partial pressures <= 2 kPa skipped.
Objective: least squares on ln(p_calc/p_exp) for p(NH3) and p(CO2).
Run: python validation/regress_co2.py"""
import os, sys, math
import numpy as np
from scipy.optimize import least_squares
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, '..', 'dwsim')]
import SWEQ_Solver as S
import validate_co2 as V

ROWS = V.load()
FIT = [r for r in ROWS if round(r['T_K']) in (313, 353)]
HOLD = [r for r in ROWS if round(r['T_K']) == 333]
NAMES = ['lam_NH3_HCO3', 'lam_NH3_NH2COO', 'lam_CO2_HCO3']

def apply(x):
    for k, v in zip(NAMES, x):
        S.PITZER_PARAMS[k][0] = float(v)

def residuals(x, rows):
    apply(x)
    out = []
    for r in rows:
        pp = S.calculate_equilibrium(r['T_K'], {'NH3': r['m_NH3'], 'H2S': 0.0, 'CO2': r['m_CO2']})['pp']
        for k, e in (('NH3', r['p_NH3_kPa']), ('CO2', r['p_CO2_kPa'])):
            if e > 2: out.append(math.log(max(pp[k] * V.PSI_TO_KPA, 1e-9) / e))
    return np.array(out)

if __name__ == '__main__':
    x0 = np.zeros(len(NAMES))
    print("Before:"); print("  fit  ", V.summary(V.evaluate(FIT)).strip()); print("  hold ", V.summary(V.evaluate(HOLD)).strip())
    sol = least_squares(residuals, x0, args=(FIT,), x_scale='jac', max_nfev=300)
    apply(sol.x)
    print("After:"); print("  fit  ", V.summary(V.evaluate(FIT)).strip()); print("  hold ", V.summary(V.evaluate(HOLD)).strip())
    J = sol.jac; s2 = (sol.fun**2).sum() / max(len(sol.fun) - len(sol.x), 1)
    se = np.sqrt(np.diag(np.linalg.pinv(J.T @ J) * s2))
    for k, v, e in zip(NAMES, sol.x, se): print("  %-15s = %+.5f ± %.5f" % (k, v, e))
