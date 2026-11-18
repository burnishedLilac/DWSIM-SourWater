"""Sanity checks for SWEQ against textbook limiting values.
Run: python tests/test_limits.py   or   pytest tests/   (no DWSIM needed)
These are limiting-case checks, NOT a validation against sour-water VLE data
(API SWEQ tables, Rumpf 1999, Kurz 1995) - that is still pending."""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'dwsim'))
import SWEQ_Solver as S

def eq(T, nh3=0.0, h2s=0.0, co2=0.0):
    return S.calculate_equilibrium(T, {'NH3': nh3, 'H2S': h2s, 'CO2': co2})

checks = []
def check(name, got, ref, tol):
    ok = abs(got - ref) <= tol
    checks.append(ok)
    print("%-4s %-48s got=%9.4f ref=%9.4f tol=%.3f" % ("OK" if ok else "FAIL", name, got, ref, tol))

atm = lambda psi: psi / S.CONST['ATM_TO_PSI']
r = eq(298.15, h2s=0.01)
check("H2S volatility 25C (atm*kg/mol, 1/0.10)", atm(r['pp']['H2S'])/r['liq']['H2S'], 10.0, 0.5)
r = eq(298.15, co2=0.01)
check("CO2 volatility 25C (atm*kg/mol, 1/0.034)", atm(r['pp']['CO2'])/r['liq']['CO2'], 29.4, 1.5)
check("pH of 0.01 m CO2 water at 25C", r['ph'], 4.18, 0.15)
r = eq(298.15, nh3=0.01)
check("pH of 0.01 m NH3 at 25C", r['ph'], 10.62, 0.15)
r = eq(373.15)
check("Pure water bubble P at 100C (atm)", r['P_bubble']/S.CONST['ATM_TO_PA'], 1.0, 0.02)
r = eq(298.15)
check("Pure water pH at 25C", r['ph'], 7.0, 0.05)
r = eq(298.15, h2s=0.01)
check("pH of 0.01 m H2S at 25C (pK1~7.0)", r['ph'], 4.50, 0.1)

# Flash mass balance
n = {'NH3': 0.5, 'H2S': 0.3, 'CO2': 0.05, 'H2O': 55.5}
V, y, L, res = S.pt_flash(353.15, 0.7*101325.0, n)
bal = max(abs(n[k] - L[k] - (V*y[k] if y else 0)) for k in n)
check("Flash 80C/0.7atm: component balance error (mol)", bal, 0.0, 1e-6)
check("Flash 80C/0.7atm: P_bubble of liquid = P_op (atm)", res['P_bubble']/101325.0, 0.7, 1e-4)
print("V/F = %.5f  y = %s" % (V/sum(n.values()), {k: round(v, 4) for k, v in (y or {}).items()}))

# Above the dew point (96.8 % water at 1 atm, 100 C): all vapor, exact balance
n2 = {'NH3': 0.02, 'H2S': 0.01, 'CO2': 0.002, 'H2O': 0.968}
V2, y2, L2, _ = S.pt_flash(373.15, 101325.0, n2)
check("Flash 100C/1atm dilute sour water: V/F (all vapor)", V2 / sum(n2.values()), 1.0, 1e-9)

def test_all_limits():
    assert all(checks), "some limiting-case checks failed (see output)"

if __name__ == "__main__":
    print("\n%d/%d checks passed" % (sum(checks), len(checks)))
    sys.exit(0 if all(checks) else 1)
