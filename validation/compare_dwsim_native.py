"""Compare SWEQ and DWSIM's built-in Sour Water package (bubble points) and, where available,
Rumpf et al. (1999) measurements. Run: python validation/compare_dwsim_native.py"""
import csv, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'dwsim'))
import SWEQ_Solver as S

RUMPF = {'353.14K m_NH3=6.010 m_H2S=2.539': (1.568, 0.537, 0.607),
         '393.14K m_NH3=3.210 m_H2S=1.056': (4.113, 0.984, 1.186)}   # p_total, p_NH3, p_H2S (bar), Table 1
f = 0.0689476
def load(name):
    return [r for r in csv.DictReader(l for l in open(os.path.join(HERE, 'data', name)) if not l.startswith('#'))]
rows = load('dwsim_native_sourwater.csv')
print("%-6s %-8s | %-27s | %-27s | %-27s" % ("T(C)", "case", "measured P / pNH3 / pH2S", "DWSIM native", "SWEQ " + S.__version__))
for r in rows:
    T = float(r['T_C']) + 273.15; xw = float(r['x_H2O']); kw = xw * S.CONST['MW']['H2O'] / 1000
    q = S.calculate_equilibrium(T, {'NH3': float(r['x_NH3']) / kw, 'H2S': float(r['x_H2S']) / kw, 'CO2': 0.0})
    P = float(r['P_bubble_bar'])
    nat = (P, P * float(r['y_NH3']), P * float(r['y_H2S']))
    swq = (q['P_bubble'] / 1e5, q['pp']['NH3'] * f, q['pp']['H2S'] * f)
    ref = RUMPF.get(r['rumpf_ref'])
    fmt = lambda t: "%.4f / %.4f / %.4f" % t
    print("%-6s %-8s | %-27s | %-27s | %-27s" % (r['T_C'], "Rumpf" if ref else "dilute",
          fmt(ref) if ref else "-", fmt(nat), fmt(swq)))
    if ref:
        dev = lambda t: " ".join("%+.0f%%" % (100 * (a / b - 1)) for a, b in zip(t, ref))
        print("%-15s | %-27s | %-27s | %-27s" % ("", "deviation vs measured:", dev(nat), dev(swq)))

print("\nWithout electrolyte chemistry (DWSIM NRTL), same 80 C Rumpf point:")
for r in load('dwsim_nrtl.csv'):
    P = float(r['P_bubble_bar']); ref = RUMPF[r['rumpf_ref']]
    nrtl = (P, P * float(r['y_NH3']), P * float(r['y_H2S']))
    print("  NRTL %.3f / %.3f / %.3f bar  ->  %s of measured" % (nrtl + (" ".join("%.1fx" % (a / b) for a, b in zip(nrtl, ref)),)))
