"""
SOUR WATER EQUILIBRIUM SOLVER (SWEQ) - v8.2.0
--------------------------------------------------------------------------------
Thermodynamics: Edwards (1978) + PITZER MODEL (Rumpf et al., 1999)
Activity: Pitzer specific interactions for NH3/H2S + Davies fallback for CO2
Henry's Law: Continuous Physical Scaling + Pitzer Activity Coupling
Water Activity: Osmotic Mole Fraction Correction
Stability: Successive Substitution Damping for Non-Linear Convergence
Integration: Optimized for DWSIM Python Script Unit Operation

Changelog v8.2.0:
- FIX: Henry constants for H2S/CO2 were solubilities (mol/kg/atm) used as
  volatilities (atm*kg/mol); partial pressures were ~100x (H2S) and ~900x (CO2) low.
- FIX: Consistent activity convention: a(H+) = 10^-pH; each equilibrium now uses
  the activity coefficient of each species once (no double gamma^2 on H2S/CO2/Kw,
  added missing gamma(NH4+)). Charge balance now in molality (m_H = a_H/gamma).
- NEW: Isothermal PT-flash (Rachford-Rice-style bisection on vapor flow) at the
  operating pressure. oms1 = liquid, oms2 = vapor (as documented).
- CHG: Errors are no longer silently swallowed; they are raised to the DWSIM log.
- FIX: NH4+ and CO2/HCO3- constants replaced (Bates & Pinching; Plummer & Busenberg);
  previous values were ~0.46 pK units off at 25 C.
- CHG: Report is written to disk but no longer auto-opened (OPEN_REPORT flag).

Changelog v8.0.0:
- MAJOR UPGRADE: Replaced global Davies model with Pitzer species-specific model.
- Removed arbitrary Davies limits (Model is now valid up to ~10 molal).
- Added water activity depression for concentrated brines.

Author: Alexander Francisco Cescon
Release: Enterprise Extreme Concentration Edition (GPLv3 Only).
--------------------------------------------------------------------------------
"""

import math

__version__ = "8.2.0"
import sys
import os

# --- 1. CONSTANTS ---

CONST = {
    'R_GAS': 8.314,
    'PSI_TO_PA': 6894.76,
    'ATM_TO_PA': 101325.0,
    'ATM_TO_PSI': 14.6959,
    'MW': {'NH3': 17.031, 'H2S': 34.08, 'CO2': 44.01, 'H2O': 18.015},
    # Henry volatility at 298.15 K, atm*kg/mol (= 1 / solubility constant, Sander compilation)
    # NH3 ~ 1/59, H2S ~ 1/0.10, CO2 ~ 1/0.034
    'H_REF': {'NH3': 1.0/59.0, 'H2S': 1.0/0.10, 'CO2': 1.0/0.034},
    'H_SCALE': {'NH3': 4100.0, 'H2S': 2100.0, 'CO2': 2400.0}
}

SOLVER = {
    'TOL_CHARGE': 1e-9, 'TOL_IONIC': 1e-6,
    'MAX_ITER_PH': 350, 'MAX_ITER_ION': 100, # Increased for Pitzer stability
    'MIN_VAL': 1e-18,
    'CO2_DAVIES_WARN': 0.5, # We still track CO2 limits
    'MAX_ITER_FLASH': 200, 'TOL_FLASH': 1e-7
}

OPEN_REPORT = False  # True = open the datasheet after each run (Windows)

EDWARDS_PARAMS = {
    'NH3':   [1.587,   11160.0,   0.0,         0.0,        0.0],
    'H2S':   [-293.88, 683858.0, -6.27125e8,  2.555e11,   -3.91757e13],
    'HS':    [-220.07, 258273.0, -1.8396e8,   6.809e10,   -1.0267e13],
    'CO2':   [-241.79, 536256.0, -4.8123e8,   1.94e11,    -2.96445e13],
    'HCO3':  [-294.74, 364385.0, -2.841e8,    1.23323e11, -2.0759e13],
    'Kw':    [39.5554, -177822.0, 1.843e8,    -0.8541e11, 1.4292e13]
}

# Pitzer Parameters (Standard Natural Base)
PITZER_PARAMS = {
    # Gases (Rumpf converted)
    'beta0_NH3_NH3':     [-0.01979,   9.864,      0.0],
    'beta0_H2S_H2S':     [-0.26156,   69.751,     0.0],
    
    # Ions NH4-HS (Standard Pitzer 1-1 Electrolyte at 25C reference)
    'beta0_NH4_HS':      [0.052, 0.0, 0.0], 
    'beta1_NH4_HS':      [0.180, 0.0, 0.0], 
    'Cphi_NH4_HS':       [-0.001, 0.0, 0.0]
}

COMP_MAP = {
    'NH3': ['ammonia', 'nh3', 'amonia'],
    'H2S': ['hydrogen sulfide', 'h2s', 'sulfeto'],
    'CO2': ['carbon dioxide', 'co2', 'dioxido'],
    'H2O': ['water', 'h2o', 'agua']
}

try:
    import clr
    from System import Array, Double, DateTime
except ImportError:
    pass

# --- 2. THERMODYNAMICS & PITZER ENGINE ---

def calc_ln_k(coeffs, T_R):
    A, B, C, D, E = coeffs
    return A + B/T_R + C/(T_R**2) + D/(T_R**3) + E/(T_R**4)

def calc_davies_A(T_K):
    t = T_K - 273.15
    return 0.4913 + 6.08e-4 * t + 5.95e-6 * t**2

def get_activity_coefficients(I, A):
    if I < 1e-11: return 1.0, 1.0
    sqI = math.sqrt(I)
    f = -A * (sqI / (1.0 + sqI) - 0.3 * I)
    f = max(min(f, 2.0), -2.0)  # Davies is meaningless beyond I ~ 0.5; keep numerics bounded
    return 10**f, 10**(f * 4)

def calc_pitzer_p(q, T_K):
    return q[0] + q[1]/T_K + (q[2]*math.log(T_K) if q[2] != 0 else 0.0)

def get_pitzer_gammas(T_K, m_NH3, m_H2S, m_NH4, m_HS):
    """Calculates activity coefficients using standard Pitzer formulation."""
    A_phi = 0.3915 + 0.00021*(T_K - 298.15) + 2.2e-6*((T_K - 298.15)**2)
    
    b0_A_A = calc_pitzer_p(PITZER_PARAMS['beta0_NH3_NH3'], T_K)
    b0_S_S = calc_pitzer_p(PITZER_PARAMS['beta0_H2S_H2S'], T_K)
    
    b0_C_A = calc_pitzer_p(PITZER_PARAMS['beta0_NH4_HS'], T_K)
    b1_C_A = calc_pitzer_p(PITZER_PARAMS['beta1_NH4_HS'], T_K)
    cphi_C_A = calc_pitzer_p(PITZER_PARAMS['Cphi_NH4_HS'], T_K)
    
    I = 0.5 * (m_NH4 + m_HS)
    if I < 1e-10: return 1.0, 1.0, 1.0
        
    sqI = math.sqrt(I)
    f_I = -A_phi * ( (sqI / (1.0 + 1.2*sqI)) + (2.0/1.2)*math.log(1.0 + 1.2*sqI) )
    
    alpha = 2.0
    x = alpha * sqI
    g_x = 2.0 * (1.0 - (1.0 + x)*math.exp(-x)) / (x**2) if x > 1e-8 else 1.0
    
    # Standard Pitzer B_gamma term
    B_gamma = 2.0 * b0_C_A + (2.0 * b1_C_A / (alpha**2 * I)) * (1.0 - (1.0 + alpha*sqI - 0.5*(alpha**2)*I)*math.exp(-alpha*sqI))
    
    # C_gamma term
    C_gamma = 1.5 * cphi_C_A
    
    m = (m_NH4 + m_HS) / 2.0 
    
    # Mean ionic activity coefficient (Standard Pitzer 1-1)
    ln_gamma_pm = f_I + m * B_gamma + (m**2) * C_gamma
    
    # Neutral species (Simplified)
    ln_gamma_NH3 = 2.0 * b0_A_A * m_NH3 
    ln_gamma_H2S = 2.0 * b0_S_S * m_H2S
    
    # Safety caps
    ln_gamma_pm = max(min(ln_gamma_pm, 10.0), -10.0)
    ln_gamma_NH3 = max(min(ln_gamma_NH3, 10.0), -10.0)
    ln_gamma_H2S = max(min(ln_gamma_H2S, 10.0), -10.0)
    
    return math.exp(ln_gamma_NH3), math.exp(ln_gamma_H2S), math.exp(ln_gamma_pm)

def solve_charge_balance(T_K, ph, K, molals, A_davies):
    """Speciation at fixed pH. Convention: a_H = 10^-pH (activity).
    Concentrations are molalities; each equilibrium uses each species' gamma once.
    gamma(+1/-1 ions of NH4/HS) = Pitzer mean gamma; other ions = Davies."""
    a_h = max(10**(-ph), SOLVER['MIN_VAL'])
    m_NH3_total, m_H2S_total, m_CO2_total = molals

    g1_davies, g2_davies = 1.0, 1.0
    g_NH3, g_H2S, g_pm = 1.0, 1.0, 1.0
    nh3_aq, h2s_aq = m_NH3_total, m_H2S_total
    nh4_ion, hs_ion = 0.0, 0.0
    ionic_strength = 0.0
    prev_I = -1.0
    w = 0.5  # damping

    for _ in range(SOLVER['MAX_ITER_ION']):
        if abs(ionic_strength - prev_I) < SOLVER['TOL_IONIC']: break
        prev_I = ionic_strength

        g_NH3_new, g_H2S_new, g_pm_new = get_pitzer_gammas(T_K, nh3_aq, h2s_aq, nh4_ion, hs_ion)
        g1_dav_new, g2_dav_new = get_activity_coefficients(ionic_strength, A_davies)
        g_NH3 = w * g_NH3 + (1-w) * g_NH3_new
        g_H2S = w * g_H2S + (1-w) * g_H2S_new
        g_pm = w * g_pm + (1-w) * g_pm_new
        g1_davies = w * g1_davies + (1-w) * g1_dav_new
        g2_davies = w * g2_davies + (1-w) * g2_dav_new

        # AMMONIA: NH4+ <=> NH3 + H+ ; Ka = a_NH3 a_H / a_NH4
        ka_nh4 = 1.0 / K['NH3']
        r_nh4 = (g_NH3 * a_h) / (ka_nh4 * g_pm)          # m_NH4 / m_NH3
        nh3_aq = m_NH3_total / (1.0 + r_nh4)
        nh4_ion = m_NH3_total - nh3_aq

        # SULFIDE: H2S <=> H+ + HS- ; HS- <=> H+ + S--
        r1 = K['H2S'] * g_H2S / (g_pm * a_h)              # m_HS / m_H2S
        r2 = K['HS'] * g_pm / (g2_davies * a_h)           # m_S / m_HS
        h2s_aq = m_H2S_total / (1.0 + r1 + r1*r2)
        hs_ion = h2s_aq * r1
        s_ion = hs_ion * r2

        # CARBONATE (Davies): CO2 + H2O <=> H+ + HCO3- ; HCO3- <=> H+ + CO3--
        c1 = K['CO2'] / (g1_davies * a_h)
        c2 = K['HCO3'] * g1_davies / (g2_davies * a_h)
        co2_aq = m_CO2_total / (1.0 + c1 + c1*c2)
        hco3_ion = co2_aq * c1
        co3_ion = hco3_ion * c2

        # WATER: Kw = a_H a_OH
        oh_ion = K['Kw'] / (g1_davies * a_h)
        h_ion = a_h / g1_davies

        ionic_strength = 0.5 * (h_ion + oh_ion + nh4_ion + hs_ion + hco3_ion + 4.0*s_ion + 4.0*co3_ion)

    error = (h_ion + nh4_ion) - (oh_ion + hs_ion + 2.0*s_ion + hco3_ion + 2.0*co3_ion)
    spec = {'NH3': nh3_aq, 'NH4': nh4_ion, 'H2S': h2s_aq, 'HS': hs_ion, 'S': s_ion,
            'CO2': co2_aq, 'HCO3': hco3_ion, 'CO3': co3_ion, 'H': h_ion, 'OH': oh_ion,
            'I': ionic_strength, 'g_NH3': g_NH3, 'g_H2S': g_H2S, 'g_pm': g_pm, 'g_dav': g1_davies}
    return error, spec

def literature_K(T_K):
    """Replacements for Edwards-form constants that were ~0.46 pK off at 25 C.
    NH4+ : Bates & Pinching (1949)   pKa = 0.09018 + 2729.92/T
    CO2/HCO3- : Plummer & Busenberg (1982)"""
    lT = math.log10(T_K)
    pKa_nh4 = 0.09018 + 2729.92 / T_K
    logK1 = -356.3094 - 0.06091964*T_K + 21834.37/T_K + 126.8339*lT - 1684915.0/T_K**2
    logK2 = -107.8871 - 0.03252849*T_K + 5151.79/T_K + 38.92561*lT - 563713.9/T_K**2
    return {'NH3': 10**pKa_nh4, 'CO2': 10**logK1, 'HCO3': 10**logK2}

def calculate_equilibrium(T_K, m_in):
    if T_K < 273.15: raise ValueError("SWEQ: temperature below 0 C")
    T_R = T_K * 1.8
    K_vals = {k: math.exp(calc_ln_k(v, T_R)) for k, v in EDWARDS_PARAMS.items()}
    K_vals.update(literature_K(T_K))  # K['NH3'] = 1/Ka(NH4+)
    A_davies = calc_davies_A(T_K)
    
    low, high = 0.0, 14.0
    ph_guess = 7.0
    spec = {}
    
    for _ in range(SOLVER['MAX_ITER_PH']):
        ph_guess = (low + high) / 2.0
        err, spec = solve_charge_balance(T_K, ph_guess, K_vals, (m_in['NH3'], m_in['H2S'], m_in['CO2']), A_davies)
        if abs(err) < SOLVER['TOL_CHARGE']: break
        if err > 0: low = ph_guess
        else: high = ph_guess
    
    nh3, h2s, co2 = spec['NH3'], spec['H2S'], spec['CO2']
    inv_T_term = (1.0/298.15 - 1.0/T_K)
    
    # Henry's Law with Pitzer Activities
    p_nh3 = (CONST['H_REF']['NH3'] * nh3 * spec['g_NH3'] * math.exp(CONST['H_SCALE']['NH3'] * inv_T_term)) * CONST['ATM_TO_PSI']
    p_h2s = (CONST['H_REF']['H2S'] * h2s * spec['g_H2S'] * math.exp(CONST['H_SCALE']['H2S'] * inv_T_term)) * CONST['ATM_TO_PSI']
    p_co2 = (CONST['H_REF']['CO2'] * co2 * 1.0 * math.exp(CONST['H_SCALE']['CO2'] * inv_T_term)) * CONST['ATM_TO_PSI']
    
    # Water Activity (Osmotic Depression)
    # Ideal (phi = 1) osmotic approximation
    sum_moles = sum([v for k,v in spec.items() if k not in ['I', 'g_NH3', 'g_H2S', 'g_pm', 'g_dav']])
    a_w = math.exp(-0.018015 * sum_moles)
    
    p_w_pure = (10**(5.20389 - 1733.926/(T_K - 39.485))) * 14.5038
    p_w = p_w_pure * a_w
    
    p_tot_pa = (p_nh3 + p_h2s + p_co2 + p_w) * CONST['PSI_TO_PA']
    
    return {'ph': ph_guess, 'P_bubble': p_tot_pa, 'pp': {'NH3': p_nh3, 'H2S': p_h2s, 'CO2': p_co2, 'H2O': p_w}, 'liq': spec, 'aw': a_w}

def calculate_density(T_K, P_Pa, spec, kg_water, mass_total_kg):
    t_c = T_K - 273.15
    rho_w = 1000.0 * (1 - (t_c + 288.9414)/(508929.2 * (t_c + 68.12963)) * (t_c - 3.9863)**2)
    v_inf = {'NH4': 18.0, 'HS': 20.0, 'NH3': 24.5, 'H2S': 35.0, 'CO2': 34.0, 'HCO3': 25.0, 'CO3': 20.0, 'S': 20.0}
    vol_water_L = kg_water / (rho_w / 1000.0)
    vol_solutes_mL = (
        spec['NH4']*v_inf['NH4'] + spec['HS']*v_inf['HS'] + spec['NH3']*v_inf['NH3'] + spec['H2S']*v_inf['H2S'] +
        spec['CO2']*v_inf['CO2'] + spec['HCO3']*v_inf['HCO3'] + spec['CO3']*v_inf['CO3'] + spec['S']*v_inf['S']
    ) * kg_water
    vol_solutes_L = vol_solutes_mL / 1000.0
    total_vol_L = vol_water_L + vol_solutes_L
    rho_kg_L = mass_total_kg / total_vol_L
    return (rho_kg_L * 1000.0) * (1 + 4.5e-10 * (P_Pa - 101325.0))

def generate_report(res, T, P_op, flows_h, total_h, rho):
    l = res['liq']; pp = res['pp']
    W = 80; HR = "=" * W; SR = "-" * W
    is_f = res['P_bubble'] > (P_op * 1.01)
    
    rho_kg_L = rho / 1000.0
    h2s_mgL = (l['H2S'] + l['HS'] + l['S']) * CONST['MW']['H2S'] * rho_kg_L * 1000.0
    nh3_mgL = (l['NH3'] + l['NH4']) * CONST['MW']['NH3'] * rho_kg_L * 1000.0

    lines = [HR, " SWEQ - SOUR WATER EQUILIBRIUM SOLVER v8.2.0 ".center(W), " Edwards / Pitzer (NH4HS) + Davies (carbonate) ".center(W), HR]
    try: lines.append(" Date: %s " % DateTime.Now.ToString("yyyy-MM-dd HH:mm").center(W))
    except: pass
    lines.append("\n 1. EXECUTIVE SUMMARY & SAFETY CHECK ".ljust(W))
    lines.append(SR)
    lines.append(" STATUS: " + ("!!! FLASHING DETECTED !!!" if is_f else "STABLE LIQUID"))
    lines.append("")
    lines.append(" %-25s | %-15s | %-10s | %-20s" % ("PARAMETER", "VALUE", "UNIT", "COMMENT"))
    lines.append(" %s-+-%s-+-%s-+-%s" % ("-"*25, "-"*15, "-"*10, "-"*20))
    lines.append(" %-25s | %15.2f | %-10s | Inlet Condition" % ("Temperature", T-273.15, "C"))
    lines.append(" %-25s | %15.3f | %-10s | Operating P" % ("System Pressure", P_op/CONST['ATM_TO_PA'], "atm"))
    lines.append(" %-25s | %15.3f | %-10s | Equilibrium P" % ("Bubble Pressure", res['P_bubble']/CONST['ATM_TO_PA'], "atm"))
    lines.append(" %-25s | %15.4f | %-10s | Buffering Status" % ("Calculated pH", res['ph'], "-"))
    lines.append(" %-25s | %15.4f | %-10s | Driver" % ("Ionic Strength (I)", l['I'], "molal"))
    lines.append(" %-25s | %15.4f | %-10s | Osmotic Correction" % ("Water Activity", res['aw'], "a_w"))
    
    lines.append("\n 2. STREAM COMPOSITION & FLOWS ".ljust(W))
    lines.append(SR)
    lines.append(" %-25s | %15s | %-10s " % ("COMPONENT", "FLOW (kg/h)", "MASS %"))
    lines.append(" %s-+-%s-+-%s" % ("-"*25, "-"*15, "-"*10))
    lines.append(" %-25s | %15.2f | %10.2f " % ("Water (H2O)", flows_h['H2O'], (flows_h['H2O']/total_h)*100))
    for k in ['NH3', 'H2S', 'CO2']:
        lines.append(" %-25s | %15.2f | %10.2f " % (k, flows_h[k], (flows_h[k]/total_h)*100))
        
    lines.append("\n 3. CHEMICAL SPECIATION & PITZER ACTIVITIES ".ljust(W))
    lines.append(SR)
    lines.append(" Concentrations in molality (mol/kg H2O)")
    lines.append("")
    lines.append("  AMMONIA SYSTEM   ::  NH3(aq): %7.4f (g=%.2f) <==> NH4+: %7.4f (g=%.2f)" % (l['NH3'], l['g_NH3'], l['NH4'], l['g_pm']))
    lines.append("  SULFIDE SYSTEM   ::  H2S(aq): %7.4f (g=%.2f) <==> HS-:  %7.4f (g=%.2f)" % (l['H2S'], l['g_H2S'], l['HS'], l['g_pm']))
    lines.append("  CARBONATE SYSTEM ::  CO2(aq): %7.4f          <==> HCO3: %7.4f (Davies)" % (l['CO2'], l['HCO3']))
    
    lines.append("\n 4. ENVIRONMENTAL & VAPOR PREDICTIONS ".ljust(W))
    lines.append(SR)
    lines.append(" %-38s | %s" % ("ENVIRONMENTAL (Liquid Effluent)", "VAPOR PHASE (Equilibrium Partial P)"))
    lines.append(" %s-+-%s" % ("-"*38, "-"*39))
    lines.append(" H2S Total: %15.1f mg/L      | p(H2S): %10.4f psia" % (h2s_mgL, pp['H2S']))
    lines.append(" NH3 Total: %15.1f mg/L      | p(NH3): %10.4f psia" % (nh3_mgL, pp['NH3']))
    lines.append(" Density:   %15.2f kg/m3     | p(CO2): %10.4f psia" % (rho, pp['CO2']))
    lines.append("                                        | p(H2O): %10.4f psia" % (pp['H2O']))

    if res.get('V_frac') is not None:
        lines.append("\n Flash @ operating P: vapor fraction (molar) = %.6f" % res['V_frac'])
    if (l['CO2'] + l['HCO3']) > SOLVER['CO2_DAVIES_WARN']:
        lines.append("\n > Note: High CO2 detected. Carbonate interactions using Davies fallback.")
        
    lines.append(HR + "\n End of Report")
    return "\n".join(lines)

# --- 5. PT FLASH ---

KEYS = ['NH3', 'H2S', 'CO2', 'H2O']

def bubble_from_moles(T_K, n):
    """n: mol of each key in liquid. Returns equilibrium result (P_bubble, pp in psia)."""
    kg_w = n['H2O'] * CONST['MW']['H2O'] / 1000.0
    molals = {k: n[k] / kg_w for k in ['NH3', 'H2S', 'CO2']}
    return calculate_equilibrium(T_K, molals)

def _rachford_rice(z, K):
    g = lambda b: sum(z[k]*(K[k]-1.0)/(1.0+b*(K[k]-1.0)) for k in KEYS)
    if g(0.0) <= 0: return 0.0
    if g(1.0) >= 0: return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if g(mid) > 0: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

def pt_flash(T_K, P_Pa, n_feed):
    """Isothermal flash: Rachford-Rice with K_i = y_i/x_i updated by successive
    substitution from the electrolyte model (K_i = p_i / (P x_i)).
    Mass balance is exact by construction. Returns (V, y, L_moles, res_liquid)."""
    F = sum(n_feed.values())
    z = {k: n_feed[k]/F for k in KEYS}
    res = bubble_from_moles(T_K, n_feed)
    if res['P_bubble'] <= P_Pa:
        return 0.0, None, dict(n_feed), res
    psi_P = P_Pa / CONST['PSI_TO_PA']
    def k_values(res, x):
        return {k: max((res['pp'][k]/psi_P) / max(x[k], 1e-30), 1e-12) for k in KEYS}
    K = k_values(res, z)
    beta = 0.0
    for _ in range(SOLVER['MAX_ITER_FLASH']):
        beta = _rachford_rice(z, K)
        if beta >= 1.0:
            raise ValueError("SWEQ: feed fully vaporizes at this T/P (outside model range)")
        x = {k: z[k]/(1.0+beta*(K[k]-1.0)) for k in KEYS}
        sx = sum(x.values()); x = {k: v/sx for k, v in x.items()}
        L = {k: F*(1.0-beta)*x[k] for k in KEYS}
        res = bubble_from_moles(T_K, L)
        K_new = k_values(res, x)
        d = max(abs(math.log(K_new[k]/K[k])) for k in KEYS)
        K = {k: math.exp(0.5*math.log(K[k]) + 0.5*math.log(K_new[k])) for k in KEYS}
        if d < 1e-8: break
    V = F*beta
    y = {k: K[k]*x[k] for k in KEYS}
    sy = sum(y.values()); y = {k: v/sy for k, v in y.items()}
    L = {k: n_feed[k] - V*y[k] for k in KEYS}
    return V, y, L, bubble_from_moles(T_K, L)

# --- 6. MAIN BRIDGE (DWSIM) ---

def Main():
    T = ims1.GetTemperature()
    P_op = ims1.GetPressure()
    F_mol = ims1.GetMolarFlow()
    comp = ims1.GetOverallComposition()
    ids = list(ims1.ComponentIds)

    if F_mol < SOLVER['MIN_VAL']: return

    n = {k: 0.0 for k in KEYS}
    idx = {}
    for i, name in enumerate(ids):
        nc = name.lower()
        for k, tags in COMP_MAP.items():
            if any(t in nc for t in tags):
                n[k] = F_mol * comp[i]; idx[k] = i
                break
    if n['H2O'] < SOLVER['MIN_VAL']:
        raise ValueError("SWEQ: no water found in inlet stream")
    unmapped = [ids[i] for i in range(len(ids)) if i not in idx.values() and comp[i] > 1e-6]
    if unmapped:
        raise ValueError("SWEQ: unsupported components in feed: %s" % ", ".join(unmapped))

    V, y, L, res = pt_flash(T, P_op, n)
    Lt = sum(L.values())

    x_full = [0.0] * len(ids); y_full = [0.0] * len(ids)
    for k, i in idx.items():
        x_full[i] = L[k] / Lt
        if y: y_full[i] = y[k]

    oms1.SetTemperature(T); oms1.SetPressure(P_op)
    oms1.SetOverallComposition(Array[Double](x_full))
    oms1.SetMolarFlow(Lt)

    oms2.SetTemperature(T); oms2.SetPressure(P_op)
    oms2.SetOverallComposition(Array[Double](y_full if y else x_full))
    oms2.SetMolarFlow(V if V > 0 else 0.0)

    oms1.Calculate(); oms2.Calculate()

    res['V_frac'] = V / F_mol
    m_kg_s = {k: v * CONST['MW'][k] / 1000.0 for k, v in L.items()}
    kg_w = m_kg_s['H2O']
    rho = calculate_density(T, P_op, res['liq'], kg_w, sum(m_kg_s.values()))
    report_str = generate_report(res, T, P_op, {k: m*3600 for k, m in m_kg_s.items()}, sum(m_kg_s.values())*3600, rho)

    for folder in [os.path.join(os.path.expanduser('~'), 'Desktop'), 'C:\\Temp']:
        try:
            if not os.path.exists(folder): os.makedirs(folder)
            path = os.path.join(folder, 'SWEQ_Datasheet.txt')
            with open(path, 'w') as fh: fh.write(report_str)
            if OPEN_REPORT and os.name == 'nt': os.startfile(path)
            break
        except (IOError, OSError):
            continue

if 'ims1' in globals():
    Main()
