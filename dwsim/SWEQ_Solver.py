"""
SOUR WATER EQUILIBRIUM SOLVER (SWEQ) - v8.4.0
--------------------------------------------------------------------------------
Thermodynamics: Edwards (1978) + PITZER MODEL (Rumpf et al., 1999)
Activity: Pitzer specific interactions for NH3/H2S + Davies fallback for CO2
Henry's Law: Continuous Physical Scaling + Pitzer Activity Coupling
Water Activity: Osmotic Mole Fraction Correction
Stability: Successive Substitution Damping for Non-Linear Convergence
Integration: Optimized for DWSIM Python Script Unit Operation

Changelog v8.4.0:
- NEW: Carbamate (NH2COO-) chemistry, Edwards et al. (1978) constant; validated against Kurz et al. (1995).
- FIX: CO2 Henry constant from Carroll, Slupsky & Mather (1991).
- FIX: Davies equation replaced by the Pitzer Debye-Hueckel term for all ions (Davies gave gamma = 2-6
  for HCO3- at high ionic strength); ionic strength now includes every ion.
- CHG: Ion-molecule lambda placed on the anion (Kurz et al. convention); identical for NH4HS solutions.
- NEW: lambda(NH3-NH2COO-) regressed to Kurz et al. (1995).

Changelog v8.3.1:
- FIX: Outlet 2 (vapor) is optional when nothing vaporizes; clear message when it is needed.
- FIX: Feed above its dew point no longer raises an error; all flow goes to the vapor outlet (oms2).

Changelog v8.3.0:
- FIX: H2S Henry constant temperature dependence from IUPAC SDS-32 correlation (was ~1.7x too
  volatile at 393 K).
- NEW: Ion-molecule Pitzer terms (NH3/H2S with NH4HS) and NH3 ternary term.
- CHG: beta1(NH4HS), lambda(NH3-NH4HS), lambda(H2S-NH4HS) regressed to Rumpf et al. (1999).

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

__version__ = "8.4.0"
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
    # v8.3.0: beta1 and the lambda terms regressed to Rumpf et al. (1999) 313 + 393 K
    # (validation/regress.py); beta0 and Cphi kept from previous versions.
    'beta1_NH4_HS':      [0.0178, 0.0, 0.0],
    'Cphi_NH4_HS':       [-0.001, 0.0, 0.0],
    # Ion-molecule interactions (lambda of NH3 / H2S with NH4+ and HS-, taken equal for both ions)
    'lam_NH3_salt':      [-0.0181, 0.0, 0.0],
    # Ternary NH3-NH3-NH3 (Rumpf et al. convention, sum over unordered triples)
    'tau_NH3_NH3_NH3':   [0.0, 0.0, 0.0],
    # NH3-CO2 system: molecule-anion interactions (lambda), regressed to Kurz et al. (1995)
    'lam_NH3_HCO3':      [0.0, 0.0, 0.0],
    'lam_NH3_NH2COO':    [0.0528, 0.0, 0.0],   # v8.4.0, Kurz 1995 313+353 K (validation/regress_co2.py)
    'lam_CO2_HCO3':      [0.0, 0.0, 0.0],
    'lam_H2S_salt':      [-0.0184, 0.0, 0.0]
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

def get_pitzer_gammas(T_K, m_NH3, m_H2S, m_NH4, m_HS, I, m_HCO3=0.0, m_CARB=0.0, m_CO2=0.0):
    """Activity coefficients (Pitzer). I is the TOTAL ionic strength of the solution.
    - Debye-Hueckel term F (Pitzer form, b = 1.2) applies to every ion: ln g = z^2 F + ...
    - NH4+/HS- binary terms (beta0, beta1, Cphi) apply between NH4+ and HS- only.
    - Ion-molecule terms follow Kurz et al. (1995): lambda(molecule, NH4+) = 0 and the whole
      salt effect is put on the anion, lambda(molecule, HS-) = 2*lam ('lam_*_salt' = B(G,NH4HS)/2).
      For NH4HS solutions this is identical to splitting it equally; it keeps NH4+ paired with
      HCO3-/NH2COO- from inheriting the NH4HS salting effect. No data yet for HCO3-/CO3--/NH2COO-.
    - NH3/CO2 with HCO3- and NH2COO-: lambda terms regressed to Kurz et al. (1995).
    Returns dict with g_NH3, g_H2S, g_NH4, g_HS, g_HCO3, g_CARB, g_CO2, g_pm (NH4HS mean) and F."""
    A_phi = 0.3915 + 0.00021*(T_K - 298.15) + 2.2e-6*((T_K - 298.15)**2)

    b0_A_A = calc_pitzer_p(PITZER_PARAMS['beta0_NH3_NH3'], T_K)
    b0_S_S = calc_pitzer_p(PITZER_PARAMS['beta0_H2S_H2S'], T_K)
    b0_C_A = calc_pitzer_p(PITZER_PARAMS['beta0_NH4_HS'], T_K)
    b1_C_A = calc_pitzer_p(PITZER_PARAMS['beta1_NH4_HS'], T_K)
    cphi_C_A = calc_pitzer_p(PITZER_PARAMS['Cphi_NH4_HS'], T_K)
    lam_N = calc_pitzer_p(PITZER_PARAMS['lam_NH3_salt'], T_K)
    lam_S = calc_pitzer_p(PITZER_PARAMS['lam_H2S_salt'], T_K)
    tau_N = calc_pitzer_p(PITZER_PARAMS['tau_NH3_NH3_NH3'], T_K)

    if I > 1e-12:
        sqI = math.sqrt(I)
        F = -A_phi * ((sqI / (1.0 + 1.2*sqI)) + (2.0/1.2)*math.log(1.0 + 1.2*sqI))
        alpha = 2.0
        B_gamma = 2.0 * b0_C_A + (2.0 * b1_C_A / (alpha**2 * I)) * (1.0 - (1.0 + alpha*sqI - 0.5*(alpha**2)*I)*math.exp(-alpha*sqI))
    else:
        F, B_gamma = 0.0, 2.0 * b0_C_A + b1_C_A
    C_gamma = 1.5 * cphi_C_A

    ln_NH4 = F + m_HS * B_gamma + m_HS * m_NH4 * C_gamma
    ln_HS = F + m_NH4 * B_gamma + m_HS * m_NH4 * C_gamma + 4.0 * (lam_N * m_NH3 + lam_S * m_H2S)
    lNa = calc_pitzer_p(PITZER_PARAMS['lam_NH3_HCO3'], T_K)
    lNc = calc_pitzer_p(PITZER_PARAMS['lam_NH3_NH2COO'], T_K)
    lCa = calc_pitzer_p(PITZER_PARAMS['lam_CO2_HCO3'], T_K)
    ln_NH3 = 2.0 * b0_A_A * m_NH3 + 3.0 * tau_N * m_NH3**2 + 4.0 * lam_N * m_HS + 2.0 * (lNa * m_HCO3 + lNc * m_CARB)
    ln_H2S = 2.0 * b0_S_S * m_H2S + 4.0 * lam_S * m_HS
    ln_CO2 = 2.0 * lCa * m_HCO3
    ln_HCO3 = F + 2.0 * (lNa * m_NH3 + lCa * m_CO2)
    ln_CARB = F + 2.0 * lNc * m_NH3

    cap = lambda v: max(min(v, 10.0), -10.0)
    ln_NH4, ln_HS, ln_NH3, ln_H2S = cap(ln_NH4), cap(ln_HS), cap(ln_NH3), cap(ln_H2S)
    ln_CO2, ln_HCO3, ln_CARB = cap(ln_CO2), cap(ln_HCO3), cap(ln_CARB)
    return {'g_NH3': math.exp(ln_NH3), 'g_H2S': math.exp(ln_H2S), 'g_NH4': math.exp(ln_NH4),
            'g_HS': math.exp(ln_HS), 'g_HCO3': math.exp(ln_HCO3), 'g_CARB': math.exp(ln_CARB),
            'g_CO2': math.exp(ln_CO2), 'g_pm': math.exp(0.5 * (ln_NH4 + ln_HS)), 'F': F}

def solve_charge_balance(T_K, ph, K, molals, A_davies=None):
    """Speciation at fixed pH. Convention: a_H = 10^-pH (activity); balances in molality;
    each equilibrium uses each species' activity coefficient once (see get_pitzer_gammas)."""
    a_h = max(10**(-ph), SOLVER['MIN_VAL'])
    m_NH3_total, m_H2S_total, m_CO2_total = molals

    G = {'g_NH3': 1.0, 'g_H2S': 1.0, 'g_NH4': 1.0, 'g_HS': 1.0, 'g_HCO3': 1.0, 'g_CARB': 1.0,
         'g_CO2': 1.0, 'g_pm': 1.0, 'F': 0.0}
    nh3_aq, h2s_aq = m_NH3_total, m_H2S_total
    co2_aq, hco3_ion, carb_ion = m_CO2_total, 0.0, 0.0
    nh4_ion, hs_ion = 0.0, 0.0
    a_w = 1.0
    ionic_strength = 0.0
    prev_I = -1.0
    w = 0.5  # damping

    for _ in range(SOLVER['MAX_ITER_ION']):
        if abs(ionic_strength - prev_I) < SOLVER['TOL_IONIC']: break
        prev_I = ionic_strength

        G_new = get_pitzer_gammas(T_K, nh3_aq, h2s_aq, nh4_ion, hs_ion, ionic_strength,
                                  hco3_ion, carb_ion, co2_aq)
        G = {k: w * G[k] + (1 - w) * G_new[k] for k in G}
        g_NH3, g_H2S, g_NH4, g_HS = G['g_NH3'], G['g_H2S'], G['g_NH4'], G['g_HS']
        g1 = math.exp(G['F'])          # other 1:1 ions (H+, OH-)
        g_HCO3, g_CARB, g_CO2 = G['g_HCO3'], G['g_CARB'], G['g_CO2']
        g2 = math.exp(4.0 * G['F'])    # 2:1 ions (S--, CO3--)

        # SULFIDE: H2S <=> H+ + HS- ; HS- <=> H+ + S--
        r1 = K['H2S'] * g_H2S / (g_HS * a_h)              # m_HS / m_H2S
        r2 = K['HS'] * g_HS / (g2 * a_h)                  # m_S / m_HS
        h2s_aq = m_H2S_total / (1.0 + r1 + r1*r2)
        hs_ion = h2s_aq * r1
        s_ion = hs_ion * r2

        # AMMONIA: NH4+ <=> NH3 + H+ ; Ka = a_NH3 a_H / a_NH4
        r_nh4 = (g_NH3 * a_h) / ((1.0 / K['NH3']) * g_NH4)   # m_NH4 / m_NH3

        # CARBONATE: CO2 + H2O <=> H+ + HCO3- ; HCO3- <=> H+ + CO3--
        c1 = K['CO2'] * g_CO2 / (g_HCO3 * a_h)
        c2 = K['HCO3'] * g_HCO3 / (g2 * a_h)
        # CARBAMATE: NH2COO- + H2O <=> NH3 + HCO3- ; K5 = a_NH3 a_HCO3 / (a_w a_NH2COO)
        # -> m_carb = kc * m_NH3 * m_HCO3
        kc = g_NH3 * g_HCO3 / (K['CARB'] * a_w * g_CARB)
        # Coupled NH3 / CO2 balances -> quadratic in free NH3 (n):
        #   n(1+r) + q n CO2tot/(B + q n) = NH3tot,  q = kc*c1, B = 1 + c1 + c1*c2
        B = 1.0 + c1 + c1*c2
        q = kc * c1
        a2 = (1.0 + r_nh4) * q
        a1 = (1.0 + r_nh4) * B + q * (m_CO2_total - m_NH3_total)
        a0 = -m_NH3_total * B
        if a2 > 1e-300 and m_NH3_total > 0:
            nh3_aq = (-a1 + math.sqrt(a1*a1 - 4.0*a2*a0)) / (2.0*a2)
        else:
            nh3_aq = m_NH3_total / (1.0 + r_nh4)
        co2_aq = m_CO2_total / (B + q * nh3_aq)
        hco3_ion = co2_aq * c1
        co3_ion = hco3_ion * c2
        carb_ion = kc * nh3_aq * hco3_ion
        nh4_ion = nh3_aq * r_nh4

        # WATER: Kw = a_H a_OH
        oh_ion = K['Kw'] / (g1 * a_h)
        h_ion = a_h / g1

        ionic_strength = 0.5 * (h_ion + oh_ion + nh4_ion + hs_ion + hco3_ion + carb_ion + 4.0*s_ion + 4.0*co3_ion)
        a_w = math.exp(-0.018015 * (nh3_aq + nh4_ion + h2s_aq + hs_ion + s_ion + co2_aq + hco3_ion
                                    + co3_ion + carb_ion + h_ion + oh_ion))

    error = (h_ion + nh4_ion) - (oh_ion + hs_ion + 2.0*s_ion + hco3_ion + 2.0*co3_ion + carb_ion)
    spec = {'NH3': nh3_aq, 'NH4': nh4_ion, 'H2S': h2s_aq, 'HS': hs_ion, 'S': s_ion,
            'CO2': co2_aq, 'HCO3': hco3_ion, 'CO3': co3_ion, 'NH2COO': carb_ion, 'H': h_ion, 'OH': oh_ion,
            'I': ionic_strength, 'g_NH3': g_NH3, 'g_H2S': g_H2S, 'g_CO2': G['g_CO2'], 'g_pm': G['g_pm'], 'g_dav': g1}
    return error, spec

def literature_K(T_K):
    """Replacements for Edwards-form constants that were ~0.46 pK off at 25 C.
    NH4+ : Bates & Pinching (1949)   pKa = 0.09018 + 2729.92/T
    CO2/HCO3- : Plummer & Busenberg (1982)"""
    lT = math.log10(T_K)
    pKa_nh4 = 0.09018 + 2729.92 / T_K
    logK1 = -356.3094 - 0.06091964*T_K + 21834.37/T_K + 126.8339*lT - 1684915.0/T_K**2
    logK2 = -107.8871 - 0.03252849*T_K + 5151.79/T_K + 38.92561*lT - 563713.9/T_K**2
    # Carbamate hydrolysis NH2COO- + H2O <=> NH3 + HCO3-, Edwards et al. (1978) as tabulated by later
    # authors (ln K = -2900/T + 8.6, molality scale). Sign set so that carbamate is less stable at
    # higher temperature; checked against Kurz et al. (1995) data (validation/validate_co2.py).
    lnK_carb = CARB_LNK[0] / T_K + CARB_LNK[1]
    return {'NH3': 10**pKa_nh4, 'CO2': 10**logK1, 'HCO3': 10**logK2, 'CARB': math.exp(lnK_carb)}

CARB_LNK = [-2900.0, 8.6]

def henry_co2(T_K):
    """CO2 Henry volatility, atm*kg/mol, from Carroll, Slupsky & Mather (1991), J. Phys. Chem. Ref.
    Data 20, 1201: ln(H/MPa) = -6.8346 + 1.2817e4/T - 3.7668e6/T^2 + 2.997e8/T^3 (mole-fraction
    basis, 273-433 K), converted to molality."""
    H_MPa = math.exp(-6.8346 + 1.2817e4 / T_K - 3.7668e6 / T_K**2 + 2.997e8 / T_K**3)
    return H_MPa / 0.101325 * CONST['MW']['H2O'] / 1000.0

def henry_h2s(T_K):
    """H2S Henry volatility, atm*kg/mol, from the IUPAC-evaluated solubility of H2S in water at
    p(H2S) = 1 atm: ln x = -24.912 + 3477.1/T + 0.3993 ln T + 0.015700 T
    (Fogg & Young, IUPAC Solubility Data Series vol. 32, 1988). Replaces a constant-enthalpy
    van 't Hoff form that was ~1.7x too volatile at 393 K."""
    x = math.exp(-24.912 + 3477.1 / T_K + 0.3993 * math.log(T_K) + 0.015700 * T_K)
    m = x / ((1.0 - x) * CONST['MW']['H2O'] / 1000.0)
    return 1.0 / m

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
    p_h2s = (henry_h2s(T_K) * h2s * spec['g_H2S']) * CONST['ATM_TO_PSI']
    p_co2 = (henry_co2(T_K) * co2 * spec['g_CO2']) * CONST['ATM_TO_PSI']
    
    # Water Activity (Osmotic Depression)
    # Ideal (phi = 1) osmotic approximation
    sum_moles = sum([v for k,v in spec.items() if not (k == 'I' or k.startswith('g_'))])
    a_w = math.exp(-0.018015 * sum_moles)
    
    p_w_pure = (10**(5.20389 - 1733.926/(T_K - 39.485))) * 14.5038
    p_w = p_w_pure * a_w
    
    p_tot_pa = (p_nh3 + p_h2s + p_co2 + p_w) * CONST['PSI_TO_PA']
    
    return {'ph': ph_guess, 'P_bubble': p_tot_pa, 'pp': {'NH3': p_nh3, 'H2S': p_h2s, 'CO2': p_co2, 'H2O': p_w}, 'liq': spec, 'aw': a_w}

def calculate_density(T_K, P_Pa, spec, kg_water, mass_total_kg):
    t_c = T_K - 273.15
    rho_w = 1000.0 * (1 - (t_c + 288.9414)/(508929.2 * (t_c + 68.12963)) * (t_c - 3.9863)**2)
    v_inf = {'NH4': 18.0, 'HS': 20.0, 'NH3': 24.5, 'H2S': 35.0, 'CO2': 34.0, 'HCO3': 25.0, 'CO3': 20.0, 'S': 20.0, 'NH2COO': 30.0}
    vol_water_L = kg_water / (rho_w / 1000.0)
    vol_solutes_mL = (
        spec['NH4']*v_inf['NH4'] + spec['HS']*v_inf['HS'] + spec['NH3']*v_inf['NH3'] + spec['H2S']*v_inf['H2S'] +
        spec['CO2']*v_inf['CO2'] + spec['HCO3']*v_inf['HCO3'] + spec['CO3']*v_inf['CO3'] + spec['S']*v_inf['S'] + spec['NH2COO']*v_inf['NH2COO']
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

    lines = [HR, " SWEQ - SOUR WATER EQUILIBRIUM SOLVER v8.4.0 ".center(W), " Edwards / Pitzer, NH3-H2S-CO2-H2O with carbamate ".center(W), HR]
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
    lines.append("  CARBONATE SYSTEM ::  CO2(aq): %7.4f          <==> HCO3: %7.4f  CO3: %7.4f" % (l['CO2'], l['HCO3'], l['CO3']))
    lines.append("  CARBAMATE        ::  NH2COO-: %7.4f" % l['NH2COO'])
    
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
        beta = min(_rachford_rice(z, K), 1.0 - 1e-9)  # keep a trace of liquid while iterating
        x = {k: z[k]/(1.0+beta*(K[k]-1.0)) for k in KEYS}
        sx = sum(x.values()); x = {k: v/sx for k, v in x.items()}
        L = {k: F*(1.0-beta)*x[k] for k in KEYS}
        res = bubble_from_moles(T_K, L)
        K_new = k_values(res, x)
        d = max(abs(math.log(K_new[k]/K[k])) for k in KEYS)
        K = {k: math.exp(0.5*math.log(K[k]) + 0.5*math.log(K_new[k])) for k in KEYS}
        if d < 1e-8: break
    if beta > 1.0 - 1e-6:
        # Above the dew point: single vapor phase with the feed composition
        return F, dict(z), {k: 0.0 for k in KEYS}, None
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
        x_full[i] = L[k] / Lt if Lt > 0 else n[k] / F_mol
        if y: y_full[i] = y[k]

    liq = globals().get('oms1')
    vap = globals().get('oms2')
    if liq is None:
        raise ValueError("SWEQ: connect a material stream to OUTLET 1 (liquid) of the Python Script")
    if vap is None and V > 0:
        raise ValueError("SWEQ: %.1f %% of the feed vaporizes at these conditions; connect a material "
                         "stream to OUTLET 2 (vapor) of the Python Script" % (100.0 * V / F_mol))

    liq.SetTemperature(T); liq.SetPressure(P_op)
    liq.SetOverallComposition(Array[Double](x_full))
    liq.SetMolarFlow(Lt)
    liq.Calculate()

    if vap is not None:
        vap.SetTemperature(T); vap.SetPressure(P_op)
        vap.SetOverallComposition(Array[Double](y_full if y else x_full))
        vap.SetMolarFlow(V if V > 0 else 0.0)
        vap.Calculate()

    if res is None:  # all vapor: no liquid datasheet to report
        return
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
