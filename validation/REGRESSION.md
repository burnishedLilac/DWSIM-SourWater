# Pitzer parameter regression

Script: `validation/regress.py`. Objective: least squares on ln(p_calc/p_exp) for p(NH3) and p(H2S);
points with p ≤ 0.01 bar are skipped (below or near the detection limit in Rumpf 1999).

## Attempt 2 — Rumpf et al. 1999 (adopted in v8.3.0)

**Physics fix first.** Before fitting, the H2S Henry constant was replaced by the IUPAC SDS-32
correlation. The previous constant-enthalpy form was ~1.27× too volatile at 353 K and ~1.7× at 393 K,
which showed up as a p(H2S) error growing with temperature (Rumpf AAD 81.6 % → 34.6 % with this fix alone).

**Split.** Fit: Rumpf 313 K + 393 K. Hold-out: Rumpf 353 K (interpolation in T). External: Oratovskii 1964, Leyko 1959.

AAD, % — p(NH3) / p(H2S):

| Variant | fit 313+393 K | hold-out 353 K | Oratovskii | Leyko |
|---|---|---|---|---|
| defaults after Henry fix (no regression) | 11.9 / 29.4 | 19.5 / 43.6 | 26.3 / 16.6 | 25.3 / 30.6 |
| **β1, λ_H2S, λ_NH3 (adopted)** | **5.8 / 14.7** | **17.4 / 18.6** | 35.6 / 27.7 | **12.7 / 10.9** |
| β0, β1, λ_H2S, λ_NH3 | 6.4 / 13.0 | 19.1 / 16.0 | 38.2 / 30.3 | 18.8 / 33.0 |
| β1, λ_H2S(T), λ_NH3 | 5.9 / 14.5 | 17.6 / 17.9 | 35.5 / 27.8 | 12.5 / 10.9 |
| β0(T), β1, λ_H2S(T), λ_NH3 | 6.2 / 12.7 | 20.2 / 16.5 | 38.5 / 29.3 | 19.8 / 34.2 |
| λ_H2S, λ_NH3 | 6.2 / 18.6 | 20.1 / 28.5 | 38.0 / 23.1 | 18.7 / 22.3 |

Adopted (no temperature dependence; β0 = 0.052 and Cφ = −0.001 unchanged):
β1(NH4HS) = 0.0178, λ(H2S–NH4HS) = −0.0184, λ(NH3–NH4HS) = −0.0181.
The 3-parameter set was chosen because it is the simplest model that improves both the hold-out
isotherm and the independent Leyko data; extra temperature terms did not improve the hold-out.

**Trade-off.** Agreement with Oratovskii 1964 gets worse (both pressures under-predicted by
~25–35 %). Those pressures are estimates made by the IUPAC compiler, and they cannot be reconciled
with Rumpf 1999 within this model. Rumpf is treated as the reference.

**Why not use the parameters published by Rumpf et al. (Table 2)?** They were fitted together with
their own equilibrium constants, Henry constants and a virial-EoS vapor phase. Inserted into SWEQ's
framework they made p(H2S) worse (AAD ~67 % vs 35 %), so they are not transferable.

## Attempt 1 — Oratovskii 1964 fit, Leyko hold-out (not adopted)

All 2–7-parameter variants fitted Oratovskii to 6–14 % but failed on Leyko (>100 % in p(NH3)); the
hold-out was an extrapolation in T and concentration with back-calculated compositions.

## Attempt 3 — NH3–CO2–H2O, Kurz et al. 1995 (adopted in v8.4.0)

**Model fixes first.** With carbamate added but Davies activity coefficients, p(CO2) was over-predicted ~8×:
Davies gave γ(HCO3−) of 2–6 at I = 2–7 mol/kg while γ(NH4+) was 0.5. Replacing Davies by the Pitzer
Debye–Hückel term for all ions, and placing the ion–molecule λ on the anion (Kurz convention), brought the
errors to p(NH3) 22 % / p(CO2) 15.5 % without any fitting. The remaining p(NH3) error grew with CO2 loading
(−3 % without CO2, up to −50 % near loading 0.5): a missing salting-out of NH3 by the ammonium salts.

**Split.** Fit: 313 + 353 K. Hold-out: 333 K. Two-phase points only; partial pressures ≤ 2 kPa skipped.

| Variant | hold-out p(NH3) | hold-out p(CO2) | hold-out p_total | parameters |
|---|---|---|---|---|
| no regression | 22.8 % | 16.7 % | 12.1 % | — |
| **λ(NH3–NH2COO−) (adopted)** | **12.8 %** | **18.1 %** | **7.1 %** | +0.0528 ± 0.0076 |
| λ(NH3–HCO3−), λ(NH3–NH2COO−) | 13.2 % | 12.5 % | 6.7 % | +0.141 ± 0.046, +0.042 ± 0.006 |
| + λ(CO2–HCO3−) | 13.1 % | 14.6 % | 8.0 % | λ(CO2–HCO3−) not identifiable (± 50) |

The one-parameter variant was adopted: it is well determined, and the second parameter does not improve the
hold-out consistently. Independent check: Kurz et al. report β(NH3, NH2COO−) = −0.03933 + 25.263/T, i.e. +0.037 at
333 K in their own framework — same sign and magnitude as the value regressed here.
