# SWEQ — Sour Water Equilibrium Solver for DWSIM

[![tests](https://github.com/burnishedLilac/DWSIM-SourWater/actions/workflows/tests.yml/badge.svg)](https://github.com/burnishedLilac/DWSIM-SourWater/actions/workflows/tests.yml)

**Version:** 8.4.0 · **License:** GPLv3 · **Runs in:** DWSIM Python Script Unit Operation (IronPython/Python.NET) or plain Python 3

SWEQ calculates vapor–liquid equilibrium and ionic speciation of sour water
($NH_3 - H_2S - CO_2 - H_2O$) and performs an isothermal flash inside a DWSIM flowsheet.

## Repository layout

```
dwsim/SWEQ_Solver.py        single-file solver — paste this into DWSIM
tests/test_limits.py        limiting-case checks (Henry constants, dilute pH, water vapor pressure, flash balance)
validation/validate.py      comparison against published experimental VLE data
validation/data/*.csv       experimental data with full source citation in each header
validation/RESULTS.md       latest validation report (generated)
validation/compare_dwsim_native.py   SWEQ vs DWSIM built-in Sour Water package
validation/validate_co2.py   NH3–CO2–H2O validation (Kurz 1995)
validation/validate_ternary.py   NH3–H2S–CO2–H2O independent test (Wilson 1985, Badger & Silver 1938, van Krevelen 1949)
CHANGELOG.md · ROADMAP.md
```

## Model

| Part | Implementation |
|---|---|
| Chemical equilibria | NH4+ (Bates & Pinching 1949), CO2/HCO3− (Plummer & Busenberg 1982), carbamate NH2COO− (Edwards et al. 1978), H2S/HS−/Kw (Edwards-form correlations) |
| Activity, NH4+/HS− and neutral NH3/H2S | Pitzer (Debye–Hückel + β0, β1, Cφ); neutral self-interaction from Rumpf et al. |
| Activity, all ions | Pitzer Debye–Hückel term with total ionic strength; NH4+/HS− binary terms |
| Activity, ion–molecule | λ(NH3–NH4HS), λ(H2S–NH4HS), β1(NH4HS) regressed to Rumpf et al. 1999; λ(NH3–NH2COO−) regressed to Kurz et al. 1995 (`validation/REGRESSION.md`) |
| Vapor | Henry's law with activity coupling (ideal gas); H2S Henry from IUPAC SDS-32; CO2 Henry from Carroll et al. 1991; water by Antoine × ideal osmotic activity |
| Flash | PT-flash, Rachford–Rice with K-values from the electrolyte model |

Convention: pH is the H+ **activity**; balances are in molality; each equilibrium uses each species' activity coefficient once.

## Usage in DWSIM

1. Add a **Python Script** Unit Operation.
2. Connect one inlet (`ims1`) and two outlets: liquid `oms1`, vapor `oms2`.
3. Paste `dwsim/SWEQ_Solver.py` into the script editor and run.
4. Feed must contain only Water, Ammonia, Hydrogen Sulfide and Carbon Dioxide. The flash runs at the inlet T and P;
   a datasheet (`SWEQ_Datasheet.txt`) is written to the Desktop (or `C:\Temp`). Errors are shown in the DWSIM log.

## Testing and validation

```bash
python tests/test_limits.py      # or: pytest tests/
python validation/validate.py    # regenerates validation/RESULTS.md
```

Current results (v8.4.0), NH3–H2S–H2O (293–393 K), average absolute deviation of partial pressures:

| Data set | n | p(H2S) | p(NH3) | role |
|---|---|---|---|---|
| Rumpf et al. 1999, 353 K isotherm | 23 | 18.6 % | 17.4 % | hold-out |
| Rumpf et al. 1999, all isotherms | 64 | 16.1 % | 10.7 % | 313/393 K used in fit |
| Leyko 1959 | 14 | 10.9 % | 12.7 % | independent |
| Oratovskii et al. 1964 | 28 | 27.7 % | 35.6 % | independent (compiler-estimated pressures) |
| *v8.2.0 on Rumpf 1999, for reference* | 64 | 81.6 % | 15.1 % | |

NH3–CO2–H2O (Kurz et al. 1995, 313–353 K, two-phase points), average absolute deviation:

| Data set | p(NH3) | p(CO2) | p_total | role |
|---|---|---|---|---|
| Kurz 1995, 333 K isotherm | 12.8 % | 18.1 % | 7.1 % | hold-out |
| Kurz 1995, all isotherms | 11.0 % | 18.6 % | 8.0 % | 313/353 K used in fit |
| *v8.3.1 (no carbamate), for reference* | 62.5 % | 1242 % | 305 % | |

`python validation/validate_co2.py` reproduces these numbers (`validation/RESULTS_CO2.md`).

NH3–H2S–CO2–H2O, independent test (no parameter fitted to these data; IUPAC SDS-32 compilation):

| Data set | n | p(NH3) | p(CO2) | p(H2S) |
|---|---|---|---|---|
| Wilson, Gillespie & Owens 1985 (GPA), 366 K | 3 | 9.2 % | 20.6 % | 17.7 % |
| Wilson et al. 1985, 422 K (extrapolation) | 4 | 11.1 % | 16.4 % | 22.0 % |
| Badger & Silver 1938, 293 K | 12–13 | 17.2 % | 29.7 % | 53.8 % |
| van Krevelen et al. 1949, 293–333 K | 27 | – | – | 23.4 % |

The older data (1938, 1949) are in mol/dm³ (density assumed 1.0 kg/L) and were obtained by chemical analysis; at
293 K p(H2S) is over-predicted. `python validation/validate_ternary.py` reproduces these numbers.

**Comparison with DWSIM's built-in Sour Water package** (bubble points, `validation/compare_dwsim_native.py`),
deviation from Rumpf 1999 measurements, p_total / p(NH3) / p(H2S):

| Condition | DWSIM built-in | SWEQ |
|---|---|---|
| 80 °C, 6.0 m NH3, 2.5 m H2S (SWEQ hold-out) | +1 % / −38 % / +39 % | −5 % / −16 % / +3 % |
| 120 °C, 3.2 m NH3, 1.1 m H2S (in SWEQ fit) | +12 % / −10 % / +58 % | −2 % / +1 % / +1 % |

| 80 °C, same point, DWSIM **NRTL** (no electrolyte chemistry) | +364 % / +652 % / +366 % | |

Only two measured points so far. The built-in Sour Water package is part of DWSIM's Patreon edition and over-predicts
p(H2S) in all five conditions tested; a non-electrolyte model such as NRTL ignores NH3 + H2S ⇌ NH4⁺ + HS⁻ and
over-predicts all pressures 4–8×.

The data sets disagree with each other at high concentration (see `validation/REGRESSION.md`).
Treat results above ~6 mol/kg NH3 or above 393 K as extrapolation. Every point is listed in `validation/RESULTS.md`.

## Known limitations

* β0 and Cφ of NH4HS are uncited legacy values; regressed parameters have no temperature dependence.
* No solid phase: NH4HCO3 precipitation (Kurz et al. observed it at 313 K and high CO2 loading) is not modelled.
* The three-gas system combines parameters fitted to the two sub-systems; it is tested (not fitted) against the
  limited ternary data above. p(H2S) at 293 K is over-predicted by ~30–50 %.
* p(CO2) is under-predicted by ~15 % on average; no CO2–ion interaction parameter was identifiable from the data.
* Carbonate uses Davies; no Pitzer CO2 network yet.
* Ideal vapor phase; isothermal only (no energy balance); no NH4HS solid.

## References

* Edwards, Maurer, Newman & Prausnitz (1978), AIChE J. 24, 966.
* Bates & Pinching (1949), J. Res. Natl. Bur. Stand. 42, 419.
* Plummer & Busenberg (1982), Geochim. Cosmochim. Acta 46, 1011.
* Rumpf, Pérez-Salado Kamps, Sing & Maurer (1999), Fluid Phase Equilib. 158–160, 923 (reference data set).
* Kurz, Rumpf & Maurer (1995), Fluid Phase Equilib. 104, 261 (NH3–CO2–H2O reference data).
* Carroll, Slupsky & Mather (1991), J. Phys. Chem. Ref. Data 20, 1201 (CO2 Henry constant).
* Fogg & Young (eds.) (1988), IUPAC Solubility Data Series vol. 32, *Hydrogen Sulfide, Deuterium Sulfide and Hydrogen Selenide*.
