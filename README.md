# SWEQ — Sour Water Equilibrium Solver for DWSIM

[![tests](https://github.com/burnishedLilac/DWSIM-SourWater/actions/workflows/tests.yml/badge.svg)](https://github.com/burnishedLilac/DWSIM-SourWater/actions/workflows/tests.yml)

**Version:** 8.2.0 · **License:** GPLv3 · **Runs in:** DWSIM Python Script Unit Operation (IronPython/Python.NET) or plain Python 3

SWEQ calculates vapor–liquid equilibrium and ionic speciation of sour water
($NH_3 - H_2S - CO_2 - H_2O$) and performs an isothermal flash inside a DWSIM flowsheet.

## Repository layout

```
dwsim/SWEQ_Solver.py        single-file solver — paste this into DWSIM
tests/test_limits.py        limiting-case checks (Henry constants, dilute pH, water vapor pressure, flash balance)
validation/validate.py      comparison against published experimental VLE data
validation/data/*.csv       experimental data with full source citation in each header
validation/RESULTS.md       latest validation report (generated)
CHANGELOG.md · ROADMAP.md
```

## Model

| Part | Implementation |
|---|---|
| Chemical equilibria | NH4+ (Bates & Pinching 1949), CO2/HCO3− (Plummer & Busenberg 1982), H2S/HS−/Kw (Edwards-form correlations) |
| Activity, NH4+/HS− and neutral NH3/H2S | Pitzer (Debye–Hückel + β0, β1, Cφ); neutral self-interaction from Rumpf et al. |
| Activity, carbonate and S2− | Davies (fallback, valid to I ≈ 0.5 mol/kg) |
| Vapor | Henry's law with activity coupling (ideal gas); water by Antoine × ideal osmotic activity |
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

Current results (v8.2.0, NH3–H2S–H2O, 293–363 K), average absolute deviation of partial pressures:

| Data set | n | p(H2S) | p(NH3) |
|---|---|---|---|
| Oratovskii et al. 1964 (wt%, exact conversion) | 28 | 21.5 % | 26.3 % |
| Leyko 1959 (mol/L, density assumed 1.0) | 14 | 41.0 % | 25.3 % |
| *v8.1.0 on Oratovskii, for reference* | 28 | 99.6 % | 30.5 % |

**Read this before using the numbers:** the model is reasonable at moderate concentration and drifts at
high concentration (H2S over-predicted, NH3 under-predicted). Treat results above ~3 mol/kg NH3 as
qualitative until the Pitzer parameters are regressed (see ROADMAP). See `validation/RESULTS.md` for every point.

## Known limitations

* NH4+/HS− Pitzer parameters are constant (no temperature dependence) and have no cited source; ion–molecule
  parameters (NH3–NH4HS, H2S–NH4HS) are missing.
* Carbonate uses Davies; no Pitzer CO2 network yet.
* Ideal vapor phase; isothermal only (no energy balance); no NH4HS solid.

## References

* Edwards, Maurer, Newman & Prausnitz (1978), AIChE J. 24, 966.
* Bates & Pinching (1949), J. Res. Natl. Bur. Stand. 42, 419.
* Plummer & Busenberg (1982), Geochim. Cosmochim. Acta 46, 1011.
* Rumpf, Pérez-Salado Kamps, Sing & Maurer (1999), Fluid Phase Equilib. 158–160, 923.
* Fogg & Young (eds.) (1988), IUPAC Solubility Data Series vol. 32, *Hydrogen Sulfide, Deuterium Sulfide and Hydrogen Selenide*.
