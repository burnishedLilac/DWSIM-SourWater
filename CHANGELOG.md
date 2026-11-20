# Changelog

## [8.4.0] — 2026-11-20

### Added
- Carbamate (NH2COO−) chemistry: NH2COO− + H2O ⇌ NH3 + HCO3−, Edwards et al. (1978) constant. The sign of its
  temperature dependence (lost in the secondary source) was fixed by thermodynamics and confirmed with the data.
- Kurz, Rumpf & Maurer (1995) NH3–CO2–H2O data (63 points, 313–353 K) and `validation/validate_co2.py`.
- λ(NH3–NH2COO−) = 0.0528 regressed to Kurz 313 + 353 K; 333 K kept as hold-out (`validation/regress_co2.py`).

### Fixed
- CO2 Henry constant now from Carroll, Slupsky & Mather (1991).
- Davies equation replaced by the Pitzer Debye–Hückel term for all ions: Davies gave γ(HCO3−) = 2–6 at high ionic
  strength and over-predicted p(CO2) ~8×. Ionic strength now includes every ion.
- Ion–molecule λ terms placed on the anion (Kurz et al. convention), so NH4+ paired with HCO3−/NH2COO− no longer
  inherits the NH4HS salting effect. NH3–H2S results are unchanged.

### Validation
- New independent NH3–H2S–CO2–H2O test (`validation/validate_ternary.py`): Wilson et al. 1985 (GPA), Badger & Silver
  1938, van Krevelen et al. 1949. At 366 K (stripper conditions): p(NH3) 9 %, p(CO2) 21 %, p(H2S) 18 %.
- Kurz 1995 (all VLE points): p(NH3) 62.5 % → 11.0 %, p(CO2) 1242 % → 18.6 %, p_total 305 % → 8.0 %.
- Rumpf 1999 (NH3–H2S): unchanged (p(H2S) 16.1 %, p(NH3) 10.7 %).

## [8.3.1] — 2026-11-18

### Added
- First run inside DWSIM confirmed. Bubble-point comparison with DWSIM's built-in Sour Water package (`validation/compare_dwsim_native.py`, `validation/data/dwsim_native_sourwater.csv`).

### Fixed
- A feed above its dew point raised "feed fully vaporizes" and stopped the flowsheet; it now goes entirely to the vapor outlet (`oms2`).
- Outlet 2 (vapor) is optional when nothing vaporizes; a clear message asks for it when it is needed (was `NameError: oms2`).

## [8.3.0] — 2026-11-18

### Fixed
- H2S Henry constant: temperature dependence from the IUPAC SDS-32 correlation; the previous form was ~1.7× too volatile at 393 K.
- v8.2.0 release was cut before the repository restructure was merged; this release contains it.

### Added
- Experimental data of Rumpf et al. (1999), 64 points, 313–393 K, as the reference validation set.
- Ion–molecule Pitzer terms λ(NH3–NH4HS), λ(H2S–NH4HS) and the NH3 ternary term.
- `validation/regress.py` (fit/hold-out regression) and `validation/REGRESSION.md`.

### Changed
- β1(NH4HS), λ(NH3–NH4HS), λ(H2S–NH4HS) regressed to Rumpf 1999 (313 + 393 K); 353 K kept as hold-out.
  Rumpf p(H2S) AAD 81.6 % → 16.1 %, p(NH3) 15.1 % → 10.7 %.

## [8.2.0] — 2026-11-18

### Fixed
- Henry constants for H2S and CO2 were solubilities (mol·kg⁻¹·atm⁻¹) used as volatilities; partial pressures were ~100× (H2S) and ~900× (CO2) too low.
- Activity coefficients applied consistently (pH = H+ activity, γ used once per species; γ(NH4+) added).
- NH4+ and CO2/HCO3− equilibrium constants were ~0.46 pK units off at 25 °C; replaced by Bates & Pinching and Plummer & Busenberg.
- Outlet streams: `oms1` is liquid and `oms2` vapor, as documented.
- Errors are no longer silently swallowed.

### Added
- Isothermal PT-flash (Rachford–Rice).
- `tests/test_limits.py` (limiting cases) and `validation/` (comparison with IUPAC-compiled experimental data).
- GitHub Actions workflow running tests and validation.

### Changed
- Repository reorganized: solver moved to `dwsim/SWEQ_Solver.py`.
- Datasheet no longer opens automatically (`OPEN_REPORT` flag).

## [8.1.0] and earlier
- Edwards + Pitzer (NH4+/HS−) model with Davies fallback for carbonate; bubble-point calculation.
