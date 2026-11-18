# Changelog

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
