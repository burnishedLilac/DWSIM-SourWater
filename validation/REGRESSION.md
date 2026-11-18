# Pitzer parameter regression — first attempt (not adopted)

Script: `validation/regress.py`. Fit set: Oratovskii et al. 1964 (343–363 K, NH3 ≤ ~7 mol/kg).
Hold-out: Leyko 1959 (293–323 K, NH3 up to ~19 mol/kg, density assumed). AAD of partial pressures, %.

| Variant | fit p(H2S) | fit p(NH3) | hold-out p(H2S) | hold-out p(NH3) |
|---|---|---|---|---|
| v8.2.0 defaults (no regression) | 21.5 | 26.3 | 41.0 | 25.3 |
| 7 params, β0(T), β1, λ_NH3(T), λ_H2S(T) | 5.8 | 8.6 | 234 | 201 |
| 4 params, β0(T), λ_NH3, λ_H2S | 7.7 | 7.8 | 397 | 135 |
| 4 params, β0, β1, λ_NH3, λ_H2S (no T) | 10.1 | 9.3 | 29.5 | 129 |
| 3 params, β0, λ_NH3, λ_H2S | 13.6 | 8.7 | 17.1 | 120 |
| 2 params, λ_NH3, λ_H2S | 12.9 | 9.8 | 26.3 | 108 |

## Conclusion

Every variant fits the fit set well and fails the hold-out, especially p(NH3). The hold-out is an
extrapolation in both temperature and concentration, and its NH3 concentrations were back-calculated
by the IUPAC compiler, so it cannot discriminate between "overfitting" and "bad hold-out data".
**No regressed parameters are adopted**; the solver keeps the v8.2.0 defaults
(ion–molecule λ terms are implemented but set to 0).

Next: add Rumpf et al. 1999 (NH3–H2S–H2O, 313–393 K, measured by the same group that proposed the
model) and re-run with a temperature-stratified hold-out.
