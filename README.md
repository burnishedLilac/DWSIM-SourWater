# SWEQ - Sour Water Equilibrium Solver

**Version:** 8.2.0  
**Integration:** DWSIM Python Script Unit Operation  
**License:** GPLv3  

SWEQ is a Python-based thermodynamic solver for calculating the vapor-liquid equilibrium (VLE) and ionic speciation of sour water systems ($NH_3 - H_2S - CO_2 - H_2O$). It is designed to be integrated directly into the DWSIM process simulator.

## Features

* **Thermodynamic Framework:** Based on the weak electrolyte model proposed by Edwards, Newman & Prausnitz (1978).
* **High Concentration Support:** Replaces the standard Davies equation with the Pitzer Activity Coefficient Model (calibrated via Rumpf et al., 1999). This allows the solver to converge in highly concentrated brines (ionic strengths > 6.0 molal) typical of stripper bottoms and atmospheric overheads.
* **Hybrid Activity Model:** * Uses Pitzer parameters for major $NH_4^+/HS^-$ interactions and "salting-out" effects on neutral gases ($\gamma_{NH_3}$, $\gamma_{H_2S}$).
  * Uses the Davies equation as a fallback for minor carbonate species ($HCO_3^-, CO_3^{2-}$).
* **Water activity:** ideal osmotic approximation ($\phi = 1$); full Pitzer $\phi$ is on the roadmap.
* **PT-Flash:** isothermal flash at the inlet stream pressure (Rachford-Rice, K-values from the electrolyte model).
* **Numerical Stability:** The charge-balance bisection loop includes a successive substitution damping factor ($w = 0.5$) to prevent non-linear oscillation during the calculation of activity coefficients.
* **Automated Output:** Generates a formatted text report (`SWEQ_Datasheet.txt`) containing phase status, calculated pH, bubble pressure, environmental metrics (mg/L), and specific activity coefficients.

## Usage in DWSIM

1. Add a **Python Script** Unit Operation to your DWSIM flowsheet.
2. Connect one inlet stream (`ims1`) and two outlet streams (liquid `oms1`, vapor `oms2`).
3. Paste the contents of `SWEQ_Solver.py` into the script editor.
4. Run the flowsheet. The script maps standard components (Ammonia, Hydrogen Sulfide, Carbon Dioxide, and Water), performs the flash calculation, updates the streams, and exports the datasheet to your Desktop or `C:\Temp\`. Errors are reported in the DWSIM log.

## Validation

`python tests/validate.py` runs limiting-case checks (Henry constants, pH of dilute NH3/H2S/CO2, water vapor pressure, flash mass balance). Validation against multicomponent sour-water VLE data (API SWEQ, Rumpf 1999, Kurz 1995) is **still pending**; until then, treat results at high concentration as qualitative.

## Known limitations

* NH4+/HS- Pitzer parameters are constant (no temperature dependence) and have no cited source; no NH3-salt or H2S-salt interaction parameters.
* Carbonate species use the Davies equation (not valid above I ~ 0.5 molal).
* Ideal vapor phase; no enthalpy balance (isothermal only).

## References

* Edwards, T. J., Maurer, G., Newman, J., & Prausnitz, J. M. (1978). *Vapor‐liquid equilibria in multicomponent aqueous solutions of volatile weak electrolytes*. AIChE Journal.
* Bates, R. G. & Pinching, G. D. (1949). J. Res. Natl. Bur. Stand. 42, 419 (NH4+ dissociation).
* Plummer, L. N. & Busenberg, E. (1982). Geochim. Cosmochim. Acta 46, 1011 (CO2/HCO3- constants).
* Rumpf, B., Pérez-Salado Kamps, Á., Sing, R., & Maurer, G. (1999). *Simultaneous solubility of ammonia and hydrogen sulfide in water at temperatures from 313 K to 393 K*. Fluid Phase Equilibria.
