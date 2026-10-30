# NEIVATranslation: MCE-dependent decks from NEIVAv1.0

NEIVAv1.0 (Shahid et al., GMD 17, 7679, 2024; please cite) Integrated emission factors -> ASP decks as a function of modified combustion efficiency.

    python3 build_neiva_decks.py --type "temperate forest" --mce 0.80 0.975 0.025 --template Williams
    python3 coverage.py                        # route/mass coverage of the fire type
    python3 plot_mce_enhancement.py <runs_root> out.png

Files: `neiva_data.py` (CSV/SQL reader), `neiva_map.py` (compound -> deck species), `mce_fit.py` (EF vs MCE), `build_neiva_decks.py`, `coverage.py`, `plot_mce_enhancement.py`.

Method
* CO, CO2 from the carbon balance: dC = dCO/(1-MCE), dCO2 = dCO MCE/(1-MCE); dCO fixed at the template value.
* Other species: carbon-normalised EF y (g/gC) per study column, linear fit vs MCE (log-linear for OA, BC) if n>=10, span>=0.05, p<0.10 and leave-one-out beats the mean; otherwise the mean y (no MCE dependence). OA is not significant for temperate forest and uses the mean.
* Lumped NEIVA rows are apportioned to member compounds. Routes: explicit deck species; baseline lumps (carbon-conserving); other identified NMOG -> IVOC bins by log10 C* (mass-conserving; C* >= 1e6 ug/m3 all land in IVOC9); halogenated/inorganic skipped.
* Template's unidentified-NMOG pool (gas IVOC3-7) is kept as a constant multiple of NEIVA identified NMOG (`--no-unid-scale` to keep it per unit dCO instead).
* Species without NEIVA data (H2, CH4, PAN, PPN, HCN, H2O2, HNO3, 1-propanol) keep the template value. Salts NH4Cl/NH4NO3/KCl/K2SO4 from the NEIVA ions; Nt scales with aerosol mass.
* Temperate-forest data cover MCE 0.781-0.968; the 0.975 deck is an extrapolation (BC especially). Each deck has NEIVA_PROVENANCE.txt and NEIVA_SPECIES_MAPPING.csv.
