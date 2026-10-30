# MOZART-T1 and GEOS-Chem translations (ASP hybrid decks)

* `t1_kpp.py`, `mozt1.py`, `build_t1_decks.py`, `verify_t1.py`: WRF-Chem KPP T1_MOZCART (tag MZ148_T1mozcart_20170323) -> ASP.
  `python3 build_t1_decks.py` writes `InorganicInputDecks/<Case>_MOZT1/`; `python3 verify_t1.py` checks every translated rate
  against the KPP Fortran functions (T 220-320 K, P 200-1013 mb).
* `hybrid.py`: generic hybrid builder (same policy as the MCM decks: a baseline reaction is kept iff one reactant is not in the
  new mechanism).  `aspfmt.py`: reaction-line formatting.
* This T1 is the 2017 WRF-Chem version (142 species, 344 reactions), not the later CAM-chem T1.1 / T1.2.
* Rates are structural (no T or P frozen).  T1 aerosol-uptake reactions (usr_*_aer) use ASP's first-order heterogeneous flag 3
  (gamma, MW) so they follow ASP's own particle surface area; T1's fixed Dg = 0.1 cm2/s and prescribed sulfate/OC/BC areas are not used.
* T1 species without a baseline row use approximate molecular weights (`build_t1_decks.MW`; lumped MOZART species are approximate).
* T1's VBS SOA yield function is not used; ASP's SOA counters are re-injected on first-generation VOC reactions.
