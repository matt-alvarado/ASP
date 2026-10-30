"""Lumping of the paper compounds (keys of CB6r5Translation/species_map.MAP) onto GEOS-Chem 14.8 and MOZART-T1 species.

Mole-for-mole onto the species whose OH reactivity / structure is closest (same convention as the CRACMM / SAPRC decks); carbon is not
conserved where a lump is smaller than the compound.  Entries: key -> ([(species, coefficient)], note).  Keys not listed here
(CO, CO2, NO, NO2, O3, H2O2, HONO, SO2, HNO3, HCl, NH3, H2, CH4, N2O5, PAN, HCN, N2O, H2SO4) are direct.
"""
def _m(d):
    return {k: ([(s, c) for s, c in (v[0] if isinstance(v[0], list) else [(v[0], 1.0)])], v[1]) for k, v in d.items()}

GC = _m({
 'HCHO': ('CH2O', 'explicit'), 'MEOH': ('MOH', 'explicit'), 'FACD': ('HCOOH', 'explicit'), 'ETH': ('C2H4', 'explicit'),
 'ETHA': ('C2H6', 'explicit'), 'ETHY': ('C2H2', 'explicit'), 'ALD2': ('ALD2', 'explicit'), 'ETOH': ('EOH', 'explicit'),
 'AACD': ('ACTA', 'explicit'), 'GLY': ('GLYX', 'explicit'), 'GLYD': ('GLYC', 'explicit (glycolaldehyde)'),
 'MeFormate': ('C2H6', 'kOH 2.3e-13 ~ ethane 2.4e-13'), 'PRPA': ('C3H8', 'explicit'), 'Propene': ('PRPE', 'explicit (>=C3 alkenes)'),
 'Acrolein': ('ACR', 'explicit'), 'Propanal': ('RCHO', 'GEOS-Chem RCHO = C3 aldehyde'), 'Acetone': ('ACET', 'explicit'),
 'MGLY': ('MGLY', 'explicit'), 'Hydroxyacetone': ('HAC', 'explicit'), 'nPropanol': ('ROH', 'GEOS-Chem ROH = C3 alcohol'),
 'C3H6O': ('PRPE', 'Timbavati lumped alkene treated as propene'), 'Butadiene': ('C4H6', 'explicit'),
 '1-Butene': ('PRPE', '>=C3 alkenes'), 'iButene': ('PRPE', '>=C3 alkenes'), '2-Butene': ('PRPE', '>=C3 alkenes'),
 'C4H8': ('PRPE', '>=C3 alkenes'), 'Butane': ('ALK4', 'GEOS-Chem ALK4 = >=C4 alkanes'), 'Butanal': ('RCHO', 'higher aldehyde'),
 'MEK': ('MEK', 'explicit'), 'nButanol': ('ROH', 'C3+ alcohols'), 'Furan': ('FURA', 'explicit'),
 'Biacetyl': ([('MGLY', 1.33)], 'biacetyl: 1.33 MGLY (kOH-weighted)'), 'ISOP': ('ISOP', 'explicit'),
 'Pentane': ('ALK4', '>=C4 alkanes'), '1-Pentene': ('PRPE', '>=C3 alkenes'), '2-Pentene': ('PRPE', '>=C3 alkenes'),
 'Cyclopentanol': ('ROH', 'alcohol'), 'Pentanone': ('MEK', 'ketone >= C4'), 'Hexane': ('ALK4', '>=C4 alkanes'),
 '1-Hexene': ('PRPE', '>=C3 alkenes'), 'Heptane': ('ALK4', '>=C4 alkanes'), 'Heptanone': ('MEK', 'ketone >= C4'),
 'Octene': ('PRPE', '>=C3 alkenes'), 'Octanone': ('MEK', 'ketone >= C4'), 'MeAcetate': ('C2H6', 'kOH 3.4e-13 ~ ethane'),
 'Benzene': ('BENZ', 'explicit'), 'TOL': ('TOLU', 'explicit'), 'Xylene': ('XYLE', 'explicit'), 'Ethylbenzene': ('EBZ', 'explicit'),
 'Styrene': ('STYR', 'explicit'), 'Phenol': ('PHEN', 'explicit'), 'Benzaldehyde': ('BALD', 'explicit'), 'PAH': ('NAP', 'naphthalene'),
 'APIN': ('MTPA', 'alpha-pinene and other cyclic monoterpenes'),
 'TERP': ([('MTPA', 0.5), ('LIMO', 0.5)], 'monoterpenes split 50/50 between MTPA and LIMO'),
 'ALD_hi': ('RCHO', 'C3-C7 aldehydes'), 'RACM_ALD': ('RCHO', 'higher aldehydes'), 'RACM_HC5': ('ALK4', 'RACM2 HC5'),
 'RACM_HC8': ('ALK4', 'RACM2 HC8 (GEOS-Chem has no larger alkane lump)'), 'RACM_KET': ('MEK', 'RACM2 KET'),
 'RACM_OLI': ('PRPE', 'internal olefins'), 'RACM_OLT': ('PRPE', 'terminal olefins'), 'RACM_ROH': ('ROH', 'C3+ alcohols'),
 'PPN': ('PPN', 'explicit'), 'CH3NO3': ('MENO3', 'methyl nitrate (ASP placeholder species CH3NO3 stays empty)'),
})

T1 = _m({
 'HCHO': ('CH2O', 'explicit'), 'MEOH': ('CH3OH', 'explicit'), 'FACD': ('HCOOH', 'explicit'), 'ETH': ('C2H4', 'explicit'),
 'ETHA': ('C2H6', 'explicit'), 'ETHY': ('C2H2', 'explicit'), 'ALD2': ('CH3CHO', 'explicit'), 'ETOH': ('C2H5OH', 'explicit'),
 'AACD': ('CH3COOH', 'explicit'), 'GLY': ('GLYOXAL', 'explicit'), 'GLYD': ('GLYALD', 'explicit'),
 'MeFormate': ('C2H6', 'kOH 2.3e-13 ~ ethane 2.4e-13'), 'PRPA': ('C3H8', 'explicit'), 'Propene': ('C3H6', 'explicit'),
 'Acrolein': ('MACR', 'no acrolein in T1: unsaturated carbonyl MACR (kOH 2.9e-11 vs 2.0e-11)'),
 'Propanal': ('CH3CHO', 'no propanal in T1: acetaldehyde (kOH 1.5e-11 vs 2.0e-11)'), 'Acetone': ('CH3COCH3', 'explicit'),
 'MGLY': ('CH3COCHO', 'explicit'), 'Hydroxyacetone': ('HYAC', 'explicit'), 'nPropanol': ('C2H5OH', 'no C3+ alcohol in T1: ethanol'),
 'C3H6O': ('C3H6', 'Timbavati lumped alkene treated as propene'), 'Butadiene': ('BIGENE', 'no butadiene in T1: BIGENE'),
 '1-Butene': ('BIGENE', 'BIGENE = C4+ alkenes'), 'iButene': ('BIGENE', 'BIGENE'), '2-Butene': ('BIGENE', 'BIGENE'),
 'C4H8': ('BIGENE', 'BIGENE'), 'Butane': ('BIGALK', 'BIGALK = C4+ alkanes'), 'Butanal': ('CH3CHO', 'no higher aldehyde in T1: acetaldehyde'),
 'MEK': ('MEK', 'explicit'), 'nButanol': ('C2H5OH', 'no C3+ alcohol in T1: ethanol'), 'Furan': ('BIGENE', 'no furan in T1: BIGENE (kOH 4e-11)'),
 'Biacetyl': ([('CH3COCHO', 1.33)], 'biacetyl: 1.33 MGLY (kOH-weighted)'), 'ISOP': ('ISOP', 'explicit'),
 'Pentane': ('BIGALK', 'BIGALK'), '1-Pentene': ('BIGENE', 'BIGENE'), '2-Pentene': ('BIGENE', 'BIGENE'),
 'Cyclopentanol': ('C2H5OH', 'alcohol'), 'Pentanone': ('MEK', 'ketone'), 'Hexane': ('BIGALK', 'BIGALK'), '1-Hexene': ('BIGENE', 'BIGENE'),
 'Heptane': ('BIGALK', 'BIGALK'), 'Heptanone': ('MEK', 'ketone'), 'Octene': ('BIGENE', 'BIGENE'), 'Octanone': ('MEK', 'ketone'),
 'MeAcetate': ('C2H6', 'kOH 3.4e-13 ~ ethane'), 'Benzene': ('BENZENE', 'explicit'), 'TOL': ('TOLUENE', 'explicit'),
 'Xylene': ('XYLENES', 'explicit'), 'Ethylbenzene': ('TOLUENE', 'no ethylbenzene in T1: toluene (kOH 7.0e-12 vs 5.6e-12)'),
 'Styrene': ('BIGENE', 'no styrene in T1: reactivity dominated by the C=C bond, BIGENE'), 'Phenol': ('PHENOL', 'explicit'),
 'Benzaldehyde': ('BZALD', 'explicit'), 'PAH': ('TOLUENE', 'no PAH in T1: toluene (approximation, small amounts)'),
 'APIN': ('APIN', 'explicit'), 'TERP': ([('APIN', 0.5), ('LIMON', 0.5)], 'monoterpenes split 50/50 between APIN and LIMON'),
 'ALD_hi': ('CH3CHO', 'C3-C7 aldehydes -> acetaldehyde'), 'RACM_ALD': ('CH3CHO', 'higher aldehydes -> acetaldehyde'),
 'RACM_HC5': ('BIGALK', 'RACM2 HC5'), 'RACM_HC8': ('BIGALK', 'RACM2 HC8'), 'RACM_KET': ('MEK', 'RACM2 KET'),
 'RACM_OLI': ('BIGENE', 'internal olefins'), 'RACM_OLT': ('C3H6', 'terminal olefins as propene'), 'RACM_ROH': ('C2H5OH', 'C3+ alcohols'),
 'PPN': ('PAN', 'no PPN in T1: PAN'), 'CH3NO3': ('ONITR', 'organic nitrates (ASP placeholder species CH3NO3 stays empty)'),
})

def check():
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'CB6r5Translation'))
    import species_map as sm
    direct = {'CO','CO2','NO','NO2','O3','H2O2','HONO','SO2','HNO3','HCl','NH3','H2','CH4','N2O5','PAN','HCN','N2O','H2SO4'}
    for nm, d in (('GC', GC), ('T1', T1)):
        miss = [k for k in sm.MAP if k not in d and k not in direct]
        print(nm, 'unmapped keys:', miss)
if __name__ == '__main__': check()
