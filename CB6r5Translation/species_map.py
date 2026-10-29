"""Mapping of the compounds in the papers' initial / background tables to CB6r5 model species.

Every entry:  key -> (components, n_carbon, basis)
  components : [(CB6r5 species, ppb of that species per ppb of compound)]
  n_carbon   : carbon atoms in the compound (None where not meaningful)
  basis      : 'doc'    - stated in the CB6 documentation / Ramboll CB6 mapping memo (see README)
               'rule'   - follows from the CB6 structural lumping rules (PAR=1 C, OLE=2 C terminal C=C,
                          IOLE=4 C internal C=C, KET=1 C, ALDX=2 C, TOL=7 C, XYLMN=8 C, TERP=10 C)
               'assump' - a judgement call that has no explicit statement in the documentation
                          (documented in the README so it can be changed)
               'direct' - explicit CB6r5 species

Carbon represented by a CB6r5 species (for the carbon-balance column of the mapping tables).
"""
from collections import OrderedDict

CARBON = {
    'PAR': 1, 'OLE': 2, 'IOLE': 4, 'KET': 1, 'ALDX': 2, 'ALD2': 2, 'TOL': 7, 'XYLMN': 8,
    'TERP': 10, 'APIN': 10, 'ISOP': 5, 'BENZENE': 6, 'NAPH': 10, 'CRES': 7,
    'ETH': 2, 'ETHA': 2, 'ETHY': 2, 'ETOH': 2, 'AACD': 2, 'FACD': 1, 'FORM': 1, 'MEOH': 1,
    'GLY': 2, 'GLYD': 2, 'MGLY': 3, 'ACET': 3, 'PRPA': 3, 'ISPD': 4, 'CO': 1, 'CH4': 1,
    'CO2': 1, 'HCN': 1, 'PAN': 2, 'PANX': 3, 'NTR1': 1,
}

MAP = OrderedDict()


def m(key, comp, nC=None, basis='rule', note=''):
    MAP[key] = (comp, nC, basis, note)


# ---- inorganic / direct ------------------------------------------------------------------
for k in ('CO', 'CO2', 'NO', 'NO2', 'O3', 'H2O2', 'HONO', 'SO2', 'HNO3', 'HCl', 'NH3', 'H2',
          'CH4', 'N2O5', 'PAN', 'HCN', 'N2O'):
    m(k, [(k, 1.0)], None, 'direct')
m('H2SO4', [('H2SO4', 1.0)], None, 'direct', 'CMAQ name SULF; ASP condensation/aqueous decks use H2SO4')
m('CO2', [('CO2', 1.0)], None, 'direct', 'inert tracer (not in CB6r5)')
m('HCN', [('HCN', 1.0)], None, 'direct', 'inert tracer (not in CB6r5)')
m('N2O', [('N2O', 1.0)], None, 'direct', 'inert tracer (not in CB6r5)')

# ---- C1 - C2 explicit --------------------------------------------------------------------
m('HCHO', [('FORM', 1.0)], 1, 'direct')
m('MEOH', [('MEOH', 1.0)], 1, 'direct')
m('FACD', [('FACD', 1.0)], 1, 'direct')
m('ETH', [('ETH', 1.0)], 2, 'direct')
m('ETHA', [('ETHA', 1.0)], 2, 'direct')
m('ETHY', [('ETHY', 1.0)], 2, 'direct')
m('ALD2', [('ALD2', 1.0)], 2, 'direct')
m('ETOH', [('ETOH', 1.0)], 2, 'direct')
m('AACD', [('AACD', 1.0)], 2, 'direct')
m('GLY', [('GLY', 1.0)], 2, 'direct')
m('GLYD', [('GLYD', 1.0)], 2, 'direct', 'hydroxyacetaldehyde / glycolaldehyde (OHCH2CHO, HOCH2CHO)')
m('MeFormate', [('PAR', 1.0)], 2, 'rule',
  'methyl formate: ester carbon is unreactive (UNR, not a CB6r5 model species) -> 1 PAR; cf. methyl acetate = 2 PAR + UNR')

# ---- C3 ----------------------------------------------------------------------------------
m('PRPA', [('PRPA', 1.0)], 3, 'direct')
m('Propene', [('OLE', 1.0), ('PAR', 1.0)], 3, 'doc', 'propene = OLE + PAR')
m('Acrolein', [('ISPD', 0.75)], 3, 'assump',
  'CMAQ BUTADIENE13/ACROLEIN are no-feedback toxics tracers, so acrolein is represented by ISPD (the CB6 unsaturated '
  'carbonyl lump, whose photolysis already uses acrolein cross sections), scaled to conserve carbon (3/4)')
m('Propanal', [('ALDX', 1.0), ('PAR', 1.0)], 3, 'doc', 'propionaldehyde = ALDX + PAR')
m('Acetone', [('ACET', 1.0)], 3, 'direct')
m('MGLY', [('MGLY', 1.0)], 3, 'direct')
m('Hydroxyacetone', [('KET', 1.0), ('PAR', 2.0)], 3, 'assump',
  'CH3C(O)CH2OH (ACETOL, C3H6O2): ketone group + 2 PAR; no explicit mapping found')
m('nPropanol', [('PAR', 3.0)], 3, 'assump', 'C3+ alcohols other than methanol/ethanol: one PAR per carbon')
m('C3H6O', [('OLE', 1.0), ('PAR', 1.0)], 3, 'assump',
  'Timbavati Table S3 "C3H6O" is classified OLEL (lumped alkene) by the paper; mapped as propene')

# ---- C4 ----------------------------------------------------------------------------------
m('Butadiene', [('OLE', 2.0)], 4, 'assump',
  'two terminal C=C bonds -> 2 OLE (CB treats each double bond separately; 4 C conserved)')
m('1-Butene', [('OLE', 1.0), ('PAR', 2.0)], 4, 'rule')
m('iButene', [('OLE', 1.0), ('PAR', 2.0)], 4, 'rule', 'isobutene: terminal olefin + 2 PAR (branching ignored)')
m('2-Butene', [('IOLE', 1.0)], 4, 'doc', 'cis/trans-2-butene = IOLE')
m('C4H8', [('OLE', 1.0), ('PAR', 2.0)], 4, 'assump', 'Timbavati "C4H8" (unspeciated butenes): treated as 1-butene (as the baseline deck did)')
m('Butane', [('PAR', 4.0)], 4, 'doc', 'butane = 4 PAR (n- and iso-)')
m('Butanal', [('ALDX', 1.0), ('PAR', 2.0)], 4, 'rule')
m('MEK', [('KET', 1.0), ('PAR', 3.0)], 4, 'doc', 'butanone = KET + 3 PAR')
m('nButanol', [('PAR', 4.0)], 4, 'assump', 'C3+ alcohols: one PAR per carbon')
m('Furan', [('OLE', 2.0)], 4, 'rule',
  'furan ring = 2 C=C -> 2 OLE (cf. 2-butylfuran = 2 OLE + 4 PAR, pyrrole = 2 OLE in the Ramboll memo)')
m('Biacetyl', [('MGLY', 1.5)], 4, 'doc', '2,3-butanedione = 1.5 MGLY (CB6 documentation)')

# ---- C5 - C8 -----------------------------------------------------------------------------
m('ISOP', [('ISOP', 1.0)], 5, 'direct')
m('Pentane', [('PAR', 5.0)], 5, 'rule', 'n-, i-pentane, 2-methylbutane')
m('1-Pentene', [('OLE', 1.0), ('PAR', 3.0)], 5, 'rule', '1-pentene, 3-methyl-1-butene, 2-methyl-1-butene')
m('2-Pentene', [('IOLE', 1.0), ('PAR', 1.0)], 5, 'doc', '2-pentene isomers = IOLE + PAR; also 2-methyl-2-butene, cyclopentene')
m('Cyclopentanol', [('PAR', 5.0)], 5, 'assump', 'alcohol: one PAR per carbon')
m('Pentanone', [('KET', 1.0), ('PAR', 4.0)], 5, 'rule')
m('Hexane', [('PAR', 6.0)], 6, 'rule', 'n-hexane, isohexanes')
m('1-Hexene', [('OLE', 1.0), ('PAR', 4.0)], 6, 'rule', '1-hexene, 2-methyl-1-pentene, 4-methyl-1-pentene')
m('Heptane', [('PAR', 7.0)], 7, 'rule')
m('Heptanone', [('KET', 1.0), ('PAR', 6.0)], 7, 'rule')
m('Octene', [('OLE', 1.0), ('PAR', 6.0)], 8, 'assump', 'octenes assumed terminal (1-octene)')
m('Octanone', [('KET', 1.0), ('PAR', 7.0)], 8, 'rule')
m('MeAcetate', [('PAR', 2.0)], 3, 'doc', 'methyl acetate = 2 PAR + UNR')

# ---- aromatics / terpenes ----------------------------------------------------------------
m('Benzene', [('BENZENE', 1.0)], 6, 'direct')
m('TOL', [('TOL', 1.0)], 7, 'direct')
m('Xylene', [('XYLMN', 1.0)], 8, 'direct', 'xylene isomers')
m('Ethylbenzene', [('TOL', 1.0), ('PAR', 1.0)], 8, 'doc', 'ethylbenzene = TOL + PAR')
m('Styrene', [('TOL', 1.0), ('OLE', 1.0), ('PAR', -1.0)], 8, 'assump',
  'Ramboll memo: styrenes = TOL, OLE and PAR (coefficients not given); TOL + OLE - PAR conserves 8 C')
m('Phenol', [('CRES', 1.0)], 6, 'assump', 'phenol -> CRES (cresol and higher phenols); per-mole, carbon not conserved')
m('Benzaldehyde', [('TOL', 1.0)], 7, 'assump', 'aromatic ring + one carbon: TOL (7 C); no explicit mapping found')
m('PAH', [('NAPH', 1.0)], None, 'assump', 'PAHs represented as naphthalene (the only PAH in CB6r5)')
m('APIN', [('APIN', 1.0)], 10, 'direct', 'alpha-pinene (CMAQ explicit species)')
m('TERP', [('TERP', 1.0)], 10, 'direct', 'limonene and other monoterpenes; CB6 TERP = all monoterpenes')

# ---- aldehyde / alkane / olefin lumps used by the Williams (RACM2-style) deck ------------
m('ALD_hi', [('ALDX', 1.0), ('PAR', 3.0)], 5, 'assump', 'higher aldehydes (C3-C7): average C5 = ALDX + 3 PAR')
m('RACM_ALD', [('ALDX', 1.0), ('PAR', 3.0)], 5, 'assump',
  'Williams "ALD" (acetaldehyde is listed separately, so C3+): same as C3-C7 aldehydes')
m('RACM_HC5', [('PAR', 5.0)], 5, 'assump', 'RACM2 HC5 (pentanes/hexanes) as 5 PAR')
m('RACM_HC8', [('PAR', 8.0)], 8, 'assump', 'RACM2 HC8 as 8 PAR')
m('RACM_KET', [('KET', 1.0), ('PAR', 3.0)], 4, 'assump', 'RACM2 KET (MEK and larger) as MEK')
m('RACM_OLI', [('IOLE', 1.0)], 4, 'assump', 'RACM2 OLI (internal olefins) as 2-butene')
m('RACM_OLT', [('OLE', 1.0), ('PAR', 1.0)], 3, 'assump', 'RACM2 OLT (terminal olefins) as propene')
m('RACM_ROH', [('PAR', 3.0)], 3, 'assump', 'RACM2 ROH (C3+ alcohols) as propanol')

# ---- nitrates ----------------------------------------------------------------------------
m('PPN', [('PANX', 1.0)], 3, 'direct', 'peroxypropionyl nitrate -> PANX (higher PANs)')
m('CH3NO3', [('NTR1', 1.0)], 1, 'direct', 'methyl nitrate -> NTR1 (monofunctional organic nitrates)')

# =========================================================================================
# per-case source tables  (source name -> key)
# =========================================================================================
# JGR Alvarado & Prinn (2009) auxiliary tables S1 (Otavi), S2 (Alaska), S3 (Timbavati)
JGR_NAMES = {
    'CO': 'CO', 'CO2': 'CO2', 'CH4': 'CH4', 'O3': 'O3', 'NO2': 'NO2', 'NO': 'NO', 'HNO3': 'HNO3',
    'H2O2': 'H2O2', 'H2': 'H2', 'HONO': 'HONO', 'SO2': 'SO2', 'NH3': 'NH3', 'N2O': 'N2O',
    'N2O5': 'N2O5', 'CH3C(O)O2NO2': 'PAN',
    'C2H6': 'ETHA', 'C2H4': 'ETH', 'C2H2': 'ETHY', 'C3H8': 'PRPA', 'C3H6': 'Propene',
    'CH3OH': 'MEOH', 'CH2O': 'HCHO', 'HCHO': 'HCHO', 'OHCH2CHO': 'GLYD', 'CH3CHO': 'ALD2',
    'CH3C(O)CH3': 'Acetone', 'HC(O)OH': 'FACD', 'HCOOH': 'FACD', 'CH3C(O)OH': 'AACD',
    'CH3COOH': 'AACD', 'CH2=C(CH3)CH=CH2': 'ISOP', 'C5H8': 'ISOP', 'C6H5CH3': 'TOL',
    'CH3C(O)C2H5': 'MEK', 'CH3COC2H5': 'MEK', 'C6H5OH': 'Phenol', 'CH3C(O)CH2OH': 'Hydroxyacetone',
    'C3H6O2': 'Hydroxyacetone', 'C4H8': 'C4H8', 'C3H6O': 'C3H6O',
    # Timbavati (Table S3) - verbose names
    'i-butane': 'Butane', 'n-butane': 'Butane', 'trans-2-butene': '2-Butene', 'cis-2-butene': '2-Butene',
    'i-pentane': 'Pentane', 'n-pentane': 'Pentane', '2-methyl-butane': 'Pentane',
    'butadiene': 'Butadiene', '3-methyl-1-butene': '1-Pentene', '2-methyl-1-butene': '1-Pentene',
    '1-pentene': '1-Pentene', 'trans-2-pentene': '2-Pentene', 'c-2-pentene': '2-Pentene',
    '2-methyl-2-butene': '2-Pentene', 'cyclopentene': '2-Pentene',
    'n-hexane': 'Hexane', 'isohexanes': 'Hexane', '2-methyl-1-pentene': '1-Hexene',
    '4-methyl-1-pentene': '1-Hexene', '1-hexene': '1-Hexene', 'heptane': 'Heptane',
    '1-butene': '1-Butene', 'i-butene': 'iButene', 'octenes': 'Octene', 'terpenes': 'TERP',
    'xylenes': 'Xylene', 'ethylbenzene': 'Ethylbenzene', 'styrene': 'Styrene', 'PAHs': 'PAH',
    'ethanol': 'ETOH', '1-propanol': 'nPropanol', 'butanols': 'nButanol', 'cyclopentanol': 'Cyclopentanol',
    'acrolein': 'Acrolein', 'C3-C7_aldehydes': 'ALD_hi', '2,3-butanedione': 'Biacetyl',
    'pentanones': 'Pentanone', 'heptanones': 'Heptanone', 'octanones': 'Octanone',
    'benzaldehyde': 'Benzaldehyde', 'methyl_formate': 'MeFormate', 'methyl_acetate': 'MeAcetate',
}

# Williams: names in the ACP-2015 supplement deck (GasPhaseChems.in, identical to Table S3)
WILLIAMS_NAMES = {
    'CO': 'CO', 'CO2': 'CO2', 'NO': 'NO', 'NO2': 'NO2', 'O3': 'O3', 'H2O2': 'H2O2', 'HONO': 'HONO',
    'SO2': 'SO2', 'HNO3': 'HNO3', 'HCl': 'HCl', 'NH3': 'NH3', 'H2': 'H2', 'CH4': 'CH4',
    'HCHO': 'HCHO', 'MEOH': 'MEOH', 'Formic_Acid': 'FACD', 'HCN': 'HCN',
    'ETHE': 'ETH', 'Ethane': 'ETHA', 'Acetaldehyde': 'ALD2', 'ETOH': 'ETOH', 'Acetic_Acid': 'AACD',
    'GLY': 'GLY', 'Acetylene': 'ETHY', 'Propene': 'Propene', 'Propane': 'PRPA', 'ACR': 'Acrolein',
    'Propanal': 'Propanal', 'Acetone': 'Acetone', 'NPROPOL': 'nPropanol', 'MGLY': 'MGLY',
    'Butadiene': 'Butadiene', '1-Butene': '1-Butene', 'cis-2-Butene': '2-Butene',
    'trans-2-Butene': '2-Butene', 'i-Butene': 'iButene', 'n-Butane': 'Butane', 'i-Butane': 'Butane',
    'Butanal': 'Butanal', 'MEK': 'MEK', 'NBUTOL': 'nButanol', 'ISOP': 'ISOP', 'FURAN': 'Furan',
    'ALD': 'RACM_ALD', 'API': 'APIN', 'BALD': 'Benzaldehyde', 'BEN': 'Benzene', 'DIEN': 'Butadiene',
    'HC5': 'RACM_HC5', 'KET': 'RACM_KET', 'LIM': 'TERP', 'OLI': 'RACM_OLI', 'OLT': 'RACM_OLT',
    'PHEN': 'Phenol', 'ROH': 'RACM_ROH', 'TOL': 'TOL', 'XYM': 'Xylene', 'XYP': 'Xylene', 'XYO': 'Xylene',
    'HOCH2CHO': 'GLYD', 'ACETOL': 'Hydroxyacetone', 'BIACET': 'Biacetyl', 'PAN': 'PAN', 'PPN': 'PPN',
    'CH3NO3': 'CH3NO3',
}
# names that are carried over unchanged from the existing deck (ASP-specific, not gas-phase chemistry)
WILLIAMS_PASSTHRU = {'Trace', 'O2', 'N2', 'OffGas', 'Wall_O3', 'Wall_NOx'}


def apply(entries):
    """entries: [(source_name, key, init, bg)] -> OrderedDict species -> [init, bg], plus detail rows"""
    out = OrderedDict()
    rows = []
    for src, key, ini, bg in entries:
        comp, nC, basis, note = MAP[key]
        for sp, c in comp:
            v = out.setdefault(sp, [0.0, 0.0])
            v[0] += c * ini
            v[1] += c * bg
        rows.append((src, key, ini, bg, comp, nC, basis, note))
    return out, rows


def carbon_represented(comp):
    tot = 0.0
    for sp, c in comp:
        if sp not in CARBON:
            return None
        tot += c * CARBON[sp]
    return tot
