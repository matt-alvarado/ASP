"""Species data and paper-species mapping for the CRACMM2/3 and SAPRC07TC/TIC decks."""
import json
import os
import re
from collections import OrderedDict

import mech2asp as m2a

HERE = os.path.dirname(os.path.abspath(__file__))

# ASP-specific / inorganic species whose molar mass is not in the mechanism tables
MW_EXTRA = {'H2SO4': 98.08, 'HCl': 36.46, 'CH4': 16.04, 'H2': 2.016, 'CO2': 44.01, 'NH3': 17.03, 'HCN': 27.03,
            'N2O': 44.01, 'CH3NO3': 77.04, 'RO2_T': 12.0, 'RCO3_T': 12.0, 'CO': 28.01, 'O1D': 16.0, 'O': 16.0,
            'O3': 48.0, 'NO': 30.01, 'NO2': 46.01, 'NO3': 62.0, 'N2O5': 108.0, 'HONO': 47.01, 'HNO3': 63.01,
            'SO2': 64.06, 'H2O2': 34.01, 'OH': 17.01, 'HO2': 33.01, 'HNO4': 79.02, 'PAN': 121.05, 'Trace': 12.0}


def read_mw(mech):
    """-> OrderedDict name -> (MW, description) using the mechanism's own species tables (CMAQ names)"""
    cfg = m2a.CFG[mech]
    out = OrderedDict()
    if cfg['family'] == 'cracmm':
        for ln in open(cfg['gc'], errors='ignore'):
            mm = re.match(r"^\s*'([^']+)'\s*,\s*([\d.]+)\s*,.*?!([^,]*),", ln)
            if mm:
                out[mm.group(1)] = (float(mm.group(2)), mm.group(3).strip())
    else:
        for ln in open(cfg['sp_table'], errors='ignore'):
            f = [x.strip() for x in ln.strip().strip('|').split('|')]
            if len(f) >= 4:
                try:
                    out[f[0]] = (float(f[3]), f[1])
                except ValueError:
                    pass
    return out


def inv_rename(mech):
    return {v: k for k, v in m2a.CFG[mech]['rename'].items()}


def mw_of(mech, name, table):
    if name in MW_EXTRA:
        return MW_EXTRA[name]
    src = inv_rename(mech).get(name, name)
    if src in table:
        return table[src][0]
    raise KeyError('no molar mass for %s (%s)' % (name, mech))


ACYL = re.compile(r'(?i)acyl|acetylperoxy|propionylperoxy|enoylperoxy|enoyl\)?peroxy|benzoyl|peroxyacyl|methylprop-2-enoyl|acylperoxy')


def peroxy_flag(mech, name, table):
    """1 = counted in RO2_T, 2 = also in RCO3_T (acyl peroxy), 0 = not a peroxy radical"""
    src = inv_rename(mech).get(name, name)
    if name in ('RO2_T', 'RCO3_T'):
        return 0
    desc = table.get(src, (0, ''))[1]
    fam = m2a.CFG[mech]['family']
    if fam == 'saprc':
        if name in ('RO2C', 'RO2XC', 'MEO2', 'TERPNRO2'):
            return 1
        if name in ('MECO3', 'RCO3', 'BZCO3', 'MACO3', 'IMACO3', 'HCOCO3', 'NIT1NO3OOA'):
            return 2
        if re.search(r'(?i)peroxy radical', desc) and not re.search(r'(?i)tracker|tracer|operator', desc):
            return 2 if ACYL.search(desc) else 1
        return 0
    if name in ('PAN', 'PPN', 'MPAN'):
        return 0
    if name in ('ACO3', 'RCO3', 'MACP', 'BALP'):
        return 2
    if re.search(r'(?i)peroxy', desc) and not re.search(r'(?i)peroxyacetic|peroxynitr|hydroperox|nitrate$', desc) \
            and not name.startswith('V'):
        return 2 if ACYL.search(desc) else 1
    return 0


# ---------------------------------------------------------------------------------------------
# paper-species mapping.  key = species_map.MAP key.  Mole-for-mole onto the lump that CMAQ's own
# definition assigns by kOH class / structure (documented assumption; carbon is not conserved).
# ---------------------------------------------------------------------------------------------
def _d(*pairs):
    return list(pairs)


DIRECT = ['CO', 'CO2', 'NO', 'NO2', 'O3', 'H2O2', 'HONO', 'SO2', 'HNO3', 'HCl', 'NH3', 'H2', 'CH4', 'N2O5', 'PAN',
          'HCN', 'N2O', 'H2SO4']

SAPRC_MAP = {
    'HCHO': ('HCHO', 'explicit'), 'MEOH': ('MEOH', 'explicit'), 'FACD': ('HCOOH', 'explicit'),
    'ETH': ('ETHENE', 'explicit'), 'ETHA': ('ALK1', 'kOH 2.4e-13 cm3/molec/s = 3.6e2 ppm-1 min-1 (ALK1: 2-5e2)'),
    'ETHY': ('ACETYLENE', 'explicit'), 'ALD2': ('CCHO', 'explicit'), 'ETOH': ('ETOH', 'explicit'),
    'AACD': ('CCOOH', 'explicit'), 'GLY': ('GLY', 'explicit'), 'GLYD': ('HOCCHO', 'explicit (glycolaldehyde)'),
    'MeFormate': ('ALK1', 'kOH 2.3e-13 (3.4e2 ppm-1 min-1): ALK1'),
    'PRPA': ('ALK2', 'kOH 1.1e-12 (1.6e3): ALK2'), 'Propene': ('PROPENE', 'explicit'),
    'Acrolein': ('ACROLEIN', 'explicit'), 'Propanal': ('RCHO', 'SAPRC RCHO = propionaldehyde and higher aldehydes'),
    'Acetone': ('ACETONE', 'explicit'), 'MGLY': ('MGLY', 'explicit'),
    'nPropanol': ('ALK4', 'kOH 5.3e-12 (7.8e3 ppm-1 min-1): ALK4'),
    'C3H6O': ('PROPENE', 'Timbavati "C3H6O" is a lumped alkene in the paper: treated as propene'),
    'Butadiene': ('BUTADIENE13', 'explicit'), '1-Butene': ('OLE1', 'kOH 3.1e-11 (4.6e4 ppm-1 min-1 < 7e4): OLE1'),
    'iButene': ('OLE2', 'kOH 5.1e-11 (7.5e4 > 7e4): OLE2'),
    '2-Butene': ('OLE2', 'kOH 5.6e-11 (8.3e4): OLE2'), 'C4H8': ('OLE1', 'unspeciated butenes treated as 1-butene'),
    'Butane': ('ALK3', 'kOH 2.1-2.4e-12 (3.1-3.5e3): ALK3'), 'Butanal': ('RCHO', 'higher aldehyde'),
    'MEK': ('MEK', 'explicit'), 'nButanol': ('ALK5', 'kOH 8.5e-12 (1.3e4 > 1e4): ALK5'),
    'Furan': ('OLE1', 'kOH 4.0e-11 (5.9e4): OLE1 (no explicit furan in SAPRC07)'), 'Biacetyl': ('BACL', 'explicit'),
    'ISOP': ('ISOPRENE', 'explicit'), 'Pentane': ('ALK4', 'kOH 3.7-3.8e-12 (5.6e3): ALK4'),
    '1-Pentene': ('OLE1', 'kOH 3.1e-11: OLE1'), '2-Pentene': ('OLE2', 'kOH 6.5e-11 (9.6e4): OLE2'),
    'Cyclopentanol': ('ALK5', 'kOH 1.1e-11 (1.6e4): ALK5'),
    'Pentanone': ('PRD2', 'kOH 4.4e-12-1e-11: PRD2 (ketones with kOH > 5e3 ppm-1 min-1)'),
    'Hexane': ('ALK4', 'kOH 5.2e-12 (7.7e3): ALK4'), '1-Hexene': ('OLE1', 'kOH 3.7e-11 (5.5e4): OLE1'),
    'Heptane': ('ALK4', 'kOH 6.8e-12 (1.0e4): on the ALK4/ALK5 boundary, ALK4 chosen'),
    'Heptanone': ('PRD2', 'kOH ~1e-11: PRD2'), 'Octene': ('OLE1', 'kOH 3.3e-11: OLE1'),
    'Octanone': ('PRD2', 'kOH ~1.2e-11: PRD2'), 'MeAcetate': ('ALK1', 'kOH 3.4e-13 (5.0e2): ALK1'),
    'Benzene': ('BENZENE', 'explicit'), 'TOL': ('TOLUENE', 'explicit'),
    'Ethylbenzene': ('ARO1', 'kOH 7.0e-12 (1.0e4 < 2e4): ARO1'),
    'Styrene': ('ARO2MN', 'kOH 5.8e-11 (8.6e4 > 2e4): ARO2MN'), 'Phenol': ('CRES', 'CRES = phenols and cresols'),
    'Benzaldehyde': ('BALD', 'explicit'), 'PAH': ('NAPHTHAL', 'PAHs represented as naphthalene'),
    'APIN': ('APIN', 'explicit'), 'TERP': ('TERP', 'SAPRC TERP = other monoterpenes'),
    'ALD_hi': ('RCHO', 'C3-C7 aldehydes'), 'RACM_ALD': ('RCHO', 'higher aldehydes'),
    'RACM_HC5': ('ALK4', 'RACM2 HC5 -> ALK4 (kOH 3.4-6.8e-12)'), 'RACM_HC8': ('ALK5', 'RACM2 HC8 -> ALK5'),
    'RACM_KET': ('PRD2', 'RACM2 KET (kOH 1.2e-11) -> PRD2'), 'RACM_OLI': ('OLE2', 'internal olefins -> OLE2'),
    'RACM_OLT': ('OLE1', 'terminal olefins -> OLE1'), 'RACM_ROH': ('ALK4', 'C3+ alcohols -> ALK4'),
    'PPN': ('PAN2', 'PAN2 = PPN and higher PANs'), 'CH3NO3': ('RNO3', 'alkyl nitrates'),
}
SAPRC_TC_OVERRIDE = {'Hydroxyacetone': ('MEK', 'no HACET in SAPRC07TC: hydroxyacetone (kOH 3e-12) -> MEK')}
SAPRC_TIC_OVERRIDE = {'Hydroxyacetone': ('HACET', 'explicit')}

CRACMM_MAP = {
    'HCHO': ('HCHO', 'explicit'), 'MEOH': ('MOH', 'explicit'), 'FACD': ('ORA1', 'explicit'),
    'ETH': ('ETE', 'explicit'), 'ETHA': ('ETH', 'explicit (ethane)'), 'ETHY': ('ACE', 'explicit'),
    'ALD2': ('ACD', 'explicit'), 'ETOH': ('EOH', 'explicit'), 'AACD': ('ORA2', 'explicit'), 'GLY': ('GLY', 'explicit'),
    'GLYD': ('HKET', 'RACM2/CRACMM HKET = hydroxy ketones incl. glycolaldehyde'),
    'MeFormate': ('HC3', 'low-kOH species -> HC3'), 'PRPA': ('HC3', 'HC3: kOH < 3.4e-12 alkanes (propane...)'),
    'Propene': ('OLT', 'terminal olefins'), 'Acrolein': ('ACRO', 'explicit'),
    'Propanal': ('ALD', 'C3+ aldehydes'), 'Acetone': ('ACT', 'explicit'), 'MGLY': ('MGLY', 'explicit'),
    'Hydroxyacetone': ('HKET', 'explicit lump'), 'nPropanol': ('ROH', 'C3+ alcohols'),
    'C3H6O': ('OLT', 'treated as propene'), 'Butadiene': ('BDE13', 'explicit'), '1-Butene': ('OLT', 'terminal olefin'),
    'iButene': ('OLT', 'terminal olefin (isobutene)'), '2-Butene': ('OLI', 'internal olefin'),
    'C4H8': ('OLT', 'treated as 1-butene'), 'Butane': ('HC3', 'kOH 2.4e-12 < 3.4e-12: HC3'),
    'Butanal': ('ALD', 'C3+ aldehydes'), 'MEK': ('MEK', 'explicit'), 'nButanol': ('ROH', 'C3+ alcohols'),
    'Furan': ('FURAN', 'explicit'), 'ISOP': ('ISO', 'explicit'), 'Pentane': ('HC5', 'kOH 3.8e-12: HC5'),
    '1-Pentene': ('OLT', 'terminal olefin'), '2-Pentene': ('OLI', 'internal olefin'),
    'Cyclopentanol': ('ROH', 'alcohol'), 'Pentanone': ('KET', 'ketones other than acetone/MEK'),
    'Hexane': ('HC5', 'kOH 5.2e-12: HC5'), '1-Hexene': ('OLT', 'terminal olefin'),
    'Heptane': ('HC5', 'kOH 6.8e-12: on the HC5/HC10 boundary, HC5 chosen'), 'Heptanone': ('KET', 'ketone'),
    'Octene': ('OLT', 'terminal olefin'), 'Octanone': ('KET', 'ketone'), 'MeAcetate': ('HC3', 'low-kOH species -> HC3'),
    'Benzene': ('BEN', 'explicit'), 'TOL': ('TOL', 'explicit'), 'Xylene': ('XYL', 'explicit lump'),
    'Ethylbenzene': ('EBZ', 'explicit'), 'Styrene': ('STY', 'explicit'), 'Phenol': ('PHEN', 'explicit'),
    'Benzaldehyde': ('BALD', 'explicit'), 'PAH': ('NAPH', 'PAHs represented as naphthalene'),
    'APIN': ('API', 'alpha-pinene and other cyclic terpenes'), 'ALD_hi': ('ALD', 'C3-C7 aldehydes'),
    'RACM_ALD': ('ALD', 'identical lump'), 'RACM_HC5': ('HC5', 'identical lump'), 'RACM_HC8': ('HC10', 'HC8 -> HC10 (CRACMM has no HC8)'),
    'RACM_KET': ('KET', 'identical lump'), 'RACM_OLI': ('OLI', 'identical lump'), 'RACM_OLT': ('OLT', 'identical lump'),
    'RACM_ROH': ('ROH', 'identical lump'), 'PPN': ('PPN', 'explicit'), 'CH3NO3': ('ONIT', 'organic nitrates'),
}


def get_map(mech):
    """-> dict key -> ([(species, coef)], basis_note)"""
    fam = m2a.CFG[mech]['family']
    out = {}
    for k in DIRECT:
        out[k] = ([(k, 1.0)], 'direct')
    if fam == 'saprc':
        base = dict(SAPRC_MAP)
        base.update(SAPRC_TIC_OVERRIDE if mech == 'SAPRC07TIC' else SAPRC_TC_OVERRIDE)
        for k, (sp, note) in base.items():
            out[k] = ([(sp, 1.0)], note)
        out['Xylene'] = ([('MXYL', 1.0 / 3), ('OXYL', 1.0 / 3), ('PXYL', 1.0 / 3)], 'xylene isomers split equally')
    else:
        for k, (sp, note) in CRACMM_MAP.items():
            out[k] = ([(sp, 1.0)], note)
        out['TERP'] = ([('API', 0.5), ('LIM', 0.5)], 'monoterpenes split 50/50 between API (cyclic, one C=C) and LIM (limonene-like)')
        out['Biacetyl'] = ([('MGLY', 1.33)], 'biacetyl: 1.33 MGLY (kOH-weighted)')
    return out


def apply(mech, entries):
    """entries [(src, key, init, bg)] -> OrderedDict species -> [init, bg], detail rows"""
    mp = get_map(mech)
    if m2a.CFG[mech]['family'] == 'saprc':
        pass
    out = OrderedDict()
    rows = []
    for src, key, ini, bg in entries:
        comp, note = mp[key]
        for sp, c in comp:
            v = out.setdefault(sp, [0.0, 0.0])
            v[0] += c * ini
            v[1] += c * bg
        rows.append((src, key, ini, bg, comp, note))
    return out, rows
