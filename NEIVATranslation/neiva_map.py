"""NEIVA compound -> baseline (RACM2-lumped, Williams) deck species.

classify(row) returns (route, target, nC_target):
    ('explicit', species, nC)   baseline species that is the compound itself (nC = carbon number of the compound)
    ('lump',     species, nC)   baseline lump; moles are carbon-conserving:  mol_lump = mol * nC_compound / nC_lump
    ('ivoc',     k, None)       unmapped identified NMOG -> IVOCk (k from log10 C*, C*(IVOCk) = 10^(k-3) ug/m3), mass-conserving
    ('skip',     reason, None)  halogenated / inorganic compounds that the mechanism does not carry
"""
import re

EXPLICIT = {
    ('CH2O', 'formaldehyde'): 'HCHO', ('CH4O', 'methanol'): 'MEOH', ('CH2O2', 'formic acid'): 'Formic_Acid',
    ('C2H4', 'ethene'): 'ETHE', ('C2H6', 'ethane'): 'Ethane', ('C2H4O', 'acetaldehyde'): 'Acetaldehyde', ('C2H6O', 'ethanol'): 'ETOH',
    ('C2H4O2', 'acetic acid'): 'Acetic_Acid', ('C2H4O2', 'glycolaldehyde'): 'HOCH2CHO', ('C2H2O2', 'glyoxal'): 'GLY',
    ('C2H2', 'acetylene'): 'Acetylene', ('C3H6', 'propene'): 'Propene', ('C3H8', 'propane'): 'Propane', ('C3H4O', 'acrolein'): 'ACR',
    ('C3H6O', 'propanal'): 'Propanal', ('C3H6O', 'acetone'): 'Acetone', ('C3H4O2', 'methylglyoxal'): 'MGLY', ('C3H6O2', 'acetol'): 'ACETOL',
    ('C3H8O', '1-propanol'): 'NPROPOL', ('C4H6', '1,3-butadiene'): 'Butadiene', ('C4H8', '1-butene'): '1-Butene',
    ('C4H8', 'cis-2-butene'): 'cis-2-Butene', ('C4H8', 'trans-2-butene'): 'trans-2-Butene', ('C4H8', 'i-butene'): 'i-Butene',
    ('C4H10', 'butane'): 'n-Butane', ('C4H10', 'i-butane'): 'i-Butane', ('C4H8O', 'n-butanal'): 'Butanal',
    ('C4H8O', 'methyl ethyl ketone'): 'MEK', ('C4H10O', '1-butanol'): 'NBUTOL', ('C5H8', 'isoprene'): 'ISOP', ('C4H4O', 'furan'): 'FURAN',
    ('C6H6', 'benzene'): 'BEN', ('C6H6O', 'phenol'): 'PHEN', ('C7H8', 'toluene'): 'TOL', ('C7H6O', 'benzaldehyde'): 'BALD',
    ('C4H6O2', '2,3-butadione'): 'BIACET', ('C8H10', 'o-xylene'): 'XYO', ('CH3NO3', 'methyl nitrate'): 'CH3NO3',
    ('C10H16', 'alpha-pinene'): 'API', ('C10H16', 'd-limonene'): 'LIM',
}
# baseline lumps: species -> (MW, carbon number represented)
LUMPS = {'HC5': 5, 'HC8': 8, 'OLT': 3, 'OLI': 5, 'DIEN': 4, 'ALD': 3, 'KET': 5, 'ROH': 3, 'XYM': 8, 'PHEN': 6, 'CSL': 7, 'TOL': 7,
         'API': 10, 'LIM': 10, 'FURAN': 4}
PINENE_LIKE = ('pinene', 'carene', 'camphene', 'tricyclene', 'thujene', 'sabinene', 'fenchene', 'pyronene')
IVOC_MW = {1: 524., 2: 479., 3: 434., 4: 389., 5: 344., 6: 299., 7: 254., 8: 208., 9: 163.}


def ncarbon(formula):
    m = re.match(r'^C(\d*)', str(formula))
    if not m or not str(formula).startswith('C') or str(formula).startswith(('Cl', 'Co', 'Ca', 'Cs', 'Cu', 'Cd', 'Cr')):
        return 0
    return int(m.group(1) or 1)


def _terminal_alkene(smile):
    s = str(smile)
    return s.startswith('C=') or s.endswith('=C') or '(C=' in s or '(=C)' in s or 'C(=C)' in s


def classify(formula, compound, voc_type, smile, kOH, cstar):
    f = str(formula); c = str(compound).strip().lower(); v = str(voc_type); nC = ncarbon(f)
    if (f, c) in EXPLICIT:
        sp = EXPLICIT[(f, c)]
        return ('explicit', sp, nC)
    if nC == 0:
        return ('skip', 'no carbon', None)
    els = set(re.findall(r'[A-Z][a-z]?', f))
    if 'halogen' in v or els & {'Cl', 'Br', 'F', 'I'}:
        return ('skip', 'halogenated', None)
    if f == 'C10H16':
        if any(k in c for k in PINENE_LIKE):
            return ('lump', 'API', nC)
        if 'unknown' in c or 'isomers' in c:
            return ('lump', 'LIM', nC)
        return ('lump', 'LIM', nC)
    if 'N-containing' in v or els & {'N', 'S', 'P', 'Hg', 'Si'}:
        return ('ivoc', None, None)
    cyc = 'Cyclic' in v
    if v.startswith('Aromatic'):
        if 'polycyclic' in v:
            return ('ivoc', None, None)
        if 'phenol' in v:
            return ('lump', 'CSL' if nC >= 7 else 'PHEN', nC)
        if 'alkene' in v or nC >= 8:
            return ('lump', 'XYM', nC)
        return ('lump', 'TOL', nC)
    if 'furan' in v and nC <= 8:
        return ('lump', 'FURAN', nC)
    simple = v.replace('Aliphatic HC:', '')
    if simple == 'alkane' or simple == 'Cyclic;alkane':
        if nC <= 2:
            return ('ivoc', None, None)
        k = kOH if (kOH is not None) else 5e-12
        return ('lump', 'HC5' if k < 6.8e-12 else 'HC8', nC)
    if simple == 'alkene' or simple == 'Cyclic;alkene':
        nd = smile.count('=') if smile else 1
        if nd >= 2 or ('alkyne' in v):
            return ('lump', 'DIEN', nC)
        return ('lump', 'OLT' if (_terminal_alkene(smile) and not cyc) else 'OLI', nC)
    if simple in ('alkene;alkyne', 'alkyne'):
        return ('lump', 'DIEN', nC)
    if simple == 'aldehyde' or simple == 'alkene;aldehyde':
        return ('lump', 'ALD', nC)
    if simple == 'ketone' or simple == 'alkene;ketone':
        return ('lump', 'KET', nC)
    if simple in ('alcohol',):
        return ('lump', 'ROH', nC)
    return ('ivoc', None, None)


def ivoc_bin(cstar):
    """log10 C* -> IVOC index (C*(IVOCk) = 10^(k-3) ug/m3, k = 1..9)"""
    if cstar is None or cstar != cstar:
        return 9
    return int(min(9, max(1, round(cstar) + 3)))
