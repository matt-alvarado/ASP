"""CMAQ gas->particle uptake (tier 1) and first-order particle-phase organic reactions (tier 2) for the translated decks.

apply(case_dir, mech, acid=True) edits a finished deck in place (call it after the deck builder; it is idempotent, the added
lines sit between the markers `CMAQ_HET_BEGIN` / `CMAQ_HET_END` in GasChemicalMechanism.in):

  tier 1  GasChemicalMechanism.in   uptake reactions   GLY / MGLY (flag 3, gamma 2.9e-3 and 0.09 x 2.9e-3, Liggio et al. 2005 as in CMAQ),
                                    IEPOX / IMAE / IHMML / IPX / INALD with the acidity-dependent flag 18 (common_asp >= 'rate flag 18').
          GasPhaseChems.in, OrgPhaseChems.in, OrganicDissolution.in, Hydrophilic*  nonvolatile product species (like ASOAT / AXYL3)
  tier 2  OrgPhaseReactions.in      oligomerisation (9.48816e-6 /s), organic-nitrate hydrolysis (9.259e-5 /s), CRACMM ATRPN / AHONIT
                                    -> AHOM + HNO3 (9.26e-5 /s) and, for CRACMM3, AHOM photolysis (0.01 x jNO2)

Also (where the CMAQ mech has them): NTR2 hydrolysis (flag 19), N2O5 uptake (flag 20, Davis et al. gamma), NO3 (gamma 1e-3) and HO2
(gamma 0.2, no products as in CMAQ) uptake.  NO2 uptake (flag 21, 1e-4 m/s x fine surface area, 0.5 HONO + 0.5 HNO3).  Not represented: chlorine / bromine / iodine and SO2 heterogeneous chemistry.
CMAQ's CB6r5 / AERO6 IEPOX goes through the intermediates IEPOXP (+ ASO4J); here the epoxide goes directly to the product
(sulfate is not consumed), and the AERO7i dimers are folded into the tetrol / organosulfate channels.
"""
import glob
import os
import re

import cmaq_ae
import cmaq_oa

MECHS = '/home/claude/usepa/cmaq/CCTM/src/MECHS/'
NML = {'CB6R5': 'cb6r5_ae7_aq', 'SAPRC07TC': 'saprc07tc_ae6_aq', 'SAPRC07TIC': 'saprc07tic_ae7i_aq', 'CRACMM2': 'cracmm2', 'CRACMM3': 'cracmm3'}
GAMMA_GLY = 2.9e-3
BEGIN, END = '!! CMAQ_HET_BEGIN', '!! CMAQ_HET_END'
ASOA_CSTAR, ASOA_ENTH = 1e-10, 1.0e3

# (gas reactant, [(product species, channel)], kind, scale, extra gas product)   kind None -> flag 3 with gamma
GLY = [('GLY', 'AGLY', 1.0), ('MGLY', 'AGLY', 0.09)]
GLY3 = [('GLY', 'AGLYOLIG', 1.0), ('MGLY', 'AGLYOLIG', 0.09)]
UPTAKE = {
    'CB6R5':      dict(const=GLY, acid=[('EPOX', [('AISO3', 0)], 1, 1.0, None)]),
    'SAPRC07TC':  dict(const=[], acid=[('IEPOX', [('AISO3', 0)], 1, 1.0, None)]),
    'SAPRC07TIC': dict(const=GLY, acid=[('IEPOX', [('AIETET', 1), ('AIEOS', 2)], 3, 1.0, None),
                                        ('IMAE', [('AIMGA', 1), ('AIMOS', 2)], 2, 1.0, None),
                                        ('IHMML', [('AIMGA', 1), ('AIMOS', 2)], 2, 1.0, None)]),
    'CRACMM2':    dict(const=GLY, acid=[('IEPOX', [('AISO3NOS', 1), ('AISO3OS', 2)], 1, 1.0, None),
                                        ('IPX', [('AISO4', 0)], 1, 2.0, None),
                                        ('INALD', [('AISO5', 0)], 1, 0.5, 'HNO3')]),
}
UPTAKE['CRACMM3'] = dict(const=GLY3, acid=UPTAKE['CRACMM2']['acid'])
# inorganic / NTR2 pathways: (reactant, product text, flag line), where used in the CMAQ mech definition
EXTRA = {
    'CB6R5':      [('NTR2', 'HNO3', '19 ; 1.4 ;'), ('N2O5', '2 HNO3', '20 ; 1 ; 108.01 ;')],
    'SAPRC07TC':  [('N2O5', '2 HNO3', '20 ; 1 ; 108.01 ;')],
    'SAPRC07TIC': [('N2O5', '2 HNO3', '20 ; 1 ; 108.01 ;'), ('NO3', 'HNO3', '3 ; 1e-3 ; 62.0 ; 0 ;')],
    'CRACMM2':    [('N2O5', '2 HNO3', '20 ; 1 ; 108.01 ;'), ('NO3', 'HNO3', '3 ; 1e-3 ; 62.0 ; 0 ;'), ('HO2', 'OffGas', '3 ; 0.2 ; 33.0 ; 0 ;')],
}
for _m in list(EXTRA):
    EXTRA[_m].append(('NO2', '0.5 HONO + 0.5 HNO3', '21 ; 1e-4 ;'))
EXTRA['CRACMM3'] = EXTRA['CRACMM2']
GAS_MW = {'GLY': 58.0, 'MGLY': 72.0, 'EPOX': 118.1, 'IEPOX': 118.1, 'IMAE': 102.0, 'IHMML': 102.0, 'IPX': 150.1, 'INALD': 195.1}

# CMAQ particle species -> ASP species (semivolatile V/A pairs are one species named like the vapor)
CB_MAP = {'AAVB2': 'SVAVB2', 'AAVB3': 'SVAVB3', 'AAVB4': 'SVAVB4', 'AISO1': 'SVISO1', 'AISO2': 'SVISO2', 'ASQT': 'SVSQT',
          'AMTNO3': 'MTNO3', 'AISOPNN': 'ISOPNN',
          'AXYL1': 'SVXYL1', 'AXYL2': 'SVXYL2', 'ATOL1': 'SVTOL1', 'ATOL2': 'SVTOL2', 'ABNZ1': 'SVBNZ1', 'ABNZ2': 'SVBNZ2',
          'ATRP1': 'SVTRP1', 'ATRP2': 'SVTRP2', 'APAH1': 'SVPAH1', 'APAH2': 'SVPAH2', 'AALK1': 'SVALK1', 'AALK2': 'SVALK2',
          'ATRPN': 'VTRPN', 'AHONIT': 'VHONIT', 'AHOM': 'VHOM'}
K_OLIG, K_HYD, K_TRPN = 9.48816e-6, 9.259e-5, 9.26e-5
OLIG = {   # reactant -> product mole coefficient per mole (from the mech_*.def files)
    'CB6R5': [('AAVB2', .907, 'AOLGA'), ('AAVB3', .925, 'AOLGA'), ('AAVB4', .943, 'AOLGA'), ('AISO1', .5, 'AOLGB'),
              ('AISO2', .5, 'AOLGB'), ('ASQT', 1.5, 'AOLGB')],
    'SAPRC07TIC': [('AISO1', .5, 'AOLGB'), ('AISO2', .5, 'AOLGB'), ('ASQT', 1.5, 'AOLGB'), ('AAVB2', .907, 'AOLGA'),
                   ('AAVB3', .925, 'AOLGA'), ('AAVB4', .943, 'AOLGA')],
    'SAPRC07TC': [('AXYL1', .8571, 'AOLGA'), ('AXYL2', 1.1429, 'AOLGA'), ('ATOL1', .8571, 'AOLGA'), ('ATOL2', 1.1429, 'AOLGA'),
                  ('ABNZ1', .7143, 'AOLGA'), ('ABNZ2', .7143, 'AOLGA'), ('ATRP1', .8, 'AOLGB'), ('ATRP2', .9, 'AOLGB'),
                  ('AISO1', .5, 'AOLGB'), ('AISO2', .5, 'AOLGB'), ('ASQT', 1.5, 'AOLGB'), ('APAH1', 1.4286, 'AOLGA'),
                  ('APAH2', 1.4286, 'AOLGA'), ('AALK1', 1.7143, 'AOLGA'), ('AALK2', 1.7143, 'AOLGA')],
}
HYD = {   # reactant -> [(coef, product)], rate
    'CB6R5': [('AMTNO3', [(1.0, 'AMTHYD')])],
    'SAPRC07TIC': [('AMTNO3', [(1.0, 'HNO3'), (1.0, 'AMTHYD')]), ('AISOPNN', [(2.0, 'HNO3'), (0.5, 'AMTHYD')])],
}
CRACMM_TRPN = [('ATRPN', [(1.0, 'AHOM'), (1.0, 'HNO3')]), ('AHONIT', [(1.0, 'AHOM'), (1.0, 'HNO3')])]


def _mw(mech):
    d = NML[mech]
    f = glob.glob(MECHS + d + '/AE_*.nml')[0]
    return {k.strip(): v for k, v in cmaq_ae.read_mw(f).items()}


def _names(path, sep=';'):
    out = []
    for ln in open(path, errors='ignore').read().replace('\r', '').split('\n'):
        s = ln.strip()
        if s and not s.startswith('!') and sep in s:
            out.append(re.split(r'[\s;]', s)[0])
    return out


def _append(path, lines, before=None):
    L = open(path, errors='ignore').read().replace('\r', '').split('\n')
    while L and not L[-1].strip():
        L.pop()
    if before:
        i = next(k for k, l in enumerate(L) if l.startswith(before))
        L[i:i] = lines
    else:
        L.extend(lines)
    open(path, 'w').write('\n'.join(L) + '\n')


def _add_species(case_dir, new, mwd):
    gas = _names(os.path.join(case_dir, 'GasPhaseChems.in'))
    org = _names(os.path.join(case_dir, 'OrgPhaseChems.in'))
    gas_new = [n for n in new if n not in gas]
    org_new = [n for n in new if n not in org]
    if gas_new:
        _append(os.path.join(case_dir, 'GasPhaseChems.in'),
                ['%s\t;1\t;0;%g;0;0;0!' % (n, mwd[n]) for n in gas_new], before='O2\t')
    if not org_new:
        return
    part = [(n, mwd[n], ASOA_CSTAR, ASOA_ENTH, True) for n in org_new]
    _append(os.path.join(case_dir, 'OrgPhaseChems.in'),
            ['!CMAQ particle-phase products of uptake / particle-phase reactions (nonvolatile, one ASP species with its gas partner)'] +
            ['%-14s;1\t;%.2f\t; %.1f ;       1.0; 1.0;' % (n, mw, cmaq_oa.RHO_ORG) for n, mw, c, h, nv in part])
    for fn in ('HydrophilicOrgChems.in', 'HydrophilicOrganicDissolution.in'):
        rows = []
        for n, mw, c, h, nv in part:
            if fn == 'HydrophilicOrgChems.in':
                rows.append('%-8s;1      ;%-9.1f ; 1.5 ;  %d ;0.04; 1.00e7 ; 1. ; 1.0;' % (n, mw, max(1, int(round(mw / 14.0)))))
            else:
                rows.append('%s <=> %s ;  0.1 ; 1.0e-16 ; 15.   ; 0.      ; 0.5  ; 0.     ; 0.5 ;' % (n, n))
        _append(os.path.join(case_dir, fn), rows)
    rows = ['! CMAQ uptake / particle-phase products: nonvolatile (C*(298) = %g ug/m3)' % ASOA_CSTAR]
    for n, mw, c, h, nv in part:
        rows.append('%s <=> %s\t; %.1f;\t1.\t; %.4e\t; %.1f\t; ! C*298 = %.3g' % (n, n, cmaq_oa.MASS_ACC, cmaq_oa.psat300_mbar(c, h, mw), h / 1e3, c))
    _append(os.path.join(case_dir, 'OrganicDissolution.in'), rows)


def _strip_block(path):
    L = open(path, errors='ignore').read().replace('\r', '').split('\n')
    if BEGIN in L:
        i, j = L.index(BEGIN), L.index(END)
        L = L[:i] + L[j + 1:]
    while L and not L[-1].strip():
        L.pop()
    return L


def spec_extra(mech):
    return EXTRA.get(mech, [])


def apply(case_dir, mech, acid=True, iepox_gamma=0.01):
    """acid=True: flag 18 (particle-acidity-dependent) for the epoxides; False: constant gamma `iepox_gamma` (flag 3)"""
    mech = mech.upper()
    mwd = _mw(mech)
    mwd.update({'AHOM': mwd.get('AHOM', 250.0)})
    gas = set(_names(os.path.join(case_dir, 'GasPhaseChems.in')))
    gm = os.path.join(case_dir, 'GasChemicalMechanism.in')
    L = _strip_block(gm)
    nrxn = sum(1 for l in L if l.strip() and not l.strip().startswith('!'))
    # index of the NO2 photolysis reaction (for the AHOM photolysis rate)
    i_no2 = None
    k = 0
    for l in L:
        s = l.strip()
        if s and not s.startswith('!'):
            k += 1
            if i_no2 is None and re.match(r'^NO2\s*=>\s*(NO\s*\+\s*O|O\s*\+\s*NO)\s*;\s*0\s*;', s):
                i_no2 = k
    spec = UPTAKE[mech]
    new, rxn, notes, skipped = [], [], [], []
    for r, prod, scale in spec['const']:
        if r not in gas:
            skipped.append(r)
            continue
        g = GAMMA_GLY * scale
        rxn.append('%s => %s ; 3 ; %.5g ; %g ; 0 ;' % (r, prod, g, GAS_MW[r]))
        new.append(prod)
        notes.append('%s -> %s: constant gamma %.3g (CMAQ GAMMA_GLY%s)' % (r, prod, g, ' x 0.09' if scale != 1 else ''))
    for r, prods, kind, scale, extra in spec['acid']:
        if r not in gas:
            skipped.append(r)
            continue
        for p, ch in prods:
            ex = (' + %s' % extra) if extra else ''
            if acid:
                rxn.append('%s => %s%s ; 18 ; %g ; %g ; %d ; %d ;' % (r, p, ex, scale, GAS_MW[r], kind, ch))
                notes.append('%s -> %s%s: acidity-dependent (flag 18, kind %d, channel %d, scale %g)' % (r, p, ex, kind, ch, scale))
            elif ch in (0, 1):      # constant gamma: one product only (the non-sulfate channel)
                rxn.append('%s => %s%s ; 3 ; %.5g ; %g ; 0 ;' % (r, p, ex, iepox_gamma * scale, GAS_MW[r]))
                notes.append('%s -> %s%s: constant gamma %.3g' % (r, p, ex, iepox_gamma * scale))
            new.append(p)
    for r, prod, tail in spec_extra(mech):
        if r not in gas:
            skipped.append(r)
            continue
        rxn.append('%s => %s ; %s' % (r, prod, tail))
        notes.append('%s -> %s: %s' % (r, prod, 'flag ' + tail.split(';')[0].strip()))
    # tier 2
    org = []
    pm2asp = CB_MAP
    for r, c, p in OLIG.get(mech, []):
        org.append((pm2asp[r], K_OLIG, [(c, p)], 'oligomerisation'))
    for r, prods in HYD.get(mech, []):
        org.append((pm2asp[r], K_HYD, prods, 'nitrate hydrolysis'))
    if mech.startswith('CRACMM'):
        for r, prods in CRACMM_TRPN:
            org.append((pm2asp[r], K_TRPN, prods, 'organic nitrate -> AHOM'))
        if mech == 'CRACMM3' and i_no2:
            org.append(('VHOM', 'PHOTO %d 0.01' % i_no2, [(0.8, 'AMTN1'), (0.2, 'ORA1')], 'AHOM photolysis, 0.01 x jNO2'))
    org_rows = []
    org_have = set(_names(os.path.join(case_dir, 'OrgPhaseChems.in')))
    for r, k, prods, why in org:
        if r not in org_have:
            skipped.append(r)
            continue
        for _, p in prods:
            if p not in gas and p not in ('HNO3', 'ORA1'):
                new.append(p)
        org_rows.append('%-8s ; %s ; %s   ! %s' % (r, ('%.6g' % k) if not isinstance(k, str) else k,
                                                  ' + '.join('%.4g %s' % (c, p) for c, p in prods), why))
    new = [n for n in dict.fromkeys(new) if n in mwd or n in gas]
    _add_species(case_dir, [n for n in new if n in mwd and n not in ('HNO3', 'ORA1')], mwd)
    out = [BEGIN, '!! CMAQ gas->particle uptake (tier 1), see cmaq_het.py / NewMechanisms_README.md',
           '!! flag 3: gamma, MW; flag 18: scale, MW, kind, channel (acid-dependent epoxide uptake)'] + \
          ['!! ' + n for n in notes] + rxn + [END]
    open(gm, 'w').write('\n'.join(L + out) + '\n')
    p = os.path.join(case_dir, 'OrgPhaseReactions.in')
    with open(p, 'w') as fh:
        fh.write('!! OrgPhaseReactions.in - first-order particle-phase organic reactions of CMAQ (%s)\n' % mech)
        fh.write('!! Reactant ; rate [1/s] (or PHOTO <gas rxn #> <scale>) ; products (mole coefficients)\n')
        fh.write('!! Needs common_asp with OrgPhaseReactions.h.  Products that are not organic-phase species go to the gas phase.\n')
        fh.write('\n'.join(org_rows) + '\n')
    return dict(reactions=rxn, particle_reactions=org_rows, skipped=sorted(set(skipped)), new_species=new, notes=notes)


if __name__ == '__main__':
    import sys
    print(apply(sys.argv[1], sys.argv[2], acid=(len(sys.argv) < 4 or sys.argv[3] != 'const')))
