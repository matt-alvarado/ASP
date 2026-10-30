"""Organic aerosol (POA / SOA) and extra heterogeneous / particle-phase chemistry for the GEOS-Chem and MOZART-T1 ASP decks.

Post-processes a finished <Case>_GCHEM / <Case>_MOZT1 deck in place (as CB6r5Translation/cmaq_het.py does for the CMAQ decks):

  apply_gc(case_dir, case)   GEOS-Chem 14.8.0 'complex SOA with semivolatile POA' (carbon_mod.F90, SOA_PARA_INIT):
        TSOA0-3 (C* 0.1/1/10/100 ug/m3 at 298 K), ASOA1-3 (C* 1/10/100) + ASOAN (non-volatile), semivolatile POA (C* 1646 / 20 at 300 K)
        and its oxidation product OPOA (C* x 0.01).  Terpene (MTPA, MTPO, LIMO) and aromatic (BENZ, TOLU, XYLE, NAP) yields with the same
        NO / HO2 branching as GEOS-Chem (BETANO = k_NO[NO] / (k_NO[NO] + k_HO2[HO2])); IVOC emissions as NAP.  The baseline counter ->
        IVOC ladder, the IVOC aging block and the baseline POA rows are removed.
        Heterogeneous uptake products (SOAGX, SOAIE, LVOCOA, IONITA, MONITA, AONITA) and IONITA / MONITA hydrolysis are added.
  apply_t1(case_dir, case)   CESM2 / CAM-chem MOZART-T1 VBS-SOA (Tilmes et al. 2019, Table 1): SOAG0-4 (C* 0.01 ... 100 ug/m3 at 300 K)
        from isoprene, monoterpenes, beta-caryophyllene, aromatics and S/IVOC; SOA photolysis J = 0.04 J(NO2); glyoxal uptake -> SOAG0.

Everything that GEOS-Chem / CAM-chem do not publish in the sources available here is a named constant below (ASSUMED_*).
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'CB6r5Translation'))
import cmaq_oa      # noqa: E402
import cmaq_het     # noqa: E402

R = 8.314462618
BEGIN, END = '!! MOZ_OA_BEGIN', '!! MOZ_OA_END'
TP = {'Williams': (288.4, 880.0), 'Otavi': (290.0, 730.0), 'Alaska': (280.6, 846.5), 'Timbavati': (288.0, 900.0)}
BASE_ORG = re.compile(r'^(IVOC\d|POA\d|LEVO|CBIO|CPD\d)$')
COUNTERS = ('TOLSOA', 'XYLSOA', 'ISOSOA', 'TERPSOA')

# ---------------------------------------------------------------- GEOS-Chem (carbon_mod.F90, v14.8.0)
GC_MW_SOA = 150.0                       # species_database.yml: TSOA3 / SOAS / ASOAN MW_g = 150; one value for all, so that ASP's molar Raoult's law equals
#                                         GEOS-Chem's mass-based absorption
GC_HEAT_VAPOR_K = 5.0e3                 # HEAT_VAPOR [K]  -> dHvap = 41.6 kJ/mol
GC_DHVAP = GC_HEAT_VAPOR_K * R
# C* [ug/m3] at the reference temperature (298 K; 300 K for POA / OPOA), product index 1..4 of SOA_PARA_INIT
GC_TERP_NAMES = ['TSOA1', 'TSOA2', 'TSOA3', 'TSOA0']          # KOM_REF(1..4) = 1/1, 1/10, 1/100, 1/0.1
GC_ARO_NAMES = ['ASOA1', 'ASOA2', 'ASOA3', 'ASOAN']           # KOM_REF(4) = 1e10 (non-volatile)
GC_PART = [('TSOA0', 0.1, 298.0, False), ('TSOA1', 1.0, 298.0, False), ('TSOA2', 10.0, 298.0, False), ('TSOA3', 100.0, 298.0, False),
           ('ASOA1', 1.0, 298.0, False), ('ASOA2', 10.0, 298.0, False), ('ASOA3', 100.0, 298.0, False), ('ASOAN', 1e-10, 298.0, True),
           ('SVPOA1', 1646.0, 300.0, False), ('SVPOA2', 20.0, 300.0, False), ('OPOA1', 16.46, 300.0, False), ('OPOA2', 0.2, 300.0, False)]
# mass yields ALPHA(branch, product, parent): high NOx (RO2 + NO), low NOx (RO2 + HO2), NO3 route
GC_TERP = {
    'MT': dict(parents=['MTPA', 'MTPO'], high=[0.0095, 0.0900, 0.0150, 0.0400], low=[0.019, 0.180, 0.030, 0.080], no3=[0.0, 0.3207, 1.0830, 0.0]),
    'LI': dict(parents=['LIMO'], high=[0.4743, 0.1174, 1.4190, 0.0], low=[0.3661, 0.3214, 0.8168, 0.0], no3=[0.0, 0.3207, 1.0830, 0.0]),
}
# aromatic RO2 (name of the RO2 in fullchem) -> parent MW, ARO2CARB, high NOx (ASOA1-3), HO2 route (ASOAN)
GC_ARO = {'BRO2': dict(mw=78.12, carb=72.0 / 78.0, high=[0.0778, 0.0, 0.7932], ho2=0.37),
          'TRO2': dict(mw=92.15, carb=84.0 / 92.0, high=[0.0315, 0.0944, 0.0800], ho2=0.30),
          'XRO2': dict(mw=106.18, carb=96.0 / 106.0, high=[0.0250, 0.0360, 0.0899], ho2=0.36),
          'NRO2': dict(mw=128.18, carb=120.0 / 128.0, high=[0.0387, 0.2956, 0.2349], ho2=0.73)}
GC_K_NO, GC_K_HO2 = '2.6e-12; 0; 350; 1;', '1.4e-12; 0; 700; 1;'      # AARO2NO, BBRO2NO ; AARO2HO2, BBRO2HO2
GC_POA_SPLIT = (0.49, 0.51)             # ALPHA(1, 1:2, PARENTPOA)
GC_POA_OH = 2.0e-11                     # POA + OH (carbon_mod.F90: RK = 2e-11 [OH])
GC_OMOC_RATIO = 2.1 / 1.4               # OM/OC of OPOA / POA: OPOA mass per POA mass oxidised (GEOS-Chem works in OC mass)
GC_NAP_MW = 128.18
# uptake products (fullchem.eqn) -> (ASP product species, MW).  IONITA / MONITA are tracked as moles of N in GEOS-Chem (MW_g 14.01): here mole for mole.
GC_PROD = {'GLYX': 'SOAGX', 'MGLY': '1.5 SOAGX', 'PYAC': '1.5 SOAGX', 'IEPOXA': 'SOAIE', 'IEPOXB': 'SOAIE', 'IEPOXD': 'SOAIE', 'ICHE': 'SOAIE',
           'HMML': '0.8 SOAIE', 'LVOC': 'LVOCOA', 'MVKN': 'IONITA', 'R4N2': 'IONITA', 'IHN1': 'IONITA', 'IHN2': 'IONITA', 'IHN3': 'IONITA',
           'IHN4': 'IONITA', 'INPD': 'IONITA', 'INPB': 'IONITA', 'IDN': 'IONITA', 'ITCN': 'IONITA', 'ITHN': 'IONITA', 'MCRHNB': 'IONITA',
           'MCRHN': 'IONITA', 'MONITS': 'MONITA', 'MONITU': 'MONITA', 'HONIT': 'MONITA', 'NPHEN': 'AONITA'}
GC_PROD_MW = {'SOAGX': 58.04, 'SOAIE': 118.15, 'LVOCOA': 154.19, 'IONITA': 147.13, 'MONITA': 215.0, 'AONITA': 189.12, 'INDIOL': 102.0}
GC_HYD = [('IONITA', 2.78e-4), ('MONITA', 2.78e-4)]        # IONITA / MONITA = INDIOL + HNO3 (fullchem.eqn)

# ---------------------------------------------------------------- MOZART-T1 / CAM-chem VBS-SOA (Tilmes et al. 2019, JAMES 11, Table 1)
T1_CSTAR = [0.01, 0.1, 1.0, 10.0, 100.0]                    # SOAG0..4 at 300 K
ASSUMED_SOAG_MW = 150.0       # NOT published in the paper (CAM-chem uses a per-species MW in chem_mech.in); the yields are molar
ASSUMED_SOAG_DHVAP = 30.0e3   # J/mol, NOT published; VBS convention
ASSUMED_IVOC_MW, ASSUMED_SVOC_MW, ASSUMED_POM_MW = 184.0, 250.0, 250.0     # NOT published
T1_YIELD = {   # (parent group, oxidant) -> molar yields to SOAG0..4
    ('ISOP', 'OH'): [0.0031, 0.0035, 0.0003, 0.0271, 0.0474], ('ISOP', 'O3'): [0, 0, 0, 0.0033, 0], ('ISOP', 'NO3'): [0, 0, 0, 0.059024, 0.025024],
    ('MTERP', 'OH'): [0.0508, 0.1149, 0.0348, 0.0554, 0.1278], ('MTERP', 'O3'): [0.0508, 0.1149, 0.0348, 0.0554, 0.1278],
    ('MTERP', 'NO3'): [0, 0, 0, 0.17493, 0.59019],
    ('BCARY', 'OH'): [0.2202, 0.2067, 0.0653, 0.1284, 0.114], ('BCARY', 'O3'): [0.2202, 0.2067, 0.0653, 0.1284, 0.114],
    ('BCARY', 'NO3'): [0, 0, 0, 0.17493, 0.59019],
    ('XYLENES', 'OH'): [0.1677, 0.0174, 0.0860, 0.0512, 0.1598], ('TOLUENE', 'OH'): [0.1364, 0.0101, 0.0763, 0.2157, 0.0738],
    ('BENZENE', 'OH'): [0.0023, 0.0008, 0.0843, 0.0443, 0.1621],
}
T1_IVOC_YIELD = [0.2381, 0.1308, 0.0348, 0.0076, 0.0113]
T1_SVOC_YIELD = [0.5931, 0.1534, 0.0459, 0.0085, 0.0128]
T1_K_IVOC_OH = '1.34e-11; 0; 0; 1;'
T1_GROUP = {'APIN': 'MTERP', 'BPIN': 'MTERP', 'LIMON': 'MTERP', 'MYRC': 'MTERP', 'ISOP': 'ISOP', 'BCARY': 'BCARY', 'XYLENES': 'XYLENES',
            'TOLUENE': 'TOLUENE', 'BENZENE': 'BENZENE'}
T1_J_SOA = 0.04                  # J(SOA) = 0.04 x J(NO2), particulate SOA loss, no products


# ---------------------------------------------------------------- helpers
def _read(path):
    return open(path, errors='ignore').read().replace('\r', '').split('\n')


def _write(path, L):
    while L and not L[-1].strip():
        L.pop()
    open(path, 'w').write('\n'.join(L) + '\n')


def _gas_rows(case_dir):
    """GasPhaseChems.in -> dict name -> (perox, mw, ini, bg)"""
    out = {}
    for ln in _read(os.path.join(case_dir, 'GasPhaseChems.in')):
        s = ln.strip()
        if s and not s.startswith('!') and ';' in s and not s.startswith('END'):
            f = [x.strip() for x in s.rstrip('!').split(';')]
            try:
                out[f[0]] = (int(f[2]), float(f[3]), float(f[4]), float(f[5]))
            except (ValueError, IndexError):
                pass
    return out


def _edit_gas(case_dir, drop_re, add_rows, set_ini=None):
    """drop rows matching drop_re, add rows (name, perox, mw, ini, bg) before the O2 row, optionally add to the initial value of existing rows"""
    p = os.path.join(case_dir, 'GasPhaseChems.in')
    L = []
    for ln in _read(p):
        s = ln.strip()
        nm = re.split(r'[\s;]', s)[0] if s and not s.startswith('!') else ''
        if nm and drop_re.match(nm):
            continue
        if nm and set_ini and nm in set_ini:
            f = [x.strip() for x in s.rstrip('!').split(';')]
            f[4] = '%.6g' % (float(f[4]) + set_ini[nm])
            ln = '%s\t;%s\t;%s;%s;%s;%s;%s!' % tuple(f[:7])
        L.append(ln)
    i = next(k for k, l in enumerate(L) if l.startswith('O2\t'))
    new = ['%s\t;1\t;%d;%g;%.6g;%.6g;0!' % (n, pf, mw, ini, bg) for n, pf, mw, ini, bg in add_rows]
    L[i:i] = ['!', '! organic-aerosol species added by moz_oa.py'] + new
    _write(p, L)


def _psat300(c_ref, t_ref, dh, mw):
    """C* at t_ref -> ASP Psat(300 K) [mbar], CMAQ / VBS temperature form C*(T) = C*ref (Tref/T) exp[dH/R (1/Tref - 1/T)]"""
    c300 = c_ref * (t_ref / 300.0) * math.exp(dh / R * (1.0 / t_ref - 1.0 / 300.0))
    return c300 * 1e-6 * R * 300.0 / mw / 100.0, c300


def _write_org(case_dir, part, mw, dh_of, note):
    """part: [(name, cstar_ref, t_ref, nonvol)]; replaces the baseline organic rows (IVOCk, POAk, LEVO, CBIO, CPDk) in the five organic files"""
    def keep(ln):
        s = ln.strip()
        nm = re.split(r'[\s;<]', s)[0].lstrip('!') if s else ''
        return not (not s.startswith('!') and BASE_ORG.match(nm))
    for fn in ('OrgPhaseChems.in', 'HydrophilicOrgChems.in', 'HydrophilicOrganicDissolution.in', 'OrganicDissolution.in'):
        p = os.path.join(case_dir, fn)
        L = [l for l in _read(p) if keep(l)]
        while L and not L[-1].strip():
            L.pop()
        L.append('!')
        L.append('!' + note)
        for n, c, t, nv in part:
            m = mw[n]
            dh = 1.0e3 if nv else dh_of
            ps, c300 = _psat300(c, t, dh, m)
            if fn == 'OrgPhaseChems.in':
                L.append('%-14s;1\t;%.2f\t; %.1f ;       1.0; 1.0;' % (n, m, cmaq_oa.RHO_ORG))
            elif fn == 'HydrophilicOrgChems.in':
                L.append('%-8s;1      ;%-9.1f ; 1.5 ;  %d ;0.04; 1.00e7 ; 1. ; 1.0;' % (n, m, max(1, int(round(m / 14.0)))))
            elif fn == 'HydrophilicOrganicDissolution.in':
                L.append('%s <=> %s ;  0.1 ; 1.0e-16 ; 15.   ; 0.      ; 0.5  ; 0.     ; 0.5 ;' % (n, n))
            else:
                L.append('%s <=> %s\t; %.1f;\t1.\t; %.4e\t; %.1f\t; ! C*(%g K) = %.3g ug/m3, C*(300 K) = %.3g' % (n, n, cmaq_oa.MASS_ACC, ps, dh / 1e3, t, c, c300))
        _write(p, L)


def _strip_baseline_chem(case_dir):
    """remove the counter -> IVOC chains + IVOC aging block and the counter products; returns the reaction lines"""
    p = os.path.join(case_dir, 'GasChemicalMechanism.in')
    L = _read(p)
    for i, l in enumerate(L):
        if l.startswith('!! ---------------- 2. SOA precursor counters'):
            L = L[:i]
            break
    if BEGIN in L:
        i, j = L.index(BEGIN), L.index(END)
        L = L[:i] + L[j + 1:]
    out = []
    for l in L:
        if l.strip() and not l.startswith('!'):
            for c in COUNTERS:
                l = re.sub(r'\s*\+\s*%s\b' % c, '', l)
        out.append(l)
    return out


def _add_products(L, pattern, prod_text):
    """append ' + prod_text' to the product side of every reaction line whose reaction text (before the first ';') matches `pattern`"""
    n = 0
    for i, l in enumerate(L):
        if not l.strip() or l.startswith('!'):
            continue
        head, sep, tail = l.partition(';')
        if sep and re.match(pattern, head.strip()):
            L[i] = head.rstrip() + ' + ' + prod_text + ' ;' + tail
            n += 1
    return n


def _fmt(c):
    return ('%.5g' % c)


def _modes_totals(case_dir):
    modes = cmaq_oa.read_modes(os.path.join(case_dir, 'AerosolModes.in'))
    return modes, sum(m for n, m in modes.items() if BASE_ORG.match(n))


def _ppb_to_ug(ppb, mw, T, P):
    return ppb / cmaq_oa._ppb_per_ugm3(mw, T, P)


def _env_modes(case_dir, split, notes):
    """background (environmental) aerosol: baseline organics -> the new particle species with the fractions in `split`"""
    p = os.path.join(case_dir, 'EnvAerosolModes.in')
    modes = cmaq_oa.read_modes(p)
    tot = sum(m for n, m in modes.items() if BASE_ORG.match(n))
    if tot > 0:
        cmaq_oa.rewrite_modes(p, {n for n in modes if BASE_ORG.match(n)}, {k: v * tot for k, v in split.items()})
        notes.append('background aerosol organics %.4g ug/m3 -> %s' % (tot, ', '.join(split)))


# ---------------------------------------------------------------- GEOS-Chem
def apply_gc(case_dir, case, mw_soa=GC_MW_SOA):
    T, P = TP[case]
    gas = _gas_rows(case_dir)
    L = _strip_baseline_chem(case_dir)
    notes = []
    # ---- initial organics: baseline POA / IVOC / particle organics -> SVPOA1/2 (gas + particle), NAP (IVOC)
    modes, mp_ug = _modes_totals(case_dir)
    svoc_gas_ug = nap_ppb = 0.0
    for n, (pf, mw, ini, bg) in gas.items():
        m = re.match(r'^IVOC(\d)$', n)
        if m and ini > 0:
            if int(m.group(1)) <= 6:
                svoc_gas_ug += _ppb_to_ug(ini, mw, T, P)
            else:
                nap_ppb += ini * mw / GC_NAP_MW
                notes.append('gas %s %.4g ppb (MW %.0f) -> NAP %.4g ppb (IVOC surrogate)' % (n, ini, mw, ini * mw / GC_NAP_MW))
    part_new = {'SVPOA1': mp_ug * GC_POA_SPLIT[0], 'SVPOA2': mp_ug * GC_POA_SPLIT[1]}
    gas_ini = {'SVPOA1': svoc_gas_ug * GC_POA_SPLIT[0] * cmaq_oa._ppb_per_ugm3(mw_soa, T, P),
               'SVPOA2': svoc_gas_ug * GC_POA_SPLIT[1] * cmaq_oa._ppb_per_ugm3(mw_soa, T, P)}
    notes.append('particle organics %.4g ug/m3 -> SVPOA1 %.4g + SVPOA2 %.4g ug/m3; gas SVOC (IVOC bins <= 6) %.4g ug/m3' %
                 (mp_ug, part_new['SVPOA1'], part_new['SVPOA2'], svoc_gas_ug))
    cmaq_oa.rewrite_modes(os.path.join(case_dir, 'AerosolModes.in'), {n for n in modes if BASE_ORG.match(n)}, part_new)
    _env_modes(case_dir, {'SVPOA1': GC_POA_SPLIT[0], 'SVPOA2': GC_POA_SPLIT[1]}, notes)
    # ---- species files
    mw = {n: mw_soa for n, c, t, nv in GC_PART}
    _write_org(case_dir, [(n, c, t, nv) for n, c, t, nv in GC_PART], mw, GC_DHVAP,
               'GEOS-Chem 14.8.0 complex SOA with semivolatile POA (carbon_mod.F90 SOA_PARA_INIT): C*, one MW (mass-based absorption), dHvap = 5000 K x R')
    add = [(n, 0, mw_soa, 0.0, 0.0) for n, c, t, nv in GC_PART] + [('SRO2MT', 0, mw_soa, 0.0, 0.0), ('SRO2LI', 0, mw_soa, 0.0, 0.0)]
    _edit_gas(case_dir, re.compile(r'^(IVOC\d(_R1O2)?|POA\d|LEVO|CBIO|CPD\d|TOLSOA|XYLSOA|ISOSOA|TERPSOA)$'), add,
              set_ini={'NAP': nap_ppb} if nap_ppb else None)
    # gas-phase initial SVPOA: written into the rows just added
    p = os.path.join(case_dir, 'GasPhaseChems.in')
    G = _read(p)
    for i, l in enumerate(G):
        for n, v in gas_ini.items():
            if l.startswith(n + '\t'):
                G[i] = '%s\t;1\t;0;%g;%.6g;0;0!' % (n, mw_soa, v)
    _write(p, G)
    # ---- reactions
    rxn = ['!! GEOS-Chem 14.8.0 complex SOA (carbon_mod.F90): yields are mass yields x MW(parent)/MW(SOA); NO / HO2 branching as BETANO',
           '!! terpene OH / O3 oxidation -> dummy RO2 (SRO2MT = MTPA + MTPO, SRO2LI = LIMO) that reacts with NO (high-NOx alphas) or HO2 (low-NOx alphas)']
    for key, d in GC_TERP.items():
        sro2 = 'SRO2' + key
        mwp = 136.0
        for par in d['parents']:
            n1 = _add_products(L, r'^%s \+ (OH|O3) =>' % par, sro2)
            notes.append('%s + OH / O3: +%s (%d reactions)' % (par, sro2, n1))
            y = ' + '.join('%s %s' % (_fmt(a * mwp / mw_soa), nm) for a, nm in zip(d['no3'], GC_TERP_NAMES) if a > 0)
            n2 = _add_products(L, r'^%s \+ NO3 =>' % par, y)
            notes.append('%s + NO3: +%s' % (par, y))
        for lab, al, k in (('NO', d['high'], GC_K_NO), ('HO2', d['low'], GC_K_HO2)):
            y = ' + '.join('%s %s' % (_fmt(a * mwp / mw_soa), nm) for a, nm in zip(al, GC_TERP_NAMES) if a > 0)
            rxn.append('%s + %s => %s + %s ; %s' % (sro2, lab, lab, y, k))
    rxn.append('!! aromatics: the SOA yields ride on the fullchem RO2 + NO / RO2 + HO2 loss reactions (GEOS-Chem: DARO2 = loss of LBRO2N / LBRO2H ...)')
    for ro2, d in GC_ARO.items():
        f = d['mw'] * d['carb'] / mw_soa
        hi = ' + '.join('%s %s' % (_fmt(a * f), nm) for a, nm in zip(d['high'], GC_ARO_NAMES) if a > 0)
        n1 = _add_products(L, r'^%s \+ NO =>' % ro2, hi)
        n2 = _add_products(L, r'^%s \+ HO2 =>' % ro2, '%s ASOAN' % _fmt(d['ho2'] * f))
        notes.append('%s + NO: +%s; + HO2: +%s ASOAN (%d / %d reactions)' % (ro2, hi, _fmt(d['ho2'] * f), n1, n2))
    rxn.append('!! semivolatile POA oxidation (gas phase POG + OH, k = 2e-11) -> OPOA; OPOA mass = POA mass x (OM/OC 2.1 / 1.4)')
    for a, b in (('SVPOA1', 'OPOA1'), ('SVPOA2', 'OPOA2')):
        rxn.append('%s + OH => %s %s + OH ; %.4g; 0; 0; 1;' % (a, _fmt(GC_OMOC_RATIO), b, GC_POA_OH))
    # ---- uptake products (replace the dropped 'M' product of the flag-3 uptake lines)
    prods = []
    for i, l in enumerate(L):
        m = re.match(r'^(\w+) => M ; 3;', l)
        if m and m.group(1) in GC_PROD:
            L[i] = l.replace('=> M ;', '=> %s ;' % GC_PROD[m.group(1)], 1)
            prods.append(GC_PROD[m.group(1)].split()[-1])
    prods = list(dict.fromkeys(prods))
    notes.append('uptake products added: ' + ', '.join(prods))
    L += ['', BEGIN] + rxn + [END]
    _write(os.path.join(case_dir, 'GasChemicalMechanism.in'), L)
    newsp = [n for n in prods if n in GC_PROD_MW] + ['INDIOL']
    cmaq_het._add_species(case_dir, [n for n in dict.fromkeys(newsp)], GC_PROD_MW)
    # ---- particle-phase reactions
    rows = ['!! OrgPhaseReactions.in - GEOS-Chem 14.8.0 particle-phase reactions',
            '!! Reactant ; rate [1/s] ; products (mole coefficients).  Needs common_asp with OrgPhaseReactions.h.']
    for r, k in GC_HYD:
        rows.append('%-8s ; %.4g ; 1 INDIOL + 1 HNO3   ! organic nitrate hydrolysis (fullchem.eqn)' % (r, k))
    open(os.path.join(case_dir, 'OrgPhaseReactions.in'), 'w').write('\n'.join(rows) + '\n')
    return notes


# ---------------------------------------------------------------- MOZART-T1
def apply_t1(case_dir, case, mw_soag=ASSUMED_SOAG_MW, dh=ASSUMED_SOAG_DHVAP):
    T, P = TP[case]
    gas = _gas_rows(case_dir)
    L = _strip_baseline_chem(case_dir)
    notes = []
    # ---- initial organics: particle organic -> POM (non-volatile); measured gas IVOC bins -> SVOC / IVOC species
    modes, mp_ug = _modes_totals(case_dir)
    svoc_ug = ivoc_ug = 0.0
    measured = False
    for n, (pf, mw, ini, bg) in gas.items():
        m = re.match(r'^IVOC(\d)$', n)
        if m and ini > 0:
            measured = True
            u = _ppb_to_ug(ini, mw, T, P)
            if int(m.group(1)) <= 6:
                svoc_ug += u
            else:
                ivoc_ug += u
    if not measured:
        svoc_ug = 0.6 * mp_ug                                  # SVOC = additional 60 % of POM (Tilmes et al. 2019)
        nm_ug = 0.0
        for n, (pf, mw, ini, bg) in gas.items():
            if n in ('CO', 'CO2', 'CH4', 'H2', 'Trace', 'O3', 'NO', 'NO2', 'HNO3', 'SO2', 'NH3', 'HCl', 'H2O2', 'HONO', 'PAN', 'N2O5', 'N2O', 'HCN', 'H2SO4', 'NO3', 'N2', 'O2', 'OffGas', 'H2O', 'IVOC', 'SVOC', 'POM'):
                continue
            if ini > 0:
                nm_ug += _ppb_to_ug(ini, mw, T, P)
        ivoc_ug = 0.2 * nm_ug                                  # IVOC = additional 20 % of NMVOC
        notes.append('no measured IVOC bins: SVOC = 0.6 x POM = %.4g ug/m3, IVOC = 0.2 x NMVOC = %.4g ug/m3 (Tilmes et al. 2019)' % (svoc_ug, ivoc_ug))
    else:
        notes.append('measured IVOC bins: gas bins <= 6 -> SVOC %.4g ug/m3, bins >= 7 -> IVOC %.4g ug/m3' % (svoc_ug, ivoc_ug))
    gas_ini = {'SVOC': svoc_ug * cmaq_oa._ppb_per_ugm3(ASSUMED_SVOC_MW, T, P), 'IVOC': ivoc_ug * cmaq_oa._ppb_per_ugm3(ASSUMED_IVOC_MW, T, P)}
    part_new = {'POM': mp_ug}
    notes.append('particle organics %.4g ug/m3 -> POM (non-volatile, MAM4 POM)' % mp_ug)
    cmaq_oa.rewrite_modes(os.path.join(case_dir, 'AerosolModes.in'), {n for n in modes if BASE_ORG.match(n)}, part_new)
    _env_modes(case_dir, {'POM': 1.0}, notes)
    names = ['SOAG%d' % k for k in range(5)]
    mw = {n: mw_soag for n in names}
    mw['POM'] = ASSUMED_POM_MW
    part = [(n, c, 300.0, False) for n, c in zip(names, T1_CSTAR)] + [('POM', 1e-10, 300.0, True)]
    _write_org(case_dir, part, mw, dh, 'MOZART-T1 / CESM2 VBS-SOA (Tilmes et al. 2019): SOAG0-4 C* = 0.01 ... 100 ug/m3 at 300 K; MW and dHvap ASSUMED, POM non-volatile')
    add = [(n, 0, mw_soag, 0.0, 0.0) for n in names] + [('POM', 0, ASSUMED_POM_MW, 0.0, 0.0),
          ('SVOC', 0, ASSUMED_SVOC_MW, gas_ini['SVOC'], 0.0), ('IVOC', 0, ASSUMED_IVOC_MW, gas_ini['IVOC'], 0.0)]
    _edit_gas(case_dir, re.compile(r'^(IVOC\d(_R1O2)?|POA\d|LEVO|CBIO|CPD\d|TOLSOA|XYLSOA|ISOSOA|TERPSOA)$'), add)
    # ---- reactions
    rxn = ['!! MOZART-T1 / CESM2 VBS-SOA (Tilmes et al. 2019, Table 1): molar yields to SOAG0-4 on the parent + oxidant reactions']
    for par, grp in T1_GROUP.items():
        for ox in ('OH', 'O3', 'NO3'):
            y = T1_YIELD.get((grp, ox))
            if not y:
                continue
            txt = ' + '.join('%s %s' % (_fmt(a), nm) for a, nm in zip(y, names) if a > 0)
            n1 = _add_products(L, r'^%s \+ %s =>' % (par, ox), txt)
            notes.append('%s + %s: +SOAG yields (%d reactions)' % (par, ox, n1))
    for sp, y in (('IVOC', T1_IVOC_YIELD), ('SVOC', T1_SVOC_YIELD)):
        rxn.append('%s + OH => OH + %s ; %s' % (sp, ' + '.join('%s %s' % (_fmt(a), nm) for a, nm in zip(y, names)), T1_K_IVOC_OH))
    # glyoxal uptake -> SOAG0 (Knote et al. 2014)
    for i, l in enumerate(L):
        m = re.match(r'^GLYOXAL => (CO2|M)\s*;\s*3;', l)
        if m:
            L[i] = l.replace('=> %s' % m.group(1), '=> SOAG0', 1)
            notes.append('GLYOXAL uptake: product CO2 -> SOAG0 (1 mol per mol, assumed)')
    L += ['', BEGIN] + rxn + [END]
    _write(os.path.join(case_dir, 'GasChemicalMechanism.in'), L)
    # ---- SOA photolysis J = 0.04 J(NO2)
    k = 0
    i_no2 = None
    for l in L:
        s = l.strip()
        if s and not s.startswith('!'):
            k += 1
            if i_no2 is None and re.match(r'^NO2\s*=>\s*(NO\s*\+\s*O|O\s*\+\s*NO)\s*;\s*0\s*;', s):
                i_no2 = k
    rows = ['!! OrgPhaseReactions.in - MOZART-T1 VBS-SOA: particulate SOA photolysis J = 0.04 x J(NO2), no products (Tilmes et al. 2019)',
            '!! Reactant ; rate (PHOTO <gas reaction #> <scale>) ; products.  Needs common_asp with OrgPhaseReactions.h.']
    if i_no2:
        for n in names:
            rows.append('%-8s ; PHOTO %d %.3g ;   ! SOA photolysis' % (n, i_no2, T1_J_SOA))
    else:
        notes.append('WARNING: NO2 photolysis reaction not found; no SOA photolysis')
    open(os.path.join(case_dir, 'OrgPhaseReactions.in'), 'w').write('\n'.join(rows) + '\n')
    return notes


if __name__ == '__main__':
    d, case, kind = sys.argv[1], sys.argv[2], sys.argv[3]
    for n in (apply_gc if kind == 'gc' else apply_t1)(d, case):
        print(n)
