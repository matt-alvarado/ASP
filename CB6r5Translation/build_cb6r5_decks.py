#!/usr/bin/env python3
"""Build InorganicInputDecks/<Case>_CB6r5/ for the four cases.

  * GasChemicalMechanism.in  - CB6r5 core gas chemistry (translated from CMAQ mech_cb6r5_ae7_aq.def),
                               then the SOA-precursor -> IVOC chains and IVOC aging block of the case's old deck
  * GasPhaseChems.in         - CB6r5 species; initial / background values from the papers' original compound
                               tables, mapped with species_map.py
  * PhotoRates.in            - 51-column table re-derived from the case's TUV table (CB6r5 column set)
  * everything else          - copied unchanged from the old deck
Also writes SPECIES_MAPPING.csv in each new deck directory and CB6r5_README.md next to this script.
"""
import csv
import json
import math
import os
import re
import shutil
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cmaq2asp as c2a                      # noqa: E402
import species_map as sm                    # noqa: E402

DECKS = os.path.join(HERE, '..', 'InorganicInputDecks')
DEF = '/home/claude/usepa/cmaq/CCTM/src/MECHS/cb6r5_ae7_aq/mech_cb6r5_ae7_aq.def'
ATT = [d for d in os.popen('ls -d /mnt/user-data/attachments/7be0c5a0*').read().split()][0]
TABLES = {'Otavi': 'jgrd15228-sup-0002-ts01.txt', 'Alaska': 'jgrd15228-sup-0003-ts02.txt',
          'Timbavati': 'jgrd15228-sup-0004-ts03.txt'}
CASES = ['Williams', 'Otavi', 'Alaska', 'Timbavati']

MW = json.load(open(os.path.join(HERE, 'cb6r5_molwt.json')))
MW.update({'SULF': 98.08, 'H2SO4': 98.08, 'HCl': 36.46, 'CH4': 16.04, 'H2': 2.016, 'CO2': 44.01,
           'NH3': 17.03, 'HCN': 27.03, 'N2O': 44.01, 'CH3NO3': 77.04, 'RO2_T': 12.0, 'RCO3_T': 12.0,
           'CO': 28.01, 'O1D': 16.0, 'O': 16.0, 'PRPA': 44.1, 'ETHA': 30.07, 'ETOH': 46.07, 'NH3x': 17.03,
           'PAN': 121.05, 'HNO3': 63.01, 'N2O5': 108.0, 'HONO': 47.01, 'NO3': 62.0, 'SO2': 64.06,
           'H2O2': 34.01, 'O3': 48.0, 'NO': 30.01, 'NO2': 46.01})
MW['HCl'] = 36.46
PEROXY = {'RO2': 1, 'C2O3': 2, 'CXO3': 2, 'OPO3': 2}
# placeholders required by hard-wired ASP output code even though CB6r5 has no such species
PLACEHOLDERS = ['CH3NO3', 'RO2_T', 'RCO3_T']
# gases that the aqueous / dissolution decks refer to
INTERFACE = ['HCl', 'NH3', 'HNO3', 'H2SO4']
SOA_ROWS_KEEP = ['TOLSOA', 'XYLSOA', 'ISOSOA', 'TERPSOA']
SOA_ROWS_DROP = ['HC5SOA', 'HC8SOA', 'OLTSOA', 'OLISOA']


def nz(x):
    return abs(x) > 0


# ----------------------------------------------------------------------------------------
def read_case_TP(deck):
    T = P = None
    for ln in open(os.path.join(deck, 'ASPInputDeck.in'), errors='ignore'):
        ln = ln.replace('\r', '')
        m = re.match(r'^\s*([0-9.]+)\s*!!\s*Temperature', ln)
        if m:
            T = float(m.group(1))
        m = re.match(r'^\s*([0-9.]+)\s*!!\s*Pressure', ln)
        if m:
            P = float(m.group(1))
    assert T and P, deck
    return T, P


def read_old_gas(deck):
    """old GasPhaseChems.in -> (header_lines, OrderedDict name -> [evolve, peroxy, mw, init, bg, dep], raw rows)"""
    hdr, rows, raw = [], OrderedDict(), OrderedDict()
    inhdr = True
    for ln in open(os.path.join(deck, 'GasPhaseChems.in'), errors='ignore'):
        ln = ln.replace('\r', '').rstrip('\n')
        if inhdr:
            hdr.append(ln)
            if ln.startswith('!! Put CO first'):
                inhdr = False
            continue
        s = ln.strip()
        if not s or s.startswith('!') or s.startswith('END-OF-FILE'):
            continue
        body = s.split('!')[0]
        f = [x.strip() for x in body.split(';')]
        if len(f) < 7:
            continue
        rows[f[0]] = [float(f[1]), float(f[2]), float(f[3]), float(f[4]), float(f[5]), float(f[6])]
        raw[f[0]] = ln
    return hdr, rows, raw


def read_jgr_table(path):
    out = []
    for i, ln in enumerate(open(path, errors='ignore')):
        if i == 0 or not ln.strip():
            continue
        f = ln.split()
        name = f[0]
        bg = float(f[1])
        init = float(f[3])
        out.append((name, init, bg))
    return out


def case_entries(case, old_rows):
    """-> list of (source name, mapping key, init, bg)"""
    ent = []
    if case == 'Williams':
        for name, v in old_rows.items():
            if name in sm.WILLIAMS_NAMES:
                ent.append((name, sm.WILLIAMS_NAMES[name], v[3], v[4]))
        return ent
    f = [x for x in os.listdir(ATT) if x.endswith(TABLES[case])][0]
    for name, init, bg in read_jgr_table(os.path.join(ATT, f)):
        key = sm.JGR_NAMES.get(name.strip())
        if key is None:
            raise KeyError('%s: no mapping for table species %r' % (case, name))
        ent.append((name.strip(), key, init, bg))
    return ent


# ----------------------------------------------------------------------------------------
def photo_table(case_dir, tr, out_path):
    """derive CB6r5 PhotoRates.in from the case's TUV table."""
    src = open(os.path.join(case_dir, 'PhotoRates.in'), errors='ignore').read().replace('\r', '').split('\n')
    hdr_site = [ln for ln in src if ln.startswith('!!') and ('Site' in ln or 'PBL' in ln or 'lat' in ln or
                                                           'aerosol' in ln or 'Command' in ln)]
    data = []
    for ln in src:
        f = ln.split()
        if len(f) == 53 and re.match(r'^[0-9.]+$', f[0]):
            data.append(f)
    assert len(data) == 96, (case_dir, len(data))
    cols = list(tr.photo.cols.items())
    out = []
    out.append('!! PhotoRates.in  (CB6r5 column set)')
    out.append('!! ASP')
    out.append('!! Derived by derive_cb6r5_photorates (build_cb6r5_decks.py) from the 51-column TUV 5.4 table of the')
    out.append('!! original deck of this case; column k below is a linear combination of the original columns')
    out.append('!! (weights in the list).  Columns beyond the last used one are zero.')
    out.append('!!')
    out.extend(hdr_site)
    out.append('!!')
    out.append('!! Columns: time (UT, hrs), solar zenith angle (deg), then 51 photolysis rate coefficients (s-1)')
    for (lab, scale), idx in cols:
        combo, desc = c2a.PHOTO_COMBOS[lab]
        w = ' + '.join('%g x old%d' % (wt * scale, oc) for oc, wt in combo)
        out.append('!  %2d = %s  [CMAQ %s%s]   = %s' % (idx, desc, lab, '' if scale == 1 else ' x %g' % scale, w))
    for k in range(len(cols) + 1, 52):
        out.append('!  %2d = unused (0)' % k)
    out.append('!Columns: time, sza, photo-reactions')
    for f in data:
        old = [float(x) for x in f[2:]]
        new = [0.0] * 51
        for (lab, scale), idx in cols:
            combo, _ = c2a.PHOTO_COMBOS[lab]
            new[idx - 1] = scale * sum(wt * old[oc - 1] for oc, wt in combo)
        out.append('%10s %10s ' % (f[0], f[1]) + ' '.join('%10.3e' % v for v in new))
    out.append('END-OF-FILE-012343210-ELIF-FO-DNE')
    open(out_path, 'w').write('\n'.join(out) + '\n')


# ----------------------------------------------------------------------------------------
def mech_file(case, case_dir, tr, T, P, out_path):
    old = open(os.path.join(case_dir, 'GasChemicalMechanism.in'), errors='ignore').read().replace('\r', '').split('\n')
    # SOA -> IVOC chains (only the four counters that CB6r5 still produces) + IVOC aging block
    a = next(i for i, l in enumerate(old) if re.match(r'^TOLSOA\s*\+\s*NO\s', l))
    b = next(i for i, l in enumerate(old) if re.match(r'^TERPSOA\s*=>', l))
    chain = [l for l in old[a:b + 1] if l.strip() and not l.startswith('HC') and not l.startswith('OL')]
    aging = old[b + 1:]
    aging = [l.replace('Acetic_Acid', 'AACD') for l in aging]
    while aging and not aging[-1].strip():
        aging.pop()
    L = []
    L.append('!! GasChemicalMechanism.in - Carbon Bond 6 revision 5 (CB6r5), %s case' % case)
    L.append('!! ASP')
    L.append('!! Translated from the CMAQ v5.4+ mechanism definition CB6r5-AE7-AQ (mech_cb6r5_ae7_aq.def) by')
    L.append('!! CB6r5Translation/cmaq2asp.py - see CB6r5Translation/CB6r5_README.md for every choice and warning.')
    L.append('!! Gas-phase core only (the %d R<n> reactions of the .def file, minus R216 ECH4+OH, which is the CH4 tracer' % len(tr.all and [r for r in tr.all if c2a.CORE.match(r.label)]))
    L.append('!! reaction; CH4 + OH is R<n> above).  Chlorine, DMS, mercury, toxics-diagnostic, SA, BL/BP, POA and')
    L.append('!! aerosol heterogeneous groups are NOT included (gas chemistry only, per project scope).')
    L.append('!! Case conditions: T = %.1f K, P = %.1f mb (used only to choose the IUPAC/JPL fall-off forms).' % (T, P))
    L.append('!! Rate-flag conventions are those of the original ASP deck (see the file header of the old decks).')
    L.append('!! REQUIRES common_asp >= f3f1107: R46 uses 2-body rate flag 5 and R24 3-body flag 2, both mis-coded in earlier versions.')
    L.append('!!')
    L.append('!! Format: "!! <label>: <CMAQ statement>" precedes each translated reaction.')
    L.append('!! WARNING lines name rate forms ASP cannot reproduce exactly.')
    L.append('')
    notes_by = {}
    for n in tr.notes:
        notes_by.setdefault(n['label'], []).append(n)
    for lab, text, line, nb, params, comment in tr.lines:
        L.append('!! %s: %s' % (lab, text))
        for n in notes_by.get(lab, []):
            if n['level'] in ('approx', 'unsupported', 'expanded', 'soa'):
                tag = {'approx': 'APPROXIMATION', 'unsupported': 'WARNING', 'expanded': 'NOTE', 'soa': 'NOTE'}[n['level']]
                L.append('!!   %s: %s' % (tag, n['text']))
        L.append(line)
    for lab, ns in notes_by.items():
        if lab not in {x[0] for x in tr.lines}:
            for n in ns:
                L.append('!! NOT TRANSLATED %s: %s  [%s]' % (lab, n['text'], n['cmaq']))
    L.append('')
    L.append('!! ------------------------------------------------------------------------------------')
    L.append('!! SOA precursor counters -> IVOC volatility bins (unchanged from the %s deck of the' % case)
    L.append('!! baseline mechanism; the counters TOLSOA, XYLSOA, ISOSOA, TERPSOA are produced by the')
    L.append('!! CB6r5 aromatic / isoprene / terpene oxidation reactions above).  HC5SOA, HC8SOA, OLTSOA and')
    L.append('!! OLISOA chains of the baseline are not re-attached (see README, open decisions).')
    L.append('!! ------------------------------------------------------------------------------------')
    L.extend(chain)
    L.append('!')
    L.extend(aging)
    L.append('!')
    open(out_path, 'w').write('\n'.join(L) + '\n')


# ----------------------------------------------------------------------------------------
def build_case(case):
    old_dir = os.path.join(DECKS, case)
    new_dir = os.path.join(DECKS, case + '_CB6r5')
    if os.path.isdir(new_dir):
        shutil.rmtree(new_dir)
    shutil.copytree(old_dir, new_dir)
    T, P = read_case_TP(old_dir)
    tr = c2a.Translation(DEF, T, P, soa_counters='asp')
    hdr, old_rows, old_raw = read_old_gas(old_dir)
    ent = case_entries(case, old_rows)
    mapped, rows = sm.apply(ent)

    # carry-overs that are not in the papers' gas tables
    carry = []
    for nm in ('CO2', 'Trace', 'O3'):
        if nm not in mapped and nm in old_rows:
            v = old_rows[nm]
            mapped[nm] = [v[3], v[4]]
            carry.append((nm, v[3], v[4]))
    # species the mechanism needs
    all_sp = OrderedDict()
    for nm in mapped:
        all_sp[nm] = 1
    for nm in tr.species:
        all_sp.setdefault(nm, 1)
    for nm in INTERFACE + PLACEHOLDERS:
        all_sp.setdefault(nm, 1)
    for nm in ('CO', 'CO2', 'Trace', 'NO', 'NO2', 'O3', 'PAN', 'NO3', 'N2O5', 'HONO', 'NH3', 'HNO3'):
        all_sp.setdefault(nm, 1)
    # IVOC chain species of the old deck
    ivoc = [n for n in old_rows if re.match(r'^(IVOC|POA|LEVO|CBIO|CPD)', n) or n in SOA_ROWS_KEEP]
    # ---- GasPhaseChems.in --------------------------------------------------------------
    out = list(hdr)
    out[0] = '!! GasPhaseChems.in - %s case, Carbon Bond 6 revision 5 (CB6r5)' % case
    L = out
    # CO init is the second value of the "2; <CO>!" line, which is kept from the old header
    L.append('CO\t;1\t;0;28;%g;%g;0!' % tuple(mapped.get('CO', [0, 0])))
    first = ['CO']
    species_rows = []
    skip = set(first) | set(ivoc) | {'O2', 'N2', 'OffGas'}
    for nm in all_sp:
        if nm in skip:
            continue
        if nm in ('O1D', 'O'):
            pass
        ini, bg = mapped.get(nm, [0.0, 0.0])
        mw = MW.get(nm) or MW.get(nm.upper())
        if mw is None:
            if nm in ('Trace',):
                mw = old_rows['Trace'][2]
            else:
                raise KeyError('no molecular weight for ' + nm)
        if nm in ('Trace', 'CO2') and nm in old_rows:
            mw = old_rows[nm][2]
        pf = PEROXY.get(nm, 0)
        evolve = 1
        species_rows.append('%s\t;%d\t;%d;%g;%.6g;%.6g;0!' % (nm, evolve, pf, mw, ini, bg))
    L.extend(species_rows)
    L.append('!')
    L.append('!VBS / organic-phase compounds referenced by the organic dissolution decks (unchanged from the %s deck)' % case)
    L.append('!MW from Grieshop et al.; ACP; 2009')
    for nm in ivoc:
        L.append(old_raw[nm].replace('\t', '\t'))
    L.append('!')
    L.append('O2\t;0\t;0;32.;2.095e+08;2.095e+08;0!')
    L.append('N2\t;0\t;0;28.;7.808e+08;7.808e+08;0!')
    L.append('OffGas\t;0\t;0;12.;1000;1000;0!')
    L.append('END-OF-FILE-012343210-ELIF-FO-DNE')
    open(os.path.join(new_dir, 'GasPhaseChems.in'), 'w').write('\n'.join(L) + '\n')

    mech_file(case, old_dir, tr, T, P, os.path.join(new_dir, 'GasChemicalMechanism.in'))
    photo_table(old_dir, tr, os.path.join(new_dir, 'PhotoRates.in'))

    with open(os.path.join(new_dir, 'SPECIES_MAPPING.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['paper_species', 'mapping_key', 'initial_ppb', 'background_ppb', 'CB6r5_representation',
                    'carbon_in_compound', 'carbon_in_CB6r5', 'basis', 'note'])
        for src, key, ini, bg, comp, nC, basis, note in rows:
            rep = ' + '.join(('%g ' % cf if cf != 1 else '') + sp for sp, cf in comp)
            w.writerow([src, key, ini, bg, rep, nC if nC is not None else '', sm.carbon_represented(comp) or '',
                        basis, note])
    return dict(case=case, T=T, P=P, tr=tr, rows=rows, mapped=mapped, carry=carry, nspecies=len(all_sp),
                nphoto=len(tr.photo.cols), ivoc=ivoc)


if __name__ == '__main__':
    res = OrderedDict()
    for cs in (sys.argv[1:] or CASES):
        r = build_case(cs)
        res[cs] = r
        print('%-10s T=%.1f P=%.1f  %d reactions, %d species, %d photolysis columns; carried over: %s' %
              (cs, r['T'], r['P'], len(r['tr'].lines), r['nspecies'], r['nphoto'], r['carry']))
