#!/usr/bin/env python3
"""Build InorganicInputDecks/WilliamsFlats/ and WilliamsFlats_<MECH>/ for the 8 mechanisms.

Source: the SAM-ASP-TUV v1.0.0 deck set (Zenodo 10.5281/zenodo.13385240, njune74/SAM-ASP-TUV
SRC_TUV/ASP/InorganicInputDecks), which was used for the Williams Flats fire (WA, 3 Aug 2019,
FIREX-AQ; SAM case SND1E1S1). Its case-specific files are kept in zenodo_samasp_decks/.

  python3 build_williamsflats_decks.py            build all 9 decks
  python3 build_williamsflats_decks.py --check    only the validation (reproduce Williams_<MECH> from Williams)

Baseline deck (WilliamsFlats/) = the Williams/ baseline deck with these case changes:
  GasPhaseChems.in   SAM deck emission factors (flag 0, g/kg) -> parcel initial mixing ratios with ASP's own
                     flag-0 formula at an initial CO of 10,000 ppb (flag 2, as all ASP case decks); backgrounds
                     from the SAM deck. Volatility bins: the SAM deck's two 4-bin sets (IVOC1-4, IVOC5-8, same
                     C* and MW in each set) go to the baseline 9-bin set by identical C* and MW:
                     IVOC1/5 -> IVOC1, IVOC2/6 -> IVOC3, IVOC3/7 -> IVOC5, IVOC4/8 -> IVOC7.
                     OXYARO and OTHERO (lumped SOA precursors) are added.
  GasChemicalMechanism.in  + OXYARO/OTHERO + OH reactions of the SAM deck (OH regenerated), products remapped.
  AerosolModes.in    the SAM deck's smoke mode ("Williams Flats Mod"), organics remapped as above; the
                     baseline's 40-bin sectional set-up is kept.
  ASPInputDeck.in    T, P, RH: surface values of the SAM SND1 sounding (937.3 hPa); photolysis start
                     19:36 UTC (SND1 day0 215.817); mixing parameter 15 s (SAM deck); rest as Williams.
  PhotoRates.in      TUV table for 47.94 N, 118.62 W, 3 Aug 2019 (PhotoRates_TUV531.in; PROVISIONAL, see there).
  everything else    Williams baseline (EnvAerosolModes = the SAM deck's background, same composition).

Mechanism decks (WilliamsFlats_<MECH>/) = Williams_<MECH>/ with:
  * gas species values from the deck's own species-mapping CSV applied to the WilliamsFlats baseline values;
  * organic aerosol (gas and particle) moved with the same rules the mechanism builders used for Williams
    (cmaq_oa.BIN_MAP for CRACMM, cmaq_ae.BIN_MAP for CB6r5/SAPRC07 (C* >= 1e4 dropped), moz_oa for GEOS-Chem
    and MOZART-T1, unchanged IVOC bins for MCM331);
  * OXYARO/OTHERO added with products mapped by the same rules;
  * PhotoRates.in: each mechanism column rebuilt from the WilliamsFlats TUV table as described in the
    Williams_<MECH> header (copy, sum, or TUV column x CSQY ratio; the ratio is taken from the Williams pair
    as a function of solar zenith angle, held at its smallest-SZA value below 53.9 deg).
The rules are first checked by rebuilding the Williams_<MECH> values from Williams/ (--check).
"""
import csv
import math
import os
import re
import shutil
import sys
from collections import OrderedDict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
DECKS = os.path.join(ROOT, 'InorganicInputDecks')
sys.path[:0] = [os.path.join(ROOT, 'CB6r5Translation'), os.path.join(ROOT, 'MozGCTranslation')]
import cmaq_oa  # noqa: E402
import cmaq_ae  # noqa: E402
import moz_oa   # noqa: E402

SRC = os.path.join(HERE, 'zenodo_samasp_decks')
PHOTO = os.path.join(HERE, 'PhotoRates_TUV531.in')
CASE = 'WilliamsFlats'
T, P, RH = 297.0, 937.3, 0.249          # SND1 sounding, interpolated to the surface (937.31 hPa)
WILLIAMS_TP = (288.4, 880.0)
PHOTO_START = 19.6                     # 19:36 UTC = SND1 day0 215.817
MIX = 15.0                             # mixing parameter of the SAM deck (s)
CO_INIT = 10000.0
IVOC_REMAP = {'IVOC1': 'IVOC1', 'IVOC2': 'IVOC3', 'IVOC3': 'IVOC5', 'IVOC4': 'IVOC7',
              'IVOC5': 'IVOC1', 'IVOC6': 'IVOC3', 'IVOC7': 'IVOC5', 'IVOC8': 'IVOC7'}
KEEP_TEMPLATE = ('Trace', 'OffGas', 'O2', 'N2')
MECHS = ['CB6r5', 'CRACMM2', 'CRACMM3', 'GCHEM', 'MCM331', 'MOZT1', 'SAPRC07TC', 'SAPRC07TIC']
FAMILY = {'CB6r5': 'ae7', 'SAPRC07TC': 'ae7', 'SAPRC07TIC': 'ae7', 'CRACMM2': 'cracmm', 'CRACMM3': 'cracmm',
          'GCHEM': 'gc', 'MOZT1': 'moz', 'MCM331': 'base'}
MAPFILE = {'GCHEM': 'GCHEM_SPECIES_MAPPING.csv', 'MCM331': 'MCM_SPECIES_MAPPING.csv',
           'MOZT1': 'MOZT1_SPECIES_MAPPING.csv'}
ORG_FILES = ['OrgPhaseChems.in', 'OrganicDissolution.in', 'HydrophilicOrgChems.in', 'HydrophilicOrganicDissolution.in']
MODE_RE = re.compile(r'^[0-9.eE+]+\s*;\s*[0-9.eE+-]+\s*;')
ROW_DATA = re.compile(r'^\s*[0-9.]+\s+[0-9.]+\s')


# ------------------------------------------------------------------------------------- file helpers
EOL = {}      # path -> line ending of the file as read (the Williams templates are CRLF); kept on write


def rd(path):
    raw = open(path, 'rb').read().decode('latin-1')
    EOL[os.path.abspath(path)] = '\r\n' if '\r\n' in raw else '\n'
    return raw.replace('\r', '').split('\n')


def wr(path, L):
    eol = EOL.get(os.path.abspath(path), '\n')
    open(path, 'wb').write((eol.join(L).rstrip('\r\n') + eol).encode('latin-1'))


def gas_rows(L):
    """GasPhaseChems lines -> OrderedDict name -> (line index, [fields])  (the flag line is skipped)"""
    out, first = OrderedDict(), True
    for i, ln in enumerate(L):
        s = ln.strip()
        if not s or s.startswith('!') or s.startswith('END-OF-FILE'):
            continue
        if first:
            first = False
            continue
        f = [x.strip() for x in s.split('!')[0].split(';')]
        if len(f) >= 7:
            out[f[0]] = (i, f)
    return out


def num(x):
    return float(x.replace('d', 'e').replace('D', 'e'))


def set_gas(L, rows, name, init, bg):
    i, f = rows[name]
    body, _, tail = L[i].partition('!')
    parts = body.split(';')
    same = lambda a, b: abs(num(a) - b) <= 1e-6 * max(abs(b), 1e-12)   # noqa: E731
    if same(parts[4].strip(), init) and same(parts[5].strip(), bg):
        return
    parts[4] = '%.6g' % init
    parts[5] = '%.6g' % bg
    L[i] = ';'.join(parts) + '!' + tail


def insert_before_o2(L, new_lines):
    i = next(k for k, l in enumerate(L) if re.match(r'^O2\s*;', l.strip()))
    L[i:i] = new_lines


def gas_line(name, mw, init, bg, perox=0):
    return '%s\t;1\t;%d;%g;%.6g;%.6g;0.!' % (name, perox, mw, init, bg)


# ------------------------------------------------------------------------------------- SAM deck -> baseline values
def sam_values():
    """-> (OrderedDict name -> [mw, init, bg] in the 9-bin naming, notes)"""
    L = rd(os.path.join(SRC, 'GasPhaseChems.in'))
    rows = gas_rows(L)
    flag = [l for l in L if l.strip() and not l.strip().startswith('!')][0]
    assert flag.split(';')[0].strip() == '0', flag
    ef = {n: (num(f[3]), num(f[4]), num(f[5])) for n, (i, f) in rows.items()}
    mw_co, ef_co, bg_co = ef['CO']
    out = OrderedDict()
    for n, (mw, e, bg) in ef.items():
        init = CO_INIT if n == 'CO' else (e / mw) / (ef_co / mw_co) * (CO_INIT - bg_co) + bg
        if re.match(r'^IVOC[1-9]$', n):
            if n == 'IVOC9':
                assert e == 0
                continue
            m = IVOC_REMAP[n]
            prev = out.get(m, [mw, 0.0, 0.0])
            assert abs(prev[0] - mw) < 1e-6, (n, m)
            out[m] = [mw, prev[1] + (e / mw) / (ef_co / mw_co) * (CO_INIT - bg_co), 0.0]
            continue
        out[n] = [mw, init, bg]
    return out


def sam_mode():
    """first mode of the SAM deck's AerosolModes.in -> (header, name, [(species, value)]) with remapped IVOC names"""
    L = [l for l in rd(os.path.join(SRC, 'AerosolModes.in')) if l.strip() and not l.strip().startswith('!')]
    i = next(k for k, l in enumerate(L) if MODE_RE.match(l.strip()))
    hdr = ' ; '.join(x.strip() for x in L[i].split('!')[0].split(';') if x.strip()) + ';'
    name = L[i + 1].split('!')[0].strip().rstrip(';')
    sp = OrderedDict()
    for l in L[i + 2:]:
        s = l.split('!')[0].strip()
        if MODE_RE.match(s) or ';' not in s:
            break
        n, v = [x.strip() for x in s.split(';')[:2]]
        n = IVOC_REMAP.get(n, n)
        sp[n] = sp.get(n, 0.0) + float(v)
    return hdr, name, list(sp.items())


def sam_oxy_reactions():
    out = []
    for l in rd(os.path.join(SRC, 'GasChemicalMechanism.in')):
        if re.match(r'^(OXYARO|OTHERO)\s*\+\s*OH', l):
            lhs, rest = l.split('=>')
            prod, _, rate = rest.partition(';')
            terms = []
            for t in prod.split('+'):
                t = t.strip()
                m = re.match(r'^([0-9.eE+-]+)\s+(\S+)$', t)
                c, s = (float(m.group(1)), m.group(2)) if m else (1.0, t)
                terms.append((c, s))
            out.append((lhs.strip(), terms, rate.strip()))
    assert len(out) == 2, out
    return out


# ------------------------------------------------------------------------------------- aerosol modes
def replace_first_mode(path, hdr, name, species):
    L = rd(path)
    i = next(k for k, l in enumerate(L) if MODE_RE.match(l.strip()) and not l.strip().startswith('!'))
    j = i + 1
    while j < len(L) and L[j].strip() and not L[j].strip().startswith('!') and not (j > i + 1 and MODE_RE.match(L[j].strip())):
        j += 1
    L[i:j] = [hdr, name + ';'] + ['%s; %.5g;' % (n, v) for n, v in species if v > 0]
    wr(path, L)


def first_mode_species(path):
    L = [l for l in rd(path) if l.strip() and not l.strip().startswith('!')]
    i = next(k for k, l in enumerate(L) if MODE_RE.match(l.strip()))
    out = OrderedDict()
    for l in L[i + 2:]:
        s = l.split('!')[0].strip()
        if MODE_RE.match(s) or ';' not in s:
            break
        n, v = [x.strip() for x in s.split(';')[:2]]
        out[n] = float(v)
    return out


# ------------------------------------------------------------------------------------- organic-aerosol rules
def ivoc_k(n):
    m = re.match(r'^IVOC([1-9])$', n)
    return int(m.group(1)) if m else None


def oa_gas(fam, ivoc, mw_of, T_, P_):
    """gas IVOC bins {k: (ppb, mw)} -> {species: ppb} (GC: NAP is an addition to the mapped NAP)"""
    out = {}
    if fam == 'base':
        return {'IVOC%d' % k: v for k, (v, mw) in ivoc.items()}
    if fam in ('cracmm', 'ae7'):
        bm = cmaq_oa.BIN_MAP if fam == 'cracmm' else cmaq_ae.BIN_MAP
        for k, (ppb, mw) in ivoc.items():
            if k in bm and ppb > 0:
                s = bm[k]
                out[s] = out.get(s, 0.0) + ppb * mw / mw_of[s]
        return out
    ug_lo = sum(cmaq_oa._ppb_per_ugm3(mw, T_, P_) ** -1 * ppb for k, (ppb, mw) in ivoc.items() if k <= 6)
    if fam == 'gc':
        f = cmaq_oa._ppb_per_ugm3(moz_oa.GC_MW_SOA, T_, P_)
        out['SVPOA1'] = ug_lo * moz_oa.GC_POA_SPLIT[0] * f
        out['SVPOA2'] = ug_lo * moz_oa.GC_POA_SPLIT[1] * f
        out['NAP+'] = sum(ppb * mw / moz_oa.GC_NAP_MW for k, (ppb, mw) in ivoc.items() if k >= 7)
        return out
    if fam == 'moz':
        ug_hi = sum(cmaq_oa._ppb_per_ugm3(mw, T_, P_) ** -1 * ppb for k, (ppb, mw) in ivoc.items() if k >= 7)
        out['SVOC'] = ug_lo * cmaq_oa._ppb_per_ugm3(moz_oa.ASSUMED_SVOC_MW, T_, P_)
        out['IVOC'] = ug_hi * cmaq_oa._ppb_per_ugm3(moz_oa.ASSUMED_IVOC_MW, T_, P_)
        return out
    raise ValueError(fam)


def oa_particle(fam, species):
    """first-mode species list -> new list (organics replaced)"""
    keep = [(n, v) for n, v in species if ivoc_k(n) is None]
    org = [(ivoc_k(n), v) for n, v in species if ivoc_k(n) is not None]
    tot = sum(v for k, v in org)
    if fam == 'base':
        return keep + [('IVOC%d' % k, v) for k, v in org]
    if fam in ('cracmm', 'ae7'):
        bm = cmaq_oa.BIN_MAP if fam == 'cracmm' else cmaq_ae.BIN_MAP
        new = OrderedDict()
        for k, v in org:
            if k in bm:
                new[bm[k]] = new.get(bm[k], 0.0) + v
        return keep + list(new.items())
    if fam == 'gc':
        return keep + [('SVPOA1', tot * moz_oa.GC_POA_SPLIT[0]), ('SVPOA2', tot * moz_oa.GC_POA_SPLIT[1])]
    if fam == 'moz':
        return keep + [('POM', tot)]
    raise ValueError(fam)


def oa_products(fam, terms, mw_of):
    """reaction product terms with IVOCk -> mechanism species (mass-conserving molar stoichiometry)"""
    out = OrderedDict()

    def add(s, c):
        out[s] = out.get(s, 0.0) + c
    for c, s in terms:
        k = ivoc_k(s)
        if k is None:
            add(s, c)
            continue
        mw = mw_of['IVOC%d' % k]
        if fam == 'base':
            add(s, c)
        elif fam in ('cracmm', 'ae7'):
            bm = cmaq_oa.BIN_MAP if fam == 'cracmm' else cmaq_ae.BIN_MAP
            if k in bm:
                add(bm[k], c * mw / mw_of[bm[k]])
        elif fam == 'gc':
            if k <= 6:
                add('SVPOA1', c * mw / moz_oa.GC_MW_SOA * moz_oa.GC_POA_SPLIT[0])
                add('SVPOA2', c * mw / moz_oa.GC_MW_SOA * moz_oa.GC_POA_SPLIT[1])
            else:
                add('NAP', c * mw / moz_oa.GC_NAP_MW)
        elif fam == 'moz':
            add('SVOC' if k <= 6 else 'IVOC', c * mw / (moz_oa.ASSUMED_SVOC_MW if k <= 6 else moz_oa.ASSUMED_IVOC_MW))
    return out


# ------------------------------------------------------------------------------------- species-mapping CSV
def read_map(mech):
    p = os.path.join(DECKS, 'Williams_%s' % mech, MAPFILE.get(mech, 'SPECIES_MAPPING.csv'))
    rows = list(csv.reader(open(p)))
    out = []
    for r in rows[1:]:
        src = r[0].strip()
        if mech == 'MCM331' and r[4].strip() == 'not used':     # baseline species with no MCM counterpart
            out.append((src, []))
            continue
        rep = r[1].strip() if mech == 'MCM331' else r[4].strip()
        terms = []
        for t in rep.split(' + '):
            m = re.match(r'^\s*([0-9.eE+-]+)\s+(\S+)\s*$', t)
            terms.append((float(m.group(1)), m.group(2)) if m else (1.0, t.strip()))
        out.append((src, terms))
    return rows, out


def apply_map(mapping, base):
    """base: name -> (init, bg) -> mech species -> [init, bg]"""
    out = OrderedDict()
    for src, terms in mapping:
        if src not in base or src in KEEP_TEMPLATE:
            continue
        ini, bg = base[src]
        if not terms and (ini > 0 or bg > 0):
            print('WARNING: %s (%.4g ppb) has no counterpart in this mechanism and is dropped' % (src, ini))
        for c, s in terms:
            v = out.setdefault(s, [0.0, 0.0])
            v[0] += c * ini
            v[1] += c * bg
    return out


# ------------------------------------------------------------------------------------- photolysis
def photo_rows(path):
    return [[float(x) for x in l.split()] for l in rd(path) if ROW_DATA.match(l) and len(l.split()) >= 10]


def photo_defs(path):
    """'!  j = ...' header lines -> {j: ('zero',) | ('lin', [(w, col)]) | ('ratio', col)}"""
    d = {}
    for l in rd(path):
        m = re.match(r'^!\s+(\d+) = (.*)$', l)
        if not m:
            continue
        j, txt = int(m.group(1)), m.group(2)
        if 'unused' in txt:
            d[j] = ('zero',)
        elif 'x old' in txt:
            d[j] = ('lin', [(float(a), int(b)) for a, b in re.findall(r'([0-9.eE+-]+) x old(\d+)', txt)])
        elif 'CSQY ratio' in txt:
            d[j] = ('ratio', int(re.search(r'TUV col (\d+)\)', txt).group(1)))
        elif '<- TUV col' in txt:
            cols = [int(x) for x in re.search(r'<- TUV col ([0-9 +]+)', txt).group(1).split('+')]
            sc = re.search(r'CMAQ \S+ x ([0-9.eE+-]+)\s*<-', txt)          # e.g. "CMAQ ACET_06 x 0.5  <- TUV col 26"
            f = float(sc.group(1)) if sc else 1.0
            d[j] = ('lin', [(f, c) for c in cols])
        else:
            raise ValueError(l)
    return d


def build_photo(mech_template, base_williams, base_new, out_path, site_header):
    defs = photo_defs(mech_template)
    Wm, Wb, Nb = photo_rows(mech_template), photo_rows(base_williams), photo_rows(base_new)
    assert len(Wm) == len(Wb) == len(Nb) == 96
    ncol = len(Wm[0]) - 2
    new = [[r[0], r[1]] + [0.0] * ncol for r in Nb]
    for j in range(1, ncol + 1):
        kind = defs.get(j, ('zero',))
        if kind[0] == 'zero':
            continue
        if kind[0] == 'lin':
            for i in range(96):
                new[i][j + 1] = sum(w * Nb[i][c + 1] for w, c in kind[1])
            # check against the Williams pair
            for i in range(96):
                ref = sum(w * Wb[i][c + 1] for w, c in kind[1])
                assert abs(ref - Wm[i][j + 1]) <= 2e-3 * max(abs(ref), 1e-30) + 1e-33, (mech_template, j, i, ref, Wm[i][j + 1])
            continue
        c = kind[1]
        pts = sorted((Wb[i][1], Wm[i][j + 1] / Wb[i][c + 1]) for i in range(96) if Wb[i][c + 1] > 0 and Wb[i][1] < 90)
        sz, fr = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
        for i in range(96):
            new[i][j + 1] = Nb[i][c + 1] * float(np.interp(Nb[i][1], sz, fr))
    L = rd(mech_template)
    out, in_site = [], False
    for l in L:
        if ROW_DATA.match(l) and len(l.split()) >= 10:
            continue
        if l.startswith('END-OF-FILE'):
            continue
        if l.startswith('!! Command:') or l.startswith('!! Site:'):
            if not in_site:
                out.extend(site_header)
                in_site = True
            continue
        if in_site and l.startswith('!!   '):
            continue
        out.append(l)
    k = next(i for i, l in enumerate(out) if l.startswith('!! Columns:'))
    rows = ['%10.4f %10.3f ' % (r[0], r[1]) + ' '.join('%10.3e' % v for v in r[2:]) for r in new]
    out = out[:k + 1 + sum(1 for l in out[k + 1:] if l.startswith('!'))] + rows + ['END-OF-FILE-012343210-ELIF-FO-DNE---']
    # keep the column list (comment lines after "Columns:") before the data
    wr(out_path, out)


def site_header_lines():
    L = rd(PHOTO)
    keep = [l for l in L if l.startswith('!! Command:')]
    i = next(k for k, l in enumerate(L) if 'PROVISIONAL' in l)
    while L[i].strip() != '!!':
        keep.append(L[i])
        i += 1
    keep += [l for l in L if l.startswith('!! Site:') or l.startswith('!!   ')]
    return keep


# ------------------------------------------------------------------------------------- baseline deck
def build_baseline(vals):
    dst = os.path.join(DECKS, CASE)
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(os.path.join(DECKS, 'Williams'), dst)
    # gas
    p = os.path.join(dst, 'GasPhaseChems.in')
    L = rd(p)
    rows = gas_rows(L)
    kept = []
    for n in rows:
        if n in KEEP_TEMPLATE:
            continue
        if ivoc_k(n) is not None:
            v = vals.get(n, [0, 0.0, 0.0])
            set_gas(L, rows, n, v[1], 0.0)
        elif n in vals:
            set_gas(L, rows, n, vals[n][1], vals[n][2])
        else:
            kept.append(n)
    insert_before_o2(L, ['!', '! lumped SOA precursors of the SAM-ASP-TUV Williams Flats deck (OXYARO, OTHERO + OH -> IVOC, OH regenerated)',
                         gas_line('OXYARO', vals['OXYARO'][0], vals['OXYARO'][1], vals['OXYARO'][2]),
                         gas_line('OTHERO', vals['OTHERO'][0], vals['OTHERO'][1], vals['OTHERO'][2])])
    for i, l in enumerate(L):
        if l.startswith('!! NOTE: Values below are for'):
            L[i] = ('!! NOTE: Values below are for the Williams Flats Fire (WA, 3 Aug 2019, FIREX-AQ): emission factors of the '
                    'SAM-ASP-TUV v1.0.0 deck (Zenodo 10.5281/zenodo.13385240) at an initial CO of 10,000 ppb; see '
                    'WilliamsFlatsTranslation/')
    wr(p, L)
    # mechanism
    p = os.path.join(dst, 'GasChemicalMechanism.in')
    L = rd(p)
    L += ['!', '!! Williams Flats (SAM-ASP-TUV deck): lumped SOA precursors -> IVOC (products moved to the 9-bin set by C*)']
    for lhs, terms, rate in sam_oxy_reactions():
        L.append('%s => %s ; %s' % (lhs, ' + '.join('%.4g %s' % (c, IVOC_REMAP.get(s, s)) if s != 'OH' else 'OH'
                                                    for c, s in terms), rate))
    wr(p, L)
    # aerosol
    hdr, name, sp = sam_mode()
    replace_first_mode(os.path.join(dst, 'AerosolModes.in'), hdr, name, sp)
    # main deck
    p = os.path.join(dst, 'ASPInputDeck.in')
    L = rd(p)
    sub = [(r'!!\s*Temperature', '%.1f' % T), (r'!!\s*Pressure', '%.1f' % P),
           (r'!!\s*Initial Relative Humidity', '%.3f' % RH), (r'!!\s*Background Relative Humidity', '%.3f' % RH),
           (r'!!\s*Starting Time for Photolysis', '%.4f' % PHOTO_START), (r'!!\s*Mixing parameter', '%.1f' % MIX)]
    for pat, v in sub:
        i = next(k for k, l in enumerate(L) if re.search(pat, l) and not l.strip().startswith('!'))
        L[i] = re.sub(r'^\s*\S+', v, L[i], count=1)
    wr(p, L)
    shutil.copy(PHOTO, os.path.join(dst, 'PhotoRates.in'))
    return kept


# ------------------------------------------------------------------------------------- mechanism decks
def template_values(mech):
    L = rd(os.path.join(DECKS, 'Williams_%s' % mech, 'GasPhaseChems.in'))
    return {n: (num(f[3]), num(f[4]), num(f[5])) for n, (i, f) in gas_rows(L).items()}


def base_values(path):
    L = rd(path)
    return {n: (num(f[3]), num(f[4]), num(f[5])) for n, (i, f) in gas_rows(L).items()}


def mech_gas(mech, base, T_, P_):
    """baseline values {name: (mw, init, bg)} -> {mech species: (init, bg)} (mapping + organic aerosol)"""
    fam = FAMILY[mech]
    rows, mapping = read_map(mech)
    tv = template_values(mech)
    mw_of = {n: v[0] for n, v in tv.items()}
    mw_of.update({n: v[0] for n, v in base.items() if ivoc_k(n)})
    if fam == 'ae7':                       # AERO7 species a SAPRC07 deck may lack (added from CB6r5 in build_mech)
        cb = template_values('CB6r5')
        for s in cmaq_ae.BIN_MAP.values():
            mw_of.setdefault(s, cb[s][0])
    out = apply_map(mapping, {n: (v[1], v[2]) for n, v in base.items()})
    ivoc = {ivoc_k(n): (v[1], v[0]) for n, v in base.items() if ivoc_k(n) and v[1] > 0}
    oa = oa_gas(fam, ivoc, mw_of, T_, P_)
    for s, v in oa.items():
        if s == 'NAP+':
            out.setdefault('NAP', [0.0, 0.0])[0] += v
        else:
            out[s] = [v, 0.0]
    return out, mw_of


def check(mech):
    """rebuild Williams_<mech> gas values and first-mode particles from Williams/"""
    base = base_values(os.path.join(DECKS, 'Williams', 'GasPhaseChems.in'))
    got, _ = mech_gas(mech, base, *WILLIAMS_TP)
    tv = template_values(mech)
    bad = []
    for s, (ini, bg) in got.items():
        if s not in tv:
            bad.append('%s: not in template' % s)
            continue
        for a, b, what in ((ini, tv[s][1], 'init'), (bg, tv[s][2], 'bg')):
            if abs(a - b) > 2e-4 * max(abs(a), abs(b)) + 1e-6:
                bad.append('%s %s: rebuilt %.6g, template %.6g' % (s, what, a, b))
    unexplained = [s for s, v in tv.items() if s not in got and v[1] > 1e-6 and abs(v[1] - v[2]) > 1e-6
                   and s not in KEEP_TEMPLATE]
    sp = list(first_mode_species(os.path.join(DECKS, 'Williams', 'AerosolModes.in')).items())
    want = first_mode_species(os.path.join(DECKS, 'Williams_%s' % mech, 'AerosolModes.in'))
    have = OrderedDict(oa_particle(FAMILY[mech], sp))
    for n in set(want) | set(have):
        a, b = have.get(n, 0.0), want.get(n, 0.0)
        if abs(a - b) > 2e-3 * max(a, b) + 1e-6:
            bad.append('AerosolModes %s: rebuilt %.5g, template %.5g' % (n, a, b))
    return bad, unexplained


def copy_rows(species, src_mech, dst):
    """add a species' rows from another mechanism deck (GasPhaseChems + the four organic-property files)"""
    src = os.path.join(DECKS, 'Williams_%s' % src_mech)
    gL = rd(os.path.join(src, 'GasPhaseChems.in'))
    i, f = gas_rows(gL)[species]
    p = os.path.join(dst, 'GasPhaseChems.in')
    L = rd(p)
    insert_before_o2(L, ['! %s added from Williams_%s (needed for the Williams Flats organic aerosol)' % (species, src_mech), gL[i]])
    wr(p, L)
    for fn in ORG_FILES:
        sl = [l for l in rd(os.path.join(src, fn)) if re.match(r'^%s\s*(;|<=>)' % species, l)]
        assert sl, (species, fn)
        p = os.path.join(dst, fn)
        L = rd(p)
        L += sl
        wr(p, L)


def build_mech(mech, base_dir, notes):
    fam = FAMILY[mech]
    dst = os.path.join(DECKS, '%s_%s' % (CASE, mech))
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(os.path.join(DECKS, 'Williams_%s' % mech), dst)
    base = base_values(os.path.join(base_dir, 'GasPhaseChems.in'))
    got, mw_of = mech_gas(mech, base, T, P)
    # species the Williams deck did not need
    tv = template_values(mech)
    for s in list(got):
        if s not in tv:
            if fam == 'ae7' and s in cmaq_ae.BIN_MAP.values():
                copy_rows(s, 'CB6r5', dst)
                mw_of[s] = template_values('CB6r5')[s][0]
                notes.append('%s: %s added from Williams_CB6r5 (AERO7 species; partitioning only, as the other POA species of this deck)' % (mech, s))
            else:
                raise KeyError('%s: %s not in the Williams_%s deck' % (mech, s, mech))
    p = os.path.join(dst, 'GasPhaseChems.in')
    L = rd(p)
    rows = gas_rows(L)
    for s, (ini, bg) in got.items():
        set_gas(L, rows, s, ini, bg)
    # OXYARO / OTHERO
    insert_before_o2(L, ['!', '! lumped SOA precursors of the SAM-ASP-TUV Williams Flats deck (OXYARO, OTHERO + OH, OH regenerated)'] +
                     [gas_line(n, base[n][0], base[n][1], base[n][2]) for n in ('OXYARO', 'OTHERO')])
    for i, l in enumerate(L):
        if 'Williams' in l and l.startswith('!!') and 'Values below' in l:
            L[i] = '!! NOTE: Values below are for the Williams Flats Fire (WA, 3 Aug 2019); see WilliamsFlatsTranslation/'
    wr(p, L)
    p = os.path.join(dst, 'GasChemicalMechanism.in')
    L = rd(p)
    L += ['!', '!! Williams Flats (SAM-ASP-TUV deck): lumped SOA precursors, products mapped to this mechanism\'s organic aerosol species']
    for lhs, terms, rate in sam_oxy_reactions():
        t9 = [(c, IVOC_REMAP.get(s, s)) for c, s in terms]
        prod = oa_products(fam, t9, mw_of)
        L.append('%s => %s ; %s' % (lhs, ' + '.join(('%.4g %s' % (c, s)) if s != 'OH' else 'OH' for s, c in prod.items()), rate))
    wr(p, L)
    # aerosol: organics that map to a species with no particle phase in this deck (CRACMM VROCP4OXY2, C* = 1e4) are dropped
    hdr, name, sp = sam_mode()
    org_names = {l.split(';')[0].strip() for l in rd(os.path.join(dst, 'OrgPhaseChems.in'))
                 if l.strip() and not l.strip().startswith('!') and ';' in l}
    part = oa_particle(fam, sp)
    keep_part = [(n, v) for n, v in part if n in dict(sp) and ivoc_k(n) is None or n in org_names]
    for n, v in part:
        if (n, v) not in keep_part:
            notes.append('%s: particle %s (mass fraction %.4g of the smoke mode) dropped: no particle phase in this deck '
                         '(C* >= 1e4); it is in the gas phase' % (mech, n, v))
    replace_first_mode(os.path.join(dst, 'AerosolModes.in'), hdr, name, keep_part)
    # main deck, photolysis
    shutil.copy(os.path.join(base_dir, 'ASPInputDeck.in'), os.path.join(dst, 'ASPInputDeck.in'))
    tp = os.path.join(DECKS, 'Williams_%s' % mech, 'PhotoRates.in')
    if open(tp).read() == open(os.path.join(DECKS, 'Williams', 'PhotoRates.in')).read():
        shutil.copy(os.path.join(base_dir, 'PhotoRates.in'), os.path.join(dst, 'PhotoRates.in'))
    else:
        build_photo(tp, os.path.join(DECKS, 'Williams', 'PhotoRates.in'), os.path.join(base_dir, 'PhotoRates.in'),
                    os.path.join(dst, 'PhotoRates.in'), site_header_lines())
    # mapping table with this case's values
    rows, mapping = read_map(mech)
    fn = MAPFILE.get(mech, 'SPECIES_MAPPING.csv')
    with open(os.path.join(dst, fn), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(rows[0])
        for r in rows[1:]:
            if r[0] in base:
                r = list(r)
                r[2], r[3] = '%.6g' % base[r[0]][1], '%.6g' % base[r[0]][2]
            w.writerow(r)
    if os.path.exists(os.path.join(dst, 'CMAQ_OA_MAPPING.txt')):
        os.remove(os.path.join(dst, 'CMAQ_OA_MAPPING.txt'))
    oa = ['Organic-aerosol mapping for Williams Flats (%s, rules of the %s builder)' % (mech, fam)]
    oa += ['gas %s %.4g ppb' % (s, v[0]) for s, v in got.items() if s in set(cmaq_oa.BIN_MAP.values()) | set(cmaq_ae.BIN_MAP.values())
           | {'SVPOA1', 'SVPOA2', 'SVOC', 'IVOC', 'NAP'} or ivoc_k(s)]
    oa += ['particle %s %.4g (mass fraction weight)' % (n, v) for n, v in oa_particle(fam, sp) if n not in dict(sp) or ivoc_k(n)]
    wr(os.path.join(dst, 'OA_MAPPING_WilliamsFlats.txt'), oa)


# ------------------------------------------------------------------------------------- main
def main():
    ok = True
    for m in MECHS:
        bad, unexplained = check(m)
        print('check %-11s %s%s' % (m, 'OK' if not bad else '%d differences' % len(bad),
                                    ('; template-only nonzero species: ' + ' '.join(unexplained)) if unexplained else ''))
        for b in bad[:15]:
            print('    ', b)
        ok &= not bad
    if '--check' in sys.argv:
        return 0 if ok else 1
    if not ok:
        print('validation failed; decks not written')
        return 1
    vals = sam_values()
    kept = build_baseline(vals)
    print('%s: written (species kept from Williams: %s)' % (CASE, ' '.join(kept) or 'none'))
    notes = []
    for m in MECHS:
        build_mech(m, os.path.join(DECKS, CASE), notes)
        print('%s_%s: written' % (CASE, m))
    for n in notes:
        print('  ' + n)
    return 0


if __name__ == '__main__':
    sys.exit(main())
