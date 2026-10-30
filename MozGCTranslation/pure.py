"""Deck builder for a new gas mechanism WITHOUT any retained baseline lumped gas chemistry.

  GasChemicalMechanism.in : header, the translated mechanism, the SOA-counter -> IVOC chains for the counters that appear in the
                            mechanism, and the IVOC aging block (both copied from the case's baseline deck)
  GasPhaseChems.in        : mechanism species + ASP interface / placeholder species; the paper compounds of the case are lumped onto
                            mechanism species through `lump` (key -> [(species, coefficient)], keys of CB6r5Translation/species_map.py);
                            IVOC / POA rows from the baseline deck
  PhotoRates.in and all other files: copied from the baseline deck.
"""
import collections, csv, os, re, shutil, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'CB6r5Translation'))
import build_cb6r5_decks as b6
import species_map as sm
from mech_species import MW_EXTRA
import build_mech_decks as bm

DECKS, CASES = b6.DECKS, b6.CASES
NEUTRAL = {'M', 'H2O', 'O2', 'N2', 'RO2_T', 'RCO3_T', 'hv'}
INTERFACE = ['HCl', 'NH3', 'HNO3', 'H2SO4']
PLACEHOLDERS = ['CH3NO3', 'RO2_T', 'RCO3_T']
DIRECT = ['CO', 'CO2', 'NO', 'NO2', 'O3', 'H2O2', 'HONO', 'SO2', 'HNO3', 'HCl', 'NH3', 'H2', 'CH4', 'N2O5', 'PAN', 'HCN', 'N2O', 'H2SO4']


def side_names(head):
    out = []
    for side in head.split('=>'):
        for t in side.split('+'):
            t = t.strip()
            m = re.match(r'^[0-9.]+\s+(\S+)$', t)
            out.append(m.group(1) if m else t)
    return out


def build_case(case, cfg, mech_lines, header):
    """cfg: suffix, tag, lump{key->([(sp,coef)], note)}, species_info{sp->(mw, peroxy)}, counters, acetic, mapping_csv"""
    old_dir = os.path.join(DECKS, case); new_dir = os.path.join(DECKS, case + cfg['suffix'])
    if os.path.isdir(new_dir): shutil.rmtree(new_dir)
    shutil.copytree(old_dir, new_dir)
    hdr, old_rows, old_raw = b6.read_old_gas(old_dir)
    lump = dict(cfg['lump'])
    for k in DIRECT: lump.setdefault(k, ([(cfg.get('direct_rename', {}).get(k, k), 1.0)], 'direct'))
    info = cfg['species_info']
    # --- lump the case's compounds onto mechanism species
    ent = b6.case_entries(case, old_rows)
    mapped = collections.OrderedDict(); rows = []
    for src, key, ini, bg in ent:
        comp, note = lump[key]
        for sp, c in comp:
            v = mapped.setdefault(sp, [0.0, 0.0]); v[0] += c * ini; v[1] += c * bg
        rows.append((src, key, ini, bg, comp, note))
    carry = []
    for nm in ('CO2', 'Trace', 'O3'):
        if nm not in mapped and nm in old_rows:
            v = old_rows[nm]; mapped[nm] = [v[3], v[4]]; carry.append(nm)
    used = collections.OrderedDict()
    for l in mech_lines:
        for n in side_names(l.split(' ; ')[0]):
            if n not in NEUTRAL and n != 'PROD' and n not in cfg['counters']: used.setdefault(n, 1)
    for n in mapped: used.setdefault(n, 1)
    for n in INTERFACE + PLACEHOLDERS + ['CO', 'CO2', 'Trace', 'NO', 'NO2', 'O3', 'PAN', 'NO3', 'N2O5', 'HONO']: used.setdefault(n, 1)
    ivoc = [n for n in old_rows if re.match(r'^(IVOC|POA|LEVO|CBIO|CPD)', n) or n in cfg['counters']]
    skip = {'CO', 'O2', 'N2', 'OffGas'} | set(ivoc)
    for n in cfg['counters']:
        assert n in old_rows, n
    # --- GasPhaseChems.in
    L = list(hdr); L[0] = '!! GasPhaseChems.in - %s case, %s' % (case, cfg['tag'])
    L.append('CO\t;1\t;0;28;%g;%g;0!' % tuple(mapped.get('CO', [0, 0])))
    nperox = 0; nsp = 0
    for n in used:
        if n in skip: continue
        ini, bg = mapped.get(n, [0.0, 0.0])
        if n in ('Trace', 'CO2') and n in old_rows: mw, pf = old_rows[n][2], 0
        elif n in info: mw, pf = info[n]
        elif n in MW_EXTRA: mw, pf = MW_EXTRA[n], 0
        else: raise SystemExit('no molar mass for species %s' % n)
        if n in ('H2SO4', 'HCl', 'NH3', 'HNO3', 'SO2', 'NO', 'NO2', 'O3', 'CO2', 'N2O', 'H2', 'CH4', 'HCN', 'PAN', 'N2O5', 'HONO', 'H2O2') and n in MW_EXTRA:
            mw = MW_EXTRA[n]
        nperox += pf > 0; nsp += 1
        L.append('%s\t;1\t;%d;%g;%.6g;%.6g;0!' % (n, pf, mw, ini, bg))
    L.append('!')
    L.append('!VBS / organic-phase compounds referenced by the organic dissolution decks (unchanged from the %s deck)' % case)
    L.append('!MW from Grieshop et al.; ACP; 2009')
    for n in ivoc: L.append(old_raw[n])
    L.append('!')
    for nm in ('O2', 'N2', 'OffGas'): L.append(old_raw[nm].rstrip())
    L.append('END-OF-FILE-012343210-ELIF-FO-DNE')
    open(os.path.join(new_dir, 'GasPhaseChems.in'), 'w').write('\n'.join(L) + '\n')
    # --- GasChemicalMechanism.in
    old_lines = open(os.path.join(old_dir, 'GasChemicalMechanism.in'), errors='ignore').read().replace('\r', '').split('\n')
    chain, aging = bm.old_chain(old_lines, set(cfg['counters']))
    aging = [l.replace('Acetic_Acid', cfg['acetic']) for l in aging]
    M = [header(len(mech_lines)), '!! ---------------- 1. %s ----------------' % cfg['tag']] + mech_lines
    M += ['', '!! ---------------- 2. SOA precursor counters -> IVOC volatility bins (unchanged from the %s baseline deck) ----' % case,
          '!! Counters present in this mechanism: %s' % ', '.join(cfg['counters'])] + chain + ['!'] + aging + ['!']
    open(os.path.join(new_dir, 'GasChemicalMechanism.in'), 'w').write('\n'.join(M) + '\n')
    with open(os.path.join(new_dir, cfg['mapping_csv']), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['paper_species', 'mapping_key', 'initial_ppb', 'background_ppb', '%s_representation' % cfg['tag'], 'basis_note'])
        for src, key, ini, bg, comp, note in rows:
            w.writerow([src, key, ini, bg, ' + '.join(('%.4g ' % cf if abs(cf - 1) > 1e-9 else '') + sp for sp, cf in comp), note])
    return dict(case=case, n_species=nsp, n_perox=nperox, rows=rows, mapped=mapped)
