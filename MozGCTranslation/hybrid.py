"""Generic hybrid-deck builder: new gas mechanism (already translated to ASP lines) + retained baseline lumped chemistry.

Policy (same as the MCM v3.3.1 decks): a baseline reaction is kept iff at least one reactant, after renaming baseline compounds
to the new mechanism's names, is not a species of the new mechanism.  SOA counters, IVOC/POA chains and the aging block are kept.
PhotoRates.in and the remaining deck files are copied from the baseline deck unchanged.
"""
import collections, csv, os, re, shutil, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'CB6r5Translation')); sys.path.insert(0, os.path.join(HERE, '..', 'MCMTranslation'))
import build_cb6r5_decks as b6
from baseline_mech import read_baseline, fmt_side
DECKS, CASES = b6.DECKS, b6.CASES
NEUTRAL = {'M', 'H2O', 'O2', 'N2', 'RO2_T', 'RCO3_T', 'hv'}
SPECIAL_KEEP = ['CO', 'CO2', 'Trace', 'NO', 'NO2', 'O3', 'CH3NO3', 'RO2_T', 'RCO3_T', 'NH3', 'H2SO4', 'HCl', 'HNO3', 'SO2']

def side_names(head):
    out = []
    for side in head.split('=>'):
        for t in side.split('+'):
            t = t.strip()
            m = re.match(r'^[0-9.]+\s+(\S+)$', t)
            out.append(m.group(1) if m else t)
    return out

def build_case(case, cfg, mech_lines, header):
    """cfg: suffix, rename{base->new}, ic_scale{base->factor}, species_info{new->(mw, peroxy_flag)}, tag, mapping_csv"""
    old_dir = os.path.join(DECKS, case); new_dir = os.path.join(DECKS, case + cfg['suffix'])
    if os.path.isdir(new_dir): shutil.rmtree(new_dir)
    shutil.copytree(old_dir, new_dir)
    hdr, old_rows, old_raw = b6.read_old_gas(old_dir)
    base = read_baseline(os.path.join(old_dir, 'GasChemicalMechanism.in'))
    R = cfg['rename']; ren = lambda n: R.get(n, n)
    mech_set = set(cfg['species_set'])
    keep, drop = [], []
    for r in base:
        names = [ren(n) for c, n in r['lhs']]
        (drop if all(n in mech_set or n in NEUTRAL for n in names) else keep).append(r)
    def renamed_line(r):
        return '%s => %s ;%s' % (fmt_side([(c, ren(n)) for c, n in r['lhs']]), fmt_side([(c, ren(n)) for c, n in r['rhs']]), r['tail'])
    keep_lines = [renamed_line(r) for r in keep]
    used = collections.OrderedDict()
    def add(n):
        if n not in NEUTRAL and n != 'PROD': used.setdefault(n, 1)
    for l in mech_lines:
        for n in side_names(l.split(' ; ')[0]): add(n)
    for r in keep:
        for c, n in r['lhs'] + r['rhs']: add(ren(n))
    rows0 = collections.OrderedDict()
    for n, v in old_rows.items():
        rn = ren(n); f = cfg['ic_scale'].get(n, 1.0)
        if rn in rows0:
            rows0[rn][3] += v[3] * f; rows0[rn][4] += v[4] * f
        else:
            rows0[rn] = [v[0], v[1], v[2], v[3] * f, v[4] * f, v[5]]
    for n, v in rows0.items():
        if n in NEUTRAL or n in used: continue
        if v[3] != 0 or v[4] != 0 or n in SPECIAL_KEEP or re.match(r'^(IVOC|POA|LEVO|CBIO|CPD)', n): used.setdefault(n, 1)
    for n in SPECIAL_KEEP: used.setdefault(n, 1)
    info = cfg['species_info']
    rows = collections.OrderedDict()
    for n in used:
        if n == 'CO': continue
        b = rows0.get(n)
        if n in info:
            mw, pf = info[n]
            if b is not None and n in ('H2SO4', 'HCl', 'SO2', 'NO', 'NO2', 'O3', 'CO2', 'HNO3', 'NH3'): mw = b[2]
            rows[n] = (b[0] if b else 1.0, pf, b[2] if (b and n not in cfg.get('force_mw', ())) else mw, b[3] if b else 0.0, b[4] if b else 0.0, b[5] if b else 0.0)
        elif b is not None:
            rows[n] = tuple(b)
        else:
            raise SystemExit('%s: species %s has no row (no baseline row, no species_info)' % (case, n))
    L = list(hdr); L[0] = '!! GasPhaseChems.in - %s case, %s + baseline lumped species' % (case, cfg['tag'])
    c0 = rows0.get('CO', [1, 0, 28, 0, 0, 0])
    L.append('CO\t;1\t;0;28;%g;%g;0!' % (c0[3], c0[4]))
    nperox = 0
    for n, (ev, pf, mw, ini, bg, dep) in rows.items():
        if n in ('CO', 'O2', 'N2', 'OffGas'): continue
        nperox += pf > 0
        L.append('%s\t;%d\t;%d;%g;%.8g;%.8g;%g!' % (n, int(ev), int(pf), mw, ini, bg, dep))
    L.append('!')
    for nm in ('O2', 'N2', 'OffGas'): L.append(old_raw[nm].rstrip())
    L.append('END-OF-FILE-012343210-ELIF-FO-DNE')
    open(os.path.join(new_dir, 'GasPhaseChems.in'), 'w').write('\n'.join(L) + '\n')
    M = [header(len(mech_lines), len(keep), len(drop)), '!! ---------------- 1. %s ----------------' % cfg['tag']]
    M += mech_lines + ['', '!! ---------------- 2. baseline reactions kept (reactant not in %s; compound names renamed) ------' % cfg['tag']]
    M += keep_lines + ['']
    open(os.path.join(new_dir, 'GasChemicalMechanism.in'), 'w').write('\n'.join(M) + '\n')
    with open(os.path.join(new_dir, cfg['mapping_csv']), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['baseline_name', 'model_name', 'ic_scale', 'initial_ppb', 'background_ppb', 'treatment'])
        for n, v in old_rows.items():
            rn = ren(n); f = cfg['ic_scale'].get(n, 1.0)
            t = ('renamed to %s species' % cfg['tag']) if rn != n else (('%s species (same name)' % cfg['tag']) if rn in mech_set else
                 ('lumped species kept, baseline chemistry' if rn in used else 'not used'))
            w.writerow([n, rn, f, v[3] * f, v[4] * f, t])
    return dict(case=case, n_keep=len(keep), n_drop=len(drop), n_species=len(rows), n_perox=nperox, keep=keep, drop=drop)
