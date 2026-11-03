#!/usr/bin/env python3
"""Build InorganicInputDecks/<Case>_<MECH>/ for CRACMM2, CRACMM3, SAPRC07TC, SAPRC07TIC x {Williams, Otavi, Alaska, Timbavati}.

  * GasChemicalMechanism.in - the translated gas-phase core of the CMAQ mechanism (mech2asp.py), then the SOA-counter -> IVOC
                               chains and the IVOC aging block of the case's baseline deck
  * GasPhaseChems.in        - mechanism species; initial / background values from the papers' compound tables mapped with mech_species.py
  * PhotoRates.in           - 51-column table built from the case's TUV table (mech_photo.py)
  * SPECIES_MAPPING.csv     - every paper species and the model species that represent it
  * everything else         - copied from the baseline deck
Also sets "Run Length in minutes" to 300 in every deck (baseline, CB6r5 and new).
"""
import csv
import glob
import os
import re
import shutil
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_cb6r5_decks as b6                    # noqa: E402  (case tables, old-deck readers)
import mech2asp as m2a                            # noqa: E402
import mech_photo as mph                          # noqa: E402
import mech_species as ms                         # noqa: E402
import species_map as sm                          # noqa: E402
import cmaq_oa                                    # noqa: E402
import cmaq_ae                                    # noqa: E402

DECKS = b6.DECKS
CASES = b6.CASES
INTERFACE = ['HCl', 'NH3', 'HNO3', 'H2SO4']
PLACEHOLDERS = ['CH3NO3', 'RO2_T', 'RCO3_T']
RUN_MIN = 300


def set_runlength(deck_dir, minutes=RUN_MIN):
    p = os.path.join(deck_dir, 'ASPInputDeck.in')
    raw = open(p, 'rb').read().decode('latin-1')
    new, n = re.subn(r'(?m)^(\s*)[0-9.]+(\s*!!\s*Run Length in minutes)', r'\g<1>%d\2' % minutes, raw)
    assert n == 1, p
    open(p, 'wb').write(new.encode('latin-1'))


def old_chain(old_lines, keep):
    """SOA-counter -> IVOC lines of the baseline deck for the counters in `keep`, and the IVOC aging block"""
    a = next(i for i, l in enumerate(old_lines) if re.match(r'^HC5SOA\s*\+\s*NO\s', l))
    b = next(i for i, l in enumerate(old_lines) if re.match(r'^TERPSOA\s*=>', l))
    chain = []
    for l in old_lines[a:b + 1]:
        mm = re.match(r'^(\w+)', l)
        if mm and mm.group(1) in keep:
            chain.append(l)
    aging = old_lines[b + 1:]
    while aging and not aging[-1].strip():
        aging.pop()
    return chain, aging


def mech_file(cfg, case, tr, T, P, old_lines, out_path):
    cmaq = cfg.get('oa') in ('cmaq', 'ae')
    chain, aging = ([], []) if cmaq else old_chain(old_lines, set(cfg['ivoc_counters']))
    aging = [l.replace('Acetic_Acid', cfg['acetic']) for l in aging]
    nm = cfg['label']
    L = []
    L.append('!! GasChemicalMechanism.in - %s, %s case' % (nm, case))
    L.append('!! ASP')
    L.append('!! Translated from the CMAQ mechanism definition %s by CB6r5Translation/mech2asp.py' % os.path.basename(cfg['def_']))
    L.append('!! (see CB6r5Translation/NewMechanisms_README.md for every choice and warning).')
    if cfg.get('oa') == 'ae':
        L.append('!! Gas-phase core plus CMAQ AERO organic-aerosol chemistry (AE group SOA yields -> SVAVB / SVBNZ / SVXYL / SVTOL ..., SOAALK, pcSOA).')
        L.append('!! SOA counters (ISOPRXN, TRPRXN, SESQRXN, ALKRXN, XYLNRXN ...) are replaced by their CMAQ mass-yield products; MTNO3 / ISOPNN are')
        L.append('!! partitioning vapors (see cmaq_ae.py, OrgPhaseChems.in, OrganicDissolution.in).  NOT translated: chlorine, heterogeneous (IEPOX,')
        L.append('!! glyoxal) uptake, particle-phase hydrolysis / oligomerisation (OLIG_*, HYD_*), toxics.')
    elif cmaq:
        L.append('!! Gas-phase core plus CMAQ\'s volatility-resolved organic-aerosol chemistry (ROCALK, ROCOXY, ROCARO, TRP groups; V* vapor')
        L.append('!! species are ASP gas / organic-phase species, ASOATJ -> ASOAT).  NOT translated: chlorine, isoprene heterogeneous (IEPOX/IPX/')
        L.append('!! INALD) uptake, particle-phase organic reactions (ATRPNJ/AHONITJ hydrolysis, AHOMJ photolysis) and toxics groups.  The baseline')
        L.append('!! SOA counters -> IVOC chains and the IVOC aging block are NOT used here (see OrgPhaseChems.in, OrganicDissolution.in).')
    else:
        L.append('!! Gas-phase core only; CMAQ chlorine, heterogeneous, aerosol/SVOC, POA-aging, oligomerisation and toxics groups are not')
    L.append('!! translated.  %d CMAQ statements -> %d ASP reactions.' % (len(tr.core), len(tr.lines)))
    L.append('!! REQUIRES common_asp >= c506969 (1-body rate flags 18-21 for the CMAQ heterogeneous chemistry and OrgPhaseReactions.in 53d416e; 1-body rate flag 13, 3-body flag 4 and the earlier flag fixes 40c1760 are older)')
    L.append('!! (2-body flag 5 = 508ba8d, 3-body flag 2 = f3f1107).')
    L.append('!! Rate-flag conventions: see the header of the baseline decks; new flags:')
    L.append('!!   1-body flag 13:  k = kf / Keq,  kf = fall-off with explicit Troe N')
    L.append('!!       13; A0; n0; E0; Ai; ni; Ei; F; N; KeqA; KeqE;   k0 = A0 (T/300)^-n0 exp(E0/T) [M],  ki likewise,  Keq = KeqA exp(KeqE/T)')
    L.append('!!   3-body flag 4:   same fall-off form, no Keq:  4; A0; n0; E0; Ai; ni; Ei; F; N;')
    L.append('!!   2-body flag 11:  k = A exp(B/T) (C/T + D P[Pa] + F)/100  (CRACMM NO + HO2 -> HNO3, exact)')
    L.append('!! Case conditions: T = %.1f K, P = %.1f mb, PhotoRates.in from the %s TUV table.' % (T, P, case))
    L.append('!!')
    L.append('!! Format: "!! <label>: <CMAQ statement>" precedes each translated reaction.')
    L.append('')
    notes_by = {}
    for n in tr.notes:
        notes_by.setdefault(n['label'], []).append(n)
    emitted = {x['label'] for x in tr.lines}
    for ln in tr.lines:
        L.append('!! %s: %s' % (ln['label'], ln['text']))
        for n in notes_by.get(ln['label'], []):
            tag = {'approx': 'APPROXIMATION', 'unsupported': 'WARNING', 'soa': 'NOTE', 'info': 'NOTE'}.get(n['level'])
            if tag:
                L.append('!!   %s: %s' % (tag, n['text']))
        L.append(ln['line'])
    L.append('')
    L.append('!! ------------------------------------------------------------------------------------')
    L.append('!! CMAQ statements that were NOT translated (core groups only):')
    for lab, ns in notes_by.items():
        if lab in emitted:
            continue
        for n in ns:
            if n['level'] in ('dropped', 'unsupported', 'zero', 'skipped', 'merged'):
                L.append('!! NOT TRANSLATED %s (%s): %s' % (lab, n['level'], n['text']))
                L.append('!!     [%s]' % n['cmaq'])
    L.append('')
    if not cmaq:
        L.append('!! ------------------------------------------------------------------------------------')
        L.append('!! SOA precursor counters -> IVOC volatility bins (unchanged from the %s baseline deck).' % case)
        L.append('!! Counters present in this mechanism: %s' % ', '.join(cfg['ivoc_counters']))
        L.append('!! ------------------------------------------------------------------------------------')
        L.extend(chain)
        L.append('!')
        L.extend(aging)
        L.append('!')
    open(out_path, 'w').write('\n'.join(L) + '\n')


def photo_file(cfg, case, case_dir, tr, T, out_path):
    tuv = mph.read_tuv(case_dir)
    pb = mph.PhotoBuilder(cfg['csqy'], tuv, T)
    src = open(os.path.join(case_dir, 'PhotoRates.in'), errors='ignore').read().replace('\r', '').split('\n')
    hdr_site = [ln for ln in src if ln.startswith('!!') and ('Site' in ln or 'PBL' in ln or 'lat' in ln or
                                                           'aerosol' in ln or 'Command' in ln)]
    cols = list(tr.photo.cols.items())
    vals, descs = [], []
    for (lab, sc), idx in cols:
        v, d = pb.column(lab, sc)
        vals.append(v)
        descs.append(d)
    out = ['!! PhotoRates.in  (%s column set)' % cfg['label'], '!! ASP',
           '!! Built by build_mech_decks.py / mech_photo.py from the 51-column TUV 5.4 table of the baseline deck of this case.',
           '!! "TUV col k" = same process in the TUV table; "CSQY ratio" = CMAQ cross-section x quantum-yield ratio relative to a',
           '!! TUV process (see mech_photo.py).  Columns beyond the last used one are zero.', '!!']
    out.extend(hdr_site)
    out.append('!!')
    out.append('!! Columns: time (UT, hrs), solar zenith angle (deg), then 51 photolysis rate coefficients (s-1)')
    for ((lab, sc), idx), d in zip(cols, descs):
        out.append('!  %2d = CMAQ %s%s  <- %s' % (idx, lab, '' if sc == 1 else ' x %g' % sc, d))
    for k in range(len(cols) + 1, 52):
        out.append('!  %2d = unused (0)' % k)
    out.append('!Columns: time, sza, photo-reactions')
    for i in range(tuv.shape[0]):
        new = [0.0] * 51
        for ((lab, sc), idx), v in zip(cols, vals):
            new[idx - 1] = v[i]
        out.append('%10s %10s ' % ('%g' % tuv[i, 0], '%g' % tuv[i, 1]) + ' '.join('%10.3e' % x for x in new))
    out.append('END-OF-FILE-012343210-ELIF-FO-DNE')
    open(out_path, 'w').write('\n'.join(out) + '\n')


def build_case(case, mech):
    cfg = m2a.CFG[mech]
    old_dir = os.path.join(DECKS, case)
    new_dir = os.path.join(DECKS, '%s_%s' % (case, cfg['suffix']))
    if os.path.isdir(new_dir):
        shutil.rmtree(new_dir)
    shutil.copytree(old_dir, new_dir)
    T, P = b6.read_case_TP(old_dir)
    tr = m2a.Translation(mech, T, P)
    table = ms.read_mw(mech)
    hdr, old_rows, old_raw = b6.read_old_gas(old_dir)
    ent = b6.case_entries(case, old_rows)
    mapped, rows = ms.apply(mech, ent)

    carry = []
    for nm in ('CO2', 'Trace', 'O3'):
        if nm not in mapped and nm in old_rows:
            v = old_rows[nm]
            mapped[nm] = [v[3], v[4]]
            carry.append((nm, v[3], v[4]))
    all_sp = OrderedDict()
    for nm in mapped:
        all_sp[nm] = 1
    for nm in tr.species:
        all_sp.setdefault(nm, 1)
    for nm in INTERFACE + PLACEHOLDERS + ['CO', 'CO2', 'Trace', 'NO', 'NO2', 'O3', 'PAN', 'NO3', 'N2O5', 'HONO']:
        all_sp.setdefault(nm, 1)
    cmaq = cfg.get('oa') in ('cmaq', 'ae')
    ae = tr.ae
    if cmaq:
        ivoc = [n for n in old_rows if re.match(r'^(POA|LEVO|CBIO|CPD)', n)]
        modes = cmaq_oa.read_modes(os.path.join(old_dir, 'AerosolModes.in'))
        if ae is not None:
            g_new, p_new, replaced, oa_notes = ae.emissions(old_rows, modes, T, P)
            for sp in list(p_new) + [x[0] for x in ae.partition(all_sp)] + ['PCVOC']:
                all_sp.setdefault(sp, 1)
        else:
            g_new, p_new, replaced, oa_notes = cmaq_oa.emissions(case, old_rows, modes, table, T, P)
        for sp, ppb in g_new.items():
            mapped[sp] = [ppb, 0.0]
            all_sp.setdefault(sp, 1)
    else:
        ivoc = [n for n in old_rows if re.match(r'^(IVOC|POA|LEVO|CBIO|CPD)', n) or n in cfg['ivoc_counters']]
    # counters that exist in the mechanism's chains but whose rows come from the old deck
    skip = {'CO', 'O2', 'N2', 'OffGas'} | set(ivoc)
    missing_old = [n for n in cfg['ivoc_counters'] if n not in old_rows]
    assert not missing_old, missing_old

    L = list(hdr)
    L[0] = '!! GasPhaseChems.in - %s case, %s' % (case, cfg['label'])
    L.append('CO\t;1\t;0;28;%g;%g;0!' % tuple(mapped.get('CO', [0, 0])))
    n_perox = 0
    for nm in all_sp:
        if nm in skip:
            continue
        ini, bg = mapped.get(nm, [0.0, 0.0])
        if nm in ('Trace', 'CO2') and nm in old_rows:
            mw = old_rows[nm][2]
        elif ae is not None and nm in ae.mw and nm not in ms.MW_EXTRA:
            mw = ae.mw[nm]
        elif ae is not None and nm in ae.gc_mw and nm not in table and nm not in ms.MW_EXTRA:
            mw = ae.gc_mw[nm]
        else:
            mw = ms.mw_of(mech, nm, table)
        pf = ms.peroxy_flag(mech, nm, table)
        n_perox += pf > 0
        L.append('%s\t;1\t;%d;%g;%.6g;%.6g;0!' % (nm, pf, mw, ini, bg))
    L.append('!')
    L.append('!VBS / organic-phase compounds referenced by the organic dissolution decks (unchanged from the %s deck)' % case)
    L.append('!MW from Grieshop et al.; ACP; 2009')
    if cmaq:
        L.append('!(IVOC1-9 of the baseline are replaced by the CMAQ V* species above)')
    for nm in ivoc:
        L.append(old_raw[nm])
    L.append('!')
    L.append('O2\t;0\t;0;32.;2.095e+08;2.095e+08;0!')
    L.append('N2\t;0\t;0;28.;7.808e+08;7.808e+08;0!')
    L.append('OffGas\t;0\t;0;12.;1000;1000;0!')
    L.append('END-OF-FILE-012343210-ELIF-FO-DNE')
    open(os.path.join(new_dir, 'GasPhaseChems.in'), 'w').write('\n'.join(L) + '\n')

    old_lines = open(os.path.join(old_dir, 'GasChemicalMechanism.in'), errors='ignore').read().replace('\r', '').split('\n')
    mech_file(cfg, case, tr, T, P, old_lines, os.path.join(new_dir, 'GasChemicalMechanism.in'))
    photo_file(cfg, case, old_dir, tr, T, os.path.join(new_dir, 'PhotoRates.in'))
    set_runlength(new_dir)
    if cmaq:
        part = ae.partition(list(all_sp)) if ae is not None else cmaq_oa.partitioning_species(list(all_sp), table)
        cmaq_oa.write_org_files(new_dir, part)
        cmaq_oa.rewrite_modes(os.path.join(new_dir, 'AerosolModes.in'), replaced, p_new)
        with open(os.path.join(new_dir, 'CMAQ_OA_MAPPING.txt'), 'w') as fh:
            fh.write('CMAQ %s organic-aerosol emissions mapping for %s (%s)\n' % ('AERO6/7' if ae is not None else 'CRACMM', case, cfg['label']))
            fh.write('\n'.join(oa_notes) + '\n\nPartitioning species (name, MW, C*298 ug/m3, dHvap J/mol, nonvolatile):\n')
            fh.write('\n'.join('%s %.2f %.4g %.0f %s' % x for x in part) + '\n')
            if ae is not None:
                fh.write('\nSOA counters replaced by mass-yield products (counter -> species: coefficient per counter mole):\n')
                for c in sorted(ae.used_ctr):
                    fh.write('%s -> %s\n' % (c, ' + '.join('%.4g %s' % (y, sp) for sp, y in ae.ctr[c])))

    with open(os.path.join(new_dir, 'SPECIES_MAPPING.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['paper_species', 'mapping_key', 'initial_ppb', 'background_ppb', '%s_representation' % cfg['label'], 'basis_note'])
        for src, key, ini, bg, comp, note in rows:
            rep = ' + '.join(('%.4g ' % cf if abs(cf - 1) > 1e-9 else '') + sp for sp, cf in comp)
            w.writerow([src, key, ini, bg, rep, note])
    return dict(case=case, mech=mech, T=T, P=P, tr=tr, rows=rows, mapped=mapped, carry=carry, nspecies=len(all_sp),
                nphoto=len(tr.photo.cols), nperoxy=n_perox, ivoc=ivoc)


def set_all_run_lengths():
    for d in sorted(glob.glob(os.path.join(DECKS, '*', 'ASPInputDeck.in'))):
        set_runlength(os.path.dirname(d))


if __name__ == '__main__':
    args = sys.argv[1:]
    mechs = [a for a in args if a in m2a.CFG] or list(m2a.CFG)
    cases = [a for a in args if a in CASES] or CASES
    for mech in mechs:
        for cs in cases:
            r = build_case(cs, mech)
            print('%-10s %-11s %d reactions, %d species (%d peroxy), %d photolysis columns' %
                  (cs, mech, len(r['tr'].lines), r['nspecies'], r['nperoxy'], r['nphoto']))
    set_all_run_lengths()
