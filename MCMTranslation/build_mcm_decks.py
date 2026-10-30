#!/usr/bin/env python3
"""Build InorganicInputDecks/<Case>_MCM331/ (hybrid MCM v3.3.1 + retained baseline lumped chemistry) for the four cases.

  GasChemicalMechanism.in : MCM v3.3.1 gas-phase mechanism (mcm2asp.py; every rate structural, T/P independent)
                            + baseline reactions whose reactants include a species the MCM does not contain
                            + SOA-precursor counters on the MCM first-generation VOC reactions
  GasPhaseChems.in        : all MCM species + retained baseline species; initial / background values carried over from the
                            baseline deck (baseline compound names renamed to MCM names where the MCM treats them explicitly)
  everything else (including PhotoRates.in, the 51-column TUV table) is copied unchanged from the baseline deck.
"""
import csv, glob, os, re, shutil, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'CB6r5Translation'))
import build_cb6r5_decks as b6                       # noqa: E402
from mcm2asp import *                                # noqa: E402
from baseline_mech import read_baseline, ren, renamed_line   # noqa: E402

DECKS = b6.DECKS
CASES = b6.CASES
OMIT_MCM = {33: 'HNO3 -> NA: crude first-order uptake to nitrate aerosol; ASP computes gas-aerosol exchange explicitly',
            34: 'N2O5 -> 2 NA: crude first-order uptake to nitrate aerosol; ASP computes gas-aerosol exchange explicitly'}
NEUTRAL = {'M', 'H2O', 'O2', 'N2', 'RO2_T', 'RCO3_T', 'hv'}
SPECIAL_KEEP = ['CO', 'CO2', 'Trace', 'NO', 'NO2', 'O3', 'CH3NO3', 'RO2_T', 'RCO3_T', 'NH3', 'H2SO4', 'HCl', 'HNO3', 'SO2']


def header_text(n_mcm, n_keep, n_drop):
    return """!! GasChemicalMechanism.in - hybrid MCM v3.3.1 + baseline lumped chemistry
!! ASP
!! Translated from the Master Chemical Mechanism v3.3.1 (KPP files, mcm.york.ac.uk) by MCMTranslation/mcm2asp.py.
!! MCM citation: Jenkin et al. (1997), Saunders et al. (2003), Jenkin et al. (2003), Bloss et al. (2005), Jenkin et al. (2012, 2015).
!!
!! Contents
!!   1. MCM v3.3.1 gas-phase reactions (%d ASP reactions; MCM reaction N is written in the same order, a reaction whose rate is a
!!      sum of terms is split into one ASP reaction per term).
!!   2. Baseline (RACM2-derived, Williams et al. 2015) reactions that have at least one reactant the MCM does not contain
!!      (%d kept, %d replaced by the MCM): lumped VOC/RO2 chemistry, NH3, HCN, the SOA-counter -> IVOC chains and IVOC aging.
!!      Baseline compounds that the MCM treats explicitly are renamed to MCM names (see MCM_SPECIES_MAPPING.csv).
!!   3. SOA-precursor counters (ISOSOA, TOLSOA, XYLSOA, TERPSOA) appended to the MCM OH / O3 / NO3 reactions of isoprene, toluene,
!!      the xylenes and the terpenes, as in the baseline deck.
!!   Not translated: MCM reactions %s (see MCM_NOT_TRANSLATED.txt).
!!
!! TEMPERATURE / PRESSURE: no rate constant is evaluated at a fixed T or P.  All MCM rate expressions (Arrhenius forms, Troe
!! fall-off with M, O2/N2/H2O terms, RO2-pool terms, photolysis) are expressed through the ASP rate flags below, so this file is valid at
!! any T and P; only PhotoRates.in and GasPhaseChems.in are case specific.  verify_rates.py checks all %d translated rates against the
!! MCM Fortran expressions on a grid T = 230-320 K, P = 250-1013 mb (agreement < 0.2 %%).
!! Air composition where the MCM writes N2 / O2 explicitly: N2 = 0.79 M, O2 = 0.21 M (as in ASP flags 8 and 9).
!!
!! Rate-flag conventions used here (see also the baseline deck headers):
!!   1-body   0; col [; scale]   photolysis: PhotoRates.in column x optional scale.  MCM J(j) is mapped to the baseline TUV column of the same
!!                                process; scale = J_MCM(target)/J_MCM(analog), both MCM parameterisations at zenith 0 (1 = omitted)
!!            2; A; B; C; D      A T^B exp(C/T)^D          12; A; B; C; D   RO2_T A T^B exp(C/T)^D (MCM terms k*RO2)
!!            8                  O1D + M (hard-wired)      9                O + O2 + M (hard-wired)
!!            13; A0;n0;E0;Ai;ni;Ei;F;N;KeqA;KeqE   fall-off  k0 = A0 (T/300)^-n0 exp(E0/T) M, ki likewise, Troe F, N (KeqA=1, KeqE=0: no equilibrium)
!!            14; A; B; C; E     A T^B exp(C/T) exp(E/T^3)  (tunnelling, MCM isoprene H-shifts)               [common_asp 718cca3]
!!   2-body   A; B; C; D         A T^B exp(C/T)^D          4; A1;C1;A2;C2   k1 + k2 M          5; A1;C1;A2;C2;A3;C3   k1 + k3 M/(1 + k3 M/k2)
!!            7 / 8 / 9; A;B;C;D k*(1-1/(1+c e^{d/T})) / k/(1+c e^{d/T}) / k*(1-c e^{d/T})     10; A1;C1;A2;C2   [H2O](k1 + M k2)
!!   3-body   2; A; B; C; D      M A T^B exp(C/T)^D        4; A0;n0;E0;Ai;ni;Ei;F;N  fall-off with explicit Troe N          [common_asp 40c1760]
!!            5; A0;n0;E0;Ai;ni;Ei;Fa;Fb;Fc  as 4 with F(T) = Fa exp(-Fb/T) + exp(-T/Fc), N = 0.75-1.27 log10 F   [common_asp 710a126]
!! REQUIRES common_asp >= 710a126 (flags: photolysis scale e4888c4, 1-body 14 718cca3, 3-body 5 710a126, plus 40c1760 and earlier fixes).
!!
""" % (n_mcm, n_keep, n_drop, ', '.join(str(k) for k in sorted(OMIT_MCM) + [15631]), n_mcm)


def species_rows(case_rows, used, tr):
    sp = tr.species
    ro2 = tr.ro2
    out = collections.OrderedDict()
    for n in used:
        base = case_rows.get(n)
        if n in sp:
            mw = sp[n]['mass']
            if mw is None or n in ('H2SO4', 'HCl'):
                mw = base[2] if base else {'H2SO4': 98.08, 'HCl': 36.46}[n]
            pf = 0
            if n in ro2:
                pf = 2 if re.search(r'C\(=O\)O\[O\]|OO\[|\[O\]OC\(=O\)', sp[n]['smiles']) else 1
            if n == 'HO2': pf = 0
            ini, bg = (base[3], base[4]) if base else (0.0, 0.0)
            ev, dep = (base[0], base[5]) if base else (1.0, 0.0)
            out[n] = (ev, pf, mw, ini, bg, dep)
        else:
            out[n] = tuple(base)
    return out


def build_case(case, tr, mcm_lines, mcm_nums):
    old_dir = os.path.join(DECKS, case)
    new_dir = os.path.join(DECKS, '%s_MCM331' % case)
    if os.path.isdir(new_dir): shutil.rmtree(new_dir)
    shutil.copytree(old_dir, new_dir)
    hdr, old_rows, old_raw = b6.read_old_gas(old_dir)
    base = read_baseline(os.path.join(old_dir, 'GasChemicalMechanism.in'))
    mcm_set = set(tr.species_set)
    keep, drop = [], []
    for r in base:
        re_ = [ren(n) for c, n in r['lhs']]
        (drop if all(n in mcm_set or n in NEUTRAL for n in re_) else keep).append(r)
    keep_lines = [renamed_line(r) for r in keep]
    # species used
    used = collections.OrderedDict()
    def add(n):
        if n in NEUTRAL or n in ('PROD',): return
        used.setdefault(n, 1)
    for l in mcm_lines:
        head = l.split(' ; ')[0]
        for side in head.split('=>'):
            for t in side.split('+'):
                t = t.strip()
                m = re.match(r'^[0-9.]+\s+(\S+)$', t)
                add(m.group(1) if m else t)
    for r in keep:
        for c, n in r['lhs'] + r['rhs']: add(ren(n))
    # baseline rows with non-zero init/bg must not be lost
    renamed_rows = collections.OrderedDict()
    for n, v in old_rows.items():
        rn = ren(n)
        if rn in renamed_rows:
            renamed_rows[rn] = [v[0], v[1], v[2], renamed_rows[rn][3] + v[3], renamed_rows[rn][4] + v[4], v[5]]
        else:
            renamed_rows[rn] = list(v)
    lost = []
    for n, v in renamed_rows.items():
        if n in NEUTRAL: continue
        if n not in used:
            if v[3] != 0 or v[4] != 0 or n in SPECIAL_KEEP or re.match(r'^(IVOC|POA|LEVO|CBIO|CPD)', n):
                used.setdefault(n, 1)
    for n in SPECIAL_KEEP: used.setdefault(n, 1)
    rows = species_rows(renamed_rows, [n for n in used if n != 'CO'], tr)
    for n in used:
        if n not in tr.species and n not in renamed_rows:
            raise SystemExit('species %s used in %s retained baseline reactions has no species row' % (n, case))
    L = list(hdr)
    L[0] = '!! GasPhaseChems.in - %s case, MCM v3.3.1 + baseline lumped species' % case
    c0 = renamed_rows.get('CO', [1, 0, 28, 0, 0, 0])
    L.append('CO\t;1\t;0;28;%g;%g;0!' % (c0[3], c0[4]))
    n_perox = 0
    for n, (ev, pf, mw, ini, bg, dep) in rows.items():
        if n in ('CO', 'O2', 'N2', 'OffGas'): continue
        n_perox += pf > 0
        L.append('%s\t;%d\t;%d;%g;%.8g;%.8g;%g!' % (n, int(ev), int(pf), mw, ini, bg, dep))
    L.append('!')
    for nm in ('O2', 'N2', 'OffGas'):      # baseline rows (O2 = 0.21 M, N2 = 0.79 M, as assumed by ASP flags 8 and 9)
        L.append(old_raw[nm].rstrip())
    L.append('END-OF-FILE-012343210-ELIF-FO-DNE')
    open(os.path.join(new_dir, 'GasPhaseChems.in'), 'w').write('\n'.join(L) + '\n')

    M = [header_text(len(mcm_lines), len(keep), len(drop))]
    M.append('!! ---------------- 1. MCM v3.3.1 ----------------')
    M.extend(mcm_lines)
    M.append('')
    M.append('!! ---------------- 2. baseline reactions kept (reactant not in the MCM; compound names renamed to MCM names) ------')
    M.extend(keep_lines)
    M.append('')
    open(os.path.join(new_dir, 'GasChemicalMechanism.in'), 'w').write('\n'.join(M) + '\n')
    # mapping csv
    with open(os.path.join(new_dir, 'MCM_SPECIES_MAPPING.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['baseline_name', 'model_name', 'initial_ppb', 'background_ppb', 'treatment'])
        for n, v in old_rows.items():
            rn = ren(n)
            tmt = 'renamed to MCM species' if rn != n else ('MCM species (same name)' if rn in mcm_set else
                                                          ('lumped species kept, baseline chemistry' if rn in used else 'not used'))
            w.writerow([n, rn, v[3], v[4], tmt])
    b6_set = b6
    return dict(case=case, n_keep=len(keep), n_drop=len(drop), n_species=len(rows), n_perox=n_perox, keep=keep, drop=drop)


def main():
    tr = MCMTranslator()
    tr.species_set = set(RENAME_MCM.get(n, n) for n in tr.species)
    R = read_eqn()
    mcm_lines, nums = [], []
    for r in R:
        if r['n'] in OMIT_MCM: continue
        ls = tr.translate(r)
        mcm_lines.extend(ls); nums.extend([r['n']] * len(ls))
    print('MCM: %d reactions -> %d ASP lines; skipped %s' % (len(R), len(mcm_lines), [s[0] for s in tr.skipped]))
    with open(os.path.join(HERE, 'MCM_NOT_TRANSLATED.txt'), 'w') as fh:
        for k, v in OMIT_MCM.items(): fh.write('MCM reaction %d: %s\n' % (k, v))
        for n, eq, why in tr.skipped: fh.write('MCM reaction %d (%s): %s\n' % (n, eq, why))
    with open(os.path.join(HERE, 'MCM_reaction_index.csv'), 'w') as fh:
        fh.write('asp_line_in_block,mcm_reaction\n')
        for i, n in enumerate(nums): fh.write('%d,%d\n' % (i + 1, n))
    cases = [a for a in sys.argv[1:] if a in CASES] or CASES
    for cs in cases:
        res = build_case(cs, tr, mcm_lines, nums)
        print('%-10s kept %d / replaced %d baseline reactions; %d species (%d peroxy)' %
              (cs, res['n_keep'], res['n_drop'], res['n_species'], res['n_perox']))


if __name__ == '__main__':
    main()
