#!/usr/bin/env python3
"""Write NitrateGroups.in for ASP input decks.

NitrateGroups.in lists the gas-phase species that ASP sums for the
dTotPNs/dCO and dTotANs/dCO columns of EnhanceRatios.txt
(OutputRoutines.f90, GasEnhanceOutput / ReadNitrateGroups):

  PN  peroxy nitrates: RC(O)OONO2 (PAN, PPN, MPAN, ...) and RO2NO2 (CH3O2NO2);
      not HO2NO2 (inorganic)
  AN  alkyl nitrates and other organic nitrates (RONO2, including
      multifunctional, hydroxy, carbonyl and hydroperoxy nitrates and
      dinitrates); not nitro compounds (nitrophenols, nitrocresols, ...),
      nitriles (HCN, CH3CN) or aerosol-phase nitrates

Without NitrateGroups.in ASP uses its original positional sums (33 species
from PAN and 77 from CH3NO3), which match only the baseline mechanism's
species order, so older decks still run unchanged.

The lists are derived from each deck's own GasPhaseChems.in and
GasChemicalMechanism.in, then corrected per mechanism (CORRECTIONS below,
reviewed by hand against the mechanisms):

1. Nitrogen atoms per species: inorganic N species are known; species whose
   name has no letter N are taken as N-free unless at least two reactions
   (with every other species known) imply the same N >= 1 and none implies 0;
   the rest follow from N conservation (a reaction with one unknown species).
2. PN: a closed-shell species formed by  R + NO2 (+M) => X  that decomposes
   thermally (a non-photolysis unimolecular reaction) back to NO2.
3. Radicals (excluded): species that react with NO or HO2 (peroxy radicals),
   that decompose thermally to NO2 without being PNs (nitrooxy alkoxy radicals,
   Criegee intermediates), or that are alkoxy radicals by their formation
   (ROOH photolysis => RO + OH; RO2 + NO/NO3 => RO + NO2).
4. AN: organic, closed-shell species with N >= 1 that are not PNs.
5. MCM only (explicit chemistry): nitro compounds, i.e. products of a
   phenoxy-type radical (no reaction with NO, no thermal reaction) + NO2 that
   do not decompose back, and their products keeping the same N, are not AN.
   In the lumped mechanisms this rule also catches nitrates, so their nitro
   species are excluded by name in CORRECTIONS instead.
6. Species with no reactions at all (emitted tracers) are included only if
   listed in CORRECTIONS (methyl nitrate, CH3NO3, in every mechanism).

Usage: python3 make_nitrate_groups.py [deck_dir ...]
       (default: every deck directory next to this script with a
        GasPhaseChems.in and a GasChemicalMechanism.in, and this directory)
"""
import collections
import os
import re
import sys

INORG_N = {'NO': 1, 'NO2': 1, 'NO3': 1, 'N2O5': 2, 'HNO3': 1, 'HONO': 1,
           'HNO2': 1, 'HNO4': 1, 'HO2NO2': 1, 'PNA': 1, 'NH3': 1, 'N2O': 2,
           'N2': 2, 'N': 1, 'CLNO2': 1, 'CLNO3': 1, 'ClNO2': 1, 'ClNO3': 1,
           'BRNO2': 1, 'BRNO3': 1, 'BrNO2': 1, 'BrNO3': 1}
INORGANIC = set(INORG_N) | {
    'H2O', 'M', 'O2', 'O3', 'OH', 'HO2', 'H2O2', 'O', 'O1D', 'O3P', 'H2', 'SO2',
    'SULF', 'H2SO4', 'CO', 'CO2', 'HCL', 'HCl', 'CL', 'Cl', 'CL2', 'Cl2', 'HOCL',
    'HOCl', 'CLO', 'ClO', 'BR', 'Br', 'HBR', 'HBr', 'BRO', 'BrO', 'HOBR', 'HOBr',
    'RO2', 'RO2_T', 'RCO3_T'}
NUM = re.compile(r'^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$')

# Per-mechanism corrections to the derived lists (see the module docstring).
CORRECTIONS = {
    'baseline': dict(
        add_pn=['PAN_lump'],            # lumped PAN, no reactions in the deck
        add_an=['C42AOH', 'C42OH'],     # hydroxy nitrates with N-free names
                                        # (photolyse/oxidise to nitrates without
                                        # releasing NO2); one hand-written NO3
                                        # reaction makes C42AOH look like an RO
        drop_an=['HCN']),               # nitrile
    'CB6r5': dict(drop_an=['CRON']),    # nitrocresol (nitro compound)
    'SAPRC07TC': dict(drop_an=['NPHE']),   # nitrophenols
    'SAPRC07TIC': dict(drop_an=['NPHE']),
    'CRACMM2': dict(drop_an=['VROCN1OXY1']),  # volatility-bin lump (gets both the
    'CRACMM3': dict(drop_an=['VROCN1OXY1']),  # RO2+NO nitrate branch and RO2+HO2)
    'MOZT1': dict(),
    'GCHEM': dict(drop_an=['NPHEN',     # nitrophenols
                           'INA',       # nitrooxy alkoxy radical
                           'IONITA', 'MONITA', 'AONITA']),  # aerosol-phase nitrates
    'MCM331': dict(),
}
ALWAYS_ADD_AN = ['CH3NO3']       # methyl nitrate, also as an unreactive tracer
MECHANISMS = ('CB6r5', 'SAPRC07TC', 'SAPRC07TIC', 'CRACMM2', 'CRACMM3',
              'MOZT1', 'GCHEM', 'MCM331')


def mechanism_of(deck):
    base = os.path.basename(os.path.normpath(os.path.abspath(deck)))
    for m in MECHANISMS:
        if base.endswith('_' + m):
            return m
    return 'baseline'


def species(deck):
    names = []
    for line in open(os.path.join(deck, 'GasPhaseChems.in'), errors='replace'):
        t = line.strip()
        if not t or t.startswith('!') or ';' not in t:
            continue
        n = t.split(';')[0].strip()
        if re.match(r'^[A-Za-z0-9_\-+()]+$', n) and not NUM.match(n):
            names.append(n)
    return names


def side(s):
    out = []
    for t in re.split(r'\s+\+\s+|\s*\+\s+|\s+\+\s*', s):
        t = t.strip()
        if not t:
            continue
        parts = t.split()
        coef = float(parts[0].replace('d', 'e').replace('D', 'e')) \
            if len(parts) == 2 and NUM.match(parts[0]) else 1.0
        out.append((parts[-1], coef))
    return out


def reactions(deck):
    rx = []
    for line in open(os.path.join(deck, 'GasChemicalMechanism.in'), errors='replace'):
        t = line.strip()
        if not t or t.startswith('!') or '=>' not in t:
            continue
        eq, _, par = t.partition(';')
        lhs, _, rhs = eq.partition('=>')
        pars = [p.strip() for p in par.split(';') if p.strip()]
        # photolysis: first rate parameter 0, then the PhotoRates.in column
        photo = bool(pars) and bool(NUM.match(pars[0])) and \
            float(pars[0].replace('d', 'e')) == 0.0 and 2 <= len(pars) <= 3
        rx.append((side(lhs), side(rhs), photo))
    return rx


def classify(deck, nitro_rule):
    names = species(deck)
    rx = reactions(deck)

    # 1. nitrogen atoms per species
    N, firm = {}, {}
    for k in INORGANIC:
        N[k] = INORG_N.get(k, 0); firm[k] = True
    weak = set()
    for n in names:
        if n not in N and 'N' not in n.upper():
            N[n] = 0; weak.add(n)
    nets = []
    for lhs, rhs, photo in rx:
        net = collections.defaultdict(float)
        for s, c in lhs: net[s] -= c
        for s, c in rhs: net[s] += c
        nets.append({k: v for k, v in net.items() if abs(v) > 1e-12})
    changed = True
    while changed:
        changed = False
        for net in nets:
            unk = [s for s in net if s not in N]
            if len(unk) != 1:
                continue
            u = unk[0]
            val = -sum(c * N[s] for s, c in net.items() if s != u) / net[u]
            if abs(val - round(val)) < 1e-6 and 0 <= round(val) <= 4:
                N[u] = int(round(val)); changed = True
    votes = collections.defaultdict(list)
    for net in nets:
        w = [s for s in net if s in weak]
        if len(w) != 1 or any(s not in N for s in net):
            continue
        w = w[0]
        val = -sum(c * N[s] for s, c in net.items() if s != w) / net[w]
        if abs(val - round(val)) < 1e-6:
            votes[w].append(int(round(val)))
    for w, v in votes.items():
        pos = [x for x in v if x >= 1]
        if len(pos) >= 2 and len(set(pos)) == 1 and 0 not in v:
            N[w] = pos[0]

    # reaction patterns
    react_with = collections.defaultdict(set)
    thermal = collections.defaultdict(list)       # non-photolysis unimolecular
    assoc_no2 = set()                             # R + NO2 (+M) => X
    has_rxn = set()
    to_no_no3 = []
    alkoxy = set()
    for lhs, rhs, photo in rx:
        rs = [s for s, _ in lhs if s != 'M']
        ps = [s for s, _ in rhs if s != 'M']
        has_rxn.update(rs); has_rxn.update(ps)
        if len(rs) == 2:
            a, b = rs
            react_with[a].add(b); react_with[b].add(a)
            if 'NO2' in rs and len(ps) == 1:
                assoc_no2.add(ps[0])
            other = b if a in ('NO', 'NO3') else a if b in ('NO', 'NO3') else None
            if other:
                to_no_no3.append((other, rhs))
        if len(rs) == 1 and photo:
            prods = set(ps) - {'RO2', 'RO2_T', 'RCO3_T'}
            if len(prods) == 2 and 'OH' in prods:          # ROOH + hv => RO + OH
                alkoxy |= prods - {'OH'}
        if len(rs) == 1 and not photo:
            thermal[rs[0]].append(set(ps))
    for other, rhs in to_no_no3:                           # RO2 + NO => RO + NO2
        if not {'NO', 'HO2'} <= react_with[other]:
            continue
        cs = {s: c for s, c in rhs if s not in ('RO2', 'RO2_T', 'RCO3_T')}
        if len(cs) == 2 and 'NO2' in cs:
            y = [s for s in cs if s != 'NO2'][0]
            if abs(cs[y] - cs['NO2']) < 1e-9:
                alkoxy.add(y)

    # 5. nitro compounds (MCM)
    nitro = set()
    if nitro_rule:
        for lhs, rhs, photo in rx:
            rs = [s for s, _ in lhs if s != 'M']; ps = [s for s, _ in rhs if s != 'M']
            if len(rs) == 2 and 'NO2' in rs and len(ps) == 1:
                x = rs[0] if rs[1] == 'NO2' else rs[1]
                y = ps[0]
                if 'NO' not in react_with[x] and not thermal.get(x) and x not in INORGANIC \
                   and not any('NO2' in p for p in thermal.get(y, [])):
                    nitro.add(y)
        grow = True
        while grow:
            grow = False
            for lhs, rhs, photo in rx:
                src = [s for s, _ in lhs if s in nitro]
                if not src:
                    continue
                for p, _ in rhs:
                    if p not in nitro and p not in INORGANIC and N.get(p, 0) >= 1 \
                       and N.get(p, 0) == N.get(src[0], -1):
                        nitro.add(p); grow = True

    PN, AN = [], []
    for n in names:
        if n in INORGANIC or n not in has_rxn:
            continue
        peroxy = bool(react_with[n] & {'NO', 'HO2'})
        if not peroxy and n in assoc_no2 and any('NO2' in p for p in thermal.get(n, [])):
            PN.append(n)
            continue
        if peroxy or n in alkoxy or any('NO2' in p for p in thermal.get(n, [])):
            continue                                       # radical
        if N.get(n, 0) >= 1 and n not in nitro:
            AN.append(n)
    return names, PN, AN


def groups(deck):
    mech = mechanism_of(deck)
    names, PN, AN = classify(deck, nitro_rule=(mech == 'MCM331'))
    corr = CORRECTIONS[mech]
    nameset = set(names)
    pn = set(PN) | (set(corr.get('add_pn', [])) & nameset)
    an = (set(AN) | (set(corr.get('add_an', []) + ALWAYS_ADD_AN) & nameset)) \
        - set(corr.get('drop_an', [])) - pn
    order = {n: i for i, n in enumerate(names)}
    return mech, sorted(pn, key=order.get), sorted(an, key=order.get)


def write(deck):
    mech, pn, an = groups(deck)
    lines = [
        '!! NitrateGroups.in',
        '!! Gas-phase species summed for the dTotPNs/dCO (PN) and dTotANs/dCO (AN)',
        '!! columns of EnhanceRatios.txt. One species per line: "PN ; name" or',
        '!! "AN ; name"; names not in GasPhaseChems.in are skipped with a warning.',
        '!!   PN: peroxy nitrates, RC(O)OONO2 and RO2NO2 (not HO2NO2)',
        '!!   AN: alkyl and other organic nitrates, RONO2 (not nitro compounds,',
        '!!       nitriles or aerosol-phase nitrates)',
        '!! Without this file ASP uses its original positional sums (33 species',
        '!! from PAN, 77 from CH3NO3), valid only for the baseline mechanism.',
        '!! Generated by make_nitrate_groups.py from this deck\'s GasPhaseChems.in',
        '!! and GasChemicalMechanism.in (mechanism: %s): %d PN, %d AN.' % (mech, len(pn), len(an)),
        '!!',
    ]
    lines += ['PN ; %s' % n for n in pn]
    lines += ['AN ; %s' % n for n in an]
    with open(os.path.join(deck, 'NitrateGroups.in'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return mech, len(pn), len(an)


def decks_next_to_script():
    here = os.path.dirname(os.path.abspath(__file__))
    out = []
    for d in [here] + sorted(os.path.join(here, x) for x in os.listdir(here)):
        if os.path.isdir(d) and os.path.isfile(os.path.join(d, 'GasPhaseChems.in')) \
           and os.path.isfile(os.path.join(d, 'GasChemicalMechanism.in')):
            out.append(d)
    return out


if __name__ == '__main__':
    for d in sys.argv[1:] or decks_next_to_script():
        mech, npn, nan = write(d)
        print('%-45s %-10s PN %4d  AN %4d' % (os.path.relpath(d), mech, npn, nan))
