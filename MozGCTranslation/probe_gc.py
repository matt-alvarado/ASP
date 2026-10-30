#!/usr/bin/env python3
"""Live ASP test of the new rate flags: build a probe deck (one reaction per flag/body pattern, translated GEOS-Chem parameters, trace reactant
decaying in a large excess of its partner), run ASP, and compare the measured decay rate with the independent Python mirror."""
import os, re, sys, shutil, math, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gc_rates as g, verify_gc as v
SRC = '/home/claude/asp/InorganicInputDecks/Williams'; OUT = '/home/claude/work/probe'
T, P = 288.4, 880.0
def build():
    tr = g.GCTranslator(); pick = collections.OrderedDict()
    for r in g.read_eqn():
        if re.search(r'PHOTOL|uptk|SR_MW', r['k']): continue
        ls, _ = tr.translate(r)
        for l in ls:
            head, rest = l.split(' ; ', 1)
            reac = [x.strip() for x in head.split('=>')[0].split('+')]
            if len(set(reac)) < len(reac) and 'M' not in reac: continue
            nums = [float(t) for t in rest.strip().rstrip(';').split(';') if t.strip()]
            key = (len([x for x in reac if x != 'M']), 'M' in reac, int(nums[0]) if int(nums[0]) != 0 else 0)
            pick.setdefault(key, [])
            if len(pick[key]) < 2: pick[key].append(l)
    return pick
def main():
    pick = build(); bpp = {}; rows = []; mech = []; spec = []; exp = []
    i = 0; spec.append(('PRC', 1e-6)); spec += [('PRA0', 1e-6), ('PRB0', 1e-2), ('PRQ0', 0.0)]; mech.append('PRA0 + PRB0 => PRQ0 ; 1e-12; 0; 0; 1;')
    for key, ls in pick.items():
        for l in ls:
            i += 1; head, rest = l.split(' ; ', 1)
            reac = [x.strip() for x in head.split('=>')[0].split('+')]
            nb = [x for x in reac if x != 'M']; fixed = [x for x in nb if x in ('O2', 'N2', 'H2O')]
            free = [x for x in nb if x not in fixed]
            names = []
            for j, x in enumerate(free):
                nm = 'PRA%d' % i if j == 0 else 'PRB%d' % i; names.append(nm)
            kref = v.asp_rate(l, T, 1e-6 and (P * 100. / (1.380649e-23 * T) * 1e-6), 0.0)
            Mgas = P * 100. / (1.380649e-23 * T) * 1e-6
            for x in reac:
                if x == 'O2': kref *= .21 * Mgas
                if x == 'N2': kref *= .79 * Mgas
            if len(free) == 1 and kref > 5e-3: i -= 1; continue          # too fast to measure (O + O2 + M)
            spec.append(('PRA%d' % i, 1e-6))
            if len(free) == 2:
                bppb = min(max(1.5e-3 / kref / (Mgas * 1e-9), 1e-9), 5e3); spec.append(('PRB%d' % i, bppb)); bpp[i] = bppb
            spec.append(('PRQ%d' % i, 0.0))
            lhs = names + fixed + (['M'] if 'M' in reac else [])
            mech.append('%s => PRQ%d ; %s' % (' + '.join(lhs), i, rest.strip())); exp.append((i, l, len(free)))
    base = open(SRC + '/GasPhaseChems.in').read().split('\n')
    k = [n for n, x in enumerate(base) if x.startswith('CO\t')][0]
    new = base[:k + 1] + ['%s\t;1\t;0;100;%g;0;0!' % (n, c) for n, c in spec] + base[k + 1:]
    shutil.rmtree(OUT, ignore_errors=True); shutil.copytree(SRC, OUT + '/InorganicInputDecks'); os.makedirs(OUT + '/OutputFiles')
    d = OUT + '/InorganicInputDecks/'
    open(d + 'GasPhaseChems.in', 'w').write('\n'.join(new))
    open(d + 'GasChemicalMechanism.in', 'w').write('\n'.join(['!! probe'] + mech) + '\n')
    s = open(d + 'ASPInputDeck.in').read()
    s = re.sub(r'(?m)^(\S+)(\s+!! Temperature \(K\))', '%g\\2' % T, s)
    s = re.sub(r'(?m)^(\S+)(\s+!! Pressure    \(mb\))', '%g\\2' % P, s)
    s = re.sub(r'(?m)^(\S+)(\s+!! Run Length in minutes)', '20\\2', s)
    s = re.sub(r'(?m)^(\S+)(\s+!! Include aerosols\?)', '0\\2', s)
    open(d + 'ASPInputDeck.in', 'w').write(s)
    import pickle; pickle.dump((exp, bpp), open(OUT + '/exp.pkl', 'wb'))
    print(len(exp), 'probe reactions')
if __name__ == '__main__': main()
