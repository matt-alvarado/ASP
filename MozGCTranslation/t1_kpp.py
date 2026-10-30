"""Parser for the WRF-Chem KPP T1_MOZCART equation file."""
import re
T1DIR = '/home/claude/geoschem/wrf/chem/KPP/mechanisms/t1_mozcart/'

def _sp(s):
    out = []
    for t in s.split('+'):
        t = t.strip()
        if not t: continue
        m = re.match(r'^([0-9]*\.?[0-9]+)\s+(\S+)$', t)
        out.append((float(m.group(1)), m.group(2)) if m else (1.0, t))
    return out

def read_eqn(path=T1DIR + 't1_mozcart.eqn'):
    rx = []
    for ln in open(path):
        ln = ln.rstrip('\n')
        m = re.match(r'^\{(\d+):(\w+)\}\s*(.*?)\s*=\s*(.*?)\s*:\s*(.*?)\s*;\s*$', ln)
        if not m: continue
        n, lab, l, r, k = m.groups()
        rx.append(dict(n=int(n), lab=lab, re=_sp(l), pr=_sp(r), k=k.strip()))
    return rx

if __name__ == '__main__':
    R = read_eqn(); print(len(R), R[0], R[-1])
    sp = set(n for r in R for _, n in r['re'] + r['pr']); print(len(sp), sorted(sp))
