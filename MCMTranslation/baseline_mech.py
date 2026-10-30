import re, collections
from mcm2asp import BASE2MCM, RENAME_MCM, POOL

def split_species(s):
    out = []
    for t in s.split('+'):
        t = t.strip()
        if not t: continue
        m = re.match(r'^(-?[0-9]*\.?[0-9]+)\s+(\S+)$', t)
        out.append((float(m.group(1)), m.group(2)) if m else (1.0, t))
    return out

def read_baseline(path):
    """-> header lines, list of dict(lhs=[(c,n)], rhs=[(c,n)], tail=str, raw=str, idx=int) for every reaction line"""
    rows = []
    for i, ln in enumerate(open(path, errors='ignore').read().replace('\r', '').split('\n')):
        s = ln.strip()
        if not s or s.startswith('!') or '=>' not in s: continue
        body = s.split('!')[0] if ';' not in s.split('!')[0][-1:] and False else s
        m = re.match(r'^(.*?)=>(.*?);(.*)$', s)
        if not m: raise ValueError('cannot parse: ' + ln)
        rows.append(dict(lhs=split_species(m.group(1)), rhs=split_species(m.group(2)), tail=m.group(3), raw=ln, idx=i))
    return rows

def ren(n):
    return BASE2MCM.get(n, RENAME_MCM.get(n, n))

def fmt_side(sp):
    return ' + '.join(n if abs(c - 1) < 1e-12 else '%s %s' % (('%.6g' % c), n) for c, n in sp)

def renamed_line(r):
    l = [(c, ren(n)) for c, n in r['lhs']]
    p = collections.OrderedDict()
    pr = []
    for c, n in r['rhs']:
        pr.append((c, ren(n)))
    return '%s => %s ;%s' % (fmt_side(l), fmt_side(pr), r['tail'])
