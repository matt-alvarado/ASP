"""Parser for the GEOS-Chem fullchem.eqn (KPP) equations."""
import re
GCDIR = '/home/claude/geoschem/geos-chem/KPP/fullchem/'

def _sp(s):
    out = []
    for t in re.split(r'(?<![eEdD])\s*\+\s*(?=[A-Za-z0-9_.-])|\s+-\s+', s) if False else []:
        pass
    return out

def split_terms(s):
    """'2SO4 - SALAAL + 0.5 CH2O' -> [(2,'SO4'),(-1,'SALAAL'),(0.5,'CH2O')]"""
    s = ' ' + s.strip() + ' '
    toks = re.findall(r'([+-])\s*([0-9]*\.?[0-9]+(?:[dDeE][+-]?[0-9]+)?)?\s*([A-Za-z_][A-Za-z0-9_]*)', '+ ' + s.strip())
    out = []
    for sign, c, n in toks:
        c = float(c.replace('d', 'e').replace('D', 'e')) if c else 1.0
        out.append((c if sign == '+' else -c, n))
    return out

def read_eqn(path=GCDIR + 'fullchem.eqn'):
    txt = open(path).read().split('#EQUATIONS')[1]
    txt = re.sub(r'(?m)^\s*//.*$', '', txt)
    raw = txt
    # keep {+M}/{+O2}-style third-body annotations out; drop all brace comments
    txt = re.sub(r'\{[^}]*\}', ' ', txt)
    rx = []
    for st in txt.split(';'):
        st = ' '.join(st.split()).lstrip(', ')       # '{c1},{c2}' comment lists leave a stray comma
        if '=' not in st or ':' not in st: continue
        head, k = st.rsplit(':', 1)
        # the rate may itself contain ':'? (no) ; the equation part has no ':'
        l, r = head.split('=', 1)
        # rate expression may be missing for multi-statement; keep as is
        rx.append(dict(n=len(rx) + 1, re=split_terms(l), pr=split_terms(r), k=k.strip()))
    return rx

if __name__ == '__main__':
    import collections
    R = read_eqn(); print(len(R))
    c = collections.Counter()
    for r in R:
        for f in re.findall(r'([A-Za-z_][A-Za-z0-9_]*)\(', r['k']): c[f] += 1
        if not re.search(r'[A-Za-z_]\w*\(', r['k']): c['<const>'] += 1
    for k, v in c.most_common(): print(v, k)
    sp = set(n for r in R for _, n in r['re'] + r['pr']); print(len(sp))
