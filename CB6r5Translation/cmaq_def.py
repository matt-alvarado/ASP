"""Parser for CMAQ mechanism-definition (mech_*.def) files.

Only the pieces needed to translate gas-phase reactions to ASP are handled.
A reaction statement looks like

    <R19> HO2 + HO2 = H2O2  %3 # 2.20E-13@-600 & 1.90E-33@-980;

and may continue over several lines until the terminating ';'.
"""
import re
from dataclasses import dataclass, field

NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?'


@dataclass
class Rxn:
    label: str
    reactants: list          # [(coef, name)]
    products: list           # [(coef, name)]
    ptype: str               # '' or '%2','%3','%H' ...
    rate_raw: str            # text after '#'
    text: str                # normalised one-line statement (for documentation)
    # filled by classify_rate():
    kind: str = ''           # 'arr','falloff','photo','ref','hetero','op','pct'
    terms: list = field(default_factory=list)   # parsed & terms [(A, n, E)]
    scale: float = 1.0
    target: str = ''         # photolysis label / referenced reaction / hetero label


def _parse_side(s):
    """Parse 'A + 0.5*B + 2.0*C -1.2*D' -> [(coef,name)]; empty -> []."""
    s = s.strip()
    if not s:
        return []
    # make every sign explicit:  ' -2.5*PAR' -> ' + -2.5*PAR'
    s = re.sub(r'(?<=[\w\s*])\s-(?=\s*[\d.]*\*?\s*[A-Za-z])', ' + -', s)
    s = re.sub(r'\+\s*-\s*', '+ -', s)           # '+ - 1.0*X' -> '+ -1.0*X'
    out = []
    for tok in s.split('+'):
        tok = tok.strip()
        tok = re.sub(r'^-\s+', '-', tok)           # '- 1.0*X'  -> '-1.0*X'
        if not tok:
            continue
        m = re.match(r'^(' + NUM + r')\s*\*\s*([A-Za-z][\w]*)$', tok)
        if m:
            out.append((float(m.group(1)), m.group(2)))
            continue
        m = re.match(r'^(' + NUM + r')\s+([A-Za-z][\w]*)$', tok)
        if m:
            out.append((float(m.group(1)), m.group(2)))
            continue
        m = re.match(r'^([A-Za-z][\w]*)$', tok)
        if m:
            out.append((1.0, m.group(1)))
            continue
        raise ValueError('cannot parse term %r in %r' % (tok, s))
    return out


def read_def(path):
    """Return list[Rxn] for every statement in the REACTIONS block."""
    txt = open(path).read().replace('\r', '')
    lines = []
    for ln in txt.split('\n'):
        ln = ln.split('!')[0].rstrip()          # CMAQ comments start with '!'
        lines.append(ln)
    body = []
    inside = False
    for ln in lines:
        if ln.strip().startswith('REACTIONS'):
            inside = True
            continue
        if ln.strip().startswith('END MECH'):
            break
        if inside:
            body.append(ln)
    stmts, cur = [], ''
    for ln in body:
        if not ln.strip():
            continue
        cur += ' ' + ln.strip()
        if cur.rstrip().endswith(';'):
            stmts.append(cur.strip())
            cur = ''
    rxns = []
    for st in stmts:
        m = re.match(r'^<\s*([^>]+?)\s*>\s*(.*?)\s*#\s*(.*?)\s*;$', st)
        if not m:
            raise ValueError('cannot parse statement: ' + st)
        label, eq, rate = m.groups()
        ptype = ''
        mp = re.search(r'\s(%\w+)\s*$', eq)
        if mp:
            ptype = mp.group(1)
            eq = eq[:mp.start()]
        lhs, rhs = eq.split('=', 1)
        rxns.append(Rxn(label, _parse_side(lhs), _parse_side(rhs), ptype,
                        rate.strip(), re.sub(r'\s+', ' ', st)))
    return rxns


_TERM = re.compile(r'^\s*(' + NUM + r')\s*(?:\^\s*(' + NUM + r'))?\s*(?:@\s*(' + NUM + r'))?\s*$')


def classify_rate(r):
    """Fill r.kind / r.terms / r.scale / r.target from r.rate_raw."""
    s = r.rate_raw
    m = re.match(r'^(' + NUM + r')\s*/\s*<\s*([^>]+?)\s*>$', s)
    if m:
        r.kind, r.scale, r.target = 'photo', float(m.group(1)), m.group(2)
        return r
    m = re.match(r'^(' + NUM + r')\s*\*\s*K\s*<\s*([^>]+?)\s*>$', s)
    if m:
        r.kind, r.scale, r.target = 'ref', float(m.group(1)), m.group(2)
        return r
    m = re.match(r'^(' + NUM + r')\s*~\s*<\s*([^>]+?)\s*>$', s)
    if m:
        r.kind, r.scale, r.target = 'hetero', float(m.group(1)), m.group(2)
        return r
    m = re.match(r'^(' + NUM + r')\s*\?\s*<\s*([^>]+?)\s*>$', s)
    if m:
        r.kind, r.scale, r.target = 'op', float(m.group(1)), m.group(2)
        return r
    parts = [p for p in s.split('&')]
    terms = []
    for p in parts:
        mt = _TERM.match(p)
        if not mt:
            r.kind = 'unknown'
            return r
        A = float(mt.group(1).replace('d', 'e').replace('D', 'e'))
        n = float(mt.group(2)) if mt.group(2) else 0.0
        E = float(mt.group(3)) if mt.group(3) else 0.0
        terms.append((A, n, E))
    r.terms = terms
    if r.ptype:
        r.kind = 'pct'
    elif len(terms) == 1:
        r.kind = 'arr'
    elif len(terms) == 4:
        r.kind = 'falloff'     # k0, ki, F, n   (F and n are bare numbers: terms[2][0], terms[3][0])
    else:
        r.kind = 'unknown'
    return r


def read_classified(path):
    return [classify_rate(r) for r in read_def(path)]


if __name__ == '__main__':
    import sys, collections
    rx = read_classified(sys.argv[1])
    print(len(rx), 'statements')
    c = collections.Counter((re.match(r'[A-Za-z]+', r.label).group(0), r.kind) for r in rx)
    for k, v in sorted(c.items()):
        print(k, v)
