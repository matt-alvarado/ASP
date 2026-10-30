"""Shared ASP reaction-line formatting (same conventions as MCMTranslation/mcm2asp.py)."""
def fmt(x):
    if abs(x) < 1e-300: return '0'
    return '%.8g' % x

def line(re_, pr, flag, params):
    def side(s):
        return ' + '.join(n if abs(c - 1) < 1e-9 else '%s %s' % (fmt(c), n) for c, n in s)
    head = '%s => %s' % (side(re_), side(pr))
    ps = '; '.join(fmt(p) for p in params)
    fl = ('%s;' % flag) if flag is not None else ''
    return head + ' ; ' + fl + (' ' if fl and ps else '') + (ps + ';' if ps else '')
