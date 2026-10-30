"""Symbolic rate-expression engine for MCM v3.3.1 KPP rate strings.

Every rate is reduced to a list of Terms  c * T^B * exp(C/T) * [O2]^a [N2]^b [H2O]^c M^d RO2^e * (J | fall-off symbol)
so that the ASP deck keeps all temperature / pressure dependence in the rate flags (nothing is evaluated at a fixed T or P).
"""
import ast, math, re
from mcm_kpp import D

class Term:
    __slots__ = ('c', 'B', 'C', 'o2', 'n2', 'h2o', 'm', 'ro2', 'J', 'sym', 'e3')
    def __init__(s, c=1.0, B=0.0, C=0.0, o2=0, n2=0, h2o=0, m=0, ro2=0, J=None, sym=None, e3=0.0):
        s.c, s.B, s.C, s.o2, s.n2, s.h2o, s.m, s.ro2, s.J, s.sym, s.e3 = c, B, C, o2, n2, h2o, m, ro2, J, sym, e3
    def __mul__(s, o):
        if (s.J and o.J) or (s.sym and o.sym) or ((s.J or s.sym) and (o.J or o.sym)):
            raise ValueError('product of J/fall-off symbols')
        return Term(s.c * o.c, s.B + o.B, s.C + o.C, s.o2 + o.o2, s.n2 + o.n2, s.h2o + o.h2o, s.m + o.m,
                    s.ro2 + o.ro2, s.J or o.J, s.sym or o.sym, s.e3 + o.e3)
    def pure(s):
        return not (s.o2 or s.n2 or s.h2o or s.m or s.ro2 or s.J or s.sym or s.e3)
    def __repr__(s):
        return 'T(%g T^%g e^(%g/T) o2%d n2%d h2o%d m%d ro2%d %s %s)' % (s.c, s.B, s.C, s.o2, s.n2, s.h2o, s.m, s.ro2, s.J, s.sym)

def poly_mul(a, b):
    return [x * y for x in a for y in b]

class Eval(ast.NodeVisitor):
    def __init__(self, names):
        self.names = names            # name -> list[Term]
    def visit_Expression(self, n): return self.visit(n.body)
    def visit_Constant(self, n): return [Term(float(n.value))]
    def visit_Name(self, n):
        if n.id == 'TEMP': return [Term(1.0, 1.0)]
        if n.id == 'O2': return [Term(o2=1)]
        if n.id == 'N2': return [Term(n2=1)]
        if n.id == 'H2O': return [Term(h2o=1)]
        if n.id == 'M': return [Term(m=1)]
        if n.id == 'RO2': return [Term(ro2=1)]
        if n.id in self.names: return [Term(t.c, t.B, t.C, t.o2, t.n2, t.h2o, t.m, t.ro2, t.J, t.sym, t.e3) for t in self.names[n.id]]
        raise KeyError(n.id)
    def visit_UnaryOp(self, n):
        v = self.visit(n.operand)
        if isinstance(n.op, ast.USub): return [Term(-t.c, t.B, t.C, t.o2, t.n2, t.h2o, t.m, t.ro2, t.J, t.sym, t.e3) for t in v]
        if isinstance(n.op, ast.UAdd): return v
        raise ValueError
    def visit_BinOp(self, n):
        a, b = self.visit(n.left), self.visit(n.right)
        if isinstance(n.op, ast.Add): return a + b
        if isinstance(n.op, ast.Sub): return a + [Term(-t.c, t.B, t.C, t.o2, t.n2, t.h2o, t.m, t.ro2, t.J, t.sym, t.e3) for t in b]
        if isinstance(n.op, ast.Mult): return poly_mul(a, b)
        if isinstance(n.op, ast.Div):
            if len(b) == 1 and b[0].pure():
                return poly_mul(a, [Term(1.0 / b[0].c, -b[0].B, -b[0].C)])
            raise ValueError('division')
        if isinstance(n.op, ast.Pow):
            if len(b) == 1 and b[0].pure() and b[0].B == 0 and b[0].C == 0 and len(a) == 1 and a[0].pure():
                p = b[0].c; t = a[0]
                if t.c < 0: raise ValueError('neg base')
                return [Term(t.c ** p, t.B * p, t.C * p)]
            raise ValueError('power')
        raise ValueError
    def visit_Call(self, n):
        f = n.func.id
        if f == 'EXP':
            arg = self.visit(n.args[0])
            a = sum(t.c for t in arg if t.B == 0 and t.C == 0 and t.pure())
            iv = sum(t.c for t in arg if t.B == -1 and t.pure())
            e3 = sum(t.c for t in arg if t.B == -3 and t.pure())
            if any(not t.pure() or t.B not in (0, -1, -3) for t in arg): raise ValueError('EXP arg')
            return [Term(math.exp(a), 0.0, iv, e3=e3)]
        if f == 'J':
            return [Term(1.0, J=n.args[0].id)]
        raise ValueError(f)

def fort(s):
    s = s.replace('**', '**')
    s = re.sub(r'(?<![A-Za-z_0-9])(\d+\.?\d*)[Dd]([+-]?\d+)', r'\1E\2', s)
    s = re.sub(r'(?<![A-Za-z_0-9.])(\d+)\.(?![0-9eE])', r'\1.0', s)   # "1." -> "1.0" (valid python anyway)
    return s

def parse(s, names):
    return Eval(names).visit(ast.parse(fort(s), mode='eval'))

def load_constants():
    """Arrhenius-closed constants and fall-off parameter sets from the MCM constants module."""
    txt = open(D + 'mcm_3-3-1_kpp_constants_mcm.f90').read()
    assign = {}
    for ln in txt.split('SUBROUTINE define_constants_mcm()')[1].split('END SUBROUTINE')[0].split('\n'):
        m = re.match(r'^\s*([A-Z0-9_]+)\s*=\s*(.*?)\s*(?:!.*)?$', ln)
        if m and not m.group(1).startswith('J('): assign[m.group(1)] = m.group(2)
    names = {}
    for k in ('K14ISOM1', 'K298CH3O2', 'KAPHO2', 'KAPNO', 'KCH3O2', 'KDEC', 'KNO3AL', 'KRO2HO2', 'KRO2NO', 'KRO2NO3',
              'KROPRIM', 'KROSEC'):
        names[k] = parse(assign[k], names)
    # fall-off sets: name -> (k0 terms, ki terms, FC expression)
    fall = {}
    pairs = {'KMT01': ('K10', 'K1I', 'FC1'), 'KMT02': ('K20', 'K2I', 'FC2'), 'KMT03': ('K30', 'K3I', 'FC3'),
             'KMT04': ('K40', 'K4I', 'FC4'), 'KMT07': ('K70', 'K7I', 'FC7'), 'KMT08': ('K80', 'K8I', 'FC8'),
             'KMT09': ('K90', 'K9I', 'FC9'), 'KMT10': ('K100', 'K10I', 'FC10'), 'KMT12': ('K120', 'K12I', 'FC12'),
             'KMT13': ('K130', 'K13I', 'FC13'), 'KMT14': ('K140', 'K14I', 'FC14'), 'KMT15': ('K150', 'K15I', 'FC15'),
             'KMT16': ('K160', 'K16I', 'FC16'), 'KMT17': ('K170', 'K17I', 'FC17'),
             'KBPAN': ('KD0', 'KDI', 'FCD'), 'KBPPN': ('KPPN0', 'KPPNI', 'FCPPN'), 'KFPAN': ('KC0', 'KCI', 'FCC')}
    for nm, (k0, ki, fc) in pairs.items():
        t0 = parse(assign[k0], names); ti = parse(assign[ki], names)
        assert len(t0) == 1 and len(ti) == 1 and t0[0].m == 1 and ti[0].m == 0, nm
        fcx = assign[fc]
        if 'TEMP' in fcx:                      # only FC17: 0.17*exp(-51/T)+exp(-T/204): evaluate at 298 K (documented)
            F = 0.17 * math.exp(-51 / 298.) + math.exp(-298. / 204.)
            Tdep = (0.17, 51., 204.)          # Fc = 0.17 exp(-51/T) + exp(-T/204)
        else:
            F = float(fcx); Tdep = None
        fall[nm] = dict(k0=t0[0], ki=ti[0], F=F, Tdep=Tdep)
    return names, fall, assign

def fall_params(fs):
    """-> A0, n0, E0, Ai, ni, Ei, F, N   with  k0 = A0 (T/300)^-n0 exp(E0/T) [M],  ki = Ai (T/300)^-ni exp(Ei/T)"""
    k0, ki = fs['k0'], fs['ki']
    n0 = -k0.B; A0 = k0.c / 300.0 ** n0
    ni = -ki.B; Ai = ki.c / 300.0 ** ni
    F = fs['F']; N = 0.75 - 1.27 * math.log10(F)
    return A0, n0, k0.C, Ai, ni, ki.C, F, N
