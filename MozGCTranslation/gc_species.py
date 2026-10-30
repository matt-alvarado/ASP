"""GEOS-Chem fullchem species table: formula -> molecular weight, exclusion (halogen / aerosol-phase / dummy), peroxy classification."""
import re
from gc_kpp import GCDIR
EL = dict(H=1.008, C=12.011, N=14.007, O=15.999, S=32.06, F=18.998, Cl=35.45, Br=79.904, I=126.904)

def parse_formula(f):
    f = f.replace('-', '').replace('=', '').replace(' ', '')
    if not re.fullmatch(r'[A-Za-z0-9()]+', f) or not re.search(r'[A-Z]', f): return None
    def rec(s):
        tot = {}; i = 0
        while i < len(s):
            if s[i] == '(':
                d = 1; j = i + 1
                while d:
                    d += (s[j] == '(') - (s[j] == ')'); j += 1
                sub = rec(s[i + 1:j - 1]); m = re.match(r'\d+', s[j:]); n = int(m.group()) if m else 1
                j += len(m.group()) if m else 0
                for k, v in sub.items(): tot[k] = tot.get(k, 0) + v * n
                i = j
            else:
                m = re.match(r'(Cl|Br|[A-Z])(\d*)', s[i:])
                if not m or m.group(1) not in EL: raise ValueError(f)
                tot[m.group(1)] = tot.get(m.group(1), 0) + int(m.group(2) or 1); i += len(m.group())
        return tot
    try: return rec(f)
    except ValueError: return None

def read_species(path=GCDIR + 'fullchem.eqn'):
    txt = open(path).read(); sec = txt.split('#DEFVAR')[1].split('#DEFFIX')[0] + '\n' + txt.split('#DEFFIX')[1].split('#EQUATIONS')[0]
    out = {}
    for m in re.finditer(r'(\w+)\s*=\s*IGNORE;(?:\s*\{([^}]*)\})?', sec):
        out[m.group(1)] = ' '.join((m.group(2) or '').split())
    return out

HALO = re.compile(r'Cl|Br|I|F')
def info(name, com):
    first = com.split(';')[0].strip()
    comp = parse_formula(first)
    mw = sum(EL[e] * n for e, n in comp.items()) if comp else None
    return comp, mw

# molecular weights for lumped / formula-less species (g/mol; approximate)
MW_OVER = dict(ACR=56, BUTDI=84, BUTN=117, IEPOXAOO=149, IEPOXBOO=149, PHEN=94, R7P=118, 
    A3O2=75, ACRO2=105, ALK4=58.1, ALK4N1=147, ALK4N2=119, ALK4O2=89, ALK6=86, AROMRO2=150, AROMP4=100, AROMP5=114, AROMPN=190,
    AROMCO3=165, ATO2=89, B3O2=75, BALD=106, BENZO=93, BENZO2=125, BENZP=126, BRO2=125, BUTO2=119, BZCO3=121, BZCO3H=138, BZPAN=183,
    C4HVP1=117, C4HVP2=117, CSL=108, ETHN=107, ETHP=78, ETO=61, ETOO=77, FURA=68, GCO3=91, HACTA=76, HONIT=215, HPALD1OO=147, HPALD2OO=147,
    ICHOO=133, ICNOO=178, IDHNBOO=194, IDHNDOO1=194, IDHNDOO2=194, IDNOO=210, IHOO1=117, IHOO4=117, IHPNBOO=194, IHPNDOO=194,
    IHPOO1=149, IHPOO2=149, IHPOO3=149, INA=163, INO2B=162, INO2D=162, ISOPNOO1=178, ISOPNOO2=178, KO2=103, LIMO2=185, MACR1OO=101,
    MACRNO2=165, MCROHOO=119, MEKCO3=101, MCO3=75, MDL=48, MONITS=215, MONITU=215, MTPA=136, MTPO=136, MVKOHOO=119, NRO2=190, OLND=250,
    OLNN=250, OTHRO2=61, PIO2=185, PIP=186, PO2=91, PRN1=137, R4N1=148, R4N2=133, R4O2=89, R7N1=162, R7N2=147, R7O2=103, RCO3=89,
    RCHO=58, RCOOH=74, RNO3=227, ROH=60, TLFUO2=159, TRO2=141, XRO2=100, ZRO2=130, APINP=186, APINN=215, PINAL=168, PINPAN=231,
    PINONIC=184, PINO3H=200, C96O2H=202, C96N=231, BPINO=152, BPINN=215, BPINP=186, BPINOOH=186, BPINON=184, LIMAL=168, LIMN=215,
    LIMKET=170, LIMKB=170, LIMNB=215, LIMPAN=231, LIMO2H=186, LIMO3H=200, MYRCO=152, PIN=136, APINO2=185, PINO3=183, C96O2=199,
    BPINO2=185, BPINOO2=201, LIMKO2=201, LIMO3=199, NPHEN=139, MCT=124, MEK=72, MEKPN=163, MCRENOL=86, MVKHC=102, MVKHCB=102,
    IDC=98, LVOC=166, HMS=111, PRPE=42.1, C2H4=28.05, STYR=104, TMB=120, EBZ=106, NAP=128, TOLU=92.1, XYLE=106.2, CH4=16.04,
    O=16, O1D=16, H=1.008, N=14.007, H2=2.016, N2=28.01, O2=32.0)
def build():
    S = read_species(); T = {}
    for n, com in S.items():
        comp, mw = info(n, com)
        mw = MW_OVER.get(n, mw)
        T[n] = dict(com=com, comp=comp, mw=mw)
    return T
if __name__ == '__main__':
    T = build(); bad = [n for n, v in T.items() if v['mw'] is None]
    print(len(T), 'species; no MW:', bad)
