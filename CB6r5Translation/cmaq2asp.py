"""Translate the gas-phase reactions of a CMAQ carbon-bond mechanism
(mech_cb6r5_ae7_aq.def) into ASP's GasChemicalMechanism.in format.

What ASP can and cannot express was established by reading the ASP source
(Chemistry/ChemistryParametersInitRoutines.h, EulerianEnvironment/GasPhaseChemistry.h).
`asp_rate()` below is a line-by-line port of ResetGasPhaseReactionRates() *as written*
(including its quirks) so that every translated reaction can be checked numerically
against the CMAQ rate expression it came from.

Design rules (see CB6r5_README.md for the full discussion):
  * only the core gas-phase block 'R<n>' is translated;  CL*, DMS*, HET_*, SA*, POA*, ... are skipped
  * CMAQ SOA bookkeeping counters are dropped from product lists
  * reactions with no products get the sink product 'M' (ASP rejects an empty product list)
  * '*K<Rn>' rate references (not supported by ASP) are expanded statically and exactly
  * rate forms ASP cannot reproduce are either mapped to the closest ASP form (and reported
    with the size of the error) or replaced by a fixed value at the case T,P (and reported)
"""
import math
import re
from collections import OrderedDict

from cmaq_def import read_classified

# ----------------------------------------------------------------------------------------
# configuration
# ----------------------------------------------------------------------------------------
CORE = re.compile(r'^R\d+[a-z]*$')

# CMAQ name -> ASP name
RENAME = {'SULF': 'H2SO4',   # ASP's condensation / aqueous decks call gas sulfuric acid H2SO4
          'ECH4': 'CH4',     # (R216 is dropped, see SKIP_LABELS)
          'HCL': 'HCl'}

# CMAQ SOA bookkeeping 'counter' species: one counter is produced per oxidation event of
# the parent VOC; CMAQ then converts counters to semivolatile gases with its own yields in
# the SA*/BL* add-on reactions (not translated).  ASP's baseline decks use the same idea
# (TOLSOA, XYLSOA, ISOSOA, TERPSOA ... counters feeding the IVOC3-6 bins, see the end of
# GasChemicalMechanism.in), so by default the CB6r5 counters are re-keyed onto ASP's
# existing counters and the ASP yield chains are kept unchanged ('asp' mode).  Counters with
# no ASP analogue (benzene, PAH, SO2) are dropped.  'drop' mode removes all counters.
COUNTER_TO_ASP = {'TOLRO2': 'TOLSOA', 'XYLRO2': 'XYLSOA', 'ISOPRXN': 'ISOSOA',
                  'TRPRXN': 'TERPSOA', 'TERPNRO2': 'TERPSOA'}
DROP_ALWAYS = {'BENZRO2', 'PAHRO2', 'SULRXN'}
DROP_SPECIES = set(COUNTER_TO_ASP) | DROP_ALWAYS     # (kept for the notes below)
SOA_MODE = {'mode': 'asp'}

# CMAQ carries CH4 (fixed at ATM_CH4) and ECH4 (emitted, evolving) as two species with the
# same OH reaction.  ASP has a single evolving CH4, so only R124 is kept.
SKIP_LABELS = {'R216'}

ENV = {'M'}                 # pseudo-reactants that are not species

# boltzmann for M = P/(kT)
KB = 1.380649e-23

O2_FRAC = 0.2095            # CMAQ ATM_O2 (ppm/1e6)


def Mair(T, P_mbar):
    """molecules/cm3"""
    return P_mbar * 100.0 / (KB * T) * 1e-6


# ----------------------------------------------------------------------------------------
# CMAQ rate evaluation (from the CMAQ documentation in mech_cb6r5_ae7_aq.md)
# ----------------------------------------------------------------------------------------
def _arr(term, T):
    A, n, E = term
    return A * (T / 300.0) ** n * math.exp(-E / T)


def cmaq_rate(r, T, M, rxn_by_label, H2O=0.0):
    """k in the CMAQ convention for reaction r (cm3/molec/s etc.; no concentration factors
    for O2/H2O/M pseudo-reactants - the caller adds those)."""
    if r.kind == 'arr':
        return _arr(r.terms[0], T)
    if r.kind == 'ref':
        return r.scale * cmaq_rate(rxn_by_label[r.target], T, M, rxn_by_label, H2O)
    if r.kind == 'falloff':
        k0 = _arr(r.terms[0], T) * M
        ki = _arr(r.terms[1], T)
        F = r.terms[2][0]
        n = r.terms[3][0]
        G = 1.0 / (1.0 + (math.log10(k0 / ki) / n) ** 2)
        return (k0 / (1.0 + k0 / ki)) * F ** G
    if r.kind == 'pct':
        if r.ptype == '%3':          # k = k1 + k2 M
            return _arr(r.terms[0], T) + _arr(r.terms[1], T) * M
        if r.ptype == '%2':          # k = k0 + k3 M / (1 + k3 M / k2)   (order: k0, k2, k3)
            k0, k2, k3 = (_arr(t, T) for t in r.terms)
            return k0 + k3 * M / (1.0 + k3 * M / k2)
    raise ValueError('cannot evaluate %s (%s)' % (r.label, r.kind))


# ----------------------------------------------------------------------------------------
# port of ASP's ResetGasPhaseReactionRates  (GasPhaseChemistry.h), quirks included
#   params = the ';'-separated numbers in the deck line, exactly as typed
#   RR(L,3) = params[0], RR(L,4) = params[1], ...
# ----------------------------------------------------------------------------------------
def asp_rate(nbodies, params, T, M, C_H2O, P_mbar):
    RR = {i + 3: v for i, v in enumerate(params)}
    g = lambda i: RR.get(i, 0.0)
    f = int(g(3))
    if nbodies == 1:
        if f == 2:
            return g(4) * T ** g(5) * math.exp(g(6) / T) ** g(7)
        if f == 4:
            k0 = g(4) * (T / 300.) ** (-g(5)) * M
            ki = g(6) * (T / 300.) ** (-g(7))
            kf = (k0 * ki) / (k0 + ki) * g(8) ** (1. / (1. + math.log10(k0 / ki) ** 2))
            return kf / (g(9) * math.exp(g(10) / T))
        if f == 5:
            return g(4)
        if f == 7:
            k0 = M * 1.3e-3 * (300.0 / T) ** 3.5 * math.exp(-11000.0 / T)
            ki = 9.7e14 * (T / 300.0) ** 0.1 * math.exp(-11080.0 / T)
            ncd = 0.75 - 1.27 * math.log10(0.35)
            fd = 10. ** (math.log10(0.35) / (1. + (math.log10(k0 / ki) / ncd) ** 2))
            return fd * k0 * ki / (k0 + ki)
        if f == 8:
            return M * (0.21 * 3.2e-11 * math.exp(67. / T) + 0.79 * 2.0e-11 * math.exp(130. / T))
        if f == 9:
            n2, o2 = 0.79 * M, 0.21 * M
            return (5.6e-34 * (T / 300.) ** -2.6 * n2 + 6.0e-34 * (T / 300.) ** -2.6 * o2) * o2
        if f == 11:
            k0 = g(4) * math.exp(g(5) / T) * M
            ki = g(6) * math.exp(g(7) / T)
            ncd = 0.75 - 1.27 * math.log10(g(8))
            fd = 10. ** (math.log10(g(8)) / (1. + (math.log10(k0 / ki) / ncd) ** 2))
            return fd * k0 * ki / (k0 + ki)
        raise ValueError('1-body flag %d not ported' % f)
    if nbodies == 2:
        if f == 0:
            return g(3) * T ** g(4) * math.exp(g(5) / T) ** g(6)
        k1 = g(4) * math.exp(g(5) / T)
        k2 = g(6) * math.exp(g(7) / T)
        if f == 4:
            return k1 + k2 * M
        if f == 5:   # GasPhaseChemistry.h after common_asp 508ba8d:  k1 + k3*M/(1 + k3*M/k2)
            k3 = g(8) * math.exp(g(9) / T)
            return k1 + k3 * M / (1.0 + k3 * M / k2)
        if f == 10:
            return C_H2O * (k1 + M * k2)
        raise ValueError('2-body flag %d not ported' % f)
    if nbodies == 3:
        if f == 1:
            return g(4) * (T / 300.) ** (-g(5)) * M
        if f == 2:   # thermal third order, after common_asp 3-body flag-2 fix: M*A*T^B*exp(C/T)^D
            return M * g(4) * T ** g(5) * math.exp(g(6) / T) ** g(7)
        if f == 3:   # IUPAC
            k0 = g(4) * (T / 300.) ** (-g(5)) * M
            ki = g(6) * (T / 300.) ** (-g(7))
            ncd = 0.75 - 1.27 * math.log10(g(8))
            fd = 10. ** (math.log10(g(8)) / (1. + (math.log10(k0 / ki) / ncd) ** 2))
            return fd * k0 * ki / (k0 + ki)
        if f == 0:   # JPL, flag-less: first term is k0 itself
            k0 = g(3) * (T / 300.) ** (-g(4)) * M
            ki = g(5) * (T / 300.) ** (-g(6))
            return (k0 * ki) / (k0 + ki) * g(7) ** (1. / (1. + math.log10(k0 / ki) ** 2))
        raise ValueError('3-body flag %d not ported' % f)
    raise ValueError


# ----------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------
def fmt(x):
    return ('%.10g' % x)


def ncd(F):
    return 0.75 - 1.27 * math.log10(F)


def side_str(items):
    out = []
    for c, n in items:
        out.append(n if abs(c - 1.0) < 1e-12 else '%s %s' % (fmt(c), n))
    return ' + '.join(out)


def clean_products(prods):
    """rename, re-key/drop counters, merge duplicates, drop zeros -> [(coef,name)]"""
    d = OrderedDict()
    for c, n in prods:
        if n in DROP_ALWAYS:
            continue
        if n in COUNTER_TO_ASP:
            if SOA_MODE['mode'] == 'drop':
                continue
            n = COUNTER_TO_ASP[n]
        n = RENAME.get(n, n)
        d[n] = d.get(n, 0.0) + c
    return [(c, n) for n, c in d.items() if abs(c) > 1e-12]


def clean_reactants(reacts):
    return [RENAME.get(n, n) for _, n in reacts]


# ----------------------------------------------------------------------------------------
# photolysis columns:  CMAQ label (+ scale)  ->  PhotoRates.in column of the CB6r5 table
# The table itself is derived from the existing 51-column TUV table of each case by
# derive_cb6r5_photorates.py using PHOTO_COMBOS below.
# ----------------------------------------------------------------------------------------
# (old TUV column index, weight).  Old column numbers are those of ASP's standard
# PhotoRates.in (see PhotoRatesTUV/asp_photorates_map.py).
PHOTO_COMBOS = OrderedDict([
    ('NO2_IUPAC10',     ([(6, 1.0)],            'NO2 -> NO + O(3P)')),
    ('O3_O3P_IUPAC10',  ([(3, 1.0)],            'O3 -> O2 + O(3P)')),
    ('O3_O1D_IUPAC10',  ([(2, 1.0)],            'O3 -> O2 + O(1D)')),
    ('H2O2_IUPAC10',    ([(5, 1.0)],            'H2O2 -> 2 OH')),
    ('NO3NO2_06',       ([(8, 1.0)],            'NO3 -> NO2 + O(3P)')),
    ('NO3NO_06',        ([(7, 1.0)],            'NO3 -> NO + O2')),
    ('N2O5_IUPAC10',    ([(9, 1.0), (10, 1.0)], 'N2O5 -> NO3 + NO2 (+ NO3 + NO + O): sum of both TUV channels')),
    ('HONO_IUPAC10',    ([(12, 1.0)],           'HNO2 -> OH + NO')),
    ('HNO3_IUPAC10',    ([(13, 1.0)],           'HNO3 -> OH + NO2')),
    ('PNA_IUPAC10',     ([(14, 1.0)],           'HNO4 -> HO2 + NO2')),
    ('PAN_IUPAC10',     ([(40, 1.0), (41, 1.0)], 'PAN: sum of the NO2 and NO3 channels (CB6r5 splits 0.6/0.4 in the products)')),
    ('MEPX_IUPAC10',    ([(31, 1.0)],           'CH3OOH -> CH3O + OH')),
    ('NTR_IUPAC10',     ([(36, 1.0)],           'isopropyl nitrate -> products + NO2')),
    ('FORM_R_IUPAC13',  ([(17, 1.0)],           'CH2O -> H + HCO')),
    ('FORM_M_IUPAC13',  ([(18, 1.0)],           'CH2O -> H2 + CO')),
    ('ALD2_R_IUPAC13',  ([(19, 1.0)],           'CH3CHO -> CH3 + HCO')),
    ('ALDX_R_IUPAC13',  ([(22, 1.0)],           'C2H5CHO -> C2H5 + HCO')),
    ('GLYD_IUPAC13',    ([(23, 1.0)],           'HOCH2CHO: sum of the three TUV channels')),
    ('GLY_R_IUPAC13',   ([(44, 1.0)],           'CHOCHO -> HCO + HCO')),
    ('MGLY_IUPAC10',    ([(46, 1.0)],           'CH3COCHO -> CH3CO + HCO')),
    ('KET_IUPAC10',     ([(28, 1.0)],           'MEK -> CH3CO + C2H5')),
    ('ACET_IUPAC10',    ([(26, 1.0)],           'acetone -> CH3CO + CH3')),
    ('ISPD',            ([(24, 1.0)],           'acrolein -> products (CB6 ISPD uses acrolein cross sections)')),
    ('HPALD',           ([(51, 1.0)],           'ASP convention: 100 x J(MACR), table column 51')),
])


class PhotoColumns:
    """assign a CB6r5 PhotoRates column to every distinct (CMAQ label, scale)"""

    def __init__(self):
        self.cols = OrderedDict()      # (label, scale) -> index (1-based)

    def index(self, label, scale):
        if label not in PHOTO_COMBOS:
            raise KeyError('no TUV mapping for CMAQ photolysis label %s' % label)
        key = (label, round(scale, 12))
        if key not in self.cols:
            self.cols[key] = len(self.cols) + 1
            if len(self.cols) > 51:
                raise ValueError('more than 51 photolysis columns needed')
        return self.cols[key]


# ----------------------------------------------------------------------------------------
# the translation proper
# ----------------------------------------------------------------------------------------
class Translation:
    def __init__(self, def_path, T, P_mbar, soa_counters='asp'):
        SOA_MODE['mode'] = soa_counters
        self.T0, self.P0 = T, P_mbar
        self.all = read_classified(def_path)
        self.by_label = {r.label: r for r in self.all}
        self.core = [r for r in self.all if CORE.match(r.label) and r.label not in SKIP_LABELS]
        self.skipped = [r for r in self.all if r not in self.core]
        self.photo = PhotoColumns()
        self.lines = []          # (label, cmaq_text, asp_line, nbodies, params)
        self.notes = []          # dicts: label, level, text
        self.species = OrderedDict()
        for r in self.core:
            self._one(r)

    # -- bookkeeping -------------------------------------------------------------------
    def note(self, r, level, text):
        self.notes.append({'label': r.label, 'level': level, 'text': text, 'cmaq': r.text})

    def _emit(self, r, reacts, prods, nb, params, comment=None, add_M=False):
        rs = list(reacts)
        ps = list(prods)
        if add_M:
            rs.append('M')
        if not ps:
            ps = [(1.0, 'M')]
        elif add_M:
            ps = ps + [(1.0, 'M')]
        line = '%s => %s ; %s' % (' + '.join(rs), side_str(ps), '; '.join(fmt(p) for p in params) + ';')
        self.lines.append((r.label, r.text, line, nb, list(params), comment))
        for n in reacts:
            if n not in ENV and n not in ('H2O',):
                self.species.setdefault(n, None)
        for c, n in ps:
            if n not in ENV and n not in ('H2O',):
                self.species.setdefault(n, None)

    # -- one reaction ------------------------------------------------------------------
    def _one(self, r):
        T0, P0 = self.T0, self.P0
        reacts = clean_reactants(r.reactants)
        prods = clean_products(r.products)
        kept = [n for _, n in r.products if n in COUNTER_TO_ASP and SOA_MODE['mode'] == 'asp']
        gone = [n for _, n in r.products if n in DROP_ALWAYS or
                (n in COUNTER_TO_ASP and SOA_MODE['mode'] == 'drop')]
        if kept:
            self.note(r, 'soa', 'CMAQ SOA counter re-keyed: ' +
                      ', '.join('%s->%s' % (n, COUNTER_TO_ASP[n]) for n in kept))
        if gone:
            self.note(r, 'soa', 'CMAQ SOA counter dropped (no ASP analogue): ' + ', '.join(gone))
        body = [x for x in reacts if x not in ENV]
        has_H2O = 'H2O' in body
        # --- references to other rates: expand statically ---------------------------------
        if r.kind == 'ref':
            base, fac = r, 1.0
            chain = []
            while base.kind == 'ref':
                fac *= base.scale
                chain.append(base.target)
                base = self.by_label[base.target]
            self.note(r, 'expanded',
                      'ASP cannot reference another reaction\'s rate; expanded exactly as %s x %s'
                      % (fmt(fac), '->'.join(chain)))
            self._from_base(r, reacts, prods, base, fac)
            return
        self._from_base(r, reacts, prods, r, 1.0)

    def _from_base(self, r, reacts, prods, base, fac):
        """emit reaction r using the rate expression of `base` multiplied by `fac`"""
        T0, P0 = self.T0, self.P0
        lab = r.label
        body = [x for x in reacts if x not in ENV]
        nM = reacts.count('M')

        # ---- photolysis -----------------------------------------------------------------
        if base.kind == 'photo':
            assert len(body) == 1, lab
            col = self.photo.index(base.target, base.scale * fac)
            self._emit(r, body, prods, 1, [0.0, float(col)])
            return

        if base.kind in ('hetero', 'op'):
            self.note(r, 'unsupported', 'rate type %s<%s> not translated' % (base.kind, base.target))
            return

        # ---- special reactions (pseudo-reactants M / O2 / H2O) ---------------------------
        if lab == 'R2':                                    # O + O2 + M = O3 + M
            self._emit(r, ['O'], prods, 1, [9.0])
            self.note(r, 'approx', 'ASP hard-wires O + O2 + M (flag 9): 5.6e-34 N2 / 6.0e-34 O2 instead of CB6r5 6.0e-34 for M')
            return
        if lab == 'R10':                                   # O1D + M = O + M
            self._emit(r, ['O1D'], prods, 1, [8.0])
            self.note(r, 'approx', 'ASP hard-wires O1D + M (flag 8, IUPAC N2/O2 form) instead of CB6r5 2.23e-11 exp(115/T)')
            return
        if lab == 'R24':                                   # NO + NO + O2 = 2 NO2  (termolecular, O2 reactant)
            A, n, E = base.terms[0]
            # rate = k [NO]^2 [O2],  [O2] = O2_FRAC*M  ->  3-body flag 2:  M * A' * T^n * exp(C/T)
            self._emit(r, ['NO', 'NO'], prods, 3, [2.0, A * fac * O2_FRAC * 300.0 ** (-n), n, -E, 1.0], add_M=True)
            self.note(r, 'info', 'termolecular NO+NO+O2 -> 3-body rate flag 2 with [O2] = %.4f [M]; flag 2 is correct only '
                      'in common_asp >= the flag-2 fix commit (earlier versions mis-indexed its parameters)' % O2_FRAC)
            return
        if lab == 'R46':                                   # HNO3 + OH, %2 pressure form -> ASP flag 5
            (A0, n0, E0), (A2, n2, E2), (A3, n3, E3) = base.terms
            assert n0 == n2 == n3 == 0.0
            self._emit(r, body, prods, 2, [5.0, A0 * fac, -E0, A2 * fac, -E2, A3 * fac, -E3])
            self.note(r, 'info',
                      'CB6r5 %2 form k0 + k3[M]/(1+k3[M]/k2) = ASP 2-body rate flag 5, which is correct only in '
                      'common_asp >= 508ba8d (earlier versions evaluated k1 + k3*M*(1 + k3*M/k2), ~240x too fast)')
            return
        if has_H2O_third(body):                            # R41, R221 (and R20): pseudo-2nd order with [H2O]
            nonw = [x for x in body if x != 'H2O']
            if base.kind == 'pct' and base.ptype == '%3':
                (A1, n1, E1), (A2, n2, E2) = base.terms
                p = [10.0, A1 * fac, -E1, A2 * fac, -E2]
            else:
                A, n, E = base.terms[0]
                p = [10.0, A * fac * 300.0 ** (-n), -E, 0.0, 0.0]
                assert n == 0.0
            self._emit(r, nonw, prods, 2, p)
            return

        # ---- generic forms ----------------------------------------------------------------
        nb = len(body)
        if base.kind == 'arr':
            A, n, E = base.terms[0]
            A = A * fac
            Ap = A * 300.0 ** (-n)
            if nb == 1:
                self._emit(r, body, prods, 1, [2.0, Ap, n, -E, 1.0])
            elif nb == 2:
                assert Ap < 1.0, (lab, Ap)
                self._emit(r, body, prods, 2, [Ap, n, -E, 1.0])
            else:
                raise ValueError('%s: %d reactants not supported' % (lab, nb))
            return

        if base.kind == 'pct' and base.ptype == '%3' and nb == 2:
            (A1, n1, E1), (A2, n2, E2) = base.terms
            self._emit(r, body, prods, 2, [4.0, A1 * fac, -E1, A2 * fac, -E2])
            return

        if base.kind == 'falloff':
            (A0, n0, E0), (Ai, ni, Ei) = base.terms[0], base.terms[1]
            F, nF = base.terms[2][0], base.terms[3][0]
            iupac = abs(nF - ncd(F)) < 0.02
            jpl = abs(nF - 1.0) < 1e-9
            if nb == 2:
                assert E0 == 0.0 and Ei == 0.0, lab
                if iupac:
                    self._emit(r, body, prods, 3, [3.0, A0 * fac, -n0, Ai * fac, -ni, F], add_M=True)
                    if abs(nF - ncd(F)) > 0.003:
                        self.note(r, 'info', 'CMAQ n=%.2f vs IUPAC N(F)=%.4f used by ASP' % (nF, ncd(F)))
                elif jpl:
                    self._emit(r, body, prods, 3, [A0 * fac, -n0, Ai * fac, -ni, F], add_M=True)
                else:
                    raise ValueError('%s: falloff n=%g F=%g fits neither IUPAC nor JPL form' % (lab, nF, F))
                return
            if nb == 1:
                pure_arr = (n0 == 0.0 and ni == 0.0)
                # (a) N2O5 thermal decomposition: ASP's hard-wired flag 7 is the same IUPAC expression
                if (abs(A0 - 1.3e-3) < 1e-12 and n0 == -3.5 and E0 == 11000.0 and abs(Ai - 9.7e14) < 1 and
                        ni == 0.1 and Ei == 11080.0 and abs(F - 0.35) < 1e-9):
                    self._emit(r, body, prods, 1, [7.0])
                    return
                # (b) pure-Arrhenius IUPAC decomposition (PNA, PAN, PANX, OPAN): ASP flag 11
                if iupac and pure_arr:
                    self._emit(r, body, prods, 1, [11.0, A0 * fac, -E0, Ai * fac, -Ei, F])
                    return
                # (c) JPL-form (n=1) with power-law T dependence: ASP flag 4 with Keq = 1
                if jpl and E0 == 0.0 and Ei == 0.0:
                    self._emit(r, body, prods, 1, [4.0, A0 * fac, -n0, Ai * fac, -ni, F, 1.0, 0.0])
                    return
                self.note(r, 'unsupported', 'unimolecular fall-off form not representable in ASP')
                return
        self.note(r, 'unsupported', 'rate form %s %s not handled' % (base.kind, base.ptype))


def has_H2O_third(body):
    return 'H2O' in body and len(body) == 3


# ----------------------------------------------------------------------------------------
# numerical comparison of every translated reaction with its CMAQ source
# ----------------------------------------------------------------------------------------
TEST_POINTS = [(298.15, 1013.25), (288.4, 880.0), (280.6, 846.5), (290.0, 730.0),
               (288.0, 900.0), (260.0, 600.0), (310.0, 1013.25)]
C_H2O_TEST = 3e17          # molecules/cm3


def compare(tr):
    """returns list of dict(label, maxrel, worst_point) for non-photolysis reactions"""
    out = []
    for lab, text, line, nb, params, _ in tr.lines:
        r = tr.by_label[lab]
        if nb == 1 and int(params[0]) == 0 and len(params) == 2:
            continue                                    # photolysis
        base = r
        fac = 1.0
        while base.kind == 'ref':
            fac *= base.scale
            base = tr.by_label[base.target]
        worst, wp = 0.0, None
        for T, P in TEST_POINTS:
            M = Mair(T, P)
            ref = cmaq_rate(base, T, M, tr.by_label) * fac
            # concentration factors of CMAQ pseudo-reactants
            if lab == 'R2':
                ref = ref * M * (O2_FRAC * M)
            elif lab == 'R10':
                ref = ref * M
            elif lab == 'R24':
                ref = ref * (O2_FRAC * M)
            elif lab in ('R20', 'R41', 'R221'):
                ref = ref * C_H2O_TEST
            elif lab == 'R39':
                pass
            val = asp_rate(nb, params, T, M, C_H2O_TEST, P)
            rel = abs(val - ref) / abs(ref)
            if rel > worst:
                worst, wp = rel, (T, P)
        out.append({'label': lab, 'maxrel': worst, 'at': wp, 'asp_line': line})
    return out


if __name__ == '__main__':
    import sys
    tr = Translation(sys.argv[1], 288.4, 880.0)
    print('%d core reactions -> %d ASP lines' % (len(tr.core), len(tr.lines)))
    for c in compare(tr):
        if c['maxrel'] > 1e-6:
            print('%-6s max rel. deviation %.3g at T=%.1f P=%.1f   %s' % (c['label'], c['maxrel'], c['at'][0], c['at'][1], c['asp_line'][:90]))
    for n in tr.notes:
        if n['level'] != 'info':
            print(n['level'].upper().ljust(11), n['label'].ljust(6), n['text'][:150])
