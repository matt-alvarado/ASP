"""CMAQ CRACMM organic-aerosol chemistry for the ASP CRACMM2 / CRACMM3 decks.

CMAQ treats organic aerosol with (i) V* gas / A* particle species pairs that carry their own C*, molar mass, enthalpy of vaporization
(CCTM/src/aero/cracmm/ORG_DEFN.F, oa_list) and (ii) explicit gas-phase multigeneration chemistry of those species (ROCALK, ROCOXY,
ROCARO, TRP groups of mech_cracmm*.def).  Here the CMAQ vapor name is the ASP species name (one species, gas + hydrophobic organic
phase, as for the baseline IVOCk); ASOAT (nonvolatile SOA) is a gas species with C* = 1e-10 ug/m3.

  apply(mech, case_dir, tr_species, old_rows, T, P, mapped, table, poa)  -> list of V species rows added to `mapped`, writes
       OrgPhaseChems.in, OrganicDissolution.in, AerosolModes.in in case_dir.
Emissions (initial conditions):
  * Williams: the baseline's measured volatility distribution (gas IVOC3-7, particle IVOC3-6) is moved bin-by-bin onto the CMAQ
    oxygenated species of the same C* with the lowest O:C available (CMAQ's recommended choice for POA), mass conserved.
  * Other cases (POA1 + CBIO particle mass only): the particle mass is treated as CMAQ's PM emission and spread with CMAQ's documented
    GENERAL_POA volatility profile (CMAQ_Chem_Control_cracmm3.nml, Woody et al. 2016), mass conserved.  (CMAQ's default is nonvolatile
    APOC; the semivolatile profile is its documented recommendation.)
"""
import math
import os
import re

ORG_DEFN = '/home/claude/usepa/cmaq/CCTM/src/aero/cracmm/ORG_DEFN.F'
R = 8.314462618
# ASP IVOCk (C* = 10^(k-3) ug/m3) -> CMAQ oxygenated species of the same C*, lowest O:C available
BIN_MAP = {1: 'VROCN2OXY2', 2: 'VROCN1OXY1', 3: 'VROCP0OXY2', 4: 'VROCP1OXY1', 5: 'VROCP2OXY2', 6: 'VROCP3OXY2',
           7: 'VROCP4OXY2', 8: 'VROCP5OXY1', 9: 'VROCP6OXY1'}
# CMAQ GENERAL_POA profile: fraction of the PM emission, gas + particle species added (CMAQ_Chem_Control_cracmm3.nml)
GENERAL_POA = {'VROCN2OXY2': 0.35, 'VROCP0OXY2': 0.175 + 0.175, 'VROCP1OXY1': 0.091 + 0.009, 'VROCP2OXY2': 0.100, 'VROCP3OXY2': 0.100}
GENERAL_POA_GAS = {'VROCP0OXY2': 0.175, 'VROCP1OXY1': 0.091, 'VROCP2OXY2': 0.100, 'VROCP3OXY2': 0.100}
ASOAT_MW = 200.0
MASS_ACC = 1.0              # CMAQ: accommodation coefficient 1 (alphorg), equilibrium partitioning
RHO_ORG = 1.5


def read_oa_list():
    """ORG_DEFN.F oa_list -> dict vapor-name -> dict(pm, cstar298, enth_J, OC, OMOC, nonvol)"""
    out = {}
    for ln in open(ORG_DEFN, errors='ignore'):
        mm = re.match(r"\s*&\s*oa_type\('([^']+)'\s*,\s*'([^']*)'\s*,\s*'[^']*'\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([TF])\s*,\s*([TF])\s*,\s*([TF])\s*,\s*([TF])", ln)
        if not mm:
            continue
        pm, vap = mm.group(1).strip(), mm.group(2).strip()
        rec = dict(pm=pm, cstar=float(mm.group(4)), enth=float(mm.group(5)), OC=float(mm.group(6)), OMOC=float(mm.group(7)),
                   nonvol=mm.group(11) == 'T')
        out[vap or pm] = rec
    return out


def psat300_mbar(cstar298, enth_J, mw):
    """C*(298 K) at the CMAQ convention -> Psat(300 K) in mbar for ASP flag-1 dissolution"""
    c300 = cstar298 * (298.0 / 300.0) * math.exp(enth_J / R * (1.0 / 298.0 - 1.0 / 300.0))   # ug/m3 at 300 K
    return c300 * 1e-6 * R * 300.0 / mw / 100.0


def partitioning_species(species, table):
    """V species of `species` that have a CMAQ particle partner (plus ASOAT) -> list of (name, mw, cstar298, enth, nonvol)"""
    oa = read_oa_list()
    out = []
    for n in species:
        rec = oa.get(n)
        if n == 'ASOAT':
            out.append((n, ASOAT_MW, 1e-10, 1.0e3, True))
        elif rec is not None and n.startswith('V') and not rec['nonvol']:
            out.append((n, table[n][0], rec['cstar'], rec['enth'], False))
    return out


def _ppb_per_ugm3(mw, T, P):
    return 1.0 / (mw * (P * 100.0) / (R * T) * 1e-3)


def emissions(case, old_rows, aer_modes, table, T, P):
    """-> (gas dict V name -> ppb, particle dict V name -> ug/m3, replaced particle names, note lines)"""
    gas, part, notes = {}, {}, []
    ivoc_g = {k: old_rows['IVOC%d' % k][3] for k in range(1, 10) if 'IVOC%d' % k in old_rows}
    ivoc_mw = {k: old_rows['IVOC%d' % k][2] for k in range(1, 10) if 'IVOC%d' % k in old_rows}
    if any(v > 0 for v in ivoc_g.values()) or any(re.match(r'^IVOC\d$', n) for n in aer_modes):
        # measured volatility distribution -> same-C* CMAQ oxygenated species
        for k, ppb in ivoc_g.items():
            if ppb > 0:
                sp = BIN_MAP[k]
                gas[sp] = gas.get(sp, 0.0) + ppb * ivoc_mw[k] / table[sp][0]
                notes.append('gas IVOC%d %.4g ppb (MW %.0f) -> %s %.4g ppb (MW %.1f)' % (k, ppb, ivoc_mw[k], sp, gas[sp], table[sp][0]))
        for n, m in aer_modes.items():
            mm = re.match(r'^IVOC(\d)$', n)
            if mm:
                sp = BIN_MAP[int(mm.group(1))]
                part[sp] = part.get(sp, 0.0) + m
                notes.append('particle %s %.4g ug/m3 -> %s' % (n, m, sp))
        return gas, part, {n for n in aer_modes if re.match(r'^IVOC\d$', n)}, notes
    # POA1 + CBIO particle mass: PM emission
    tot = sum(m for n, m in aer_modes.items() if n in ('POA1', 'CBIO'))
    if True:
        for sp, fr in GENERAL_POA.items():
            fg = GENERAL_POA_GAS.get(sp, 0.0)
            if fg > 0:
                gas[sp] = tot * fg * _ppb_per_ugm3(table[sp][0], T, P)
            part[sp] = tot * (fr - fg)
            notes.append('GENERAL_POA %s: total %.3f of %.4g ug/m3 (gas %.4g ppb, particle %.4g ug/m3)' % (sp, fr, tot, gas.get(sp, 0.0), part[sp]))
    return gas, part, {n for n in aer_modes if n in ('POA1', 'CBIO')}, notes


def read_modes(path):
    """first mode's bulk composition of AerosolModes.in -> OrderedDict name -> ug/m3 (only used to find IVOC / POA names)"""
    from collections import OrderedDict
    out = OrderedDict()
    seen = False
    for ln in open(path, errors='ignore'):
        s = ln.replace('\r', '').strip()
        if s.startswith('!') or not s:
            continue
        if re.match(r'^[0-9.]+\s*;\s*[0-9.]+\s*;', s):
            if seen:
                break
            seen = True
            continue
        if seen and ';' in s and not s.endswith('Smoke;'):
            f = [x.strip() for x in s.split(';')]
            try:
                out[f[0]] = float(f[1])
            except ValueError:
                pass
    return out


def rewrite_modes(path, replaced, part_new):
    """remove the `replaced` species lines of every mode and add the CMAQ species (mass split across modes in proportion to the
    original mode totals of the replaced species)"""
    L = open(path, errors='ignore').read().replace('\r', '').split('\n')

    def spec(t):
        f = [x.strip() for x in t.split(';')]
        if len(f) < 2:
            return None
        try:
            return f[0], float(f[1])
        except ValueError:
            return None
    modes = []
    for i, s in enumerate(L):
        t = s.strip()
        if t.startswith('!') or not t:
            continue
        if re.match(r'^[0-9.]+\s*;\s*[0-9.]+\s*;', t):
            modes.append([]); continue
        if modes and spec(t):
            modes[-1].append(i)
    tot = [sum(spec(L[j])[1] for j in m if spec(L[j])[0] in replaced) for m in modes]
    s_all = sum(tot)
    last = {}
    for k, m in enumerate(modes):
        keep = [j for j in m if spec(L[j])[0] not in replaced]
        for j in m:
            if j not in keep:
                L[j] = None
        last[max(keep)] = ['%s; %.5g;' % (sp, v * (tot[k] / s_all if s_all > 0 else (1.0 if k == 0 else 0.0)))
                           for sp, v in part_new.items() if v > 0]
    res = []
    for i, s in enumerate(L):
        if s is None:
            continue
        res.append(s)
        if i in last:
            res.extend(last[i])
    open(path, 'w').write('\n'.join(res) + '\n')


def write_org_files(case_dir, partition, mass_acc=MASS_ACC):
    # OrgPhaseChems.in: drop the baseline IVOC rows, add the CMAQ species
    p = os.path.join(case_dir, 'OrgPhaseChems.in')
    L = [l for l in open(p, errors='ignore').read().replace('\r', '').split('\n') if not re.match(r'^IVOC\d', l)]
    while L and not L[-1].strip():
        L.pop()
    L.append('!')
    L.append('!CMAQ CRACMM organic-aerosol species (CCTM/src/aero/cracmm/ORG_DEFN.F): gas V* / particle A* pair = one ASP species')
    for n, mw, c, h, nv in partition:
        L.append('%-14s;1\t;%.2f\t; %.1f ;       1.0; 1.0;' % (n, mw, RHO_ORG))
    open(p, 'w').write('\n'.join(L) + '\n')
    # aqueous-phase organic lists: CMAQ has no aqueous uptake for the V*/A* species.  ASP's output routines (TotPOA/TotSOA/TotOA) index
    # AqOrgchems(J-1) in parallel with Orgchems(J), so the hydrophilic lists must keep the same species in the same order: the CMAQ species
    # are added with a negligible Henry's-law constant (as the baseline IVOCk rows had), after the baseline IVOC rows are removed.
    for fn in ('HydrophilicOrgChems.in', 'HydrophilicOrganicDissolution.in'):
        p = os.path.join(case_dir, fn)
        L = [l for l in open(p, errors='ignore').read().replace('\r', '').split('\n') if not re.match(r'^IVOC\d', l)]
        while L and not L[-1].strip():
            L.pop()
        for n, mw, c, h, nv in partition:
            if fn == 'HydrophilicOrgChems.in':
                L.append('%-8s;1      ;%-9.1f ; 1.5 ;  %d ;0.04; 1.00e7 ; 1. ; 1.0;' % (n, mw, max(1, int(round(mw / 14.0)))))
            else:
                L.append('%s <=> %s ;  0.1 ; 1.0e-16 ; 15.   ; 0.      ; 0.5  ; 0.     ; 0.5 ;' % (n, n))
        open(p, 'w').write('\n'.join(L) + '\n')
    # OrganicDissolution.in
    p = os.path.join(case_dir, 'OrganicDissolution.in')
    L = [l for l in open(p, errors='ignore').read().replace('\r', '').split('\n') if not re.match(r'^IVOC\d', l)]
    while L and not L[-1].strip():
        L.pop()
    L.append('!')
    L.append('! CMAQ CRACMM species: C*(298 K) and Delta Hvap from ORG_DEFN.F converted to Psat(300 K) in mbar (flag 1);')
    L.append('! C*(T) = C*(298) (298/T) exp[dHvap/R (1/298 - 1/T)] in CMAQ is the same functional form ASP uses.  Mass accommodation %.1f (CMAQ alphorg).' % mass_acc)
    L.append('! Reaction      ; Mass Acc ; Flag       ; Psat(300K) mbar    ; Delta Hvap (KJ/mol)  ; [C*(298K) ug/m3]')
    for n, mw, c, h, nv in partition:
        ps = psat300_mbar(c, h, mw)
        L.append('%s <=> %s\t; %.1f;\t1.\t; %.4e\t; %.1f\t; ! C*298 = %.3g' % (n, n, mass_acc, ps, h / 1e3, c))
    open(p, 'w').write('\n'.join(L) + '\n')
