"""CMAQ AERO6 / AERO7 organic-aerosol chemistry for the CB6r5, SAPRC07TC (AE6) and SAPRC07TIC (AE7i) decks.

How CMAQ forms organic aerosol with these mechanisms (CCTM/src/aero/aero{6,7}/ORG_DEFN.F, the mech_*.def, AE_ / GC_ namelists and
CMAQ_Chem_Control_*.nml):

  * SOA vapors.  Each semivolatile vapor (V) shares an absorbing organic phase with its particle partner (A) and is described by
    C*(298 K), enthalpy of vaporization, molar mass and O:C (oa_list).  Vapors come from two places:
      - explicit gas-phase reactions of the mechanism (SVAVB1-4 from aromatic / PAH / long-alkane peroxy radicals, MTNO3, ISOPNN,
        POA-ageing products VLVOO*, VSVOO*, ...);
      - 'reaction counters' (ALKRXN, XYLNRXN, TOLHRXN, ISOPRXN, TRPRXN, SESQRXN, PCSOARXN ... produced by the first oxidation step).
        ORG_DEFN converts the counter, in ug/m3, with a fixed mass yield alpha into each vapor: vapor mass = alpha x counter x MW(counter).
  * Here the counter is replaced by its products: product coefficient = alpha x MW(counter) / MW(particle species), so a counter event
    makes exactly the vapor mass CMAQ makes.  One ASP species carries both the gas and the organic-phase mass (named like the CMAQ vapor,
    or like the particle species for the nonvolatile ones: AXYL3 ...).
  * Partitioning: ideal Raoult mixing in one organic phase, C*(T) = C*(298) (298/T) exp[dH/R (1/298 - 1/T)] (same form as ASP flag 1).
  * POA (default in CMAQ is nonvolatile APOC/APNCOM; the control namelist documents a semivolatile alternative that is used here):
    ALVPO1 0.09, ASVPO1 0.045 (particle) and VSVPO1 0.045, VSVPO2 0.14, VSVPO3 0.18, VIVPO1 0.50 (gas) of the PM emission; POA ageing
    reactions (CB6r5 only) move the vapors to VLVOO* / VSVOO*.  PCVOC (pcSOA) is 6.579 x PM in CMAQ but zeroed for fire sources, as in all
    four cases here.
  * Williams: the measured volatility distribution (gas IVOC3-7, particle IVOC3-6) goes bin-by-bin onto the CMAQ POA species with the same
    C* (IVOC3 -> VSVPO1, IVOC4 -> VSVPO2, IVOC5 -> VSVPO3, IVOC6 -> VIVPO1, IVOC1-2 -> VLVPO1); CMAQ has no species for C* >= 1e4
    (IVOC7-9), that mass is dropped (ivoc_high='drop') or sent to PCVOC (ivoc_high='pcvoc').
Not represented (ASP has no particle-phase organic chemistry): oligomerisation (AAVB/AISO/ASQT -> AOLGA/B, ~29 h), MTNO3/ISOPNN
hydrolysis (3 h), IEPOX / glyoxal heterogeneous uptake, aqueous organic uptake.
"""
import math
import os
import re
from collections import OrderedDict

import cmaq_oa

MECHS = '/home/claude/usepa/cmaq/CCTM/src/MECHS/'
AERO = '/home/claude/usepa/cmaq/CCTM/src/aero/'
SPEC = {'CB6R5': ('cb6r5_ae7_aq', 'aero7'), 'SAPRC07TC': ('saprc07tc_ae6_aq', 'aero6'), 'SAPRC07TIC': ('saprc07tic_ae7i_aq', 'aero7')}
R = 8.314462618
ASOA_CSTAR, ASOA_ENTH = 1e-10, 1.0e3            # CMAQ nonvolatile species
MASS_ACC = cmaq_oa.MASS_ACC

# POA volatility profile of the CMAQ control namelist: species -> (gas fraction, particle fraction) of the PM emission
POA_GAS = OrderedDict([('VSVPO1', 0.045), ('VSVPO2', 0.14), ('VSVPO3', 0.18), ('VIVPO1', 0.50)])
POA_PART = OrderedDict([('VLVPO1', 0.09), ('VSVPO1', 0.045)])
# ASP IVOCk (C* = 10^(k-3) ug/m3) -> CMAQ POA species of the same (nearest) C*
BIN_MAP = {1: 'VLVPO1', 2: 'VLVPO1', 3: 'VSVPO1', 4: 'VSVPO2', 5: 'VSVPO3', 6: 'VIVPO1'}
# counters that are real, non-peroxy-pool species in the mechanism (kept as species, not expanded)
RO2_COUNTERS = {'TOLRO2', 'XYLRO2', 'BENZRO2', 'PAHRO2'}

_OA = re.compile(r"\s*&\s*oa_type\('([^']+)'\s*,\s*'([^']*)'\s*,\s*'([^']*)'\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,"
                 r"\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([TF])\s*,\s*([TF])\s*,\s*([TF])\s*,\s*([TF])")
_NML = re.compile(r"^\s*'([^']+)'\s*,\s*([-\d.eE+]+)\s*,")


def read_oa_list(path):
    out = OrderedDict()
    for ln in open(path, errors='ignore'):
        m = _OA.match(ln)
        if m:
            out[m.group(1).strip()] = dict(pm=m.group(1).strip(), vap=m.group(2).strip(), ctr=m.group(3).strip(), alpha=float(m.group(4)),
                                           cstar=float(m.group(5)), enth=float(m.group(6)), nonvol=m.group(12) == 'T')
    return out


def read_mw(path):
    return {m.group(1): float(m.group(2)) for m in (_NML.match(l) for l in open(path, errors='ignore')) if m}


class AE:
    def __init__(self, mech):
        d, org = SPEC[mech]
        self.mech = mech
        self.ae_mw = read_mw(MECHS + '%s/AE_%s.nml' % (d, d))
        self.gc_mw = read_mw(MECHS + '%s/GC_%s.nml' % (d, d))
        self.oa = read_oa_list(AERO + '%s/ORG_DEFN.F' % org)
        self.rows = OrderedDict((n, r) for n, r in self.oa.items() if n in self.ae_mw)     # organic species simulated by this mechanism
        # counter -> [(species, coefficient per counter mole)]
        self.ctr = OrderedDict()
        self.mw = {}
        for n, r in self.rows.items():
            sp = self._species(r)
            self.mw[sp] = self.ae_mw[n]
        for n, r in self.rows.items():
            if r['ctr'] and r['alpha'] > 0:
                sp = self._species(r)
                self.ctr.setdefault(r['ctr'], []).append((sp, r['alpha'] * self.gc_mw[r['ctr']] / self.ae_mw[n]))
        self.used_ctr = set()

    @staticmethod
    def _species(r):
        return r['vap'] if r['vap'] and not r['nonvol'] else r['pm']

    def is_counter(self, n):
        return n in self.ctr and n not in RO2_COUNTERS

    def expand(self, n, coef):
        self.used_ctr.add(n)
        return [(sp, coef * y) for sp, y in self.ctr[n]]

    def partition(self, species):
        """-> [(name, MW, C*298, dH J/mol, nonvolatile)] for the simulated organic species present in `species` (or produced via counters)"""
        want = set(species)
        for c in self.used_ctr:
            want |= {sp for sp, _ in self.ctr[c]}
        out = []
        for n, r in self.rows.items():
            sp = self._species(r)
            if sp not in want:
                continue
            if r['nonvol']:
                if r['ctr'] and r['alpha'] > 0:     # nonvolatile products made directly from counters (AXYL3 ...)
                    out.append((sp, self.ae_mw[n], ASOA_CSTAR, ASOA_ENTH, True))
                continue
            out.append((sp, self.ae_mw[n], r['cstar'], r['enth'], False))
        return out

    def emissions(self, old_rows, aer_modes, T, P, ivoc_high='drop'):
        """-> (gas dict name -> ppb, particle dict name -> ug/m3, replaced particle names, note lines)"""
        gas, part, notes = {}, {}, []
        ppb_per_ug = lambda mw: 1.0 / (mw * (P * 100.0) / (R * T) * 1e-3)
        ivoc_g = {k: old_rows['IVOC%d' % k][3] for k in range(1, 10) if 'IVOC%d' % k in old_rows}
        ivoc_mw = {k: old_rows['IVOC%d' % k][2] for k in range(1, 10) if 'IVOC%d' % k in old_rows}
        if any(v > 0 for v in ivoc_g.values()) or any(re.match(r'^IVOC\d$', n) for n in aer_modes):
            for k, ppb in ivoc_g.items():
                if ppb <= 0:
                    continue
                if k in BIN_MAP:
                    sp = BIN_MAP[k]
                    gas[sp] = gas.get(sp, 0.0) + ppb * ivoc_mw[k] / self.mw[sp]
                    notes.append('gas IVOC%d %.4g ppb (MW %.0f) -> %s (MW %.0f)' % (k, ppb, ivoc_mw[k], sp, self.mw[sp]))
                elif ivoc_high == 'pcvoc':
                    gas['PCVOC'] = gas.get('PCVOC', 0.0) + ppb * ivoc_mw[k] / self.gc_mw['PCVOC']
                    notes.append('gas IVOC%d %.4g ppb (C* >= 1e4, no CMAQ POA species) -> PCVOC (pcSOA precursor)' % (k, ppb))
                else:
                    notes.append('gas IVOC%d %.4g ppb (%.4g ug/m3, C* >= 1e4): no CMAQ species -> not carried' %
                                 (k, ppb, ppb * ivoc_mw[k] * P * 100.0 / (R * T) * 1e-3))
            for n, m in aer_modes.items():
                mm = re.match(r'^IVOC(\d)$', n)
                if mm:
                    k = int(mm.group(1))
                    sp = BIN_MAP.get(k)
                    if sp is None:
                        notes.append('particle %s %.4g ug/m3: no CMAQ species -> not carried' % (n, m))
                        continue
                    part[sp] = part.get(sp, 0.0) + m
                    notes.append('particle %s %.4g ug/m3 -> %s' % (n, m, sp))
            return gas, part, {n for n in aer_modes if re.match(r'^IVOC\d$', n)}, notes
        tot = sum(m for n, m in aer_modes.items() if n in ('POA1', 'CBIO'))
        for sp in list(POA_GAS) + [s for s in POA_PART if s not in POA_GAS]:
            fg, fp = POA_GAS.get(sp, 0.0), POA_PART.get(sp, 0.0)
            if fg > 0:
                gas[sp] = tot * fg * ppb_per_ug(self.mw[sp])
            part[sp] = tot * fp
            notes.append('POA profile %s: gas %.3f + particle %.3f of the %.4g ug/m3 PM emission' % (sp, fg, fp, tot))
        return gas, part, {n for n in aer_modes if n in ('POA1', 'CBIO')}, notes
