"""PhotoRates.in columns for the CMAQ photolysis labels of CRACMM2/3 and SAPRC07TC/TIC.

ASP's PhotoRates.in holds 51 TUV-5.4 photolysis rates per case (time series, 15 min).  For every CMAQ label
(label, scale) used by a mechanism the rate is built in one of three ways:

  analog : J = scale * sum_k w_k J_TUV,k           the TUV column(s) is the same photolysis process
  ratio  : J = scale * J_TUV,ref * S_label / S_ref  no TUV column for the process.  S_x(t) = sum_b a_b(t) F_b CS_x,b EQY_x,b,
                                                    where F_b, CS, EQY are CMAQ's 7-band flux, cross sections and effective
                                                    quantum yields (CSQY_DATA of the mechanism) at the case temperature, and
                                                    a_b(t) >= 0 are band attenuation factors fitted (non-negative least squares,
                                                    relative weights) at every time row to the analog labels of the same file.
  column : J = scale * old TUV column  (existing ASP convention, e.g. HPALD -> column 51)
"""
import re
from collections import OrderedDict

import numpy as np

import photofit as pf

ANALOG = {
    'NO2_06': [(6, 1)], 'NO2_RACM2': [(6, 1)], 'O3O3P_06': [(3, 1)], 'O3O3P_NASA06': [(3, 1)],
    'O3O1D_06': [(2, 1)], 'O3O1D_NASA06': [(2, 1)], 'H2O2': [(5, 1)], 'H2O2_RACM2': [(5, 1)],
    'NO3NO_06': [(7, 1)], 'NO3NO_RACM2': [(7, 1)], 'NO3NO2_6': [(8, 1)], 'NO3NO2_RACM2': [(8, 1)],
    'HONO_06': [(12, 1)], 'HONO_RACM2': [(12, 1)], 'HNO3': [(13, 1)], 'HNO3_RACM2': [(13, 1)],
    'HNO4_06': [(14, 1)], 'HNO4_RACM2': [(14, 1)], 'HCHOR_06': [(17, 1)], 'HCHO_RAD_JPL19': [(17, 1)],
    'HCHOM_06': [(18, 1)], 'HCHO_MOL_JPL19': [(18, 1)], 'CCHO_R': [(19, 1)], 'CH3CHO_RACM2': [(19, 1)],
    'C2CHO': [(22, 1)], 'ALD_JPL19': [(22, 1)], 'ACET_06': [(26, 1)], 'CH3COCH3A_JPL19': [(26, 1)],
    'MEK_06': [(28, 1)], 'MEK_JGR19': [(28, 1)], 'KET_JGR19': [(28, 1)], 'COOH': [(31, 1)], 'OP1_RACM2': [(31, 1)],
    'PAA': [(33, 1)], 'PAA_RACM2': [(33, 1)], 'IC3ONO2': [(36, 1)], 'ONIT_CALVERT08': [(36, 1)],
    'ONIT_RACM2': [(36, 1)], 'NOA': [(38, 1)], 'PAN': [(40, 1), (41, 1)], 'PAN1_JPL19': [(40, 1)],
    'PAN2_JPL19': [(41, 1)], 'PAN1_RACM2': [(40, 1)], 'PAN2_RACM2': [(41, 1)], 'PPN1_JPL19': [(42, 1)],
    'PPN2_JPL19': [(43, 1)], 'HKET_RACM2': [(29, 1), (30, 1)], 'HOCCHO_IUPAC': [(23, 1)],
    'GLY_07R': [(44, 1)], 'GLYHX_JPL19': [(44, 1)], 'GLYHX_RACM2': [(44, 1)], 'GLY_07M': [(45, 1)],
    'GLYF_JPL19': [(45, 1)], 'GLYF_RACM2': [(45, 1)], 'MGLY_06': [(46, 1)], 'MGLY_RACM2': [(46, 1)],
    'BACL_07': [(47, 1)], 'ACRO_09': [(24, 1)], 'MACR_06': [(25, 1)], 'MACR_RACM2': [(25, 1)],
    'MVK_06': [(27, 1)], 'MVK_JPL19': [(27, 1)], 'MVK_RACM2': [(27, 1)],
}
# labels without a TUV column:  label -> (reference label in the same CSQY file, reference TUV columns)
RATIO = {
    'CH3COCH3B_JPL19': ('CH3COCH3A_JPL19', [(26, 1)]),      # second acetone channel relative to the first
    'UALD_RACM2': ('ACRO_09', [(24, 1)]),                    # unsaturated aldehydes ~ acrolein
    'GLYH2_RACM2': ('GLYHX_RACM2', [(44, 1)]), 'GLYH2_JPL19': ('GLYHX_JPL19', [(44, 1)]),
    'BALD_RACM2': ('MEK_JGR19', [(28, 1)]), 'BALD1_CALVERT11': ('MEK_JGR19', [(28, 1)]),
    'BALD2_CALVERT11': ('MEK_JGR19', [(28, 1)]),
    'TRPN_WANG2023': ('ONIT_CALVERT08', [(36, 1)]),
    'BALD_06': ('MEK_06', [(28, 1)]), 'AFG1': ('NO2_06', [(6, 1)]),
}
COLUMN = {'HPALD': [(51, 1)]}          # ASP convention: 100 x J(MACR), TUV column 51


def read_tuv(case_dir):
    rows = []
    for ln in open(case_dir + '/PhotoRates.in', errors='ignore').read().replace('\r', '').split('\n'):
        f = ln.split()
        if len(f) == 53 and re.match(r'^[0-9.]+$', f[0]):
            rows.append([float(x) for x in f])
    a = np.array(rows)
    assert a.shape[0] == 96, a.shape
    return a


class PhotoBuilder:
    def __init__(self, csqy_path, tuv, T):
        self.flux, self.tab = pf.read_csqy(csqy_path)
        self.tuv = tuv                          # (96, 53)
        self.T = T
        # analog labels present in the file -> rows of A, J
        self.fit_labels = [l for l in ANALOG if l in self.tab]
        self.A = np.array([self.flux * pf.csqy_at(self.tab, l, T) for l in self.fit_labels])
        self.J = np.array([[sum(w * row[1 + c] for c, w in ANALOG[l]) for l in self.fit_labels] for row in tuv])
        self.a = np.array([pf.fit_weights(self.J[i], self.A) for i in range(len(tuv))])      # (96, 7)

    def S(self, label):
        return self.a @ (self.flux * pf.csqy_at(self.tab, label, self.T))

    def column(self, label, scale):
        """-> (values[96], description)"""
        t = self.tuv
        if label in ANALOG:
            v = sum(w * t[:, 1 + c] for c, w in ANALOG[label])
            d = 'TUV col ' + ' + '.join(str(c) if w == 1 else '%g x %d' % (w, c) for c, w in ANALOG[label])
        elif label in COLUMN:
            v = sum(w * t[:, 1 + c] for c, w in COLUMN[label])
            d = 'TUV col ' + ' + '.join(str(c) for c, w in COLUMN[label]) + ' (ASP convention)'
        elif label in RATIO:
            ref, rcols = RATIO[label]
            vref = sum(w * t[:, 1 + c] for c, w in rcols)
            Sx, Sr = self.S(label), self.S(ref)
            ratio = np.where(Sr > 0, Sx / np.where(Sr > 0, Sr, 1.0), 0.0)
            v = vref * ratio
            pos = ratio[(vref > 0) & (Sr > 0)]
            d = 'J(%s, TUV col %s) x band-weighted CSQY ratio %s/%s (ratio %.2g-%.2g)' % (
                ref, '+'.join(str(c) for c, w in rcols), label, ref, pos.min() if pos.size else 0, pos.max() if pos.size else 0)
        else:
            raise KeyError('no photolysis treatment for CMAQ label ' + label)
        return scale * v, d
