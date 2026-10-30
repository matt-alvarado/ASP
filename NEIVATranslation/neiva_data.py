"""NEIVAv1.0 data layer: reads the CSV / SQL dumps of the NEIVA repository (no MySQL needed).

  load_efcol_info()        -> DataFrame (efcol, fire_type, measurement_type, study, fuel_type, MCE)  from bkdb_info_efcol
  load_property()          -> DataFrame (id, formula, compound, smile, voc_type, S07, MOZT1, GEOS_chem, cstar (= log10 C*), kOH, ...)  from Property_Surrogate
  temperate / fire-type EF matrix: efs(fire_type) -> DataFrame rows = compounds, columns = study EF columns with an MCE
NEIVA: Shahid et al., GMD 17, 7679 (2024), https://github.com/NEIVA-BB-Emissions-Inventory/NEIVAv1.0  (please cite)
"""
import os, re
import numpy as np, pandas as pd

NEIVA = os.environ.get('NEIVA_DIR', '/home/claude/neiva-bb-emissions-inventory/neivav1.0')


def _split_tuples(block):
    """parse the VALUES (...),(...) list of a MySQL INSERT into python tuples"""
    out, i, n = [], 0, len(block)
    while i < n:
        if block[i] != '(':
            i += 1; continue
        i += 1; row = []
        while True:
            c = block[i]
            if c == "'":
                j = i + 1; buf = []
                while True:
                    if block[j] == '\\':
                        buf.append(block[j + 1]); j += 2; continue
                    if block[j] == "'":
                        if j + 1 < n and block[j + 1] == "'":
                            buf.append("'"); j += 2; continue
                        break
                    buf.append(block[j]); j += 1
                row.append(''.join(buf)); i = j + 1
            else:
                j = i
                while block[j] not in ',)':
                    j += 1
                tok = block[i:j]
                row.append(None if tok == 'NULL' else float(tok)); i = j
            if block[i] == ',':
                i += 1; continue
            if block[i] == ')':
                i += 1; break
        out.append(tuple(row))
    return out


def _insert_block(path, table):
    s = open(path, errors='ignore').read()
    key = 'INSERT INTO `%s` VALUES' % table
    j = s.index(key) + len(key)
    return _split_tuples(s[j:s.index(';\n', j)])


def load_efcol_info():
    rows = _insert_block(os.path.join(NEIVA, 'data', 'backend_db.sql'), 'bkdb_info_efcol')
    return pd.DataFrame([dict(efcol=r[0], fire_type=r[1], measurement_type=r[2], study=r[3], fuel_type=r[4], MCE=r[10]) for r in rows])


def load_property():
    rows = _insert_block(os.path.join(NEIVA, 'data', 'neiva_output_db.sql'), 'Property_Surrogate')
    cols = ['mm', 'formula', 'compound', 'smile', 'voc_type', 'method', 'S07', 'S18B', 'S22', 'S07T', 'MOZT1', 'GEOS_chem', 'vp_nannoolal', 'vp',
            'vp_ref', 'cstar', 'kOH', 'kOH_ref', 'ko3_exp', 'kno3_exp', 'hc', 'hc_ref', 'OCratio', 'oxidation_state', 'id']
    return pd.DataFrame(rows, columns=cols)


def load_integrated():
    return pd.read_csv(os.path.join(NEIVA, 'data', 'CSV', 'Integrated_EF.csv'), low_memory=False)


def efs(fire_type):
    """Integrated EFs (g/kg) of one fire type for the study columns that carry an MCE.  Returns (ef DataFrame, info DataFrame)"""
    ef = load_integrated()
    info = load_efcol_info()
    info = info[(info.fire_type == fire_type) & info.MCE.notna() & info.efcol.isin(ef.columns)].reset_index(drop=True)
    meta = ef[['mm', 'formula', 'compound', 'pollutant_category', 'id']].copy()
    vals = ef[list(info.efcol)].apply(pd.to_numeric, errors='coerce')
    return pd.concat([meta, vals], axis=1), info
