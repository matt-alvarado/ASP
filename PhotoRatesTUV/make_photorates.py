#!/usr/bin/env python3
"""make_photorates.py -- build an ASP PhotoRates.in table with NCAR TUV 5.4.

Usage example (Williams, AZ-like site, 21 June, polluted 1 km PBL):

    python3 make_photorates.py --tuv /path/to/V5.4 \
        --lat 35.2 --lon -112.2 --date 2010-06-21 --alt 2.0 \
        --pbl-height 1.0 --tau-aer 0.3 --tau-no2 0.02 --tau-so2 0.002 \
        --tau-o3 0.05 -o PhotoRates.in

What it does
  1. Copies the TUV 5.4 Fortran source, applies tuv_asp_pbl.patch (batch mode,
     PBL inputs read from a file, gas optical depth -> mixing ratio) and
     compiles it with gfortran.  The original TUV directory is never modified.
  2. Writes a TUV input file for one UT day (96 steps, UT 0.00 ... 23.75 h)
     and a pbl.inp file with the boundary-layer inputs, and runs TUV once.
  3. Reads the J-values TUV prints, maps them by reaction name to the 51 ASP
     columns (asp_photorates_map.py) and writes PhotoRates.in in the format
     read by ReadPhotolysisRates (Chemistry.f90).

Boundary-layer inputs (all optional; default = background atmosphere)
  --pbl-height   PBL depth above ground (km).  Rounded down to the 0.5 km
                 vertical grid of the calculation.
  --tau-aer      PBL aerosol extinction optical depth (at --wl-aer, default
                 330 nm; scaled to 330 nm as 1/wavelength, TUV's own PBL
                 assumption).  Replaces the Elterman profile inside the PBL.
  --tau-o3 / --tau-no2 / --tau-so2
                 PBL *gas absorption* optical depth at --wl-o3 / --wl-no2 /
                 --wl-so2.  Converted to a uniform PBL mixing ratio with
                 TUV's own surface-layer cross section for the wavelength
                 bin that contains the reference wavelength.
  --o3-ppb / --no2-ppb / --so2-ppb
                 alternative: give the PBL mixing ratio directly (ppb).
  Above the PBL: O3 follows the US Standard Atmosphere profile scaled to
  --o3-column; NO2 and SO2 are zero; aerosol follows Elterman (1968).
  If a PBL gas is not given, nothing is inserted for it (O3 keeps its
  standard-atmosphere value in the PBL).
"""
import argparse
import datetime as dt
import hashlib
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from asp_photorates_map import ASP_COLUMNS  # noqa: E402

NT = 96            # quarter-hour steps, UT 0.00 ... 23.75
DZ_KM = 0.5        # vertical grid spacing
ZTOP_KM = 60.0     # model top above ground (nz = 121 <= TUV kz = 125)
SENTINEL = "END-OF-FILE-012343210-ELIF-FO-DNE"


# ----------------------------------------------------------------------------
# build
# ----------------------------------------------------------------------------
def find_tuv_root(path):
    """Accept the V5.4 directory, its parent, or the distribution zip."""
    path = os.path.abspath(path)
    if os.path.isfile(path) and path.endswith(".zip"):
        import zipfile
        out = path[:-4] + "_unpacked"
        if not os.path.isdir(out):
            with zipfile.ZipFile(path) as z:
                z.extractall(out)
        path = out
    for cand in (path, os.path.join(path, "V5.4")):
        if os.path.isfile(os.path.join(cand, "TUV.f")):
            return cand
    for root, _dirs, files in os.walk(path):
        if "TUV.f" in files:
            return root
    sys.exit("Could not find TUV.f under %s" % path)


def build_tuv(tuv_root, build_dir):
    """Copy source, apply patch, compile.  Returns path of the executable."""
    patch = os.path.join(HERE, "tuv_asp_pbl.patch")
    h = hashlib.sha1()
    for fn in sorted(os.listdir(tuv_root)):
        if fn.endswith(".f") or fn == "params":
            h.update(open(os.path.join(tuv_root, fn), "rb").read())
    h.update(open(patch, "rb").read())
    tag = h.hexdigest()[:12]
    exe = os.path.join(build_dir, "tuv")
    stamp = os.path.join(build_dir, "build.stamp")
    if os.path.isfile(exe) and os.path.isfile(stamp) and open(stamp).read() == tag:
        return exe
    shutil.rmtree(build_dir, ignore_errors=True)
    os.makedirs(build_dir)
    for fn in os.listdir(tuv_root):
        if fn.endswith(".f") or fn in ("params", "Makefile"):
            shutil.copy(os.path.join(tuv_root, fn), build_dir)
    subprocess.run(["patch", "-p1", "-s", "-i", patch], cwd=build_dir, check=True)
    if shutil.which("gfortran") is None:
        sys.exit("gfortran not found on PATH")
    base = ["make", "-s", "FC=gfortran"]
    for flags in ("-O2 -w", "-O2 -w -fallow-argument-mismatch -std=legacy"):
        subprocess.run(["make", "-s", "clean"], cwd=build_dir,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        r = subprocess.run(base + ["FFLAGS=" + flags], cwd=build_dir,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if r.returncode == 0 and os.path.isfile(exe):
            break
    else:
        sys.exit("TUV build failed:\n" + r.stdout[-3000:])
    open(stamp, "w").write(tag)
    return exe


# ----------------------------------------------------------------------------
# TUV input files
# ----------------------------------------------------------------------------
def set_field(lines, name, value):
    """Replace one fixed-format field ('name =' label is 9 wide, value 11)."""
    pat = re.compile(r"(?<![A-Za-z0-9])%s =" % re.escape(name))
    for n, line in enumerate(lines[:20]):
        m = pat.search(line)
        if m:
            k = m.start()
            lines[n] = line[:k + 9] + str(value).rjust(11) + line[k + 20:]
            return
    raise KeyError(name)


def write_tuv_input(tuv_root, path, a, zstart):
    lines = open(os.path.join(tuv_root, "INPUTS", "defin2")).read().split("\n")
    nz = int(round(ZTOP_KM / DZ_KM)) + 1
    vals = dict(
        lat="%.4f" % a.lat, lon="%.4f" % a.lon, tmzone="0.0",
        iyear=a.date.year, imonth=a.date.month, iday=a.date.day,
        zstart="%.4f" % zstart, zstop="%.4f" % (zstart + ZTOP_KM), nz=nz,
        tstart="0.000", tstop="%.3f" % (0.25 * (NT - 1)), nt=NT,
        alsurf="%.3f" % a.albedo, psurf="%.1f" % a.psurf,
        o3col="%.3f" % a.o3_column, so2col="0.000", no2col="0.000",
        taucld="%.3f" % a.cloud_tau, zbase="%.3f" % a.cloud_base,
        ztop="%.3f" % a.cloud_top,
        tauaer="%.3f" % a.bg_aer550, ssaaer="%.3f" % a.bg_ssa, alpha="%.3f" % a.angstrom,
        zout="%.4f" % (zstart + a.zout), nstr=a.nstr,
    )
    for k, v in vals.items():
        set_field(lines, k, v)
    # laflux off, J-values on; select every photolysis reaction
    lines = [l.replace("laflux =           T", "laflux =           F") for l in lines]
    sec = None
    nj = 0
    for n, line in enumerate(lines):
        if line.startswith("===== Available photolysis"):
            sec = "j"
        elif line.startswith("===== Avaliable spectral"):
            sec = "s"
        elif sec == "j" and re.match(r"^[TF]\s+\d+ ", line):
            lines[n] = "T" + line[1:]
            nj += 1
    set_field(lines, "nmj", nj)
    open(path, "w").write("\n".join(lines))


def tau_to_wl330(tau, wl):
    """TUV scales PBL aerosol as aod(w) = aod330 * (w/330)**-1."""
    return tau * wl / 330.0


def write_pbl_input(path, a):
    def gas(tau, ppb, wl):
        if ppb is not None:
            return "%.6e 0 %.1f" % (ppb, wl)
        return "%.6e 1 %.1f" % (tau if tau else 0.0, wl)
    if a.pbl_height <= 0:
        open(path, "w").write("-999.\n0 0 300\n0 0 400\n0 0 300\n-1 0.6 0.9\n")
        return
    aod330 = -1.0 if a.tau_aer is None else tau_to_wl330(a.tau_aer, a.wl_aer)
    with open(path, "w") as f:
        f.write("%.4f\n" % a.pbl_height)
        f.write(gas(a.tau_o3, a.o3_ppb, a.wl_o3) + "\n")
        f.write(gas(a.tau_no2, a.no2_ppb, a.wl_no2) + "\n")
        f.write(gas(a.tau_so2, a.so2_ppb, a.wl_so2) + "\n")
        f.write("%.6e %.3f %.3f\n" % (aod330, a.ssa_uv, a.ssa_vis))


def run_tuv(exe, tuv_root, run_dir, a, zstart):
    work = os.path.join(run_dir, "work")
    shutil.rmtree(run_dir, ignore_errors=True)
    os.makedirs(os.path.join(work, "INPUTS"))
    shutil.copy(exe, os.path.join(work, "tuv"))
    shutil.copy(os.path.join(tuv_root, "helpin"), work)
    for d in ("DATAE1", "DATAJ1", "DATAS1"):
        src = os.path.join(tuv_root, d)
        if os.path.isdir(src):
            os.symlink(src, os.path.join(work, d))
    write_tuv_input(tuv_root, os.path.join(work, "INPUTS", "usrinp"), a, zstart)
    write_pbl_input(os.path.join(work, "pbl.inp"), a)
    # The patched TUV runs in batch mode, but rdinp still prompts for nothing.
    with open(os.path.join(run_dir, "tuv_stdout.log"), "w") as log:
        r = subprocess.run(["./tuv"], cwd=work, stdin=subprocess.DEVNULL,
                           stdout=log, stderr=subprocess.STDOUT)
    out = os.path.join(run_dir, "usrout.txt")
    if r.returncode != 0 or not os.path.isfile(out):
        sys.exit("TUV failed; see %s" % os.path.join(run_dir, "tuv_stdout.log"))
    return out


# ----------------------------------------------------------------------------
# read TUV output, write PhotoRates.in
# ----------------------------------------------------------------------------
def fnum(tok):
    # Fortran prints 3-digit exponents without the 'E' (1.234-100)
    tok = re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", tok)
    return float(tok)


def read_tuv_output(path):
    labels = {}
    rows = []
    z_km = None
    for line in open(path):
        m = re.match(r"^\s*(\d+) = (.*?)\s*$", line)
        if m and not rows:
            labels[m.group(2)] = int(m.group(1))
            continue
        m = re.match(r"^values at z =\s*([-\d.]+)", line)
        if m:
            z_km = float(m.group(1))
            continue
        tok = line.split()
        if len(tok) > 3 and re.match(r"^-?\d+\.\d+$", tok[0]) and re.match(r"^-?\d+\.\d+$", tok[1]):
            rows.append([fnum(t) for t in tok])
    if len(rows) != NT:
        sys.exit("Expected %d TUV time steps, found %d" % (NT, len(rows)))
    return labels, rows, z_km


def asp_rows(labels, rows):
    out = []
    for r in rows:
        j = {}
        for idx, _lab, rule in ASP_COLUMNS:
            if isinstance(rule, str):
                v = r[1 + labels[rule]]
            elif rule[0] == "scale":
                v = rule[1] * j[rule[2]]
            else:
                v = sum(r[1 + labels[x]] for x in rule)
            j[idx] = v if abs(v) >= 1.0e-99 else 0.0
        out.append((r[0], r[1], [j[i] for i, _l, _r in ASP_COLUMNS]))
    return out


def derived_pbl(logpath):
    """Mixing ratios / cross sections TUV derived from the optical depths."""
    out = []
    for line in open(logpath):
        m = re.match(r"\s*PBL (\w+)\s+tau->ppb, xs, bin:\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line)
        if m:
            out.append("PBL %s: %.4g ppb from tau (sigma %.3g cm2, TUV bin %.1f-%.1f nm)"
                       % (m.group(1), float(m.group(2)), float(m.group(3)),
                          float(m.group(4)), float(m.group(5))))
    return out


def write_photorates(path, table, z_km, a, cmdline):
    with open(path, "w", newline="\n") as f:
        f.write("!! PhotoRates.in\n!! ASP\n")
        f.write("!! Generated by make_photorates.py using TUV 5.4 (S. Madronich, NCAR)\n")
        f.write("!! Command: %s\n" % cmdline)
        f.write("!!\n")
        f.write("!! Read by subroutine ReadPhotolysisRates in Chemistry.f90.\n")
        f.write("!! 96 lines of data (every quarter hour, UT 0.0 ... 23.75), then the sentinel.\n")
        f.write("!! Columns: time (UT, hrs), solar zenith angle (deg), then 51 photolysis rate\n")
        f.write("!! coefficients (s-1) in the order listed below.\n")
        f.write("!!\n")
        f.write("!! Site: lat %.4f, lon %.4f, date %s, ground altitude %.3f km asl,\n"
                % (a.lat, a.lon, a.date.isoformat(), a.alt))
        f.write("!!   J-values at %.3f km asl; albedo %.3f; total O3 column %.0f DU;\n"
                % (z_km, a.albedo, a.o3_column))
        if a.pbl_height > 0:
            f.write("!!   PBL depth %.2f km above ground (rounded to %.1f km grid)\n"
                    % (a.pbl_height, DZ_KM))
            f.write("!!   PBL tau: aerosol %s @%.0f nm, O3 %s @%.0f nm, NO2 %s @%.0f nm, SO2 %s @%.0f nm\n"
                    % (a.tau_aer, a.wl_aer, a.tau_o3, a.wl_o3, a.tau_no2, a.wl_no2,
                       a.tau_so2, a.wl_so2))
            f.write("!!   PBL ppb (if given): O3 %s, NO2 %s, SO2 %s; aerosol SSA %.2f (<400 nm) %.2f (>=400 nm)\n"
                    % (a.o3_ppb, a.no2_ppb, a.so2_ppb, a.ssa_uv, a.ssa_vis))
            for line in getattr(a, "derived", []):
                f.write("!!   %s\n" % line)
        else:
            f.write("!!   No PBL perturbation (background atmosphere)\n")
        f.write("!!\n!! Photolysis rate coefficients, s-1\n")
        for idx, lab, _rule in ASP_COLUMNS:
            f.write("!%4d = %s\n" % (idx, lab))
        f.write("!values at z = %9.3f km\n" % z_km)
        f.write("!Columns: time, sza, photo-reactions\n")
        f.write("!time, hrs.  sza, deg." + "".join("%11d" % i for i, _l, _r in ASP_COLUMNS) + "\n")
        for ut, sza, js in table:
            f.write("%10.4f%11.3f" % (ut, sza) + "".join("%11.3e" % v for v in js) + "\n")
        f.write(SENTINEL + "\n")


# ----------------------------------------------------------------------------
def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tuv", default=os.environ.get("TUV_DIR"),
                   help="TUV 5.4 directory (the one containing TUV.f) or the zip "
                        "[env TUV_DIR]")
    p.add_argument("--lat", type=float, required=True, help="latitude, deg N")
    p.add_argument("--lon", type=float, required=True, help="longitude, deg E (west negative)")
    p.add_argument("--date", type=lambda s: dt.date.fromisoformat(s), required=True,
                   help="UT date, YYYY-MM-DD")
    p.add_argument("--alt", type=float, default=0.0, help="ground altitude, km asl [0]")
    p.add_argument("--zout", type=float, default=0.0,
                   help="height of the J-values above ground, km [0]")
    p.add_argument("--albedo", type=float, default=0.1, help="surface albedo [0.1]")
    p.add_argument("--o3-column", type=float, default=300.0, help="total O3 column, DU [300]")
    p.add_argument("--psurf", type=float, default=-999.0,
                   help="surface pressure, mbar (<0: standard atmosphere at --alt) [-999]")
    g = p.add_argument_group("boundary layer")
    g.add_argument("--pbl-height", type=float, default=0.0,
                   help="PBL depth above ground, km (0 = no PBL perturbation)")
    g.add_argument("--tau-aer", type=float, default=None,
                   help="PBL aerosol optical depth (default: keep Elterman profile)")
    g.add_argument("--wl-aer", type=float, default=330.0, help="wavelength of --tau-aer, nm [330]")
    g.add_argument("--ssa-uv", type=float, default=0.6, help="PBL aerosol SSA < 400 nm [0.6]")
    g.add_argument("--ssa-vis", type=float, default=0.9, help="PBL aerosol SSA >= 400 nm [0.9]")
    for gas, wl in (("o3", 310.0), ("no2", 400.0), ("so2", 300.0)):
        g.add_argument("--tau-" + gas, type=float, default=None,
                       help="PBL %s absorption optical depth" % gas.upper())
        g.add_argument("--wl-" + gas, type=float, default=wl,
                       help="wavelength of --tau-%s, nm [%g]" % (gas, wl))
        g.add_argument("--%s-ppb" % gas, type=float, default=None,
                       help="PBL %s mixing ratio, ppb (instead of --tau-%s)" % (gas.upper(), gas))
    h = p.add_argument_group("background atmosphere / numerics")
    h.add_argument("--bg-aer550", type=float, default=0.235,
                   help="column aerosol optical depth at 550 nm of the Elterman background "
                        "profile (0 = none) [0.235, the TUV default]")
    h.add_argument("--bg-ssa", type=float, default=0.99, help="background aerosol SSA [0.99]")
    h.add_argument("--angstrom", type=float, default=1.0, help="Angstrom exponent [1.0]")
    h.add_argument("--cloud-tau", type=float, default=0.0, help="cloud optical depth [0]")
    h.add_argument("--cloud-base", type=float, default=4.0, help="cloud base, km asl")
    h.add_argument("--cloud-top", type=float, default=5.0, help="cloud top, km asl")
    h.add_argument("--nstr", type=int, default=4, help="discrete-ordinate streams [4]")
    p.add_argument("--workdir", default="tuv_work", help="scratch/build directory [tuv_work]")
    p.add_argument("-o", "--output", default="PhotoRates.in")
    a = p.parse_args(argv)
    if not a.tuv:
        p.error("--tuv (or env TUV_DIR) is required")
    for gas in ("o3", "no2", "so2"):
        if getattr(a, "tau_" + gas) is not None and getattr(a, gas + "_ppb") is not None:
            p.error("give either --tau-%s or --%s-ppb, not both" % (gas, gas))
    if a.pbl_height > 0 and a.pbl_height < DZ_KM:
        p.error("--pbl-height must be at least %.1f km (grid spacing)" % DZ_KM)
    if a.pbl_height == 0 and any(x is not None for x in (
            a.tau_aer, a.tau_o3, a.tau_no2, a.tau_so2, a.o3_ppb, a.no2_ppb, a.so2_ppb)):
        p.error("PBL optical depths / mixing ratios need --pbl-height > 0")
    return a


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    a = parse_args(argv)
    work = os.path.abspath(a.workdir)
    os.makedirs(work, exist_ok=True)
    tuv_root = find_tuv_root(a.tuv)
    exe = build_tuv(tuv_root, os.path.join(work, "build"))
    zstart = a.alt
    out = run_tuv(exe, tuv_root, os.path.join(work, "run"), a, zstart)
    labels, rows, z_km = read_tuv_output(out)
    table = asp_rows(labels, rows)
    a.derived = derived_pbl(os.path.join(work, "run", "tuv_stdout.log"))
    for line in a.derived:
        print("  " + line)
    write_photorates(a.output, table, z_km if z_km is not None else zstart + a.zout, a,
                     "make_photorates.py " + " ".join(argv))
    print("wrote %s  (%d rows, 51 reactions; min SZA %.2f deg at UT %.2f)"
          % (a.output, len(table), min(r[1] for r in table),
             min(table, key=lambda r: r[1])[0]))


if __name__ == "__main__":
    main()
