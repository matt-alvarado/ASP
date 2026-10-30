import sys, os, math, pickle, glob, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify_gc as v
OUT = '/home/claude/work/probe'; T, P = 288.4, 880.0
data = {}
def fl(x):
    try: return float(x)
    except ValueError: return float('nan')
for f in sorted(glob.glob(OUT + '/OutputFiles/GasOut*.txt')):
    L = open(f).read().split('\n'); H = L[0]; W = 16
    hdr = [H[k:k + W].strip() for k in range(0, len(H), W) if H[k:k + W].strip()]
    rows = [x for x in L[1:] if x.strip()]
    for j, h in enumerate(hdr):
        col = [fl(r[j * W:(j + 1) * W]) for r in rows]
        if h == 'Time': data.setdefault('Time', col)
        else: data[h] = col
NT = 11
t = data['Time'][:NT]; print('times', t[:4], '...', t[-1], len(t))
exp, bpp = pickle.load(open(OUT + '/exp.pkl', 'rb')); kB = 1.380649e-23
M = P * 100. / (kB * T) * 1e-6; worst = 0
ctrl = data['PRC']
def kfit(Ai, Bi):
    a = data[Ai]; y = [math.log(a[n] / ctrl[n]) for n in range(NT)]; Bc = data[Bi]
    integ = sum(0.5 * (Bc[n] + Bc[n + 1]) * (t[n + 1] - t[n]) for n in range(len(t) - 1)) * 1e-9 * M
    return -(y[-1] - y[0]) / integ
cal = kfit('PRA0', 'PRB0') / 1e-12; print('calibration (measured / true) with a flagless 1e-12 reaction: %.4f' % cal)
for i, l, nfree in exp:
    a = data['PRA%d' % i]; tt = t[-1]
    if a[-1] <= 0: print(i, 'bad'); continue
    # ln(A/C) = -int k [B] dt   (control C has the same mixing / loss)
    y = [math.log(a[n] / ctrl[n]) for n in range(NT)]; loss = -(y[-1] - y[0])
    if nfree == 2:
        Bc = data['PRB%d' % i]; integ = sum(0.5 * (Bc[n] + Bc[n + 1]) * (t[n + 1] - t[n]) for n in range(len(t) - 1)) * 1e-9 * M
    else: integ = t[-1] - t[0]
    kmeas = loss / integ / (cal if nfree == 2 else 1.0)
    ref = v.asp_rate(l, T, M, 0.0)
    for x in [y_.strip() for y_ in l.split(' ; ')[0].split('=>')[0].split('+')]:
        if x == 'O2': ref *= .21 * M
        if x == 'N2': ref *= .79 * M
    if ref == 0: print(i, 'ref = 0 (floored)', kmeas); continue
    err = abs(kmeas - ref) / ref; worst = max(worst, err)
    print('%2d %-1s flag %-3s meas %.4e ref %.4e err %.2e  %s' % (i, 'B' if nfree == 2 else 'U', l.split(' ; ')[1].split(';')[0], kmeas, ref, err, l.split(' ; ')[0][:50]))
print('worst', worst)
