"""Digitize data markers and error bars from vector figures in a PDF (pdfplumber).
Everything works in raw PDF page coordinates; a Panel maps them to data coordinates through the axes frame."""
import pdfplumber


def is_black(o):
    c = o.get('stroking_color') or o.get('non_stroking_color')
    if c is None:
        return False
    c = tuple(c)
    if len(c) == 4:
        return c[3] > 0.6                              # CMYK black (K high)
    if len(c) == 3:
        return max(c) < 0.15
    if len(c) == 1:
        return c[0] < 0.15
    return False


class Panel:
    """frame = (xa, xb, ya, yb) raw coordinates of the axes box.
    xaxis / yaxis: ('x'|'y', v_at_low_raw, v_at_high_raw) = which raw direction carries the data axis and its value
    at the low-raw-coordinate and high-raw-coordinate ends of the frame (so rotated / flipped axes are handled)."""

    def __init__(self, page, frame, xaxis, yaxis):
        self.pg, self.f = page, frame
        self.xaxis, self.yaxis = xaxis, yaxis

    def val(self, which, rx, ry):
        d, lo, hi = which
        xa, xb, ya, yb = self.f
        if d == 'x':
            return lo + (rx - xa) / (xb - xa) * (hi - lo)
        return lo + (ry - ya) / (yb - ya) * (hi - lo)

    def to_data(self, rx, ry):
        return self.val(self.xaxis, rx, ry), self.val(self.yaxis, rx, ry)

    def inside(self, o, pad=0.8):
        xa, xb, ya, yb = self.f
        return xa - pad < o['x0'] and o['x1'] < xb + pad and ya - pad < o['top'] and o['bottom'] < yb + pad

    def markers(self, size, tol=0.25, square=True, kinds=('curves',)):
        out = []
        for k in kinds:
            for o in getattr(self.pg, k):
                w, h = o['x1'] - o['x0'], o['bottom'] - o['top']
                if self.inside(o) and is_black(o) and abs(w - size) < tol and (not square or abs(w - h) < 0.35):
                    out.append(((o['x0'] + o['x1']) / 2, (o['top'] + o['bottom']) / 2))
        return sorted(set((round(a, 2), round(b, 2)) for a, b in out))

    def bars(self, cx, cy, tol=0.5, gap=3.2, minlen=0.8):
        """error bars through (cx, cy). Bars are often drawn as two segments (either side of the marker),
        so chain every black straight line on the same row/column whose end is within `gap` of the current extent.
        returns ((raw_x_lo, raw_x_hi) | None, (raw_y_lo, raw_y_hi) | None)"""
        H, V = [], []
        for o in self.pg.lines:
            if not (self.inside(o, pad=2.0) and is_black(o)):
                continue
            w, h = o['x1'] - o['x0'], o['bottom'] - o['top']
            if h < 0.15 and w >= minlen and abs((o['top'] + o['bottom']) / 2 - cy) < tol:
                H.append((o['x0'], o['x1']))
            if w < 0.15 and h >= minlen and abs((o['x0'] + o['x1']) / 2 - cx) < tol:
                V.append((o['top'], o['bottom']))

        def chain(segs, c):
            lo = hi = c
            used = True
            n = 0
            while used:
                used = False
                for a, b in segs:
                    if (a < lo and b >= lo - gap) or (b > hi and a <= hi + gap):
                        nlo, nhi = min(lo, a), max(hi, b)
                        if (nlo, nhi) != (lo, hi):
                            lo, hi, used = nlo, nhi, True
                            n += 1
            return (lo, hi) if n else None
        return chain(H, cx), chain(V, cy)

    def points(self, size, **kw):
        """-> list of dict(x, y, xlo, xhi, ylo, yhi) in data units (bars default to the point itself)"""
        res = []
        for cx, cy in self.markers(size, **kw):
            x, y = self.to_data(cx, cy)
            bh, bv = self.bars(cx, cy)
            r = dict(x=x, y=y, xlo=x, xhi=x, ylo=y, yhi=y)
            for seg, fixed in ((bh, 'h'), (bv, 'v')):
                if not seg:
                    continue
                if fixed == 'h':
                    p1, p2 = self.to_data(seg[0], cy), self.to_data(seg[1], cy)
                else:
                    p1, p2 = self.to_data(cx, seg[0]), self.to_data(cx, seg[1])
                if abs(p1[0] - p2[0]) > abs(p1[1] - p2[1]):       # the bar runs along the data x axis
                    r['xlo'], r['xhi'] = min(p1[0], p2[0]), max(p1[0], p2[0])
                else:
                    r['ylo'], r['yhi'] = min(p1[1], p2[1]), max(p1[1], p2[1])
            # error bars are symmetric about the marker: use the shorter side (guards against bars of neighbouring
            # markers on the same line being chained into this one)
            for lo, hi, v in (('xlo', 'xhi', 'x'), ('ylo', 'yhi', 'y')):
                if r[lo] != r[v] or r[hi] != r[v]:
                    hw = min(r[v] - r[lo], r[hi] - r[v]) if (r[lo] < r[v] < r[hi]) else max(r[v] - r[lo], r[hi] - r[v])
                    r[lo], r[hi] = r[v] - hw, r[v] + hw
            res.append(r)
        return res
