"""Recolour representation1.eps from a shared RdBu stimulus ramp to neur5's
per-context hues (orangered / gray / deepskyblue), leaving geometry untouched.

Illustrator sets a colour once and reuses it across many fills, so a fill's
context can't be read from its colour alone -- it's recovered from which of the
three parallel arcs the fill sits on.
"""
import re, sys
import numpy as np
from collections import Counter
from scipy.cluster.hierarchy import fcluster, linkage
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/nuttidalab/Documents/asymmetry_design/code/python')
from utils import stim_cmap

SRC = '/home/nuttidalab/Documents/asymmetry_design/outputs/figs/ccn/poster/representation1.eps'
DST = '/home/nuttidalab/Documents/asymmetry_design/outputs/figs/ccn/poster/representation1_recolored.eps'
MIDX, MIDY = 173.0, 162.9          # panel split, in the EPS's y-down coordinates

colre = re.compile(r'^([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+) cmyk$')
ptre  = re.compile(r'(-?[\d.]+) (-?[\d.]+) (?:mo|li|cv)')

raw = open(SRC, 'rb').read()
cut = raw.find(b'%AI9_PrivateDataBegin')   # AI's private copy is left alone
head, tail = raw[:cut].decode('latin-1'), raw[cut:]
lines = head.split('\n')

cur, buf, fills = None, [], []
for n, ln in enumerate(lines):
    m = colre.match(ln)
    if m: cur = ' '.join(m.groups()); continue
    if ln == 'f':
        p = np.array([[float(a), float(b)] for a, b in ptre.findall('\n'.join(buf))])
        fills.append((n, cur, p.mean(0) if len(p) else None)); buf = []; continue
    buf.append(ln)
    if len(buf) > 400: buf = buf[-400:]

# the 20 ramp colours are the only ones reused heavily; axes and text are not
cnt  = Counter(c for _, c, _ in fills)
ramp = [c for c, k in cnt.items() if k >= 10 and c.split()[3] != '1']

def cmyk2rgb(s):
    c, m, y, k = [float(v) for v in s.split()]
    return np.array([(1-c)*(1-k), (1-m)*(1-k), (1-y)*(1-k)])

# a colour's position on RdBu is its stimulus level, kept as-is so the new ramps
# stay on the shared global norm rather than restretching per context
grid = np.linspace(0, 1, 512)
rdbu = plt.cm.RdBu(grid)[:, :3]
tpos = {c: grid[np.argmin(((rdbu - cmyk2rgb(c))**2).sum(1))] for c in ramp}

sel   = [(n, c, p) for n, c, p in fills if c in tpos and p is not None]
P     = np.array([p for _, _, p in sel])
panel = (P[:, 0] > MIDX).astype(int) + 2*(P[:, 1] > MIDY).astype(int)

assign = np.full(len(P), -1)
for pn in (0, 1, 2):                      # panel 3's 3D cloud is a raster, not fills
    s = np.where(panel == pn)[0]
    Z = linkage(P[s], 'single')
    # each arc is a connected chain of fills, and the gaps between arcs are wider
    # than the gaps within one, so grow the link distance until exactly 3 survive
    for thr in np.arange(5, 40, 0.5):
        lab = fcluster(Z, t=thr, criterion='distance')
        sz = Counter(lab.tolist())
        big = [k for k, v in sz.items() if v >= 50]
        if len(big) == 3 and sum(sz[k] for k in big) / len(lab) > 0.95: break
    else:
        raise SystemExit(f'panel {pn}: no clean 3-arc split')

    cores = {k: P[s][lab == k].mean(0) for k in big}
    # square arc sits leftmost, then circles, then triangles -- matching neur5's
    # context order, so rank the arcs by x to pin the hue assignment
    order = sorted(big, key=lambda k: cores[k][0])
    for j, i in enumerate(s):
        k = lab[j]
        if k not in big:      # a stray fragment joins whichever arc it sits nearest
            k = min(big, key=lambda b: np.hypot(*(P[i] - cores[b])))
        assign[i] = order.index(k)

cmaps = [stim_cmap(c) for c in ('orangered', 'gray', 'deepskyblue')]

def rgb2cmyk(r, g, b):
    k = 1 - max(r, g, b)
    if k >= 1 - 1e-9: return (0, 0, 0, 1)
    return ((1-r-k)/(1-k), (1-g-k)/(1-k), (1-b-k)/(1-k), k)

newcol = {}
for (n, c, _), a in zip(sel, assign):
    if a < 0: continue
    r, g, b, _ = cmaps[a](tpos[c])
    newcol[n] = ' '.join(f'{v:.6g}' for v in rgb2cmyk(r, g, b)) + ' cmyk'

# emit a colour immediately before each fill, overriding the shared one upstream
out = []
for n, ln in enumerate(lines):
    if n in newcol: out.append(newcol[n])
    out.append(ln)
open(DST, 'wb').write('\n'.join(out).encode('latin-1') + tail)

for pn in (0, 1, 2):
    for g in range(3):
        m = (panel == pn) & (assign == g)
        ts = [tpos[c] for (_, c, _), k in zip(sel, m) if k]
        print(f'panel {pn} arc {g}: {m.sum():4d} fills, t {min(ts):.2f}-{max(ts):.2f}')
print('recoloured', len(newcol), 'of', len(sel), 'ramp fills ->', DST)
