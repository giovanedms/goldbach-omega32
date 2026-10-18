"""Reproducible figures F1-F4 for main.tex.  Run:  python make_figs.py   (from figs/)
Data: ../../certificate/rigoroso/{cells_Q35_interval,cells_Q3_interval}.npz, data_Q3.json, data_Q35.json (read only);
      rlib.py of the certificate (imported read-only; ball arithmetic) for the analytic curve F3;
      F4_profiles_float.csv (float64, NOT certified; produced by F4_gen_profiles.py)."""
import sys, os, json, math
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RIG = os.path.abspath(os.path.join(HERE, "../../certificate/rigoroso"))
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#52514e"   # Q={3,5}, Q={3}, references
INK = "#222222"; GRID = "#dddddd"
plt.rcParams.update({
    "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 8.5,
    "axes.labelsize": 8.5, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": GRAY, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "axes.grid.axis": "y",
    "grid.color": GRID, "grid.linewidth": 0.5, "axes.axisbelow": True,
    "lines.linewidth": 1.6, "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
})
W = 16 / 2.54
LAB35, LAB3 = r"$Q=\{3,5\}$", r"$Q=\{3\}$"

def load(tag):
    d = np.load(f"{RIG}/cells_{tag}_interval.npz")
    return {k: d[k] for k in d.files}
C35, C3 = load("Q35"), load("Q3")
XMAXPLOT = 8.0

def binned(c, key, fn, nb=400, xmax=XMAXPLOT):
    """fn = np.min/np.max of key over cells in each of nb log10-bins; returns bin centre, value, and argext."""
    lz = 0.5 * (np.log10(c["u"]) + np.log10(c["v"]))
    m = lz <= xmax
    lz, y = lz[m], (c[key] if key == "cert" else np.exp(c[key]))[m]
    edges = np.linspace(lz.min(), xmax, nb + 1)
    ib = np.clip(np.searchsorted(edges, lz, side="right") - 1, 0, nb - 1)
    xs, ys = [], []
    for b in range(nb):
        s = ib == b
        if s.any():
            xs.append(lz[s][np.argmin(y[s]) if fn is np.min else np.argmax(y[s])]); ys.append(fn(y[s]))
    return np.array(xs), np.array(ys)

def finish(ax, ylab, xlim=(1.6, XMAXPLOT)):
    ax.set_xlim(*xlim); ax.set_xlabel(r"$\log_{10} z$"); ax.set_ylabel(ylab)
    ax.tick_params(length=3, width=0.6)

def save(fig, name):
    fig.savefig(os.path.join(HERE, name)); plt.close(fig)

# ---------------------------------------------------------------- F1
fig, ax = plt.subplots(figsize=(W, 6.2 / 2.54))
ax.axhline(0, color=GRAY, ls="--", lw=1.0, zorder=1)
ax.text(1.65, 0.12, r"$m=0$ (certificate fails below)", color=GRAY, fontsize=7.5, va="bottom")
mins = {}
for c, col, lab in [(C35, BLUE, LAB35), (C3, ORANGE, LAB3)]:
    x, y = binned(c, "cert", np.min)
    ax.plot(x, y, color=col, lw=1.8, label=lab, zorder=3)
    i = int(np.argmin(c["cert"][np.log10(c["v"]) <= XMAXPLOT]))
    zc = math.sqrt(c["u"][i] * c["v"][i]); mins[lab] = (math.log10(zc), c["cert"][i], zc)
    ax.plot(math.log10(zc), c["cert"][i], "o", ms=5, color=col, mec="white", mew=0.6, zorder=4)
(x1, y1, z1), (x2, y2, z2) = mins[LAB35], mins[LAB3]
ax.annotate(rf"min $+{y1:.3f}$ at $z\approx{z1/1e4:.2f}\cdot10^{{4}}$", (x1, y1), (4.55, 7.2), fontsize=8,
            arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.7), ha="left")
ax.annotate(rf"min $+{y2:.4f}$ at $z\approx{z2/1e4:.2f}\cdot10^{{4}}$", (x2, y2), (4.75, -0.95), fontsize=8,
            arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.7), ha="left", va="center")
ax.text(6.75, 5.7, LAB35, ha="center", fontsize=8.5); ax.text(7.3, 2.6, LAB3, ha="center", fontsize=8.5)
ax.text(2.6, 9.4, f"curves clipped above $+10$ (max $+{C35['cert'].max():.1f}$ at $z_0$)", fontsize=7.5, color=GRAY, va="center")
ax.set_ylim(-1.5, 10.2)
finish(ax, r"certified margin $m(z)$  (min per log-bin)")
ax.legend(loc="upper right", frameon=False, handlelength=1.8)
save(fig, "F1.pdf")

# ---------------------------------------------------------------- F2
fig, ax = plt.subplots(figsize=(W, 6.0 / 2.54))
KS = float(np.exp(C35["logK"].max()))
for v, txt in [(1.146, r"$K=1.146$ (Campbell, uniform)"),
               (KS, rf"$K={KS:.5f}$ (sup of the $Q=\{{3,5\}}$ cells: uniform fails)")]:
    ax.axhline(v, color=GRAY, ls="--", lw=1.0, zorder=1)
ax.text(7.95, 1.146 + 0.003, r"$K=1.146$ (Campbell, uniform)", color=GRAY, fontsize=7.5, ha="right", va="bottom")
ax.text(7.95, KS + 0.003, rf"$K={KS:.5f}$ (sup over the $Q=\{{3,5\}}$ cells)", color=GRAY, fontsize=7.5, ha="right", va="bottom")
for c, col, lab in [(C35, BLUE, LAB35), (C3, ORANGE, LAB3)]:
    x, y = binned(c, "logK", np.max)
    ax.plot(x, y, color=col, lw=1.8, label=lab, zorder=3)
ax.text(7.6, 1.0245, LAB35, ha="center", fontsize=8.5)
ax.text(7.6, 1.0475, LAB3, ha="center", fontsize=8.5)
ax.text(2.45, 1.116, "peak at $z\\approx43$", fontsize=7.5, color=INK, ha="left")
ax.set_ylim(0.995, 1.17)
finish(ax, r"per-cell constant $K_{\rm cell}(z)$  (max per log-bin)")
ax.legend(loc="center", bbox_to_anchor=(0.62, 0.39), ncol=2, frameon=False, handlelength=1.8, columnspacing=2.0)
save(fig, "F2.pdf")

# ---------------------------------------------------------------- F3 (analytic range, ball arithmetic)
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.path.insert(0, RIG)
from flint import arb
import rlib
from rlib import hi, lo, logF, logF_np, CQ

def analytic_curve(Q, Y, tmax):
    """M_lo(t), t>=log Y, with the constants of analytic.py at Z=Y and (s, alpha) fixed at the float optimum."""
    d = json.load(open(f"{RIG}/data_Q{''.join(map(str, Q))}.json"))["boundaries"][str(Y)]
    L_hi, lKex = d["L_hi"], d["logKexact_hi"]
    lK = rlib.logK_an_hi(Y, lKex)
    Ya = arb(Y); tY = Ya.log()
    B_up = arb("0.2614972128477"); sig = sum(arb(1) / q for q in Q)
    qX, c1 = (3, 1) if Q == (3,) else (5, 16)
    lcV = -arb(L_hi) + 2 * tY.log() - 2 / tY ** 2 - 1 / (Ya - 2) - 1 / (Ya - 1)
    k = 2 + lK; t0 = float(tY.mid())
    ss = np.linspace(2 * k + 3 + 1e-6, 16 - 10 / t0 - 1e-6, 20001); al = np.log(ss / 8)
    val = (16 * t0 + 0.5 * math.log(8) - math.log(qX) - L_hi - 8 * math.log(t0) + logF_np(ss, k)
           - math.log(2 * float(CQ[Q].mid())) - ss * t0 + al * ss - 8 * (0.2614972128 - float(sig.mid())) - 8 * (np.exp(al) - 1))
    j = int(np.argmax(val)); s, a = arb(float(ss[j])), arb(float(al[j]))
    lf = logF(s, arb(k))
    def M(t):
        t = arb(t); N = arb(8).sqrt() * (16 * t).exp()
        return (N.log() - arb(qX).log() + (1 - c1 / N).log() + lcV - 10 * t.log() + lf - (2 * CQ[Q]).log()
                - s * t - 0 + a * s - 8 * (B_up - sig) - 4 / t ** 2 - 8 * (a.exp() - 1))
    ts = np.linspace(float(tY.mid()), tmax, 160)
    ts[0] = float(tY.mid())
    return ts, np.array([lo(M(arb(Y).log() if i == 0 else float(t))) for i, t in enumerate(ts)])

LN10 = math.log(10); XA = 14.0; tmax = XA * LN10
fig, ax = plt.subplots(figsize=(W, 5.0 / 2.54))
ax.axhline(0, color=GRAY, ls="--", lw=1.0, zorder=1)
ts, m35 = analytic_curve((3, 5), 10**8, tmax)
ax.plot(ts / LN10, m35, color=BLUE, lw=1.8, zorder=3)
env = None; pts3 = {}
for Y in (10**8, 10**9, 10**10):
    tq, mq = analytic_curve((3,), Y, tmax)
    pts3[Y] = mq[0]
    xq = tq / LN10
    mi = np.interp(np.linspace(8, XA, 600), xq, mq, left=-np.inf)
    env = mi if env is None else np.maximum(env, mi)
xg = np.linspace(8, XA, 600)
ax.plot(xg, env, color=ORANGE, lw=1.8, zorder=3)
ax.plot(8, m35[0], "o", ms=5, color=BLUE, mec="white", mew=0.6, zorder=4)
for Y, v in pts3.items():
    ax.plot(math.log10(Y), v, "o", ms=5, color=ORANGE, mec="white", mew=0.6, zorder=4)
print("F3 check Q35@1e8: %.6f  Q3@1e8,1e9,1e10: %s" % (m35[0], ", ".join("%.6f" % v for v in pts3.values())))
ax.annotate(rf"$+{m35[0]:.3f}$ at $z=10^{{8}}$", (8, m35[0]), (8.35, 4.6), fontsize=8, ha="left",
            arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.7))
ax.annotate(rf"${pts3[10**8]:.3f}$ at $z=10^{{8}}$", (8, pts3[10**8]), (8.3, -5.4), fontsize=8, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.7))
ax.annotate(rf"$+{pts3[10**9]:.3f}$ at $10^{{9}}$", (9, pts3[10**9]), (9.4, -3.4), fontsize=8, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.7))
ax.annotate(rf"$+{pts3[10**10]:.3f}$ at $10^{{10}}$", (10, pts3[10**10]), (10.5, -1.5), fontsize=8, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.7))
i35 = np.argmin(abs(ts / LN10 - 13.0)); ax.text(12.0, m35[np.argmin(abs(ts/LN10-12.0))] + 2.6, LAB35, ha="center", fontsize=8.5)
ax.text(13.0, np.interp(13.0, xg, env) - 4.6, LAB3, ha="center", fontsize=8.5)
ax.set_ylim(-7.3, 25)
finish(ax, r"analytic margin $M_{\rm lo}$", xlim=(7.8, XA))
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], color=BLUE, lw=1.8, label=LAB35), Line2D([], [], color=ORANGE, lw=1.8, label=LAB3)],
          loc="upper left", frameon=False, handlelength=1.8)
save(fig, "F3.pdf")

# ---------------------------------------------------------------- F4 (float, not certified)
rows = [l.strip().split(",") for l in open(os.path.join(HERE, "F4_profiles_float.csv")).read().splitlines()[1:]]
P = {}
for c, z, m in rows:
    P.setdefault(c, ([], []))[0].append(math.log10(float(z))); P[c][1].append(float(m))
P = {k: (np.array(a), np.array(b)) for k, (a, b) in P.items()}
fig, ax = plt.subplots(figsize=(W, 6.0 / 2.54))
ax.axhline(0, color=GRAY, ls="--", lw=1.0, zorder=1)
for key, col, lab in [("M31_Q35", BLUE, LAB35), ("M31_Q3", ORANGE, LAB3)]:
    x, y = P[key]; ax.plot(x, y, color=col, lw=1.8, label=lab, zorder=3)
    i = int(np.argmin(y)); ax.plot(x[i], y[i], "o", ms=5, color=col, mec="white", mew=0.6, zorder=4)
    ax.annotate((lambda zz: rf"min ${y[i]:.2f}$ at $z\approx{zz/10**int(np.floor(np.log10(zz))):.2f}\cdot10^{{{int(np.floor(np.log10(zz)))}}}$")(10**x[i]), (x[i], y[i]),
                ((5.35, -1.7) if key == "M31_Q35" else (5.75, -6.6)), fontsize=8,
                ha="left", va="center", arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.7))
for key, col in [("M32_Q35", BLUE), ("M32_Q3", ORANGE)]:
    x, y = P[key]; ax.plot(x, y, color=col, lw=1.1, ls=(0, (4, 2)), zorder=2)
ax.text(3.7, 6.6, r"$\Omega(ab)\leq32$ (dashed, for comparison)", fontsize=7.5, color=GRAY, ha="left")
ax.text(1.7, -2.2, r"$\Omega(ab)\leq31$ (solid)", fontsize=8, ha="left")
ax.text(3.7, 8.9, "curves clipped above $+10$", fontsize=7.5, color=GRAY, va="center")
ax.set_ylim(-7.4, 10.2)
finish(ax, r"exploratory margin (float64)")
ax.legend(loc="lower left", frameon=False, handlelength=1.8)
save(fig, "F4.pdf")
print("ok")
