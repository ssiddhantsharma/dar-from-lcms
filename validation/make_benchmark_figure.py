"""Benchmark parity figure: this tool's DAR/CAR vs published/reference values.

Two independent public datasets:
  - van der Zon et al. 2026 (Zenodo 17637726): DFO immunoconjugate CAR, native SEC-MS.
    Reference = paper Table 2 / S6-7 / S11 averageCAR.
  - Fiala, Schuster & Heck 2025 (PRIDE PXD063887): cytotoxic ADC DAR, native vT-ESI.
    Reference = label/literature avDAR (BV ~4, EV ~3.8, T-DM1 ~3.5); the paper itself
    tabulates thermal-stability T-half, not avDAR, so these are the weaker label targets.

Numbers are this tool's measured values from validation/README.md. Run: python make_benchmark_figure.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# (label, reference, ours)
vdz = [
    ("N-suc 200mM", 2.04, 2.43), ("N-suc 400mM", 2.43, 2.08),
    ("N-suc 600mM", 2.48, 2.14), ("N-suc 1000mM", 2.06, 2.04),
    ("mal 600mM", 2.49, 2.76), ("NCS 600mM", 1.68, 2.17),
    ("TMTHSI 600mM", 0.79, 1.15), ("Ipilimumab N-suc", 1.93, 1.41),
]
adc = [
    ("Brentuximab vedotin", 4.0, 3.66), ("Enfortumab vedotin", 3.8, 3.25),
    ("Trastuzumab emtansine", 3.5, 2.70),
]

def mae(rows):
    return np.mean([abs(o - r) for _, r, o in rows])

fig, ax = plt.subplots(figsize=(6.4, 6.4))
lim = (0, 4.6)
ax.fill_between(lim, [lim[0] - 0.5, lim[1] - 0.5], [lim[0] + 0.5, lim[1] + 0.5],
                color="#9e9e9e", alpha=0.12, lw=0, label="+-0.5 band")
ax.plot(lim, lim, "--", color="#616161", lw=1, label="y = x")

for rows, color, marker, name in [(vdz, "#1f6feb", "o", "DFO immunoconjugates (van der Zon)"),
                                   (adc, "#e8710a", "s", "cytotoxic ADCs (PXD063887)")]:
    xs = [r for _, r, _ in rows]; ys = [o for _, _, o in rows]
    ax.scatter(xs, ys, c=color, marker=marker, s=55, edgecolor="white", lw=0.6, zorder=3,
               label=f"{name}, MAE {mae(rows):.2f}")
    for lab, r, o in rows:
        ax.annotate(lab, (r, o), xytext=(4, 3), textcoords="offset points",
                    fontsize=6.2, color="#333333")

ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
ax.set_xlabel("reference DAR / CAR (paper or label)")
ax.set_ylabel("this tool (average DAR / CAR)")
ax.set_title("dar-from-lcms vs reference, two public datasets", fontsize=11)
allrows = vdz + adc
ax.text(0.03, 0.97, "n = %d\noverall MAE = %.2f\n(all read slightly low:\nlow-DAR overweight)"
        % (len(allrows), mae(allrows)), transform=ax.transAxes, va="top", ha="left",
        fontsize=8, bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cfd8e3", lw=0.6))
ax.legend(loc="lower right", frameon=False, fontsize=7.5)
ax.grid(True, color="#9e9e9e", lw=0.2, alpha=0.5); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("benchmark_parity.png", dpi=300)
print("wrote benchmark_parity.png | overall MAE %.2f (n=%d)" % (mae(allrows), len(allrows)))
