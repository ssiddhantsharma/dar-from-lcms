"""Render a synthetic example figure for the README.

Uses made-up numbers only (no real data): placeholder mass M and modifier mass Δ.
Produces examples/demo.png in the same 3-panel style the tool emits, so newcomers can
see the output without any input file. Run: python examples/make_demo.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10.5,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "axes.linewidth": 0.8, "axes.edgecolor": "#333333",
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "text.color": "#222222", "axes.labelcolor": "#222222",
    "xtick.color": "#333333", "ytick.color": "#333333",
})

BASE, STEP = 14000.0, 480.0     # placeholders; not real
prim, sec, shade, band = "#1f4e79", "#9aa0a6", "#f6dfae", "#cfd8e3"
rng = np.random.default_rng(0)


def g(x, mu, amp, s):
    return amp * np.exp(-((x - mu) ** 2) / (2 * s * s))


# (a) synthetic chromatograms: void spike + protein peak + a late matrix feature
t = np.linspace(0, 20, 600)
tic = g(t, 0.3, 60, 0.05) + g(t, 8.0, 100, 0.22) + g(t, 15.5, 22, 1.6)
uv = g(t, 0.3, 25, 0.05) + g(t, 8.0, 100, 0.22)
tmin, tmax = 7.65, 8.35

# (b) synthetic charge-state envelope (conjugate dominant, unmodified minor), z = 6..12
mz = np.linspace(600, 2200, 4000)
env = np.zeros_like(mz)
for z in range(6, 13):
    env += g(mz, (BASE + STEP + z * 1.00728) / z, 100, 1.6)      # conjugate
    env += g(mz, (BASE + z * 1.00728) / z, 28, 1.6)              # unmodified
env = (env + rng.normal(0, 1.0, env.size)).clip(0)
env = env / env.max() * 100

# (c) synthetic deconvolved mass: unmodified, +1, and a small +16 adduct
m = np.arange(BASE - 400, BASE + STEP + 500, 1.0)
md = g(m, BASE, 28, 3) + g(m, BASE + STEP, 100, 3) + g(m, BASE + STEP + 16, 24, 3)
md = md / md.max() * 100
trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def area(lo, hi):
    s = (m >= lo) & (m <= hi)
    return float(trap(md[s], m[s]))


a0, a1 = area(BASE - 25, BASE + 25), area(BASE + STEP - 25, BASE + STEP + 25)
dar = a1 / (a0 + a1)

fig, (axA, axB, axC) = plt.subplots(3, 1, figsize=(7.2, 10.2))

axA.plot(t, tic, color=sec, lw=0.9, label="MS TIC")
axA.plot(t, uv, color=prim, lw=1.4, label="UV 280 nm")
axA.axvspan(tmin, tmax, color=shade, alpha=0.7, lw=0)
axA.set_xlim(0, 20); axA.set_ylim(-2, 108)
axA.set_xlabel("retention time (min)"); axA.set_ylabel("relative signal (%)")
axA.legend(loc="upper right", frameon=False, handlelength=1.4)
axA.spines[["top", "right"]].set_visible(False)

axB.plot(mz, env, color=prim, lw=0.8)
for z in range(6, 13):
    c = (BASE + STEP + z * 1.00728) / z
    axB.annotate("%d+" % z, xy=(c, env[np.argmin(abs(mz - c))]),
                 xytext=(0, 3), textcoords="offset points", ha="center",
                 fontsize=7.5, color="#666666")
axB.set_xlim(600, 2200); axB.set_ylim(0, 112)
axB.set_xlabel("m/z"); axB.set_ylabel("relative intensity (%)")
axB.set_title("raw charge-state envelope (input to deconvolution)",
              loc="left", fontsize=9.5, color="#666666")
axB.spines[["top", "right"]].set_visible(False)

ymax = 120
for lo_, hi_ in [(BASE - 25, BASE + 25), (BASE + STEP - 25, BASE + STEP + 25)]:
    axC.axvspan(lo_, hi_, color=band, alpha=0.6, lw=0)
axC.fill_between(m, md, color=prim, alpha=0.18, lw=0)
axC.plot(m, md, color=prim, lw=0.9)
for mm, lab in [(BASE, "unmodified"), (BASE + STEP, "+1")]:
    h = float(md[np.argmin(abs(m - mm))])
    axC.annotate("%s\n%.0f Da" % (lab, mm), xy=(mm, h), xytext=(mm, min(h + 15, ymax - 5)),
                 ha="center", va="bottom", fontsize=9, color="#222222",
                 arrowprops=dict(arrowstyle="-", lw=0.7, color="#999999", shrinkA=0, shrinkB=1))
am = BASE + STEP + 16
axC.annotate("+16", xy=(am, float(md[np.argmin(abs(m - am))])), xytext=(0, 3),
             textcoords="offset points", ha="center", fontsize=7, color="#999999")
axC.set_xlim(BASE - 400, BASE + STEP + 500); axC.set_ylim(0, ymax)
axC.set_xlabel("deconvolved mass (Da)"); axC.set_ylabel("relative abundance (%)")
axC.text(0.015, 0.96, "DAR = %.2f" % dar, transform=axC.transAxes,
         va="top", ha="left", fontsize=11, fontweight="bold", color=prim)
axC.spines[["top", "right"]].set_visible(False)

for ax, letter in [(axA, "a"), (axB, "b"), (axC, "c")]:
    ax.text(-0.09, 1.02, letter, transform=ax.transAxes, fontsize=13, fontweight="bold")
fig.suptitle("DAR demo (synthetic data)", x=0.5, y=0.998, fontsize=12.5,
             fontweight="bold", color="#222222")
fig.text(0.5, 0.968, "illustrative example; not real measurements", ha="center",
         fontsize=9.5, color="#555555")
fig.tight_layout(rect=(0, 0, 1, 0.955))
here = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(here, "demo.png"), dpi=150)
plt.close(fig)
print("wrote", os.path.join(here, "demo.png"))
