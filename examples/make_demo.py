"""Render the synthetic example figure for the README by calling the tool's own
render(), with made-up data only (no real measurements). Because it uses the real
render(), the demo always matches what the tool emits -- band areas, the optional
concentration panel, titles and all. Run: python examples/make_demo.py
"""
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import dar_auto as da   # noqa: E402

BASE, STEP = 14000.0, 480.0     # placeholders; not real
rng = np.random.default_rng(0)


def g(x, mu, amp, s):
    return amp * np.exp(-((x - mu) ** 2) / (2 * s * s))


# (a) chromatograms: void spike + protein peak + a late matrix feature. UV in a
#     realistic mAU scale so the demo concentration comes out sane.
t = np.linspace(0, 20, 600)
tic = np.column_stack([t, g(t, 0.3, 60, 0.05) + g(t, 8.0, 100, 0.22) + g(t, 15.5, 22, 1.6)])
uv = np.column_stack([t, g(t, 0.3, 0.7, 0.05) + g(t, 8.0, 3.0, 0.22)])   # ~3 mAU peak
tmin, tmax = 7.65, 8.35

# (b) charge-state envelope: conjugate dominant, unmodified minor, z = 6..12
mzx = np.linspace(600, 2200, 4000)
env = np.zeros_like(mzx)
for z in range(6, 13):
    env += g(mzx, (BASE + STEP + z * 1.00728) / z, 100, 1.6)     # conjugate
    env += g(mzx, (BASE + z * 1.00728) / z, 28, 1.6)             # unmodified
env = (env + rng.normal(0, 1.0, env.size)).clip(0)
mz = np.column_stack([mzx, env / env.max() * 100])

# (c) deconvolved mass: unmodified, +1, and a small +16 adduct
mm = np.arange(BASE - 400, BASE + STEP + 500, 1.0)
mdv = g(mm, BASE, 28, 3) + g(mm, BASE + STEP, 100, 3) + g(mm, BASE + STEP + 16, 24, 3)
md = np.column_stack([mm, mdv / mdv.max() * 100])

here = os.path.dirname(os.path.abspath(__file__))
os.environ.update({
    "FIG_TITLE": "DAR demo (synthetic data)",
    "FIG_SUBTITLE": "illustrative example; not real measurements",
    "CONJ_LABEL": "+1",
    "FIG_SUFFIX": "",
    # optional concentration panel, with made-up but plausible inputs
    "EPS280": "20000", "FLOW_ML_MIN": "0.4", "INJ_UL": "5",
    "DAD_UNIT": "mAU", "DILUTION": "1", "MW_DA": str(BASE + STEP),
})
res = da.render("demo", here, md, mz, tic, uv, tmin, tmax, BASE, STEP, peak_width=0.5)
for ext in ("png", "pdf"):                       # render writes dar_demo.<ext> -> demo.<ext>
    src = os.path.join(here, "dar_demo.%s" % ext)
    if os.path.exists(src):
        shutil.move(src, os.path.join(here, "demo.%s" % ext))
print("wrote %s  (DAR=%.2f, conc=%s mg/mL)"
      % (os.path.join(here, "demo.png"), res["DAR"], res.get("protein_conc_mg_ml")))
