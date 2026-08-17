"""Render the multi-state DAR/CAR example figure for the README, by calling the tool's
own render() with DAR_MAX_N set and made-up data only (no real measurements). Shows the
multi-state ladder + average DAR + dispersity that the tool emits for a native-MS-style
conjugate with several load states. Run: python examples/make_multidar_demo.py
"""
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import dar_auto as da   # noqa: E402

# placeholders only: a native-MS-scale IgG-like conjugate, chelator step ~ mal-DFO + Fe
BASE, STEP, NMAX = 148000.0, 770.0, 6
rng = np.random.default_rng(0)


def g(x, mu, amp, s):
    return amp * np.exp(-((x - mu) ** 2) / (2 * s * s))


# (c) deconvolved mass: an even-dominant CAR0..CAR6 ladder with glycoform-ish broadening
mm = np.arange(BASE - 1200, BASE + NMAX * STEP + 1200, 2.0)
amps = [14, 7, 100, 11, 58, 6, 22]                      # CAR0..CAR6, even states dominant
mdv = sum(g(mm, BASE + n * STEP, a, 130) for n, a in enumerate(amps))
md = np.column_stack([mm, mdv / mdv.max() * 100])

# (b) native charge-state envelope, ~24..32+ over m/z 4,400-6,400
mzx = np.linspace(4400, 6400, 4000)
env = sum(g(mzx, (BASE + 2.5 * STEP + z * 1.00728) / z, 100 * np.exp(-((z - 28) ** 2) / 18), 4)
          for z in range(24, 33))
mz = np.column_stack([mzx, env / env.max() * 100])

# (a) chromatograms
t = np.linspace(0, 20, 600)
tic = np.column_stack([t, g(t, 0.3, 50, 0.05) + g(t, 13.5, 100, 0.35)])
uv = np.column_stack([t, g(t, 0.3, 0.5, 0.05) + g(t, 13.5, 3.0, 0.35)])

here = os.path.dirname(os.path.abspath(__file__))
os.environ.update({
    "DAR_MAX_N": str(NMAX),
    "FIG_TITLE": "Multi-state DAR/CAR demo (synthetic data)",
    "FIG_SUBTITLE": "illustrative even-dominant CAR0-6 ladder; not real measurements",
    "FIG_SUFFIX": "",
})
res = da.render("multidar_demo", here, md, mz, tic, uv, 12.9, 14.2, BASE, STEP, peak_width=0.5)
for ext in ("png", "pdf"):                              # render writes dar_multidar_demo.<ext>
    src = os.path.join(here, "dar_multidar_demo.%s" % ext)
    if os.path.exists(src):
        shutil.move(src, os.path.join(here, "multidar_demo.%s" % ext))
print("wrote %s  (average DAR=%.2f, dispersity=%.2f)"
      % (os.path.join(here, "multidar_demo.png"), res["average_dar"], res["dispersity"]))
