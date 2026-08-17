"""Mirror (head-to-tail) demo: a conjugated sample on top vs its unmodified control below,
the standard way to show that conjugation shifted the mass. Synthetic data only. Calls the
tool's render() with MIRROR_MASS pointing at a synthetic control spectrum. The mirror layout
follows spectrum_utils' mirror() convention (Bittremieux et al., Apache-2.0), reimplemented.
Run: python examples/make_mirror_demo.py
"""
import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import dar_auto as da   # noqa: E402

BASE, STEP = 14000.0, 480.0


def g(x, mu, amp, s):
    return amp * np.exp(-((x - mu) ** 2) / (2 * s * s))


mm = np.arange(BASE - 400, BASE + STEP + 500, 1.0)
# treated: mostly +1 conjugate (plus a small +16); control: essentially all unmodified
treated_y = g(mm, BASE, 20, 3) + g(mm, BASE + STEP, 100, 3) + g(mm, BASE + STEP + 16, 24, 3)
treated = np.column_stack([mm, treated_y / treated_y.max() * 100])
control_y = g(mm, BASE, 100, 3) + g(mm, BASE + 16, 10, 3)
control = np.column_stack([mm, control_y / control_y.max() * 100])

mzx = np.linspace(600, 2200, 4000)
env = sum(g(mzx, (BASE + STEP + z * 1.00728) / z, 100, 1.6) + g(mzx, (BASE + z * 1.00728) / z, 25, 1.6)
          for z in range(6, 13))
mz = np.column_stack([mzx, env / env.max() * 100])
t = np.linspace(0, 20, 600)
tic = np.column_stack([t, g(t, 0.3, 50, 0.05) + g(t, 8.0, 100, 0.22)])
uv = np.column_stack([t, g(t, 0.3, 0.5, 0.05) + g(t, 8.0, 3.0, 0.22)])

here = os.path.dirname(os.path.abspath(__file__))
ctrl_path = os.path.join(tempfile.gettempdir(), "dar_mirror_control_mass.txt")
np.savetxt(ctrl_path, control)
os.environ.update({
    "FIG_TITLE": "Conjugated vs unmodified control (synthetic data)",
    "FIG_SUBTITLE": "head-to-tail: treated on top, control mirrored below; not real measurements",
    "MIRROR_MASS": ctrl_path, "FIG_SUFFIX": "",
})
res = da.render("mirror_demo", here, treated, mz, tic, uv, 7.65, 8.35, BASE, STEP, peak_width=0.5)
for ext in ("png", "pdf"):
    src = os.path.join(here, "dar_mirror_demo.%s" % ext)
    if os.path.exists(src):
        shutil.move(src, os.path.join(here, "mirror_demo.%s" % ext))
os.remove(ctrl_path)
print("wrote %s  (treated DAR=%.2f)" % (os.path.join(here, "mirror_demo.png"), res["DAR"]))
