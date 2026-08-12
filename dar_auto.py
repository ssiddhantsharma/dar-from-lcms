import sys, os, json, base64, zlib
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib; matplotlib.use("Agg")
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
from unidec import engine
from unidec.UniDecImporter.ImporterFactory import ImporterFactory

_b, _s = os.environ.get("BASE_MASS"), os.environ.get("MOD_MASS")
if not _b or not _s:
    sys.exit("set BASE_MASS (unmodified average mass, Da) and MOD_MASS (mass added per conjugation, Da)")
base, step = float(_b), float(_s)
mlo = float(os.environ.get("MASS_LB", base - 2000))
mhi = float(os.environ.get("MASS_UB", base + 3000))
zlo, zhi = int(os.environ.get("Z_LO", 3)), int(os.environ.get("Z_HI", 20))
void = float(os.environ.get("VOID_MIN", 2.0))
trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
_tag = lambda el: el.tag.split("}")[-1]


def chromatograms(path):
    # decode mzML <chromatogram> arrays with stdlib (base64 + optional zlib)
    out = {}
    for _, el in ET.iterparse(path, events=("end",)):
        t = _tag(el)
        if t == "chromatogram":
            cid, a = el.get("id", ""), {}
            for bda in el.iter():
                if _tag(bda) != "binaryDataArray":
                    continue
                dt, comp, kind, blob = np.float64, False, None, None
                for c in bda:
                    if _tag(c) == "cvParam":
                        nm = c.get("name", "")
                        if nm == "32-bit float": dt = np.float32
                        elif nm == "64-bit float": dt = np.float64
                        elif "zlib" in nm: comp = True
                        elif nm == "time array": kind = "t"
                        elif nm == "intensity array": kind = "i"
                    elif _tag(c) == "binary": blob = c.text
                if kind and blob:
                    raw = base64.b64decode(blob)
                    a[kind] = np.frombuffer(zlib.decompress(raw) if comp else raw, dt)
            if "t" in a and "i" in a:
                n = min(len(a["t"]), len(a["i"]))
                out[cid] = np.column_stack([a["t"][:n], a["i"][:n]])
            el.clear()
        elif t == "spectrum":
            el.clear()
    return out


def area(md, lo, hi):
    m = (md[:, 0] >= lo) & (md[:, 0] <= hi)
    return float(trap(md[m, 1], md[m, 0])) if m.any() else 0.0


def apex(md, target, w=15):
    m = (md[:, 0] >= target - w) & (md[:, 0] <= target + w)
    if not m.any():
        return None, 0.0
    s = md[m]
    k = int(np.argmax(s[:, 1]))
    return float(s[k, 0]), float(s[k, 1])


def charge_peaks(mz, mass, nmax=6, frac=0.15):
    # local maxima in the m/z envelope, one label per charge state (the tallest peak for that z)
    x, y = mz[:, 0], mz[:, 1]
    thr = frac * y.max()
    idx = [k for k in range(2, len(y) - 2) if y[k] > thr and y[k] == max(y[k - 2:k + 3])]
    best = {}
    for k in idx:
        z = int(round(mass / (x[k] - 1.00728)))
        if z > 0 and (z not in best or y[k] > y[best[z]]):
            best[z] = k
    keep = sorted(best.values(), key=lambda k: -y[k])[:nmax]
    return [(x[k], y[k], int(round(mass / (x[k] - 1.00728)))) for k in keep]


def analyze(path, out):
    name = os.path.splitext(os.path.basename(path))[0]
    tag = name.replace(" ", "_")
    ch = chromatograms(path)
    imp = ImporterFactory.create_importer(path)

    def pick(*subs):   # first chromatogram whose id contains any of subs (case-insensitive)
        for s in subs:
            for k in ch:
                if s.lower() in k.lower():
                    return ch[k]
        return None
    tic = pick("TIC")
    uv = pick("Sig=280", "280", "uv", "absorb")   # DAD/UV if present (any vendor); else None
    if tic is None:                                  # some vendors omit a TIC chromatogram in mzML
        tic = np.asarray(imp.get_tic())

    if os.environ.get("TMIN") and os.environ.get("TMAX"):
        tmin, tmax = float(os.environ["TMIN"]), float(os.environ["TMAX"])
    else:
        # widest span around the biggest TIC peak (after the void) that stays above half-max
        t, v = tic.T
        keep = t >= void; t, v = t[keep], v[keep]
        k = int(np.argmax(v)); thr = v[k] / 2
        lo = hi = k
        while lo and v[lo - 1] >= thr: lo -= 1
        while hi < len(v) - 1 and v[hi + 1] >= thr: hi += 1
        tmin, tmax = float(t[lo]), float(t[hi])

    spec = os.path.join(out, "_avg_%s.txt" % tag)
    np.savetxt(spec, np.asarray(imp.get_avg_scan(time_range=(tmin, tmax))))
    e = engine.UniDec(); e.open_file(spec)
    e.config.startz, e.config.endz = zlo, zhi
    e.config.masslb, e.config.massub = mlo, mhi
    e.config.massbins = 1.0
    e.config.minmz, e.config.maxmz = 600.0, 2200.0
    e.process_data(); e.run_unidec(silent=True)
    md = np.asarray(e.data.massdat)
    if md.ndim != 2 or not len(md):
        raise RuntimeError("empty deconvolution (massdat %s)" % (md.shape,))
    md[:, 1] *= 100.0 / md[:, 1].max()

    # DAR with masses ANCHORED at base / base+step; report robustness across integration widths
    table = []
    for w in (15, 25, 40, 60):
        a0, a1 = area(md, base - w, base + w), area(md, base + step - w, base + step + w)
        table.append({"halfwidth_Da": w, "naked_area": round(a0, 1),
                      "conj_area": round(a1, 1), "DAR": round(a1 / (a0 + a1) if a0 + a1 else 0, 3)})
    dar = table[1]["DAR"]   # headline = +/-25 Da band
    mz = np.asarray(e.data.data2)   # processed m/z envelope (the input to deconvolution)

    def uv_h(arr, a, b):
        if arr is None: return None
        m = (arr[:, 0] >= a) & (arr[:, 0] <= b)
        return round(float(arr[m, 1].max()), 3) if m.any() else 0.0
    uv_main, uv_late = uv_h(uv, tmin - .3, tmax + .3), uv_h(uv, 18, 22)

    # report a measured apex only for a species that is actually present (>=10%);
    # a trace species gets None instead of a peak-picker noise bump near the anchor.
    nm, _ = apex(md, base)
    cm, _ = apex(md, base + step)
    naked_apex = round(nm, 1) if (nm is not None and (1 - dar) >= 0.10) else None
    conj_apex = round(cm, 1) if (cm is not None and dar >= 0.10) else None
    trace = ("unmodified <10%, not separately resolved" if (1 - dar) < 0.10
             else "conjugate <10%, not separately resolved" if dar < 0.10 else None)

    # ---- figure: (a) chromatograms, (b) raw m/z envelope, (c) deconvolved mass ----
    prim, sec, shade, band = "#1f4e79", "#9aa0a6", "#f6dfae", "#cfd8e3"
    fig, (axA, axB, axC) = plt.subplots(3, 1, figsize=(7.2, 10.2))

    # (a) chromatograms
    if tic is not None:
        axA.plot(tic[:, 0], tic[:, 1] / tic[:, 1].max() * 100, color=sec, lw=0.9, label="MS TIC")
    if uv is not None:
        axA.plot(uv[:, 0], uv[:, 1] / np.abs(uv[:, 1]).max() * 100,
                 color=prim, lw=1.4, label="UV 280 nm")
    axA.axvspan(tmin, tmax, color=shade, alpha=0.7, lw=0)
    axA.set_xlim(0, max(25, tmax + 3)); axA.set_ylim(-2, 108)
    axA.set_xlabel("retention time (min)"); axA.set_ylabel("relative signal (%)")
    axA.legend(loc="upper right", frameon=False, handlelength=1.4)
    axA.spines[["top", "right"]].set_visible(False)

    # (b) raw charge-state envelope UniDec actually deconvolved
    axB.plot(mz[:, 0], mz[:, 1] / mz[:, 1].max() * 100, color=prim, lw=0.8)
    dom = base + step if dar >= 0.5 else base
    for x0, y0, z in charge_peaks(mz, dom):
        axB.annotate("%d+" % z, xy=(x0, y0 / mz[:, 1].max() * 100),
                     xytext=(0, 3), textcoords="offset points", ha="center",
                     fontsize=7.5, color="#666666")
    axB.set_xlim(mz[:, 0].min(), mz[:, 0].max()); axB.set_ylim(0, 112)
    axB.set_xlabel("m/z"); axB.set_ylabel("relative intensity (%)")
    axB.set_title("raw charge-state envelope (input to deconvolution)",
                  loc="left", fontsize=9.5, color="#666666")
    axB.spines[["top", "right"]].set_visible(False)

    # (c) deconvolved mass; DAR integration bands shaded, adducts labelled on the dominant species
    ymax = 120
    for lo_, hi_ in [(base - 25, base + 25), (base + step - 25, base + step + 25)]:
        axC.axvspan(lo_, hi_, color=band, alpha=0.6, lw=0)
    axC.fill_between(md[:, 0], md[:, 1], color=prim, alpha=0.18, lw=0)
    axC.plot(md[:, 0], md[:, 1], color=prim, lw=0.9)
    for m, lab in [(base, "unmodified"), (base + step, "+1")]:
        h = float(md[np.argmin(abs(md[:, 0] - m)), 1])
        axC.annotate("%s\n%.0f Da" % (lab, m), xy=(m, h), xytext=(m, min(h + 15, ymax - 5)),
                     ha="center", va="bottom", fontsize=9, color="#222222",
                     arrowprops=dict(arrowstyle="-", lw=0.7, color="#999999", shrinkA=0, shrinkB=1))
    for da, lab in [(16, "+16"), (178, "+178")]:
        am = dom + da
        if md[0, 0] <= am <= md[-1, 0]:
            h = float(md[np.argmin(abs(md[:, 0] - am)), 1])
            if h > 4:
                axC.annotate(lab, xy=(am, h), xytext=(0, 3), textcoords="offset points",
                             ha="center", fontsize=7, color="#999999")
    axC.set_xlim(base - 400, base + step + 500); axC.set_ylim(0, ymax)
    axC.set_xlabel("deconvolved mass (Da)"); axC.set_ylabel("relative abundance (%)")
    axC.text(0.015, 0.96, "DAR = %.2f" % dar, transform=axC.transAxes,
             va="top", ha="left", fontsize=11, fontweight="bold", color=prim)
    axC.spines[["top", "right"]].set_visible(False)

    for ax, letter in [(axA, "a"), (axB, "b"), (axC, "c")]:
        ax.text(-0.09, 1.02, letter, transform=ax.transAxes, fontsize=13, fontweight="bold")
    fig.suptitle(name, x=0.5, y=0.997, fontsize=10, color="#666666")
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    figpath = os.path.join(out, "dar_%s.png" % tag)
    fig.savefig(figpath, dpi=300)
    fig.savefig(figpath[:-4] + ".pdf")   # vector, for reports
    plt.close(fig)
    os.remove(spec)

    return {"file": os.path.basename(path), "window_min": [round(tmin, 2), round(tmax, 2)],
            "DAR": dar, "conjugated_pct": round(dar * 100, 1),
            "conj_apex": conj_apex, "naked_apex": naked_apex, "trace": trace,
            "expected": [base, round(base + step, 2)],
            "DAR_by_window": table,
            "UV280_main_h": uv_main, "UV280_late_h": uv_late,
            "UV_late_over_main": (round(uv_late / uv_main, 2) if uv_main else None)}


if __name__ == "__main__":
    out = os.environ.get("OUTDIR", "/data")
    res = []
    for p in sys.argv[1:]:
        try:
            res.append(analyze(p, out))
        except Exception as ex:
            res.append({"file": os.path.basename(p), "error": str(ex)})
    json.dump(res, open(os.path.join(out, "dar_results.json"), "w"), indent=2)
    print(json.dumps(res, indent=2))
