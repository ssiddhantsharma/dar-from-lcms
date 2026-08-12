import sys, os, json, base64, zlib
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
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


def apex(md, target, w=40):
    m = (md[:, 0] >= target - w) & (md[:, 0] <= target + w)
    return float(md[m][np.argmax(md[m, 1]), 0]) if m.any() else None


def analyze(path, out):
    name = os.path.splitext(os.path.basename(path))[0]
    tag = name.replace(" ", "_")
    ch = chromatograms(path)
    pick = lambda s: next((ch[k] for k in ch if s in k), None)
    tic, uv280 = pick("TIC"), pick("Sig=280")

    if os.environ.get("TMIN") and os.environ.get("TMAX"):
        tmin, tmax = float(os.environ["TMIN"]), float(os.environ["TMAX"])
    else:
        # reuse the TIC already parsed above; find the widest span above half-max after the void
        t, v = tic.T
        keep = t >= void; t, v = t[keep], v[keep]
        k = int(np.argmax(v)); thr = v[k] / 2
        lo = hi = k
        while lo and v[lo - 1] >= thr: lo -= 1
        while hi < len(v) - 1 and v[hi + 1] >= thr: hi += 1
        tmin, tmax = float(t[lo]), float(t[hi])

    imp = ImporterFactory.create_importer(path)
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

    def uv_h(arr, a, b):
        if arr is None: return None
        m = (arr[:, 0] >= a) & (arr[:, 0] <= b)
        return round(float(arr[m, 1].max()), 3) if m.any() else 0.0
    uv_main, uv_late = uv_h(uv280, tmin - .3, tmax + .3), uv_h(uv280, 18, 22)

    # ---- 2-panel figure ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8))
    if tic is not None:
        ax1.plot(tic[:, 0], tic[:, 1] / tic[:, 1].max() * 100, color="#aaa", lw=1, label="MS TIC")
    if uv280 is not None:
        ax1.plot(uv280[:, 0], uv280[:, 1] / np.abs(uv280[:, 1]).max() * 100,
                 color="#1f4e79", lw=1.2, label="UV 280 nm")
    ax1.axvspan(tmin, tmax, color="#f0a500", alpha=.25)
    ax1.set_xlim(0, max(25, tmax + 3)); ax1.set_ylim(-3, 105)
    ax1.set_xlabel("retention time (min)"); ax1.set_ylabel("relative signal (%)")
    ax1.set_title("%s   protein window %.1f-%.1f min" % (name, tmin, tmax),
                  fontsize=12, fontweight="bold")
    ax1.legend(fontsize=9, frameon=False); ax1.spines[["top", "right"]].set_visible(False)

    ax2.fill_between(md[:, 0], md[:, 1], color="#1f4e79", lw=0, alpha=.88)
    for m, lab in [(base, "unconjugated"), (base + step, "+1")]:
        h = md[np.argmin(abs(md[:, 0] - m)), 1]
        ax2.annotate("%s\n%.0f Da" % (lab, m), xy=(m, h), xytext=(m, min(h + 20, 112)),
                     ha="center", fontsize=9.5, fontweight="bold",
                     arrowprops=dict(arrowstyle="-", color="#888", lw=.8))
    ax2.set_xlim(base - 400, base + step + 500); ax2.set_ylim(0, 125)
    ax2.set_xlabel("deconvolved mass (Da)"); ax2.set_ylabel("relative abundance (%)")
    ax2.set_title("DAR %.2f   ·   %.0f%% conjugated" % (dar, dar * 100), fontsize=12)
    ax2.spines[["top", "right"]].set_visible(False); ax2.grid(axis="y", ls=":", alpha=.4)
    fig.tight_layout(); fig.savefig(os.path.join(out, "dar_%s.png" % tag), dpi=150); plt.close(fig)
    os.remove(spec)

    return {"file": os.path.basename(path), "window_min": [round(tmin, 2), round(tmax, 2)],
            "DAR": dar, "conjugated_pct": round(dar * 100, 1),
            "conj_apex": apex(md, base + step), "naked_apex": apex(md, base),
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
