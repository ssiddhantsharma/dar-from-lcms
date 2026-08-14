import sys, os, json, base64, zlib
import xml.etree.ElementTree as ET
import numpy as np

trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
_tag = lambda el: el.tag.split("}")[-1]


def load_config():
    b, s = os.environ.get("BASE_MASS"), os.environ.get("MOD_MASS")
    if not b or not s:
        raise SystemExit("set BASE_MASS (unmodified average mass, Da) and "
                         "MOD_MASS (mass added per conjugation, Da)")
    base, step = float(b), float(s)
    return dict(base=base, step=step,
                mlo=float(os.environ.get("MASS_LB", base - 2000)),
                mhi=float(os.environ.get("MASS_UB", base + 3000)),
                zlo=int(os.environ.get("Z_LO", 3)), zhi=int(os.environ.get("Z_HI", 20)),
                void=float(os.environ.get("VOID_MIN", 2.0)))


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


def elution_window(tic, void=2.0):
    # widest span around the biggest TIC peak (after the void) that stays above half-max
    t, v = np.asarray(tic).T
    keep = t >= void
    t, v = t[keep], v[keep]
    k = int(np.argmax(v)); thr = v[k] / 2
    lo = hi = k
    while lo and v[lo - 1] >= thr: lo -= 1
    while hi < len(v) - 1 and v[hi + 1] >= thr: hi += 1
    return float(t[lo]), float(t[hi])


def dar_from_massdat(md, base, step, widths=(15, 25, 40, 60)):
    # DAR anchored at base / base+step; robustness table across integration half-widths
    table = []
    for w in widths:
        a0, a1 = area(md, base - w, base + w), area(md, base + step - w, base + step + w)
        table.append({"halfwidth_Da": w, "naked_area": round(a0, 1),
                      "conj_area": round(a1, 1),
                      "DAR": round(a1 / (a0 + a1) if a0 + a1 else 0, 3)})
    return table[1]["DAR"], table   # headline = +/-25 Da band


_PLT = None


def _mpl():
    global _PLT
    if _PLT is None:
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
        _PLT = plt
    return _PLT


def _chromatograms_uv_tic(path, cfg, imp=None):
    """(tic, uv, tmin, tmax) from the mzML. UV/TIC come from the <chromatogram>
    arrays with stdlib; the elution window is TMIN/TMAX if set, else auto. `imp` is
    an optional UniDec importer used only to synthesise a TIC when the file has none."""
    ch = chromatograms(path)

    def pick(*subs):   # first chromatogram whose id contains any of subs (case-insensitive)
        for s in subs:
            for k in ch:
                if s.lower() in k.lower():
                    return ch[k]
        return None
    tic = pick("TIC")
    uv = pick("Sig=280", "280", "uv", "absorb")   # DAD/UV if present (any vendor); else None
    if tic is None and imp is not None:            # some vendors omit a TIC chromatogram in mzML
        tic = np.asarray(imp.get_tic())
    if os.environ.get("TMIN") and os.environ.get("TMAX"):
        tmin, tmax = float(os.environ["TMIN"]), float(os.environ["TMAX"])
    elif tic is not None:
        tmin, tmax = elution_window(tic, cfg["void"])
    else:
        tmin, tmax = 0.0, 0.0
    return tic, uv, tmin, tmax


def analyze(path, out):
    from unidec import engine
    from unidec.UniDecImporter.ImporterFactory import ImporterFactory
    cfg = load_config()
    base, step = cfg["base"], cfg["step"]

    name = os.path.splitext(os.path.basename(path))[0]
    tag = name.replace(" ", "_")
    imp = ImporterFactory.create_importer(path)
    tic, uv, tmin, tmax = _chromatograms_uv_tic(path, cfg, imp)

    spec = os.path.join(out, "_avg_%s.txt" % tag)
    np.savetxt(spec, np.asarray(imp.get_avg_scan(time_range=(tmin, tmax))))
    e = engine.UniDec(); e.open_file(spec)
    e.config.startz, e.config.endz = cfg["zlo"], cfg["zhi"]
    e.config.masslb, e.config.massub = cfg["mlo"], cfg["mhi"]
    e.config.massbins = 1.0
    e.config.minmz, e.config.maxmz = 600.0, 2200.0
    e.process_data()
    e.get_auto_peak_width()   # match peak width to the data's resolution (vs UniDec's fixed default)
    e.run_unidec(silent=True)
    md = np.asarray(e.data.massdat)
    if md.ndim != 2 or not len(md):
        raise RuntimeError("empty deconvolution (massdat %s)" % (md.shape,))
    md[:, 1] *= 100.0 / md[:, 1].max()
    mz = np.asarray(e.data.data2)   # processed m/z envelope (the input to deconvolution)
    peak_width = float(e.config.mzsig)
    os.remove(spec)
    return render(name, out, md, mz, tic, uv, tmin, tmax, base, step, peak_width)


def replot(path, out):
    """Regenerate the figure for `path` from cached UniDec output, with no
    deconvolution (so no UniDec/Docker needed). Reads the averaged-scan mass
    distribution and processed m/z envelope a previous analyze() run wrote under
    _avg_<tag>_unidecfiles/, and the chromatograms straight from the mzML."""
    cfg = load_config()
    base, step = cfg["base"], cfg["step"]
    name = os.path.splitext(os.path.basename(path))[0]
    tag = name.replace(" ", "_")
    d = os.path.join(out, "_avg_%s_unidecfiles" % tag)
    massf = os.path.join(d, "_avg_%s_mass.txt" % tag)
    mzf = os.path.join(d, "_avg_%s_input.dat" % tag)
    if not (os.path.exists(massf) and os.path.exists(mzf)):
        raise SystemExit("no cached deconvolution for %s (looked in %s)" % (name, d))
    md = np.loadtxt(massf)
    md[:, 1] *= 100.0 / md[:, 1].max()
    mz = np.loadtxt(mzf)
    tic, uv, tmin, tmax = _chromatograms_uv_tic(path, cfg)
    return render(name, out, md, mz, tic, uv, tmin, tmax, base, step)


def _uv_concentration(uv, base, tmax, void=2.0):
    """Optional UV-280 quantitation of the main protein peak, Beer-Lambert straight off
    the DAD trace. Env-driven so it stays out of the way unless you ask for it:

        EPS280       molar extinction coeff at 280 nm, M^-1 cm^-1  (from the sequence)
        FLOW_ML_MIN  LC flow rate, mL/min  (NOT reliably in the mzML -- read the method)
        INJ_UL       injection volume, uL
        PATH_CM      flow-cell path length, cm            (default 1.0 = 10 mm)
        DILUTION     sample dilution before injection      (default 1)
        MW_DA        molecular weight for the mg/mL step    (default: base mass)
        DAD_UNIT     'mAU' (default) or 'AU'   (Agilent exports mAU even when the mzML
                     labels the array 'absorbance unit' -- check a 220 nm channel: if it
                     tops out in the 100s it is mAU, not AU)

    The peak area alone (always returned) is a RELATIVE amount, needs no inputs. Absolute
    concentration needs EPS280 + FLOW_ML_MIN + INJ_UL. Method: baseline-subtract the UV
    trace (5th percentile after the void), integrate the main peak (apex -0.7/+0.9 min);
    moles through the cell = A[AU*min] * F[L/min] / (eps * path); sample conc =
    moles / injection_volume * dilution; mg/mL = molar * MW. Assumes one pure, resolved
    protein peak (co-eluting species share the 280 signal)."""
    if uv is None or not len(uv):
        return {}
    unit = os.environ.get("DAD_UNIT", "mAU")
    post = uv[uv[:, 0] >= void]
    if not len(post):
        return {}
    bl = float(np.percentile(post[:, 1], 5))
    v = np.clip(uv[:, 1] - bl, 0.0, None)
    win = (uv[:, 0] >= void) & (uv[:, 0] <= max(tmax + 2.0, void + 1.0))
    if not win.any():
        return {}
    apex = float(uv[win][np.argmax(v[win]), 0])
    m = (uv[:, 0] >= apex - 0.7) & (uv[:, 0] <= apex + 0.9)
    area = float(np.trapezoid(v[m], uv[m, 0]))
    res = {"uv280_apex_min": round(apex, 2), "uv280_peak_area": round(area, 4),
           "uv280_area_unit": "%s*min" % unit}
    eps, flow, inj = (os.environ.get(k) for k in ("EPS280", "FLOW_ML_MIN", "INJ_UL"))
    if eps and flow and inj:
        eps, flow, inj = float(eps), float(flow), float(inj)
        path = float(os.environ.get("PATH_CM", 1.0))
        dil = float(os.environ.get("DILUTION", 1.0))
        mw = float(os.environ.get("MW_DA", base))
        area_au = area / 1000.0 if unit.lower() == "mau" else area   # -> AU*min
        moles = area_au * (flow / 1000.0) / (eps * path)             # flow L/min -> mol
        conc_m = moles / (inj * 1e-6) * dil                          # mol/L
        res.update({"protein_conc_uM": round(conc_m * 1e6, 3),
                    "protein_conc_mg_ml": round(conc_m * mw, 4),
                    "conc_inputs": {"eps280": eps, "path_cm": path, "inj_ul": inj,
                                    "dilution": dil, "flow_ml_min": flow,
                                    "dad_unit": unit, "mw_da": round(mw, 2)}})
    return res


def render(name, out, md, mz, tic, uv, tmin, tmax, base, step, peak_width=None):
    """Build the 3-panel figure and result dict from an already-deconvolved mass
    distribution. Split out of analyze() so figures can be regenerated from cache
    (see replot) without re-running the deconvolution."""
    plt = _mpl()
    tag = name.replace(" ", "_")
    dar, table = dar_from_massdat(md, base, step)
    # two-state area fractions of the (0-drug + 1-drug) population; the +1 fraction
    # IS the reported DAR, so labelling it on the peak makes the number self-evident.
    naked_pct, conj_pct = round((1 - dar) * 100, 1), round(dar * 100, 1)
    # optional labelling for report figures (all default to the plain output)
    fig_title = os.environ.get("FIG_TITLE")            # bold heading; else the file name
    fig_subtitle = os.environ.get("FIG_SUBTITLE")      # grey second line under the heading
    conj_label = os.environ.get("CONJ_LABEL", "+1")    # optional name for the +1 species (the reagent)
    suffix = os.environ.get("FIG_SUFFIX", "")          # appended to the output filename

    def uv_h(arr, a, b):
        if arr is None: return None
        m = (arr[:, 0] >= a) & (arr[:, 0] <= b)
        return round(float(arr[m, 1].max()), 3) if m.any() else 0.0
    uv_main, uv_late = uv_h(uv, tmin - .3, tmax + .3), uv_h(uv, 18, 22)
    conc = _uv_concentration(uv, base, tmax)   # optional UV-280 amount/concentration

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
    if conc.get("protein_conc_mg_ml") is not None:   # UV-280 concentration, shown on its own panel
        ci = conc["conc_inputs"]
        axA.text(0.37, 0.96,
                 "protein ≈ %.3f mg/mL (%.1f µM)\nUV-280: %g mm cell, %g µL inj, "
                 "%g mL/min, dil %g, %s" % (conc["protein_conc_mg_ml"], conc["protein_conc_uM"],
                 ci["path_cm"] * 10, ci["inj_ul"], ci["flow_ml_min"], ci["dilution"], ci["dad_unit"]),
                 transform=axA.transAxes, va="top", ha="left", fontsize=7.5, color=prim,
                 bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cfd8e3", lw=0.6, alpha=0.9))

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

    ymax = 120
    for lo_, hi_ in [(base - 25, base + 25), (base + step - 25, base + step + 25)]:
        axC.axvspan(lo_, hi_, color=band, alpha=0.6, lw=0)
    axC.fill_between(md[:, 0], md[:, 1], color=prim, alpha=0.18, lw=0)
    axC.plot(md[:, 0], md[:, 1], color=prim, lw=0.9)
    for m, lab in [(base, "unmodified"), (base + step, conj_label)]:
        h = float(md[np.argmin(abs(md[:, 0] - m)), 1])
        axC.annotate("%s\n%.0f Da" % (lab, m),
                     xy=(m, h), xytext=(m, min(h + 15, ymax - 5)),
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
    a0, a1 = table[1]["naked_area"], table[1]["conj_area"]   # +-25 Da headline bands
    # compact left-column block (short lines only) so it can't reach the peak labels
    axC.text(0.015, 0.90,
             "a0 = %.0f,  a1 = %.0f\nDAR = a1/(a0+a1)\n±25 Da bands, a.u." % (a0, a1),
             transform=axC.transAxes, va="top", ha="left", fontsize=7, fontweight="bold",
             color="#333333", linespacing=1.35)
    axC.spines[["top", "right"]].set_visible(False)

    for ax, letter in [(axA, "a"), (axB, "b"), (axC, "c")]:
        ax.text(-0.09, 1.02, letter, transform=ax.transAxes, fontsize=13, fontweight="bold")
    # reserve headroom first, then place the heading into it with fig.text (not
    # suptitle) so the title and subtitle can't collide with each other or panel a
    top = 0.95 if (fig_title and fig_subtitle) else 0.965 if fig_title else 0.98
    fig.tight_layout(rect=(0, 0, 1, top))
    if fig_title:
        fig.text(0.5, 0.993, fig_title, ha="center", va="top",
                 fontsize=12.5, fontweight="bold", color="#222222")
        if fig_subtitle:
            fig.text(0.5, 0.966, fig_subtitle, ha="center", va="top",
                     fontsize=9.5, color="#666666")
    else:
        fig.text(0.5, 0.992, name, ha="center", va="top", fontsize=10, color="#666666")
    figpath = os.path.join(out, "dar_%s%s.png" % (tag, suffix))
    fig.savefig(figpath, dpi=300)
    fig.savefig(figpath[:-4] + ".pdf")   # vector, for reports
    plt.close(fig)

    res = {"file": name, "window_min": [round(tmin, 2), round(tmax, 2)],
           "DAR": dar, "conjugated_pct": conj_pct, "unmodified_pct": naked_pct,
           "peak_width_mz": (round(peak_width, 3) if peak_width else None),
           "conj_apex": conj_apex, "naked_apex": naked_apex, "trace": trace,
           "expected": [base, round(base + step, 2)],
           "DAR_by_window": table,
           "UV280_main_h": uv_main, "UV280_late_h": uv_late,
           "UV_late_over_main": (round(uv_late / uv_main, 2) if uv_main else None)}
    res.update(conc)   # optional UV-280 amount/concentration (computed above)
    if "protein_conc_mg_ml" in res:
        ci = res["conc_inputs"]
        print("[conc] %s: %.4f mg/mL (%.1f uM) | UV280 peak %.3f %s @ %.2f min | "
              "eps=%g path=%gcm inj=%guL dil=%g flow=%gmL/min %s"
              % (name, res["protein_conc_mg_ml"], res["protein_conc_uM"],
                 res["uv280_peak_area"], res["uv280_area_unit"], res["uv280_apex_min"],
                 ci["eps280"], ci["path_cm"], ci["inj_ul"], ci["dilution"],
                 ci["flow_ml_min"], ci["dad_unit"]))
    return res


if __name__ == "__main__":
    out = os.environ.get("OUTDIR", "/data")
    fn = replot if os.environ.get("REPLOT") else analyze   # REPLOT=1 reuses cached deconvolution
    res = []
    for p in sys.argv[1:]:
        try:
            res.append(fn(p, out))
        except Exception as ex:
            res.append({"file": os.path.basename(p), "error": str(ex)})
    json.dump(res, open(os.path.join(out, "dar_results.json"), "w"), indent=2)
    print(json.dumps(res, indent=2))
