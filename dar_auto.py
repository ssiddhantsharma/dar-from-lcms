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
    # regime presets: 'denaturing' = RP intact of small proteins (low charge, m/z 600-2200);
    # 'native' = native SEC-MS of large folded proteins (high charge, m/z ~2000-8000). Every
    # preset value is individually overridable by its own env var.
    native = os.environ.get("MODE", "denaturing").lower() == "native"
    mzlo0, mzhi0, zlo0, zhi0, mub0 = ((2000.0, 8000.0, 10, 45, base + 6000.0) if native
                                      else (600.0, 2200.0, 3, 20, base + 3000.0))
    return dict(base=base, step=step, native=native,
                mlo=float(os.environ.get("MASS_LB", base - 2000)),
                mhi=float(os.environ.get("MASS_UB", mub0)),
                zlo=int(os.environ.get("Z_LO", zlo0)), zhi=int(os.environ.get("Z_HI", zhi0)),
                mzlo=float(os.environ.get("MZ_LO", mzlo0)), mzhi=float(os.environ.get("MZ_HI", mzhi0)),
                massbins=float(os.environ.get("MASSBINS", 1.0)),
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


def _charge_support(mz, mass, zlo, zhi, thr=0.05, tol=1.0):
    """How many charge states independently back a deconvolved mass. For each z in [zlo,zhi]
    the species would appear at m/z = (mass + z*proton)/z; count how many of those predicted
    positions carry real signal in the processed m/z envelope. A mass supported by many charge
    states is a genuine species; one supported by only one or two is likely a harmonic or a
    deconvolution artifact, so this is a stronger 'is the peak real' check than area alone.

    Idea from FLASHDeconv's multi-charge agreement scoring (Jeong et al. 2020), adapted to
    charge-resolved (not isotope-resolved) intact data: we test charge-state presence, not the
    isotope-envelope cosine, which our resolution does not support."""
    if mz is None or not len(mz) or mass <= 0:
        return 0
    x, y = mz[:, 0], mz[:, 1]
    ymax = float(y.max())
    if ymax <= 0:
        return 0
    proton = 1.007276
    n = 0
    for z in range(int(zlo), int(zhi) + 1):
        if z <= 0:
            continue
        target = (mass + z * proton) / z
        m = (x >= target - tol) & (x <= target + tol)
        if m.any() and float(y[m].max()) >= thr * ymax:
            n += 1
    return n


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


def _state_area(md, base, step, n, satellites, w):
    """Integrated area of load state n = sum over its satellite offsets of a +-w Da band
    at base + n*step + offset. Satellites are the glycoforms and/or adducts that belong
    to the same load state (see dar_distribution)."""
    c = base + n * step
    return sum(area(md, c + g - w, c + g + w) for g in satellites)


def dar_distribution(md, base, step, nmax, w=25.0, satellites=(0.0,)):
    """Multi-state DAR/CAR from a deconvolved mass spectrum. For load state n, sum a +-w Da
    band at base + n*step + g over the satellite offsets g -- the glycoforms and/or adducts
    that belong to the same load state (e.g. IgG glycoforms 0/+162/+324, +16 oxidation,
    +178 gluconoylation). `satellites=(0.0,)` (the default) is a single band, i.e. the
    resolved-species case; passing a glycoform list makes the count glycoform-aware, which
    is the honest fix for native glycosylated IgG whose load states otherwise overlap.

    Reports per-state area and fraction, the average load, and the dispersity index
    (van der Zon et al., Anal. Chim. Acta 1395 (2026), Eqs 1-3):

        average = sum_n (n * A_n) / sum_n (A_n)                       (Eq 1)
        Mw      = sum_n (n^2 * A_n) / sum_n (n * A_n)                 (Eq 2)
        dispersity = Mw / average                                    (Eq 3)

    The two-state DAR (dar_from_massdat) is exactly the nmax=1, single-band case. Keep the
    satellite offsets below `step` so one state's envelope does not bleed into the next."""
    areas = [_state_area(md, base, step, n, satellites, w) for n in range(nmax + 1)]
    tot = sum(areas)
    if tot <= 0:
        nan = float("nan")
        return {"average_dar": nan, "dispersity": nan,
                "state_areas": [round(a, 1) for a in areas],
                "state_frac": [nan] * (nmax + 1)}
    num = sum(n * a for n, a in enumerate(areas))            # sum n*A_n (numerator of Eq 1)
    avg = num / tot
    disp = (sum(n * n * a for n, a in enumerate(areas)) / num) / avg if (num and avg) else float("nan")
    return {"average_dar": round(avg, 3), "dispersity": round(disp, 3),
            "state_areas": [round(a, 1) for a in areas],
            "state_frac": [round(a / tot, 4) for a in areas]}


def _satellites():
    """Satellite mass offsets (Da) folded into each load state, from the SATELLITES env
    (comma-separated), e.g. IgG glycoforms '0,162.05,324.11' or with adducts
    '0,16,162.05,178'. Default (0.0,) = a single band per state."""
    offs = tuple(float(x) for x in os.environ.get("SATELLITES", "").split(",") if x.strip())
    return offs or (0.0,)


def _load_mirror_mass():
    """Optional control spectrum to mirror below the sample on the deconvolved-mass panel
    (env MIRROR_MASS = path to a control's deconvolved *_mass.txt). Returns it normalised to
    100, or None if unset/missing. The mirror (treated on top, control below) is the standard
    head-to-tail comparison, following spectrum_utils' mirror() convention (reimplemented)."""
    p = os.environ.get("MIRROR_MASS")
    if not p or not os.path.exists(p):
        return None
    a = np.loadtxt(p)
    if a.ndim != 2 or not len(a):
        return None
    a[:, 1] *= 100.0 / a[:, 1].max()
    return a


def _dar_uncertainty(md, base, step, nmax, satellites, widths=(15.0, 25.0, 40.0)):
    """Reproducibility estimate for the average DAR: its spread across integration
    half-widths. Returns the population standard deviation over `widths`, or nan if any
    width gives nan. This is a lower bound on uncertainty (window choice only); true
    uncertainty also needs replicate injections."""
    vals = [dar_distribution(md, base, step, nmax, w, satellites)["average_dar"] for w in widths]
    if any(v != v for v in vals):
        return float("nan")
    mean = sum(vals) / len(vals)
    return (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5


def _mass_quality(md, base, step, nmax, satellites, w=25.0):
    """Two trust checks on the deconvolution: (1) mass_error_ppm, observed-vs-theoretical
    mass of the most abundant load state (its base glycoform); (2) captured_fraction, the
    share of total deconvolved signal that falls in the anchored state+satellite bands --
    a low value means much signal sits outside the model (bad anchors, adducts, or a noisy
    deconvolution) and the DAR should be treated with caution."""
    areas = [_state_area(md, base, step, n, satellites, w) for n in range(nmax + 1)]
    if sum(areas) <= 0 or not len(md):
        return {"mass_error_ppm": None, "captured_fraction": 0.0}
    n_dom = int(max(range(nmax + 1), key=lambda n: areas[n]))
    theo = base + n_dom * step
    obs, _ = apex(md, theo, w=w)
    ppm = round((obs - theo) / theo * 1e6, 1) if obs is not None else None
    total = area(md, md[0, 0], md[-1, 0])
    return {"mass_error_ppm": ppm,
            "captured_fraction": round(sum(areas) / total, 3) if total > 0 else 0.0}


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


# Species colour map, z-order, and the per-axis styling / label-formatter helpers below are
# adapted from the patterns in spectrum_utils (Bittremieux et al., Apache-2.0): a `colors`
# dict keyed by species with a None fallback, a `_format_ax` axis styler, and a pluggable
# label formatter. Reimplemented here for DAR/CAR load states (not peptide ion types).
COLORS = {"unmodified": "#1f4e79", "conjugate": "#1976d2", "adduct": "#9aa0a6",
          "control": "#b0403a", None: "#212121"}


def _format_ax(ax, xlabel, ylabel, italic_x=False):
    """Consistent axis styling: faint grid drawn behind the data, minor ticks, small tick
    labels, top/right spines off. Pattern from spectrum_utils (_format_ax), reimplemented."""
    import matplotlib.ticker as mticker
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.grid(True, "major", color="#9e9e9e", linewidth=0.2)
    ax.grid(True, "minor", color="#9e9e9e", linewidth=0.15)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", which="both", labelsize="small")
    ax.set_xlabel(xlabel, style="italic" if italic_x else "normal")
    ax.set_ylabel(ylabel)
    ax.spines[["top", "right"]].set_visible(False)


def _state_label(n, conj_label="+1"):
    """Label for load state n on the deconvolved-mass panel. Pluggable formatter in the
    spirit of spectrum_utils' `annot_fmt` (reimplemented for DAR load states)."""
    return "unmodified" if n == 0 else (conj_label if n == 1 else "+%d" % n)


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
    e.config.massbins = cfg["massbins"]
    e.config.minmz, e.config.maxmz = cfg["mzlo"], cfg["mzhi"]
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
    d, _, keyf = _cache_paths(out, tag)   # record the deconvolution key so process() can resume
    if os.path.isdir(d):
        json.dump({"key": _decon_key(cfg), "peak_width": peak_width}, open(keyf, "w"))
    return render(name, out, md, mz, tic, uv, tmin, tmax, base, step, peak_width)


def replot(path, out, peak_width=None):
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
    return render(name, out, md, mz, tic, uv, tmin, tmax, base, step, peak_width)


def _decon_key(cfg):
    """The parameters that actually govern the UniDec deconvolution. Everything DAR-specific
    (MOD_MASS, SATELLITES, the mirror, band integration, the whole figure) is cheap arithmetic
    on the resulting massdat, so a change to any of those does NOT invalidate the cache; a
    change to these does. BASE_MASS is folded in via the mass bounds (mlo/mhi track it)."""
    return {"zlo": cfg["zlo"], "zhi": cfg["zhi"], "mlo": cfg["mlo"], "mhi": cfg["mhi"],
            "massbins": cfg["massbins"], "mzlo": cfg["mzlo"], "mzhi": cfg["mzhi"],
            "tmin": os.environ.get("TMIN"), "tmax": os.environ.get("TMAX")}


def _cache_paths(out, tag):
    d = os.path.join(out, "_avg_%s_unidecfiles" % tag)
    return d, os.path.join(d, "_avg_%s_mass.txt" % tag), os.path.join(d, "_darcache.json")


def process(path, out):
    """Auto-resume dispatcher: reuse a cached deconvolution when its governing params are
    unchanged (the fast path -- re-anchoring, glycoform offsets and figure tweaks then cost
    nothing), and only fall back to a fresh UniDec run when the cache is missing or its
    deconvolution key differs. REPLOT=1 forces replot; NO_CACHE=1 forces a fresh run."""
    if os.environ.get("REPLOT"):
        return replot(path, out)
    if os.environ.get("NO_CACHE"):
        return analyze(path, out)
    cfg = load_config()
    tag = os.path.splitext(os.path.basename(path))[0].replace(" ", "_")
    _, massf, keyf = _cache_paths(out, tag)
    if os.path.exists(massf) and os.path.exists(keyf):
        try:
            cached = json.load(open(keyf))
            if cached.get("key") == _decon_key(cfg):
                return replot(path, out, peak_width=cached.get("peak_width"))
        except (ValueError, OSError):
            pass
    return analyze(path, out)


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
    pk_area = float(np.trapezoid(v[m], uv[m, 0]))   # not the module-level area(); local peak area
    res = {"uv280_apex_min": round(apex, 2), "uv280_peak_area": round(pk_area, 4),
           "uv280_area_unit": "%s*min" % unit}
    eps, flow, inj = (os.environ.get(k) for k in ("EPS280", "FLOW_ML_MIN", "INJ_UL"))
    if eps and flow and inj:
        eps, flow, inj = float(eps), float(flow), float(inj)
        path = float(os.environ.get("PATH_CM", 1.0))
        dil = float(os.environ.get("DILUTION", 1.0))
        mw = float(os.environ.get("MW_DA", base))
        area_au = pk_area / 1000.0 if unit.lower() == "mau" else pk_area   # -> AU*min
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
    nmax = int(os.environ.get("DAR_MAX_N", "1"))       # >1 -> multi-state DAR/CAR distribution
    sats = _satellites()                               # glycoform/adduct offsets folded per state
    dist = dar_distribution(md, base, step, nmax, satellites=sats) if nmax > 1 else None
    mirror_md = _load_mirror_mass()                    # optional control spectrum to mirror below
    # mass-accuracy + deconvolution-quality, and an uncertainty (integration-window spread)
    mq = _mass_quality(md, base, step, nmax if nmax > 1 else 1, sats if nmax > 1 else (0.0,))
    # charge-state support per load state: how many charge states back each species (real vs
    # artifact); FLASHDeconv multi-charge-agreement idea, adapted to charge-resolved data
    try:
        _c = load_config(); _zlo, _zhi = _c["zlo"], _c["zhi"]
    except SystemExit:                     # render() called directly without env: sane defaults
        _zlo, _zhi = 3, 20
    _nstates = nmax if dist is not None else 1
    _cstol = max(1.0, peak_width or 0.0)
    charge_support = [_charge_support(mz, base + n * step, _zlo, _zhi, tol=_cstol)
                      for n in range(_nstates + 1)]
    dom_idx = int(np.argmax(dist["state_areas"])) if dist is not None else (1 if dar >= 0.5 else 0)
    charge_support_dom = charge_support[min(dom_idx, len(charge_support) - 1)]
    if dist is not None:
        dar_sd = round(_dar_uncertainty(md, base, step, nmax, sats), 3)
    else:
        _dv = [r["DAR"] for r in table]
        _m = sum(_dv) / len(_dv)
        dar_sd = round((sum((x - _m) ** 2 for x in _dv) / len(_dv)) ** 0.5, 3)
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
    axA.legend(loc="upper right", frameon=False, handlelength=1.4)
    _format_ax(axA, "retention time (min)", "relative signal (%)")
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
    axB.set_title("raw charge-state envelope (input to deconvolution)",
                  loc="left", fontsize=9.5, color="#666666")
    _format_ax(axB, "m/z", "relative intensity (%)", italic_x=True)

    ymax = 120
    bands = ([(base + n * step - 25, base + n * step + 25) for n in range(nmax + 1)]
             if dist is not None else
             [(base - 25, base + 25), (base + step - 25, base + step + 25)])
    for lo_, hi_ in bands:
        axC.axvspan(lo_, hi_, color=band, alpha=0.6, lw=0)
    axC.fill_between(md[:, 0], md[:, 1], color=prim, alpha=0.18, lw=0)
    axC.plot(md[:, 0], md[:, 1], color=prim, lw=0.9)
    if mirror_md is not None:                               # head-to-tail: control mirrored below
        axC.fill_between(mirror_md[:, 0], -mirror_md[:, 1], color=COLORS["control"], alpha=0.16, lw=0)
        axC.plot(mirror_md[:, 0], -mirror_md[:, 1], color=COLORS["control"], lw=0.9)
        axC.axhline(0, color="#666666", lw=0.6)
        axC.text(0.985, 0.97, "treated", transform=axC.transAxes, ha="right", va="top",
                 fontsize=8, fontweight="bold", color=prim)
        axC.text(0.985, 0.03, "control", transform=axC.transAxes, ha="right", va="bottom",
                 fontsize=8, fontweight="bold", color=COLORS["control"])
    if dist is not None:                                    # multi-state ladder: label +0.. +nmax
        for n in range(nmax + 1):
            m = base + n * step
            h = float(md[np.argmin(abs(md[:, 0] - m)), 1])
            axC.annotate("+%d" % n, xy=(m, h), xytext=(0, 4), textcoords="offset points",
                         ha="center", fontsize=7.5, color="#444444")
        axC.set_xlim(base - 300, base + nmax * step + 300)
    else:
        for m, lab in [(base, _state_label(0)), (base + step, _state_label(1, conj_label))]:
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
        axC.set_xlim(base - 400, base + step + 500)
    axC.set_ylim(-ymax if mirror_md is not None else 0, ymax)
    pm = (" ± %.2f" % dar_sd) if dar_sd == dar_sd else ""   # append uncertainty unless nan
    # mass-accuracy + capture as a small trust line (spectrum_utils mass_errors idea, scaled
    # to our single dominant state)
    mqs = "mass err %s ppm; captured %.0f%%; %d charge states" % (
        ("%+.0f" % mq["mass_error_ppm"]) if mq["mass_error_ppm"] is not None else "n/a",
        100 * mq["captured_fraction"], charge_support_dom)
    if dist is not None:
        axC.text(0.015, 0.96, "average DAR = %.2f%s" % (dist["average_dar"], pm),
                 transform=axC.transAxes, va="top", ha="left", fontsize=11,
                 fontweight="bold", color=prim)
        axC.text(0.015, 0.90,
                 "dispersity %.2f;  states 0-%d\nsum(n*An)/sum(An); ±25 Da bands, a.u.\n%s"
                 % (dist["dispersity"], nmax, mqs),
                 transform=axC.transAxes, va="top", ha="left", fontsize=7,
                 fontweight="bold", color="#333333", linespacing=1.35)
    else:
        axC.text(0.015, 0.96, "DAR = %.2f%s" % (dar, pm), transform=axC.transAxes,
                 va="top", ha="left", fontsize=11, fontweight="bold", color=prim)
        a0, a1 = table[1]["naked_area"], table[1]["conj_area"]   # +-25 Da headline bands
        axC.text(0.015, 0.90,      # compact left-column block, short lines, clear of peak labels
                 "a0 = %.0f,  a1 = %.0f\nDAR = a1/(a0+a1)\n±25 Da bands, a.u.\n%s" % (a0, a1, mqs),
                 transform=axC.transAxes, va="top", ha="left", fontsize=7, fontweight="bold",
                 color="#333333", linespacing=1.35)
    _format_ax(axC, "deconvolved mass (Da)", "relative abundance (%)")

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
           "UV_late_over_main": (round(uv_late / uv_main, 2) if uv_main else None),
           "dar_sd": (dar_sd if dar_sd == dar_sd else None),
           "mass_error_ppm": mq["mass_error_ppm"], "captured_fraction": mq["captured_fraction"],
           "charge_support": charge_support, "charge_support_dominant": charge_support_dom}
    res.update(conc)   # optional UV-280 amount/concentration (computed above)
    if dist is not None:   # optional multi-state DAR/CAR distribution
        res.update({"average_dar": dist["average_dar"], "average_dar_sd": res["dar_sd"],
                    "dispersity": dist["dispersity"], "dar_states": nmax,
                    "satellites": list(sats), "dar_state_frac": dist["state_frac"],
                    "dar_state_areas": dist["state_areas"]})
        print("[multi-DAR] %s: average DAR %.3f ± %s, dispersity %.3f | fractions %s | satellites %s"
              % (name, dist["average_dar"], res["dar_sd"], dist["dispersity"],
                 [round(f, 3) for f in dist["state_frac"]], list(sats)))
    print("[quality] %s: mass error %s ppm, captured %.0f%% of signal in anchored bands, "
          "%d charge states support the main species"
          % (name, mq["mass_error_ppm"], 100 * mq["captured_fraction"], charge_support_dom))
    if "protein_conc_mg_ml" in res:
        ci = res["conc_inputs"]
        print("[conc] %s: %.4f mg/mL (%.1f uM) | UV280 peak %.3f %s @ %.2f min | "
              "eps=%g path=%gcm inj=%guL dil=%g flow=%gmL/min %s"
              % (name, res["protein_conc_mg_ml"], res["protein_conc_uM"],
                 res["uv280_peak_area"], res["uv280_area_unit"], res["uv280_apex_min"],
                 ci["eps280"], ci["path_cm"], ci["inj_ul"], ci["dilution"],
                 ci["flow_ml_min"], ci["dad_unit"]))
    return res


# per-sample columns a manifest row may set (each maps to the env var of the same name)
_MANIFEST_ENV = ("BASE_MASS", "MOD_MASS", "MODE", "DAR_MAX_N", "SATELLITES", "MIRROR_MASS",
                 "EPS280", "FLOW_ML_MIN", "INJ_UL", "PATH_CM", "DILUTION", "MW_DA", "DAD_UNIT",
                 "TMIN", "TMAX", "Z_LO", "Z_HI", "MZ_LO", "MZ_HI", "MASS_LB", "MASS_UB")


def run_manifest(path, out, fn):
    """Batch driver from a per-sample manifest CSV: one row per sample, each setting its own
    chemistry/params via columns named like the env vars in _MANIFEST_ENV (plus `file`). A
    plate of different conjugates thus runs in one command. Input-manifest pattern from
    SmartPeak's sequence.csv + nf-core's schema_input.json (reimplemented); a JSON schema for
    the columns lives in assets/manifest_schema.json."""
    import csv
    res = []
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        f = row.get("file") or row.get("File")
        if not f:
            continue
        saved = {k: os.environ.get(k) for k in _MANIFEST_ENV}   # set row's params, restore after
        for k in _MANIFEST_ENV:
            v = row.get(k, row.get(k.lower()))
            if v not in (None, ""):
                os.environ[k] = str(v)
        try:
            res.append(fn(f if os.path.isabs(f) else os.path.join(out, f), out))
        except Exception as ex:  # noqa: BLE001 - keep the batch going
            res.append({"file": os.path.basename(f), "error": str(ex)})
        for k, v in saved.items():
            os.environ.pop(k, None) if v is None else os.environ.__setitem__(k, v)
    return res


def plate_heatmap(res, out):
    """Batch-summary heatmap: samples x load-state (CAR/DAR) fractions. Only for multi-state
    runs (needs dar_state_frac) with >=2 samples. Heatmap idea from SmartPeak (reimplemented)."""
    rows = [r for r in res if r.get("dar_state_frac")]
    if len(rows) < 2:
        return None
    plt = _mpl()
    nmax = max(len(r["dar_state_frac"]) for r in rows)
    mat = np.array([r["dar_state_frac"] + [0.0] * (nmax - len(r["dar_state_frac"])) for r in rows])
    fig, ax = plt.subplots(figsize=(max(5.0, nmax * 0.8), max(2.5, len(rows) * 0.5)))
    im = ax.imshow(mat, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(nmax)); ax.set_xticklabels(["+%d" % n for n in range(nmax)])
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r["file"] for r in rows], fontsize=7)
    ax.set_xlabel("load state")
    ax.set_title("CAR/DAR distribution across samples", loc="left", fontsize=9.5, color="#666666")
    for i in range(len(rows)):
        for j in range(nmax):
            if mat[i, j] > 0.02:
                ax.text(j, i, "%.0f%%" % (100 * mat[i, j]), ha="center", va="center",
                        fontsize=6, color="white" if mat[i, j] > 0.5 else "#222222")
    fig.colorbar(im, ax=ax, label="fraction", fraction=0.03)
    fig.tight_layout()
    p = os.path.join(out, "dar_plate_heatmap.png")
    fig.savefig(p, dpi=200)
    plt.close(fig)
    return p


def _worker(args):
    """Module-level so multiprocessing can pickle it; env (the shared folder config) is
    inherited by each worker, so this is only used for the shared-config folder case."""
    path, out = args
    try:
        return process(path, out)
    except Exception as ex:  # noqa: BLE001
        return {"file": os.path.basename(path), "error": str(ex)}


if __name__ == "__main__":
    out = os.environ.get("OUTDIR", "/data")
    manifest = os.environ.get("MANIFEST")                  # per-sample CSV (else CLI file args)
    jobs = int(os.environ.get("JOBS", "1"))                # >1 -> parallel folder run
    if manifest:
        # per-row env mutation makes the manifest path process-global; keep it serial. process()
        # still resumes from cache per row (REPLOT=1 forces replot, NO_CACHE=1 forces fresh).
        res = run_manifest(manifest, out, process)
    elif jobs > 1 and len(sys.argv) > 2:
        from multiprocessing import Pool
        with Pool(jobs) as pool:
            res = pool.map(_worker, [(p, out) for p in sys.argv[1:]])
    else:
        res = [_worker((p, out)) for p in sys.argv[1:]]
    json.dump(res, open(os.path.join(out, "dar_results.json"), "w"), indent=2)
    print(json.dumps(res, indent=2))

    # tidy one-row-per-sample summary for plate/batch QC, plus a distribution heatmap
    import csv
    cols = ["file", "DAR", "dar_sd", "average_dar", "average_dar_sd", "dispersity",
            "mass_error_ppm", "captured_fraction", "charge_support_dominant",
            "protein_conc_mg_ml", "window_min", "error"]
    with open(os.path.join(out, "dar_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in res:
            w.writerow({k: r.get(k) for k in cols})
    hp = plate_heatmap(res, out)
    if hp:
        print("wrote", hp)
