"""Unit tests for the analytical core of dar_auto.

These cover the pure functions (no UniDec/Docker needed): mzML chromatogram parsing,
band integration, peak apex, charge-state assignment, elution-window detection, and the
DAR math. Run with `pytest`.
"""
import base64
import importlib.util
import json
import os
import pathlib
import zlib

import numpy as np
import pytest

# import dar_auto by path; it must import without env vars or the UniDec/matplotlib stack
_SPEC = importlib.util.spec_from_file_location(
    "dar_auto", pathlib.Path(__file__).resolve().parent.parent / "dar_auto.py")
da = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(da)


def test_import_needs_no_env_or_heavy_stack():
    # merely importing the module must not require BASE_MASS/MOD_MASS or UniDec
    assert hasattr(da, "analyze") and hasattr(da, "dar_from_massdat")


def test_load_config_requires_masses(monkeypatch):
    monkeypatch.delenv("BASE_MASS", raising=False)
    monkeypatch.delenv("MOD_MASS", raising=False)
    with pytest.raises(SystemExit):
        da.load_config()


def test_load_config_parses_and_derives(monkeypatch):
    monkeypatch.setenv("BASE_MASS", "10000")
    monkeypatch.setenv("MOD_MASS", "500")
    monkeypatch.delenv("MASS_LB", raising=False)
    monkeypatch.delenv("MASS_UB", raising=False)
    cfg = da.load_config()
    assert cfg["base"] == 10000.0 and cfg["step"] == 500.0
    assert cfg["mlo"] == 8000.0 and cfg["mhi"] == 13000.0   # base-2000, base+3000


def test_area_is_trapezoidal():
    md = np.array([[0.0, 0.0], [1.0, 10.0], [2.0, 0.0]])   # triangle, integral = 10
    assert abs(da.area(md, 0, 2) - 10.0) < 1e-9
    assert da.area(md, 50, 60) == 0.0                       # empty band -> 0


def test_apex_finds_local_max_and_height():
    md = np.array([[100.0, 1.0], [101.0, 5.0], [102.0, 2.0], [140.0, 9.0]])
    m, h = da.apex(md, 101, w=2)                             # window excludes the far 140 peak
    assert m == 101.0 and h == 5.0
    assert da.apex(md, 500, w=2) == (None, 0.0)             # nothing in window


def test_elution_window_brackets_the_peak():
    t = np.linspace(0, 25, 501)
    v = np.exp(-((t - 8.5) ** 2) / 0.10)                    # sharp peak at 8.5 min
    v[t < 0.5] = 5.0                                         # void spike that must be ignored
    tmin, tmax = da.elution_window(np.column_stack([t, v]), void=2.0)
    assert tmin < 8.5 < tmax and (tmax - tmin) < 3.0


def _gauss(x, mu, amp, s=1.0):
    return amp * np.exp(-((x - mu) ** 2) / (2 * s * s))


def test_dar_from_massdat_matches_area_ratio():
    base, step = 10000.0, 500.0
    m = np.arange(9000, 11000, 1.0)
    inten = _gauss(m, base, 30.0) + _gauss(m, base + step, 90.0)   # 3:1 conj:naked by area
    md = np.column_stack([m, inten])
    dar, table = da.dar_from_massdat(md, base, step)
    assert abs(dar - 0.75) < 0.03                                  # 90/(30+90)
    assert [r["halfwidth_Da"] for r in table] == [15, 25, 40, 60]


def test_dar_is_zero_when_no_conjugate():
    base, step = 10000.0, 500.0
    m = np.arange(9000, 11000, 1.0)
    inten = _gauss(m, base, 100.0)                                  # only unmodified
    dar, _ = da.dar_from_massdat(np.column_stack([m, inten]), base, step)
    assert dar < 0.02


def test_charge_peaks_recovers_charges():
    mass = 10000.0
    mz = np.arange(600, 2100, 0.5)
    inten = np.zeros_like(mz)
    for z in (6, 7, 8, 9):
        inten += _gauss(mz, (mass + z * 1.00728) / z, 100.0, s=0.7)
    peaks = da.charge_peaks(np.column_stack([mz, inten]), mass)
    zs = sorted(z for _, _, z in peaks)
    assert {6, 7, 8, 9}.issubset(set(zs))
    assert len(zs) == len(set(zs))                                  # one label per charge state


def _write_mini_mzml(path, t, i):
    def enc(a):
        return base64.b64encode(zlib.compress(np.asarray(a, np.float32).tobytes())).decode()
    xml = (
        '<?xml version="1.0"?><mzML><run><chromatogramList count="1">'
        '<chromatogram id="TIC"><binaryDataArrayList count="2">'
        '<binaryDataArray><cvParam name="32-bit float"/><cvParam name="zlib compression"/>'
        '<cvParam name="time array"/><binary>%s</binary></binaryDataArray>'
        '<binaryDataArray><cvParam name="32-bit float"/><cvParam name="zlib compression"/>'
        '<cvParam name="intensity array"/><binary>%s</binary></binaryDataArray>'
        '</binaryDataArrayList></chromatogram></chromatogramList></run></mzML>'
        % (enc(t), enc(i)))
    path.write_text(xml)


def test_dar_distribution_matches_paper_equations():
    # van der Zon et al. 2026 worked example (dispersity definition): intensities
    # 20/60/20 for CAR1/CAR2/CAR3 -> averageDAR 2.0, dispersity 1.1. Validates Eqs 1-3.
    base, step = 10000.0, 500.0
    m = np.arange(base - 200, base + 5 * step + 200, 1.0)
    inten = (_gauss(m, base + 1 * step, 20, 3) + _gauss(m, base + 2 * step, 60, 3)
             + _gauss(m, base + 3 * step, 20, 3))
    md = np.column_stack([m, inten])
    d = da.dar_distribution(md, base, step, nmax=4)
    assert abs(d["average_dar"] - 2.0) < 0.02
    assert abs(d["dispersity"] - 1.1) < 0.02
    assert abs(sum(d["state_frac"]) - 1.0) < 1e-6

    # the two-state DAR is exactly the nmax=1 case of the general average
    dar1, _ = da.dar_from_massdat(md, base, step)
    d1 = da.dar_distribution(md, base, step, nmax=1)
    assert abs(d1["average_dar"] - dar1) < 2e-3


def test_load_mirror_mass(tmp_path, monkeypatch):
    monkeypatch.delenv("MIRROR_MASS", raising=False)
    assert da._load_mirror_mass() is None                   # unset -> None
    p = tmp_path / "ctrl_mass.txt"
    np.savetxt(str(p), np.column_stack([np.arange(9000, 9100, 1.0),
                                        _gauss(np.arange(9000, 9100, 1.0), 9050, 40, 3)]))
    monkeypatch.setenv("MIRROR_MASS", str(p))
    m = da._load_mirror_mass()
    assert m is not None and abs(m[:, 1].max() - 100.0) < 1e-6   # normalised to 100
    monkeypatch.setenv("MIRROR_MASS", str(tmp_path / "nope.txt"))
    assert da._load_mirror_mass() is None                   # missing file -> None


def test_run_manifest_sets_and_restores_env(tmp_path, monkeypatch):
    # a fake analysis fn records the per-row env it sees; run_manifest must set each row's
    # params and restore the prior environment afterwards.
    for k in ("BASE_MASS", "MOD_MASS", "DAR_MAX_N"):
        monkeypatch.delenv(k, raising=False)
    man = tmp_path / "manifest.csv"
    man.write_text("file,BASE_MASS,MOD_MASS,DAR_MAX_N\n"
                   "a.mzML,9515.87,526.55,1\n"
                   "b.mzML,148057,769.7,6\n")
    seen = []

    def fake(path, out):
        seen.append((os.path.basename(path), os.environ.get("BASE_MASS"), os.environ.get("DAR_MAX_N")))
        return {"file": os.path.basename(path)}
    res = da.run_manifest(str(man), str(tmp_path), fake)
    assert [r["file"] for r in res] == ["a.mzML", "b.mzML"]
    assert seen == [("a.mzML", "9515.87", "1"), ("b.mzML", "148057", "6")]
    assert os.environ.get("BASE_MASS") is None and os.environ.get("DAR_MAX_N") is None  # restored


def test_plate_heatmap(tmp_path):
    res = [{"file": "s1", "dar_state_frac": [0.1, 0.2, 0.7]},
           {"file": "s2", "dar_state_frac": [0.6, 0.3, 0.1]}]
    p = da.plate_heatmap(res, str(tmp_path))
    assert p and os.path.exists(p)                          # heatmap written for >=2 multi-state rows
    assert da.plate_heatmap(res[:1], str(tmp_path)) is None  # <2 samples -> skip
    assert da.plate_heatmap([{"file": "x", "DAR": 0.9}], str(tmp_path)) is None  # no distribution -> skip


def test_mode_presets(monkeypatch):
    for k in ("MODE", "MZ_LO", "MZ_HI", "Z_LO", "Z_HI", "MASS_LB", "MASS_UB", "MASSBINS"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("BASE_MASS", "148000"); monkeypatch.setenv("MOD_MASS", "770")
    d = da.load_config()                                   # denaturing default
    assert (d["mzlo"], d["mzhi"], d["zlo"], d["zhi"]) == (600.0, 2200.0, 3, 20)
    assert not d["native"]
    monkeypatch.setenv("MODE", "native")                   # native preset
    n = da.load_config()
    assert n["native"] and n["mzlo"] == 2000.0 and n["mzhi"] == 8000.0 and n["zhi"] == 45
    monkeypatch.setenv("MZ_HI", "7000")                    # per-var override beats the preset
    assert da.load_config()["mzhi"] == 7000.0


def test_satellite_folding_recovers_average_when_glycoforms_differ():
    # CAR0 sits mostly on its base glycoform; CAR1 mostly on its +162/+324 glycoforms.
    # A single band (base glycoform only) under-counts CAR1 -> wrong low DAR; folding the
    # glycoform satellites recovers the true 0.5.
    base, step = 148000.0, 770.0
    m = np.arange(base - 300, base + step + 600, 1.0)

    def glyc(center, amps):
        return sum(_gauss(m, center + off, a, 3) for off, a in zip((0.0, 162.0, 324.0), amps))
    md = np.column_stack([m, glyc(base, [80, 15, 5]) + glyc(base + step, [10, 30, 60])])
    single = da.dar_distribution(md, base, step, 1)                          # base glycoform only
    glyco = da.dar_distribution(md, base, step, 1, satellites=(0.0, 162.0, 324.0))
    assert abs(single["average_dar"] - 10 / 90) < 0.03      # under-counts the +1 state
    assert abs(glyco["average_dar"] - 100 / 200) < 0.03     # glycoform-aware recovers 0.5
    assert glyco["state_areas"][1] > single["state_areas"][1]


def test_dar_uncertainty_small_for_clean_peaks():
    base, step = 10000.0, 500.0
    m = np.arange(base - 300, base + step + 300, 1.0)
    md = np.column_stack([m, _gauss(m, base, 30, 3) + _gauss(m, base + step, 90, 3)])
    sd = da._dar_uncertainty(md, base, step, 1, (0.0,))
    assert 0.0 <= sd < 0.05                                 # resolved peaks -> tiny window spread


def test_mass_error_ppm_and_captured_fraction():
    base, step = 10000.0, 500.0
    m = np.arange(base - 300, base + step + 300, 1.0)
    # dominant = conjugate, its apex shifted +5 Da from the anchor (~476 ppm at 10500)
    md = np.column_stack([m, _gauss(m, base, 10, 3) + _gauss(m, base + step + 5, 100, 3)])
    q = da._mass_quality(md, base, step, 1, (0.0,))
    assert q["mass_error_ppm"] is not None and 350 < q["mass_error_ppm"] < 600
    assert q["captured_fraction"] > 0.9                     # all signal is in the two bands


def test_captured_fraction_flags_unmodeled_signal():
    base, step = 10000.0, 500.0
    m = np.arange(base - 1000, base + step + 2500, 1.0)
    # two anchored peaks plus a large broad hump far outside the bands (e.g. deconvolution baseline)
    md = np.column_stack([m, _gauss(m, base, 30, 3) + _gauss(m, base + step, 60, 3)
                          + _gauss(m, base + step + 1500, 200, 150)])
    q = da._mass_quality(md, base, step, 1, (0.0,))
    assert q["captured_fraction"] < 0.5                     # most signal is the unmodelled hump


def test_uv_concentration(monkeypatch):
    # synthetic UV trace: flat 2 mAU baseline + triangular peak (apex 8.5, 8.0-9.0, +100),
    # so the baseline-subtracted peak area is a known 0.5*base*height = 0.5*1*100 = 50 mAU*min
    for k in ("EPS280", "FLOW_ML_MIN", "INJ_UL", "DAD_UNIT", "MW_DA", "PATH_CM", "DILUTION"):
        monkeypatch.delenv(k, raising=False)
    t = np.arange(2.0, 12.0, 0.01)
    peak = np.clip(100.0 * (1.0 - np.abs(t - 8.5) / 0.5), 0.0, None)
    uv = np.column_stack([t, 2.0 + peak])

    rel = da._uv_concentration(uv, base=14000.0, tmax=8.7)     # no env -> relative only
    assert abs(rel["uv280_peak_area"] - 50.0) < 1.0
    assert rel["uv280_area_unit"] == "mAU*min"
    assert rel["uv280_apex_min"] == 8.5
    assert "protein_conc_mg_ml" not in rel                     # absolute needs EPS280/flow/inj

    for k, v in {"EPS280": "5500", "FLOW_ML_MIN": "0.4", "INJ_UL": "5",
                 "DAD_UNIT": "mAU", "MW_DA": "14000"}.items():
        monkeypatch.setenv(k, v)
    ab = da._uv_concentration(uv, base=14000.0, tmax=8.7)
    # reproduce the documented formula from the reported area and check the code matches
    area_au = ab["uv280_peak_area"] / 1000.0                   # mAU -> AU
    conc_M = area_au * (0.4 / 1000.0) / (5500.0 * 1.0) / (5e-6)
    assert abs(ab["protein_conc_uM"] - conc_M * 1e6) < 0.5
    assert abs(ab["protein_conc_mg_ml"] - conc_M * 14000.0) < 0.01


def test_chromatograms_decodes_mzml(tmp_path):
    t = [0.0, 1.0, 2.0]
    i = [5.0, 9.0, 3.0]
    p = tmp_path / "mini.mzML"
    _write_mini_mzml(p, t, i)
    ch = da.chromatograms(str(p))
    assert "TIC" in ch
    arr = ch["TIC"]
    assert arr.shape == (3, 2)
    assert np.allclose(arr[:, 0], t) and np.allclose(arr[:, 1], i)


_DECON_UNSET = ("MASS_LB", "MASS_UB", "MODE", "MZ_LO", "MZ_HI", "Z_LO", "Z_HI",
                "MASSBINS", "TMIN", "TMAX", "SATELLITES", "REPLOT", "NO_CACHE")


def test_decon_key_ignores_dar_params_but_tracks_deconvolution(monkeypatch):
    for k in _DECON_UNSET:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("BASE_MASS", "10000"); monkeypatch.setenv("MOD_MASS", "500")
    k0 = da._decon_key(da.load_config())
    monkeypatch.setenv("MOD_MASS", "700")          # conjugate mass is post-processing only
    monkeypatch.setenv("SATELLITES", "0,162.05")   # glycoform offsets are post-processing only
    assert da._decon_key(da.load_config()) == k0   # cache stays valid
    monkeypatch.setenv("Z_LO", "5")                # charge range governs the deconvolution
    assert da._decon_key(da.load_config()) != k0   # cache now invalid


def test_process_resumes_from_cache_and_falls_back(tmp_path, monkeypatch):
    for k in _DECON_UNSET:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("BASE_MASS", "10000"); monkeypatch.setenv("MOD_MASS", "500")
    calls = []
    monkeypatch.setattr(da, "analyze", lambda p, o: calls.append("analyze") or {"file": "s"})
    monkeypatch.setattr(da, "replot",
                        lambda p, o, peak_width=None: calls.append("replot") or {"file": "s"})
    d, massf, keyf = da._cache_paths(str(tmp_path), "s")
    src = str(tmp_path / "s.mzML")
    da.process(src, str(tmp_path))                 # no cache -> analyze
    assert calls == ["analyze"]
    os.makedirs(d, exist_ok=True)
    pathlib.Path(massf).write_text("10000 100\n")
    json.dump({"key": da._decon_key(da.load_config()), "peak_width": 0.5}, open(keyf, "w"))
    da.process(src, str(tmp_path))                 # matching cache -> replot
    assert calls == ["analyze", "replot"]
    monkeypatch.setenv("Z_LO", "5")                # deconvolution param changed -> re-run
    da.process(src, str(tmp_path))
    assert calls == ["analyze", "replot", "analyze"]
    monkeypatch.delenv("Z_LO", raising=False)
    monkeypatch.setenv("NO_CACHE", "1")            # forced fresh even with a matching cache
    da.process(src, str(tmp_path))
    assert calls[-1] == "analyze"


def test_charge_support_counts_charge_states():
    mass = 10000.0
    mz = np.arange(600, 2100, 0.5)
    inten = np.zeros_like(mz)
    for z in (6, 7, 8, 9):
        inten += _gauss(mz, (mass + z * 1.00728) / z, 100.0, s=0.7)
    arr = np.column_stack([mz, inten])
    assert da._charge_support(arr, mass, 3, 20) == 4       # exactly the 4 charge states present
    assert da._charge_support(arr, mass + 137.0, 3, 20) < 4  # a wrong mass sheds support
    assert da._charge_support(None, mass, 3, 20) == 0


def test_worker_catches_errors(monkeypatch):
    monkeypatch.setattr(da, "process", lambda p, o: {"file": "ok"})
    assert da._worker(("x.mzML", "/tmp"))["file"] == "ok"

    def boom(p, o):
        raise RuntimeError("bad")
    monkeypatch.setattr(da, "process", boom)
    r = da._worker(("bad.mzML", "/tmp"))
    assert r["error"] == "bad" and r["file"] == "bad.mzML"
