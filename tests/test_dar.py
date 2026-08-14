"""Unit tests for the analytical core of dar_auto.

These cover the pure functions (no UniDec/Docker needed): mzML chromatogram parsing,
band integration, peak apex, charge-state assignment, elution-window detection, and the
DAR math. Run with `pytest`.
"""
import base64
import importlib.util
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
