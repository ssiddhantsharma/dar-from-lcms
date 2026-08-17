# Changelog

All notable changes to this project are documented here. Format based on
[Keep a Changelog](https://keepachangelog.com/); this project uses date-based notes.

## Unreleased

### Added
- **Multi-state DAR/CAR** (`DAR_MAX_N`): full load-state distribution, average DAR, and the
  dispersity index (van der Zon et al. 2026, Eqs 1-3). Two-state DAR is the `DAR_MAX_N=1` case.
- **Glycoform / adduct folding** (`SATELLITES`): sum each load state across its glycoform and
  adduct satellite offsets, the honest treatment for native glycosylated IgG.
- **Native regime preset** (`MODE=native`): native SEC-MS deconvolution settings (m/z 2000-8000,
  higher charge range), each value overridable (`MZ_LO`, `MZ_HI`, `Z_LO`, `Z_HI`, `MASSBINS`).
- **UV-280 concentration** (`EPS280`, `FLOW_ML_MIN`, `INJ_UL`, ...): Beer-Lambert protein
  concentration from the UV peak, annotated on the figure and in the JSON.
- **DAR uncertainty**: the average DAR now carries a spread across integration windows (shown as
  `± sd` on the figure and `dar_sd` in the JSON).
- **Mass-accuracy and deconvolution-quality**: `mass_error_ppm` and `captured_fraction` trust
  checks, printed and stored per sample.
- **Mirror plot** (`MIRROR_MASS`): head-to-tail overlay of a treated sample and its unmodified
  control on the deconvolved-mass panel.
- **Batch/plate mode**: an input manifest CSV (`MANIFEST`, one row per sample with its own
  chemistry) plus a `dar_summary.csv` output and a `dar_plate_heatmap.png` distribution figure.
- Figure styling refresh (consistent axis styling, species colour map, mass-error line), and
  reproducible example figures (`make_demo.py`, `make_multidar_demo.py`, `make_mirror_demo.py`).
- **Resume / cache-reuse**: `process()` auto-reuses a cached deconvolution when the parameters that
  govern it are unchanged (re-anchoring, glycoform offsets and figure tweaks then cost no UniDec
  run), and only re-deconvolves when the cache is missing or its key differs. `REPLOT=1` forces
  replot, `NO_CACHE=1` forces a fresh run.
- **Batch parallelism**: `JOBS>1` runs a shared-config folder over a process pool; the manifest
  path stays serial because it mutates process-global env per row.
- **Charge-state support** (`charge_support`, `charge_support_dominant`): how many charge states
  independently back each DAR species, a real-vs-artifact check (FLASHDeconv idea, adapted to
  charge-resolved data; see `CITATIONS.md`).
- **Thermo .raw without Wine**: the `dar` wrapper converts Thermo `.raw` with ThermoRawFileParser
  (Hulstaert et al. 2020), which runs cross-platform, so a Thermo-only run stays on the cool VM and
  never starts the QEMU/Wine VM. Other vendors still convert with ProteoWizard. `CONVERTER=pwiz`
  forces the old path.
- `CITATIONS.md`, `assets/manifest_schema.json`, `docs/output.md`, `validation/` scaffold,
  `pyproject.toml` (pip-installable analysis core), expanded unit tests.

### Notes
- Backward compatible: the default (no new env vars) is the original two-state DAR.
- `SATELLITES` band integration is exact for resolved states; for native glycosylated IgG it is an
  approximation of the reference method's glycoform-specific integration (see `docs/output.md`).

## v1.0

- Initial release: two-state DAR from intact-mass LC-MS via ProteoWizard + UniDec, containerized,
  with a per-sample figure and `dar_results.json`.
