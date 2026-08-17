# Output specification

Every run writes these into the run directory (`OUTDIR`, the mounted `/data`).

## `dar_results.json`

A JSON array, one object per input file. Fields:

| field | type | meaning |
|---|---|---|
| `file` | string | sample name |
| `window_min` | [float, float] | elution window used, min |
| `DAR` | float | two-state DAR = area(+1)/(area(0)+area(+1)) at the +-25 Da band |
| `conjugated_pct` / `unmodified_pct` | float | DAR and 1-DAR as percentages |
| `dar_sd` | float or null | uncertainty: spread of DAR across +-15/25/40 Da windows |
| `mass_error_ppm` | float or null | observed vs theoretical mass of the dominant state |
| `captured_fraction` | float | fraction of deconvolved signal inside the anchored bands (low = caution) |
| `charge_support` | list[int] | per load state (index = load number), how many charge states independently back that mass |
| `charge_support_dominant` | int | `charge_support` for the most abundant state (real vs harmonic/artifact check) |
| `DAR_by_window` | list | DAR and band areas at half-widths 15/25/40/60 Da (robustness) |
| `conj_apex` / `naked_apex` | float or null | measured apex mass of a species present at >=10% |
| `expected` | [float, float] | anchors: `BASE_MASS`, `BASE_MASS+MOD_MASS` |
| `trace` | string or null | note when a species is <10% and not separately resolved |
| `peak_width_mz` | float or null | UniDec auto peak width (analyze only) |
| `UV280_main_h`, `UV280_late_h`, `UV_late_over_main` | float or null | UV-280 main/late-feature check |

UV-280 concentration (only when `EPS280`+`FLOW_ML_MIN`+`INJ_UL` are set):

| field | meaning |
|---|---|
| `uv280_peak_area` (+ `uv280_area_unit`), `uv280_apex_min` | UV-280 main-peak area / apex |
| `protein_conc_uM`, `protein_conc_mg_ml` | absolute concentration |
| `conc_inputs` | echo of eps280 / path_cm / inj_ul / dilution / flow_ml_min / dad_unit / mw_da |

Multi-state (only when `DAR_MAX_N>1`):

| field | meaning |
|---|---|
| `average_dar` (+ `average_dar_sd`) | Sum(n*An)/Sum(An), with window-spread uncertainty |
| `dispersity` | Mw/average (1.0 = monodisperse; higher = broader) |
| `dar_states` | `DAR_MAX_N` |
| `satellites` | glycoform/adduct offsets folded into each state |
| `dar_state_frac`, `dar_state_areas` | per-state fractions and areas, index = load number |

Error rows are `{"file": ..., "error": ...}`.

## `dar_summary.csv`

One tidy row per sample for plate/batch QC: `file, DAR, dar_sd, average_dar, average_dar_sd,
dispersity, mass_error_ppm, captured_fraction, charge_support_dominant, protein_conc_mg_ml,
window_min, error`.

## Caching and parallelism

The expensive step (UniDec deconvolution) depends only on the m/z window, charge range, mass
bounds, mass-bin size and elution window. `process()` writes a small `_darcache.json` next to the
cached mass distribution and reuses it on the next run whenever those are unchanged, so
re-anchoring (`MOD_MASS`), glycoform offsets (`SATELLITES`), the mirror and any figure change
re-plot instantly with no deconvolution. `REPLOT=1` forces replot; `NO_CACHE=1` forces a fresh
run. For a shared-config folder, `JOBS=N` runs N samples in parallel (the manifest path stays
serial).

## Figures

- `dar_<sample>.png` / `.pdf` : three panels (chromatograms; charge-state envelope; deconvolved
  mass with DAR/average-DAR, uncertainty, mass-error line, and shaded bands). With `MIRROR_MASS`,
  panel c shows the control mirrored below the treated sample.
- `dar_plate_heatmap.png` : written for multi-state batch runs of >=2 samples; a samples x
  load-state fraction heatmap.

## Accuracy note

The `±25 Da` band (and its `SATELLITES` extension) is exact when load states are resolved. For
native glycosylated intact IgG the states are broad, overlapping glycoform envelopes; set
`SATELLITES` to the glycoform offsets to fold each state's envelope, but treat the result as an
approximation of the reference method's glycoform-specific integration, and always sanity-check
`captured_fraction` and the unmodified-control mirror.
