# Validation against public data

This folder checks the tool against an independent, open dataset. It is an honest check, including
where the tool's simple method does not fully transfer to a harder regime.

## Dataset

Van der Zon & Gargano (University of Amsterdam) deposited raw native SEC-HRMS data for trastuzumab
and ipilimumab DFO immunoconjugates, CC-BY-4.0:

- Data: Zenodo **10.5281/zenodo.17637726** (8 Thermo `.raw` files, ~190 MB, MD5-checked)
- Paper: van der Zon et al., *Anal. Chim. Acta* 1395 (2026) 345214,
  doi:10.1016/j.aca.2026.345214 (open access)

This is a chelator-to-antibody ratio (CAR) dataset on **native, glycosylated, intact IgG**. That is
a harder regime than the tool's core use (small, deglycosylated binders by denaturing RP, where the
load states are baseline-resolved and DAR is exact). It exercises `MODE=native` and `SATELLITES`.

Per the paper and the Zenodo record, the antibody was run across 200-1000 mM ammonium acetate, and
the reported SEC-MS CAR depends on that buffer concentration. The paper derives averageCAR by a
glycoform-specific integration of the deconvolved spectrum. See the paper for its reported values;
they are not reproduced here to avoid citing numbers second-hand.

## What the tool does here (measured, reproducible)

1. **Conversion.** ThermoRawFileParser converts all 8 `.raw` with no Wine (under the cool Rosetta
   VM it needs `MONO_ENV_OPTIONS=--interp`, which the `dar` wrapper sets). This part works.
2. **Deconvolution.** UniDec in `MODE=native` recovers the ~148 kDa envelope. The elution window
   must be set explicitly (`TMIN`/`TMAX`) because these runs end with a large salt/wash base-peak
   spike that the auto-window would otherwise pick over the protein SEC peak (~13-15.5 min).

For `Trastuzumab-N-suc-DFO_600 mM ammonium acetate` with the SEC window and literature anchors
(`BASE_MASS=148057`, `MOD_MASS=770`, glycoform `SATELLITES`), the tool reports **average_dar ~= 1.86**
(dispersity ~= 1.49), and its `captured_fraction` is low (~0.24).

## Honest limitation (why this is an approximation, not a clean match)

The native deconvolution is a continuous, glycoform-blurred envelope with **no baseline-resolved
chelator ladder**. The tool anchors bands at `BASE_MASS + n*MOD_MASS` (plus glycoform satellites);
with no resolved ladder to lock onto, `average_dar` is anchor-sensitive: sweeping plausible values
moves it from ~1.7 to ~2.4, and the coverage of total signal stays flat (~20%) with no maximum, so
there is no data-driven way to fix the anchors from these spectra. The low `captured_fraction`
correctly flags this. A faithful native CAR needs the paper's glycoform-specific integration, which
this general-purpose anchor-and-integrate tool does not implement (and should not, to stay small).

Reproduce:

```
# download the record into ./raw (8 files, ~190 MB), then:
ANALYZE_CONTEXT=colima-rosetta \
MODE=native DAR_MAX_N=6 SATELLITES=0,162.05,324.11,486.16 \
  BASE_MASS=148057 MOD_MASS=770 TMIN=13 TMAX=15.5 \
  ./dar ./raw
```

## Takeaway

The Thermo conversion path is validated on real vendor data. On native glycosylated IgG the tool
gives an approximate CAR that its own QC flags as low-confidence; its exact, validated regime is
denaturing intact-mass DAR of resolved binders. This is recorded as a known boundary, not hidden.
