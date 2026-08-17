# Validation against public data

This folder holds validation of the tool against an independent, open dataset, so the DAR/CAR
numbers can be checked by anyone.

## Dataset

Van der Zon & Gargano (University of Amsterdam) deposited raw native SEC-HRMS data for
trastuzumab and ipilimumab DFO immunoconjugates, CC-BY-4.0:

- Data: Zenodo **10.5281/zenodo.17637726** (8 Thermo `.raw` files)
- Paper (ground-truth CAR): van der Zon et al., *Anal. Chim. Acta* 1395 (2026) 345214,
  doi:10.1016/j.aca.2026.345214 (open access)

This is a **chelator**-to-antibody ratio (CAR) dataset on native, glycosylated intact IgG, a
different regime from the tool's core (small, deglycosylated binders by denaturing RP). It exercises
the `MODE=native` + `SATELLITES` (glycoform-aware) path.

## Ground truth (from the paper, at 600 mM ammonium acetate)

| conjugate | published averageCAR |
|---|---|
| trastuzumab-N-suc-DFO | 2.48 (2.04 / 2.43 / 2.48 / 2.06 across 200/400/600/1000 mM) |
| trastuzumab-mal-DFO | 2.49 |
| trastuzumab-TMTHSI-DFO | 0.79 |
| trastuzumab-NCS-DFO, ipilimumab-N-suc-DFO | in the paper's SI (Table S11) |

Method-spread reference: same sample by radiometric titration 1.68, SEC-UV 1.89, SEC-MS 1.41,
which the authors call good agreement (so +-0.2 to 0.5 is the expected inter-method spread).

## How to reproduce

The 8 files are Thermo `.raw`, so conversion uses ThermoRawFileParser (no Wine/QEMU) and the whole
run stays on the cool Rosetta VM:

```
# download the record into ./raw (8 files, ~190 MB), then:
ANALYZE_CONTEXT=colima-rosetta \
MODE=native DAR_MAX_N=6 SATELLITES=0,162.05,324.11,486.16 \
  BASE_MASS=148057 MOD_MASS=769.7 \
  ./dar ./raw
```

`BASE_MASS` and `MOD_MASS` above are *starting estimates*, not verified numbers: `BASE_MASS` is
the trastuzumab G0F/G0F glycoform (literature), and `MOD_MASS` per chelator is an estimate that
also assumes the conjugate carries iron (+55.85 Da), which the paper's SI would confirm or refute.
Rather than trust these blind, the run reads the actual load-state spacing straight off the
deconvolved mass distribution (the gap between adjacent load peaks) and anchors to the observed
unmodified glycoform, then reports the config actually used. Estimate table (chelator + Fe): mal-DFO
769.7, N-suc-DFO 700.6, NCS-DFO 810.8, TMTHSI-DFO 940.9.

## Status

**Pending the deconvolution run.** Native-IgG deconvolution needs the UniDec pipeline and a round
or two of parameter tuning, and this is an *approximation* of the paper's glycoform-specific
integration. The plan: deconvolve one file, read the observed base mass and load-state spacing to
fix the anchors from the data (not the estimates above), then run all 8 and report a table of ours
vs published CAR, honestly, with the figures and the exact config. Results will be committed here.
