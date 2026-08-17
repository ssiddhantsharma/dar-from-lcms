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

```
# download one or more .raw from the Zenodo record into ./raw, then:
MODE=native DAR_MAX_N=6 SATELLITES=0,162.05,324.11,486.16 \
  BASE_MASS=148057 MOD_MASS=769.7 \
  ./dar ./raw            # mal-DFO example: MOD_MASS = chelator 713.83 + iron 55.85
```

Per-conjugate `MOD_MASS` = chelator mass + iron (+55.85): mal-DFO 769.7, N-suc-DFO 700.6,
NCS-DFO 810.8, TMTHSI-DFO 940.9. `BASE_MASS` = trastuzumab base glycoform (confirm from sequence).

## Status

**Pending the deconvolution run** (native-IgG deconvolution needs the QEMU/UniDec pipeline and a
round or two of parameter tuning). This is an *approximation* of the paper's glycoform-specific
method, so the goal is to report how closely `average_dar` tracks the published CAR, honestly, with
the config used. Results (a table of ours vs paper, and the figures) will be committed here.
