# Citations and credits

This tool stands on other people's work. If you use it, please cite the underlying methods,
and note below exactly what was borrowed from other software (and what was changed).

## Please cite

**UniDec** (the deconvolution engine):
> Marty, M.T. et al. Bayesian Deconvolution of Mass and Ion Mobility Spectra. *Anal. Chem.*
> 2015, 87 (8), 4370-4376. doi:10.1021/acs.analchem.5b00140

**UniDec Processing Pipeline (UPP)** (biotherapeutic MS analysis):
> Phung, W. et al. UniDec Processing Pipeline for Rapid Analysis of Biotherapeutic Mass
> Spectrometry Data. *Anal. Chem.* 2023, 95 (30), 11491-11498. doi:10.1021/acs.analchem.3c02010

**ProteoWizard** (vendor raw to mzML conversion):
> Chambers, M.C. et al. A cross-platform toolkit for mass spectrometry and proteomics.
> *Nat. Biotechnol.* 2012, 30, 918-920. doi:10.1038/nbt.2377

## Methods implemented here

- **Multi-state average DAR/CAR and the dispersity index** (Eqs 1-3 in `dar_distribution`) follow
  van der Zon, A.A.M.; Weijers, B.; Vugts, D.J.; Gargano, A.F.G. *Microscale native SEC-HRMS for
  the determination of the chelator-to-antibody ratio of immunoconjugates.* Anal. Chim. Acta 2026,
  1395, 345214. doi:10.1016/j.aca.2026.345214 (open access, CC BY). We implement their averageCAR,
  weight-average, and dispersity definitions; we do **not** reproduce their tuned UniDec parameters.
- **UV-280 extinction coefficient** uses the Pace convention (5500*Trp + 1490*Tyr + 125*disulfide);
  Pace, C.N. et al. *Protein Sci.* 1995, 4, 2411-2423.

## Software patterns borrowed (reimplemented, not copied)

No source code from the projects below is copied into this repo; the listed *patterns* were
reimplemented for this tool's data (deconvolved protein mass, DAR load states). Licenses noted.

- **spectrum_utils** (Bittremieux et al., **Apache-2.0**), github.com/bittremieuxlab/spectrum_utils
  - Borrowed patterns: a `colors` dict keyed by species with a `None` fallback; a per-axis styling
    helper (`_format_ax`: grid behind data, faint minor grid, small ticks, italic m/z); a pluggable
    label formatter (their `annot_fmt`); the `mass_errors` observed-vs-theoretical visualization; and
    the `mirror()` head-to-tail (treated-on-top, control-below) convention.
  - What we changed: reimplemented for DAR/CAR load states (unmodified / +n / adduct / control)
    rather than peptide fragment ion types, over a continuous deconvolved-mass trace rather than
    centroided peaks. `mass_errors` is scaled to our single dominant state (a text annotation), not a
    per-fragment bubble plot.
- **SmartPeak** (AutoFlowResearch, **MIT**), github.com/AutoFlowResearch/SmartPeak
  - Borrowed patterns: the per-sample input manifest (`sequence.csv`) driving a batch, and reporting
    QC as tabular metrics.
  - What we changed: our manifest is a small CSV mapping columns to this tool's env vars
    (`run_manifest`), and our QC is the per-sample `mass_error_ppm` / `captured_fraction` fields.
- **nf-core** conventions, as seen in **quantms** and **mhcquant** (both **MIT**),
  github.com/bigbio/quantms and github.com/nf-core/mhcquant
  - Borrowed patterns (repo organization only): `CHANGELOG.md`, this `CITATIONS.md`, an input schema
    (`assets/manifest_schema.json`, their `assets/schema_input.json`), and an output specification
    (`docs/output.md`). No Nextflow/pipeline code is used.

## This tool

> Sharma, S. dar-from-lcms: an open, reproducible pipeline for DAR/CAR from intact-mass LC-MS.
> github.com/ssiddhantsharma/dar-from-lcms
