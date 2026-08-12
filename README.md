# dar-from-lcms

DAR (drug/modifier-to-protein ratio) from intact-mass LC-MS, via UniDec, in one command:

    BASE_MASS=<unmodified avg mass Da> MOD_MASS=<Da added per conjugation> ./dar /path/to/run_dir

Converts Agilent `.d` to mzML (ProteoWizard), deconvolves each spectrum, and writes the
unmodified and +1 masses, the conjugated fraction and DAR, plus a plot per sample. Written for a
single conjugation site (one reactive residue), where DAR is occupancy 0-1.

For how the measurement and the scripts work end to end, see [HOW-IT-WORKS.md](HOW-IT-WORKS.md).

## Setup

    brew install colima docker qemu lima-additional-guestagents   # macOS

The prebuilt UniDec/Wine binaries need a real x86_64 VM (Rosetta breaks Wine); the `dar` script
starts one automatically. On x86_64 Linux, set `COLIMA_CONTEXT=default` and skip Colima.

## Use

You must give the two masses that define the chemistry:

- `BASE_MASS` -- the unmodified protein's average mass (from its sequence).
- `MOD_MASS`  -- the mass added per conjugation (from the reagent; a thiol-maleimide adds the
  full reagent mass, since it is a Michael addition with no leaving group).

```
BASE_MASS=14000 MOD_MASS=480 ./dar ~/Downloads/run      # auto-detects the main elution peak
BASE_MASS=14000 MOD_MASS=480 TMIN=8.4 TMAX=8.7 ./dar ~/Downloads/run   # fix the window
```

Writes into the run dir:
- `dar_results.json` -- per sample: DAR, elution window, measured unmodified/conjugate masses,
  DAR across integration windows (robustness), and a UV-280 late-feature check.
- `dar_<sample>.png` -- two panels: UV-280 + MS TIC chromatograms (protein window shaded), and
  the deconvolved mass spectrum with the unmodified / +1 masses anchored at the expected values.

The masses are anchored at `BASE_MASS` and `BASE_MASS+MOD_MASS` (not free peak picking), so DAR is
read at the chemically expected positions. If the unmodified and conjugate species co-elute, UV
cannot separate them -- it confirms identity/purity, not an independent DAR; run an unmodified
control for the orthogonal check.

## Notes

- Deconvolve the elution-peak window, not the whole run, or the protein drowns in background.
- MS-area DAR assumes both species ionize alike; cross-check against UV for QC numbers.
- Common adducts (`+16` oxidation, `+178` His-tag gluconoylation) are not extra conjugation.
- Run an unmodified control to fix the base mass and confirm the mass step.

## Image notes

UniDec ships two prebuilt binaries with conflicting needs: `isogen.so` wants GLIBCXX_3.4.32
(GCC 14) and `unideclinux` wants `libhdf5_serial.so.103` (HDF5 1.10). Ubuntu 24.04 with the
explicit `libhdf5-103-1t64` runtime satisfies both (not `libhdf5-dev`, which pulls HDF5 1.14 /
`.so.310`). Install UniDec from git, not the stale PyPI release. Keep the `pip install` layer
before the `COPY`, and don't `docker builder prune`, or you evict the ~10-min UniDec build.
