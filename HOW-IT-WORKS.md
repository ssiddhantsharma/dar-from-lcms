A walkthrough of what this tool measures, how the measurement physically works, and how the two
scripts turn a folder of raw LC-MS runs into a DAR number and a figure.

Notation used below: **M** = unmodified protein mass, **Δ** = mass added per conjugation. Example
numbers (`M = 14000`, `Δ = 480`) are placeholders; you supply your own.

---

## 1. What it computes

For a protein with a **single reactive site** (e.g. one engineered cysteine plus a thiol-reactive
reagent), the "drug/modifier-to-protein ratio" (DAR) is just **occupancy**: the fraction of
molecules that carry the modifier. Each molecule is either unmodified (0) or singly modified (1),
so DAR runs from 0 to 1 and equals the conjugated fraction.

For a multi-site conjugation you read a **ladder** (0, +1, +2, …) and take the intensity-weighted
average. The default (`DAR_MAX_N=1`) is the single-site 0/1 case above; set `DAR_MAX_N` to the
highest load state and the tool reports the whole distribution, the **average DAR = Σ(n·Aₙ)/Σ(Aₙ)**,
and a **dispersity index** (how broad the load distribution is). See §5.

---

## 2. How the measurement physically works

1. **Electrospray (ESI)** puts many protons on each protein molecule at once. A single protein
   therefore does not appear as one peak; it appears as a **ladder of charge states** (+8, +9,
   +10, …), each at `m/z = (M + z·1.00728) / z`.
2. You cannot read the mass off that ladder directly, because each peak's charge `z` is unknown.
   **Deconvolution** (here, UniDec) takes the whole ladder and solves for the single neutral mass
   `M` that explains it, collapsing the charge envelope into a **zero-charge mass spectrum**
   (intensity vs. true mass in Da). That is the bottom panel of the figure.
3. When the modifier attaches it adds **Δ** Da, so every charge state shifts and the deconvolved
   spectrum gains a **second peak at M + Δ**. The area of that peak relative to the total is the DAR.

So the chain is: charge ladder → true mass → which mass is unmodified, which is +1 → ratio.

---

## 3. The two layers

The split exists for one reason: UniDec and the Agilent file reader (via Wine) only run on
**x86_64 Linux**, but a Mac is ARM. A Docker container (through Colima's QEMU VM) *is* that Linux
box, and the `dar` script drives it from outside.

**`dar`** (bash, on your machine, the orchestrator). On `./dar <folder>`:

1. Resolve the folder; select the docker context.
2. Start the x86_64 VM if it is not already up (Wine does not run under Rosetta translation).
3. Ensure both images exist: build the analysis image from the Dockerfile if missing, pull the
   ProteoWizard image if missing.
4. For each `.d` folder with no matching `.mzML` yet, run the ProteoWizard container to convert
   `.d → .mzML` (the folder is mounted into the container as `/data`). Existing mzML are skipped,
   so re-runs do not repeat the conversion.
5. Run the analysis container once over all the mzML, passing the chemistry (`BASE_MASS`,
   `MOD_MASS`) as environment variables.
6. Outputs land back in your folder (it is the mounted `/data`).

**`dar_auto.py`** (Python, inside the container, the analysis). Per file:

1. Parse the **TIC and UV-280 chromatograms** straight out of the mzML.
2. Find the **protein elution window**: the largest TIC peak after the void, widened to where it
   falls below half height.
3. **Average only the scans in that window** into one m/z spectrum (the charge envelope).
4. **UniDec deconvolves** it into a zero-charge mass spectrum.
5. Anchored at **M** (unmodified) and **M + Δ** (+1), integrate the areas and compute
   `DAR = area(+1) / (area(unmodified) + area(+1))`. Repeat at several band widths for a
   robustness table.
6. Sanity-check the UV, draw the figure, append to `dar_results.json`.

---

## 4. The two numbers it must be told

Deconvolution alone just hands you a mass spectrum; it cannot know *which* peak is unmodified and
which is +1. You supply that:

- **`BASE_MASS` (M)**, the unmodified protein's average mass, from its sequence.
- **`MOD_MASS` (Δ)**, the mass added per conjugation, from the reagent. A thiol-maleimide adds the
  **full reagent mass** (it is a Michael addition, nothing leaves), so Δ is simply the reagent's
  molecular weight.

```
BASE_MASS=14000 MOD_MASS=480 ./dar ~/Downloads/run
```

These define where to look: unmodified at `M`, conjugate at `M + Δ`.

---

## 5. How the DAR is calculated, precisely

1. Auto-pick the elution window from the TIC.
2. Average the scans in the window into one m/z spectrum.
3. UniDec → zero-charge mass spectrum.
4. Integrate the deconvolved intensity in a band around `M` and around `M + Δ`.
5. `DAR = area(M+Δ) / (area(M) + area(M+Δ))`.

The masses are **anchored** at the expected positions rather than found by free peak-picking, so
the ratio is read at the chemically correct masses even when one peak is tiny. The tool reports the
DAR at several integration band widths; if the number barely moves across them, it is not an
artifact of one arbitrary window.

**Multi-site (`DAR_MAX_N > 1`).** Integrate a band at `M + n·Δ` for n = 0…N, then

```
average DAR = Σ (n · Aₙ) / Σ (Aₙ)
Mw          = Σ (n² · Aₙ) / Σ (n · Aₙ)
dispersity  = Mw / average DAR      (1.0 = every molecule carries the same load; higher = broader)
```

(van der Zon et al., Anal. Chim. Acta 1395, 2026, Eqs 1-3; see Credits in the README). The two-state
DAR is exactly the `N = 1` case. This band integration is exact when the load states are resolved
(as for small, deglycosylated binders). For native, glycosylated intact IgG the states are broad,
overlapping glycoform envelopes, so a fixed band is approximate, the reference method resolves this
by integrating specific glycoforms (G0F/G1F).

---

## 6. Why you run an unmodified control

With only the treated run, "DAR ≈ 0.9" would rest on three assumptions:

- that `M` really is your protein (not a truncation, a wrong construct, or an unexpected mod),
- that the `+Δ` shift is really the reagent and not a coincidental adduct,
- that a peak at `M + Δ` is genuine conjugation and not deconvolution noise.

An **unmodified control**, the same protein, run identically, that never saw the reagent, turns
each assumption into a measurement:

| The control proves | How |
|---|---|
| The unmodified mass is real, not just calculated | control deconvolves to `M` |
| The `+Δ` shift is *caused by the reaction* | untreated `M` → treated `M + Δ`, the only difference being the reagent |
| The assay's noise floor (what "0% conjugated" looks like) | untreated reads a few % at the `+1` position; the treated run's high value is real signal above that |
| Whether the two species co-elute (matters for UV, see below) | compare control and sample retention times |

This is the logic of a controlled experiment: the control isolates the one variable (conjugation)
by being identical in every other respect. Without it, `+Δ` is an inference; with it, `+Δ` is a
result.

---

## 7. UV-280 concentration (optional)

Deconvolution gives *ratios* (DAR); it does not tell you *how much* protein there is, and MS
intensity is **not** proportional to concentration (ionization efficiency differs per analyte, so
you cannot read concentration off the mass spectrum without a per-analyte standard). The
quantitative signal is the **UV-280 absorbance**, which obeys Beer-Lambert (`A = ε·c·l`). When you
supply the extinction coefficient and flow parameters, the tool reports concentration from the UV
peak:

- integrate the UV-280 peak over the elution window → area `A` (AU·min);
- moles through the detector = `A · F / (ε · l)`, `F` = flow rate, `ε` = molar extinction at 280 nm
  from the sequence (`5500·Trp + 1490·Tyr + 125·disulfide`), `l` = flow-cell path length;
- sample concentration = moles / injection volume × dilution; `mg/mL = molar × MW`.

Set `EPS280`, `FLOW_ML_MIN`, `INJ_UL` (optionally `PATH_CM`, `DILUTION`, `MW_DA`, `DAD_UNIT`). **No
physical standard is needed**, the sequence-derived `ε` is the standard. Two traps: the flow rate
is usually not in the mzML (read it from the LC method), and Agilent exports **mAU** even when the
file labels the array "absorbance unit" (mis-reading it as AU is a 1000× error, sanity-check a
220 nm channel: values in the hundreds are mAU). Assumes one pure, resolved peak and a linear
detector.

---

## 8. Honest limits

- **The exact percentage assumes equal ionization** of the unmodified and modified species. This is
  reasonable for a single small modification but is not proven by the measurement; treat the number
  as "high occupancy, ≈ X%," not a certified stoichiometry. For a certified number, add an
  orthogonal method (UV with resolved peaks, or a modifier-specific assay).
- **Intact MS proves "+Δ Da added at the site," not the modifier's structure.** Another reagent of
  the same nominal mass would look identical. MS confirms the mass is consistent with your reagent
  and quantifies how much attached.
- **Adducts (e.g. +16 oxidation, +178 gluconoylation of a His-tag) are assigned by mass**, not by
  fragmentation. They are common and do not affect the DAR.
- **Deconvolve the elution window, not the whole run.** Averaging the entire run buries the protein
  in solvent/background and produces a noisy "grass" spectrum.

