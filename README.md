# Stellar-history identifiability for TOI-700 and LHS 1140

Reproducible code and derived data for testing how present-day stellar constraints identify cumulative EUV and stellar-wind histories in TOI-700 and LHS 1140.

This repository supports the manuscript:

> **Present-day stellar constraints leave multiple high-energy histories for TOI-700 and LHS 1140**

**Yuzhan Zhang**  
Independent Researcher, Japan  
ORCID: [0009-0000-3121-7972](https://orcid.org/0009-0000-3121-7972)

<!-- After the Zenodo v1.0.0 release, add the DOI badge here.
Example:
[![DOI](https://zenodo.org/badge/DOI/10.xxxx/zenodo.xxxxx.svg)](https://doi.org/10.xxxx/zenodo.xxxxx)
-->

## Overview

The atmospheric evolution of a temperate planet depends on its host star's **time-integrated high-energy and stellar-wind forcing**, not only on the star's present-day activity.

Rotation-conditioned stellar models can converge at late times even when their earlier X-ray, EUV, and wind histories differed substantially. This repository asks a specific inverse-problem question:

> **How strongly do present-day stellar constraints identify the cumulative forcing histories experienced by real habitable-zone planets?**

The analysis focuses on **TOI-700 d** and **LHS 1140 b**. It separates three distinct inference tasks:

1. **Current-state admissibility** — does a candidate stellar history reproduce the measured present state?
2. **Model adequacy** — do independent observables agree with the assumed forward activity/observation model?
3. **Historical identification** — how much uncertainty remains in the time-integrated forcing among admissible histories?

The public release is intentionally narrow. It contains only the analysis required for the paper's E5-C2 / E5-C3 / E5-C5 / E5-C6 chain rather than the full internal System-Conditioned Long-Term Habitability (SCLH) research archive.

## Main results

Within the adopted finite set of rotation-conditioned stellar histories:

- **TOI-700:** 194 rotation-compatible histories have a model-predicted current EUV spread of about **1.077**, while cumulative EUV and wind exposure span about **2.65** and **2.73**.
- The TOI-700 historical spread is not produced solely by combining neighboring mass grids. At fixed stellar-model mass, cumulative EUV still spans approximately **2.196** at 0.40 \(M_\odot\) and **2.184** at 0.45 \(M_\odot\).
- **LHS 1140:** 97 rotation-compatible histories have a model-predicted current EUV spread of about **1.0085**, while cumulative EUV and wind exposure span about **1.567** and **1.543**.
- Within the LHS 1140 rotation-matched history set, native model current \(L_X\) varies by only about **0.00382 dex**, providing little ability to rank the cumulative histories.
- The native MORS X-ray quantity and the published XMM luminosity use non-identical energy-band observation operators. Their raw discrepancy is therefore treated as an **activity/observation-mapping diagnostic**, not as a bandpass-homogenized physical rejection.

The stellar-track ensemble is treated as **deterministic finite model support**, not as an IID stellar sample and not as a population posterior.

## Scope and nonclaims

This repository does **not**:

- infer a unique atmospheric-loss history for either planet;
- assign a population prior to the MORS percentile tracks;
- claim a universal M-dwarf historical-spread factor;
- claim that present-day stellar data can never recover a unique history;
- treat the native LHS 1140 X-ray mismatch as a calibrated physical rejection;
- redistribute the original third-party MORS/Johnstone stellar-track archive.

The objective is narrower: to quantify what the adopted present-day constraints do and do not identify about the historical stellar forcing.

## Repository contents

```text
.
├── src/sclh/                 # minimal scientific utilities used by the paper
├── experiments/              # E5-C2, E5-C3, E5-C5, E5-C6
├── data/derived/             # frozen paper-supporting derived tables
├── tests/                    # focused unit tests + headline-result invariants
├── scripts/                  # verification, audit, and figure reproduction
├── figures/                  # publication figures
├── docs/                     # provenance and reproducibility guidance
├── .github/workflows/        # GitHub Actions test workflow
├── CITATION.cff
├── .zenodo.json
├── pyproject.toml
├── requirements.txt
├── SHA256SUMS.txt
└── LICENSE
```

The repository deliberately excludes unrelated SCLH experiments, internal research reports, conversation archives, obsolete manuscript drafts, submission-engineering files, and third-party raw stellar tracks.

## Quick start

Create an isolated Python environment:

```bash
python -m venv .venv
```

Activate it.

Linux/macOS:

```bash
source .venv/bin/activate
```

Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

Install the package and test dependencies:

```bash
python -m pip install --upgrade pip
pip install -e ".[dev]"
```

Run the focused test suite and release-integrity check:

```bash
pytest
python scripts/verify_release.py
```

Regenerate the five paper figures:

```bash
python scripts/make_paper_figures.py
```

## Reproduce the derived activity and information audits

The frozen E5-C2 and E5-C3 history tables are included. Therefore the E5-C5 and E5-C6 analyses can be reproduced without redistributing the third-party MORS track archive:

```bash
python scripts/reproduce_derived_audits.py
```

This reproduces the present-day activity cross-calibration and history-information diagnostics from the included history tables.

Generated outputs are written under:

```text
results/E5C5/
results/E5C6/
```

The `results/` directory is ignored by Git so reproduction does not modify the archived public-release files.

## Full raw-track reproduction

Reproducing E5-C2 and E5-C3 from the upstream stellar tracks requires the external rotation/XUV grid associated with:

C. P. Johnstone, M. Bartel & M. Güdel,  
*The active lives of stars: a complete description of rotation and XUV evolution of F, G, K, and M dwarfs*,  
**Astronomy & Astrophysics 649, A96 (2021)**  
DOI: `10.1051/0004-6361/202038407`

The associated track archive is available from Zenodo:

**DOI: `10.5281/zenodo.4266670`**

After obtaining and extracting the archive, set the track-grid path.

Linux/macOS:

```bash
export SCLH_MORS_ROOT=/path/to/TrackGrid_MstarPercentile
```

Windows PowerShell:

```powershell
$env:SCLH_MORS_ROOT = "C:\path\to\TrackGrid_MstarPercentile"
```

Then run:

```bash
python experiments/e5c2_toi700d_mors_history_inversion.py
python experiments/e5c3_lhs1140b_crosssystem_history.py
python experiments/e5c5_activity_mapping_crosscalibration.py
python experiments/e5c6_multiobservable_history_information.py
```

The original third-party stellar-track archive is **not redistributed** in this repository.

## Reproducibility layers

The release deliberately distinguishes two reproducibility levels.

### 1. Self-contained paper-level audit

Using the included derived history tables, the repository can independently reproduce:

- the headline historical-spread ratios;
- the TOI-700 fixed-mass controls;
- the LHS 1140 present-activity comparison;
- the E5-C6 diagnostic weighting/history-information audit;
- all five paper figures;
- the release checksums and regression tests.

No third-party track download is required for this layer.

### 2. Upstream raw-track regeneration

With the external Johnstone/MORS track archive installed, the repository can regenerate the TOI-700 and LHS 1140 history tables from the upstream stellar evolutionary tracks.

This separation keeps the numerical provenance explicit without silently redistributing third-party scientific assets.

## Key machine-readable outputs

The compact summary of the paper-level numerical claims is:

```text
data/derived/paper_results.json
```

The principal derived history tables are:

```text
data/derived/toi700_rotation_conditioned_histories.csv
data/derived/lhs1140_rotation_compatible_histories.csv
data/derived/lhs1140_activity_model_comparison.csv
data/derived/lhs1140_weight_sensitivity.csv
data/derived/observable_information_audit.csv
```

## Validation

The public release was validated before packaging with:

```bash
pytest
python scripts/verify_release.py
python scripts/reproduce_derived_audits.py
python scripts/make_paper_figures.py
```

Release validation status:

- focused tests: **18 passed**;
- SHA-256 release manifest: **verified**;
- E5-C5/E5-C6 derived reproduction: **passed**;
- paper-figure regeneration: **passed**.

The exact tested Python/package environment is recorded in:

```text
environment-tested.txt
```

## Provenance

The numerical scientific baseline is inherited from the frozen internal research snapshot:

```text
SCLH_v0.1.18_E5C6_dual_manuscript_snapshot.zip
```

SHA-256:

```text
75a40de65b56642ed934a2ee21923c17c46db0178084f4c27cbd5de3b6885e80
```

The public v1.0.0 release removes unrelated project material and applies one interpretation-only cleanup: the historical E5-C3 wording that treated the raw LHS 1140 X-ray interval mismatch as a hard rejection is replaced by the later bandpass-aware activity/observation-mapping interpretation.

The numerical E5-C2/E5-C3 history tables are unchanged.

See:

```text
docs/PROVENANCE.md
docs/REPRODUCIBILITY.md
```

for the detailed provenance and reproduction policy.

## Citation

The repository includes a machine-readable citation file:

```text
CITATION.cff
```

Before the manuscript receives a persistent article identifier, cite the archived software release.

After the GitHub `v1.0.0` release is archived on Zenodo:

1. add the Zenodo DOI to this README;
2. add the DOI to `CITATION.cff`;
3. add a Zenodo DOI badge near the top of this README;
4. insert the GitHub URL and Zenodo DOI into the manuscript's **Code availability** statement.

## Zenodo metadata

A draft Zenodo record description is included in:

```text
.zenodo.json
```

Recommended release tag:

```text
v1.0.0
```

Recommended Zenodo resource type:

```text
Software
```

The Zenodo creator metadata should remain:

```text
Zhang, Yuzhan
ORCID: 0009-0000-3121-7972
Affiliation: Independent Researcher, Japan
```

## License

Repository-authored code is released under the **MIT License**.

Third-party MORS/Johnstone data and any other external scientific assets remain governed by their original distribution terms and are not included in this repository.
