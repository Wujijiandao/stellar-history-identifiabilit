# Stellar-history identifiability for TOI-700 and LHS 1140

Code and derived data supporting the manuscript:

**Present-day stellar constraints leave multiple high-energy histories for TOI-700 and LHS 1140**

Author: **Yuzhan Zhang** — Independent Researcher, Japan  
ORCID: **0009-0000-3121-7972**

## Scientific scope

This repository tests how strongly present-day stellar constraints identify the **time-integrated high-energy and wind forcing histories** experienced by the temperate planets TOI-700 d and LHS 1140 b.

The public release is intentionally narrow. It contains only the analysis required for the paper's E5-C2 / E5-C3 / E5-C5 / E5-C6 chain.

It does **not** contain the full System-Conditioned Long-Term Habitability (SCLH) research archive, unrelated E0-E4 experiments, internal research reports, conversation archives, obsolete manuscript drafts, or third-party stellar-track files.

## Headline result

Within the adopted finite set of rotation-conditioned stellar histories:

- TOI-700: 194 rotation-compatible histories have a model-predicted current EUV spread of about 1.077, while cumulative EUV and wind exposure span about 2.65 and 2.73.
- LHS 1140: 97 rotation-compatible histories have a model-predicted current EUV spread of about 1.0085, while cumulative EUV and wind exposure span about 1.567 and 1.543.
- The TOI-700 cumulative-history spread persists within each fixed stellar-mass grid.
- For LHS 1140, native current model \(L_X\) varies very little across the rotation-matched set, so a present X-ray datum primarily audits the activity/observation mapping rather than ranking the latent cumulative histories.
- The model and published XMM value use non-identical X-ray observation operators; the native mismatch is therefore **not** presented as a bandpass-homogenized physical rejection.

The track ensemble is treated as **deterministic finite model support**, not as an IID stellar sample or a population posterior.

## Quick start

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell
# .\.venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
pip install -e ".[dev]"

pytest
python scripts/verify_release.py
python scripts/make_paper_figures.py
```

## Reproduce the paper's derived activity/information audits

The frozen E5-C2 and E5-C3 history tables are included, so E5-C5 and E5-C6 can be reproduced without redistributing the third-party MORS track archive:

```bash
python scripts/reproduce_derived_audits.py
```

Outputs are written under `results/E5C5/` and `results/E5C6/` and are ignored by Git.

## Full raw-track reproduction

E5-C2 and E5-C3 require the external rotation/XUV track grid associated with:

C. P. Johnstone, M. Bartel & M. Güdel, *The active lives of stars: a complete description of rotation and XUV evolution of F, G, K, and M dwarfs*, A&A 649, A96 (2021).

Track archive:
**Zenodo DOI: 10.5281/zenodo.4266670**

After obtaining and extracting the archive, set:

```bash
export SCLH_MORS_ROOT=/path/to/TrackGrid_MstarPercentile
```

or on PowerShell:

```powershell
$env:SCLH_MORS_ROOT = "C:\path\to\TrackGrid_MstarPercentile"
```

Then run:

```bash
python experiments/e5c2_toi700d_mors_history_inversion.py
python experiments/e5c3_lhs1140b_crosssystem_history.py
```

The original third-party track archive is **not redistributed** here.

## Repository layout

```text
.
├── src/sclh/                 # minimal scientific utilities used by B1
├── experiments/              # E5-C2, C3, C5, C6
├── data/derived/             # frozen paper-supporting derived tables
├── tests/                    # focused unit + headline-result invariants
├── scripts/                  # verification and figure reproduction
├── figures/                  # publication figures
├── docs/                     # provenance and release guidance
├── CITATION.cff
├── .zenodo.json
├── pyproject.toml
└── LICENSE
```

## Reproducibility layers

There are deliberately two reproducibility levels:

1. **Derived-table reproduction:** fully self-contained in this repository. It reproduces the E5-C5/E5-C6 audits and all paper figures from included derived histories.
2. **Raw-track reproduction:** requires the third-party MORS/Johnstone track archive and reproduces E5-C2/E5-C3 from those external scientific assets.

This separation avoids silently redistributing third-party data while keeping the paper's numerical provenance explicit.

## Version provenance

The numerical scientific baseline is inherited from the frozen internal snapshot:

`SCLH_v0.1.18_E5C6_dual_manuscript_snapshot.zip`

SHA-256:

`75a40de65b56642ed934a2ee21923c17c46db0178084f4c27cbd5de3b6885e80`

The public v1.0.0 release removes unrelated project material and applies one interpretation-only cleanup: the old E5-C3 "hard X-ray rejection" wording is replaced by the later bandpass-aware activity/observation-mapping audit. The numerical history tables are unchanged.

See `docs/PROVENANCE.md`.

## Citation

Before the manuscript has a persistent article identifier, cite the software release using `CITATION.cff`.

After you create the Zenodo release, add the Zenodo DOI badge and DOI to `CITATION.cff` and the manuscript's Code availability section.

## License

Repository-authored code is released under the MIT License. Third-party MORS data and any other external assets remain governed by their original terms and are not included in this repository.
