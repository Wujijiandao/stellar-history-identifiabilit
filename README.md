# Stellar-history identifiability for TOI-700 and LHS 1140

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21974562.svg)](https://doi.org/10.5281/zenodo.21974562)

Reproducible code and derived data supporting the manuscript:

**Present-day stellar constraints leave multiple high-energy histories for TOI-700 and LHS 1140**

**Yuzhan Zhang** — Independent Researcher, Beijing, China  
ORCID: **0009-0000-3121-7972**

## Scientific question

Planetary atmospheric evolution depends on time-integrated X-ray, EUV and stellar-wind forcing, while real host stars are usually constrained through measurements made at the present epoch. This repository asks an inverse-problem question:

> How strongly do present-day stellar constraints identify the cumulative forcing histories experienced by TOI-700 d and LHS 1140 b?

The stellar tracks are treated as **deterministic finite model support**, not as IID stellar samples and not as draws from a calibrated population posterior. The release therefore reports model-conditional identified-set widths rather than posterior credible intervals.

## Headline results

Within the adopted rotation-conditioned Johnstone/MORS support:

- **TOI-700:** 194 rotation-compatible histories have a model-predicted current-EUV max/min width of about **1.077**, while cumulative EUV and wind exposure have widths of about **2.65** and **2.73**.
- The TOI-700 historical spread persists at fixed stellar-model mass: cumulative-EUV widths are **2.196** at 0.40 \(M_\odot\) and **2.184** at 0.45 \(M_\odot\).
- **LHS 1140:** 97 rotation-compatible 0.20 \(M_\odot\) histories have a current-EUV width of about **1.0085**, while cumulative EUV and wind exposure have widths of **1.567** and **1.543**.
- The 0.15 \(M_\odot\) LHS 1140 model slice has **no** history that reaches the central 131-day rotation period by 12 Gyr, so support changes qualitatively across the bracketing mass grid and is not linearly interpolated.
- Native current model \(L_X\) varies by only about **0.00382 dex** across the LHS 1140 rotation-matched set; the present X-ray datum therefore mainly audits the activity/observation mapping rather than ranking cumulative histories.
- A direct native-model/XMM hard rejection is not claimed because the model and measurement use non-identical X-ray observation operators. Even using the observational 1-sigma upper edge, a bandpass-only reconciliation would require a factor of about **2.285**.

### v1.0.2: age-information leverage audit

Release v1.0.2 adds **E5-C7**, a value-of-information audit using only the frozen derived histories. It asks how much historical spread can remain if an *independent* stellar-age measurement restricts the model-conditional matching age to an interval of total width \(\Delta t\).

The diagnostic is deliberately worst-case over the unknown interval location. It is not a posterior and does not claim that the assumed age precision is observationally attainable.

Illustrative results:

- With a **100-Myr** age interval and the TOI-700 model mass fixed, the worst-case cumulative-EUV width falls to **1.119** (0.40 \(M_\odot\)) and **1.115** (0.45 \(M_\odot\)).
- The same **100-Myr** interval gives a worst-case cumulative-EUV width of **1.064** for the supported LHS 1140 0.20 \(M_\odot\) slice.
- For the **pooled** TOI-700 0.40/0.45 \(M_\odot\) support, the 100-Myr age interval still leaves a cumulative-EUV width of **1.731**. Age information alone therefore does not remove model-slice ambiguity.

This positive-control calculation shows that partial identification is **constraint-dependent**, rather than an outcome that the audit is forced to return.

## Scope and nonclaims

This repository does **not**:

- infer a unique atmospheric-loss history for either planet;
- assign a population prior to the MORS percentile tracks;
- claim a universal M-dwarf history-spread factor;
- claim that present-day stellar measurements can never recover a narrow history;
- treat the native LHS 1140 X-ray mismatch as a bandpass-homogenized physical rejection;
- redistribute the original third-party Johnstone/MORS stellar-track archive.

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
python scripts/reproduce_derived_audits.py
python scripts/make_paper_figures.py
```

The self-contained derived-table path reproduces E5-C5, E5-C6 and E5-C7 without redistributing the third-party stellar-track archive.

## Repository layout

```text
.
├── src/sclh/                 # scientific utilities
├── experiments/              # E5-C2, C3, C5, C6, C7
├── data/derived/             # frozen paper-supporting derived tables
├── tests/                    # focused tests + paper-result invariants
├── scripts/                  # verification, audit and figure reproduction
├── figures/                  # six publication figures
├── docs/                     # provenance and reproducibility guidance
├── CITATION.cff
├── .zenodo.json
├── pyproject.toml
└── LICENSE
```

## Reproduce the derived audits

The included E5-C2/E5-C3 history tables are sufficient for E5-C5, E5-C6 and E5-C7:

```bash
python scripts/reproduce_derived_audits.py
```

Outputs are written under:

```text
results/E5C5/
results/E5C6/
results/E5C7/
```

Regenerate all six paper figures with:

```bash
python scripts/make_paper_figures.py
```

## Full raw-track reproduction

Regenerating E5-C2 and E5-C3 from the upstream stellar tracks requires the external rotation/XUV grid associated with:

C. P. Johnstone, M. Bartel & M. Güdel, *The active lives of stars: a complete description of rotation and XUV evolution of F, G, K, and M dwarfs*, **A&A 649, A96 (2021)**.  
DOI: `10.1051/0004-6361/202038407`

Track archive: **Zenodo DOI `10.5281/zenodo.4266670`**

After extracting the archive, set:

```bash
export SCLH_MORS_ROOT=/path/to/TrackGrid_MstarPercentile
```

or on Windows PowerShell:

```powershell
$env:SCLH_MORS_ROOT = "C:\path\to\TrackGrid_MstarPercentile"
```

Then run:

```bash
python experiments/e5c2_toi700d_mors_history_inversion.py
python experiments/e5c3_lhs1140b_crosssystem_history.py
```

The original third-party track archive is **not redistributed** here.

## Key machine-readable outputs

```text
data/derived/paper_results.json
data/derived/toi700_rotation_conditioned_histories.csv
data/derived/lhs1140_mass_grid_rotation_support.csv
data/derived/lhs1140_rotation_compatible_histories.csv
data/derived/lhs1140_activity_model_comparison.csv
data/derived/lhs1140_weight_sensitivity.csv
data/derived/observable_information_audit.csv
data/derived/age_information_leverage.csv
data/derived/age_information_thresholds.csv
data/derived/age_information_summary.json
```

## Reproducibility layers

1. **Self-contained paper-level reproduction.** Reproduces E5-C5/E5-C6/E5-C7, headline invariants and all six figures from the included derived histories.
2. **Upstream raw-track regeneration.** Requires the external Johnstone/MORS archive and regenerates E5-C2/E5-C3.

This separation preserves provenance without silently redistributing third-party scientific assets.

## Provenance

The original E5-C2/E5-C3 numerical baseline derives from the frozen internal snapshot:

```text
SCLH_v0.1.18_E5C6_dual_manuscript_snapshot.zip
```

SHA-256:

```text
75a40de65b56642ed934a2ee21923c17c46db0178084f4c27cbd5de3b6885e80
```

Release v1.0.2 does not alter the frozen TOI-700 or LHS 1140 history tables. It adds the E5-C7 diagnostic derived from those tables, updates paper-facing interpretation, and carries forward the v1.0.1 author-affiliation correction to **Independent Researcher, Beijing, China**.

See `docs/PROVENANCE.md`, `docs/REPRODUCIBILITY.md` and `docs/AGE_INFORMATION_AUDIT.md`.

## Citation

The archival record associated with this software release is:

**DOI: `10.5281/zenodo.21974562`**

See `CITATION.cff` for machine-readable citation metadata.

## License

Repository-authored code is released under the **MIT License**. Third-party Johnstone/MORS data and other external scientific assets remain governed by their original terms and are not included in this repository.
