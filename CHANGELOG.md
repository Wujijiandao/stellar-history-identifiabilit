# Changelog

## v1.0.2 — 2026-08-17

Scientific/reproducibility update for the Nature Astronomy submission.

- Added **E5-C7 age-information leverage audit** using the frozen E5-C2/E5-C3 derived histories.
- Added a worst-case finite-support diagnostic for the residual historical width after an independent stellar-age interval is imposed.
- Added `age_information_leverage.csv`, `age_information_thresholds.csv` and `age_information_summary.json`.
- Added publication Figure 6 and its regeneration path.
- Added focused regression tests for E5-C7 and updated paper-result invariants.
- Updated the public paper summary to formalize model-conditional identified-set widths.
- Carried forward the v1.0.1 author metadata correction: **Independent Researcher, Beijing, China**.
- Updated repository/Zenodo metadata to v1.0.2 and DOI `10.5281/zenodo.21974562`.
- No E5-C2/E5-C3 frozen history values were changed.

## v1.0.1 — 2026-08-17

Metadata-only correction.

- Corrected the author affiliation from `Independent Researcher, Japan` to `Independent Researcher, Beijing, China`.
- No scientific code, derived history tables or paper results were changed.

## v1.0.0 — 2026-08-17

First public paper-specific release.

- Extracted the E5-C2/C3/C5/C6 stellar-history identifiability chain from the larger internal SCLH project.
- Added frozen derived history tables and paper-result invariants.
- Added publication-figure regeneration and focused tests.
- Added GitHub/Zenodo metadata and provenance documentation.
- Interpretation-only cleanup: native LHS 1140 X-ray interval mismatch is described as an activity/observation-mapping discrepancy alarm pending bandpass homogenization, not as a calibrated physical rejection.
- No numerical E5-C2/E5-C3 history values were changed.
