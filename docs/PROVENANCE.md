# Provenance and interpretation governance

## Scientific baseline

This public repository is a deliberately reduced extraction of the internal SCLH v0.1.18 scientific snapshot.

Upstream snapshot SHA-256:

`75a40de65b56642ed934a2ee21923c17c46db0178084f4c27cbd5de3b6885e80`

Only the B1 stellar-history inverse-problem chain is retained.

## Included numerical chain

- **E5-C2:** TOI-700 d rotation-conditioned history inversion.
- **E5-C3:** LHS 1140 b cross-system replication.
- **E5-C5:** independent fully-convective present-day rotation/activity cross-calibration.
- **E5-C6:** multi-observable history-information and observation-operator audit.

## Interpretation-only public-release patch

The historical E5-C3 implementation included a direct zero-discrepancy interval comparison between native MORS present \(L_X\) and the published LHS 1140 XMM table value.

A later audit established that these quantities do not use identical X-ray observation operators (approximately 0.1–2.4 keV for the native model quantity versus 0.2–2.4 keV for the published XMM table value). Therefore the public release retains the numerical mismatch but does **not** describe it as a bandpass-homogenized physical rejection.

This patch changes comments/summary interpretation only. It does not change the E5-C3 rotation-compatible history table or the paper's cumulative-history values.

For the same reason, the old internal E5-C3 summary JSON is not redistributed in this minimal public repository. The authoritative public headline quantities are consolidated in `data/derived/paper_results.json`.

## Third-party assets

The original Johnstone/Bartel/Güdel rotation-XUV tracks are not redistributed.

Public source:
- Zenodo DOI `10.5281/zenodo.4266670`
- associated paper: Johnstone, Bartel & Güdel (2021), A&A 649, A96.

The included `data/derived/*.csv` files are derived analysis products needed to audit the published numerical claims.

## v1.0.2 derived-only extension

E5-C7 is a derived-only value-of-information audit. It consumes the frozen E5-C2/E5-C3 history tables and does not alter or regenerate those histories. It quantifies the worst residual model-conditional cumulative-forcing width that can remain after a hypothetical independent age interval of declared total width is applied.

The E5-C7 outputs are therefore downstream diagnostics with explicit provenance to the same frozen history tables. They do not introduce a new stellar-evolution track family or a new population prior.
