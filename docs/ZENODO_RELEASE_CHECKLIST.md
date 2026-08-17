# Zenodo release checklist — v1.0.2

1. Confirm the repository root is `stellar-history-identifiability`.
2. Run `pytest`.
3. Run `python scripts/reproduce_derived_audits.py`.
4. Run `python scripts/make_paper_figures.py`.
5. Run `python scripts/verify_release.py` after refreshing `SHA256SUMS.txt`.
6. Confirm `CITATION.cff`, `.zenodo.json`, `pyproject.toml`, `README.md` and `paper_results.json` all report version `1.0.2`.
7. Confirm creator metadata:
   - Yuzhan Zhang
   - ORCID 0009-0000-3121-7972
   - Independent Researcher, Beijing, China
8. Create GitHub tag/release `v1.0.2`.
9. Update the existing Zenodo archival record as intended by the author, retaining DOI `10.5281/zenodo.21974562` if Zenodo permits that edit path.
10. Verify the public GitHub release and Zenodo record before final journal submission.
