# GitHub → Zenodo release checklist

1. Create the public GitHub repository, e.g. `stellar-history-identifiability`.
2. Upload this package as the repository root.
3. Review `README.md`, `CITATION.cff`, `.zenodo.json`, and `LICENSE`.
4. Run `pytest` and `python scripts/verify_release.py`.
5. Commit and push.
6. Connect/enable the repository in Zenodo.
7. Create GitHub release/tag `v1.0.0`.
8. Confirm the Zenodo record metadata:
   - Creator: Zhang, Yuzhan
   - ORCID: 0009-0000-3121-7972
   - Affiliation: Independent Researcher, Japan
   - Type: Software
   - License: MIT
9. Obtain the Zenodo DOI.
10. Add the DOI badge and DOI to the GitHub README/CITATION metadata.
11. Replace the manuscript Code availability placeholder with the GitHub URL + archived Zenodo DOI.
12. Do not upload third-party MORS raw tracks unless their redistribution terms are independently verified.
