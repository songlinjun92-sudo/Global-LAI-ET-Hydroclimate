# Checklist for the first public repository release

1. Replace generic contributor metadata and confirm all authors permit publication under the chosen code license.
2. Run `python -m pytest -q` and `python -m lai_et demo --root examples/demo_run` from a freshly created environment.
3. Cross-check the output of a **real data run** with the archived manuscript tables. Where numerics differ, investigate exact product versions, training cap, PFT domain and ICE definition; never silently change numbers.
4. Remove raw third-party data, personal paths, authentication tokens, and all confidential/project-only files; inspect `git status` and GitHub's public preview before publication.
5. Create the public GitHub repo and include a stable tagged release (`v1.0.0`). Consider Zenodo archive + DOI and add citation metadata.
6. Complete repository URL in the manuscript code-availability paragraph **after** upload; do not claim public availability early.
7. Cite all underlying datasets and official PML method/data papers as required by their maintainers. GitHub code licensing does not grant redistribution rights to source data.
