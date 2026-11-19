# Publishing a versioned research artifact

1. Extract the repository folder. Start with a private GitHub repository, review every file and confirm the software license, author identifiers and citation metadata. Add them only after approval. No license has been selected in this candidate.
2. Run the quick start and full reproduction. If you intentionally add metadata, regenerate SHA256SUMS with `python scripts/update_manifest.py --confirm-reviewed-change` and run verification again.
3. Create an empty repository on GitHub (do not initialize it with another README). From the extracted folder, run:

```sh
git init -b main
git add .
git commit -m "Scientific software release candidate"
git remote add origin https://github.com/YOUR-ACCOUNT/YOUR-REPOSITORY.git
git push -u origin main
```

Replace the example URL with the actual repository address. Authentication is handled by GitHub; do not commit credentials. Check that the Actions workflow actually succeeds. Inspect its uploaded evidence and compare its JSON and data files to the reference outputs. If a platform-specific stdout difference occurs, diagnose it; do not overwrite the reference before review.
4. After author approval and the final manuscript audit, tag a fixed version and create a GitHub release. Archive the identical release using a trusted archival service such as Zenodo and use its issued version-specific DOI. Do not cite an uncreated DOI or URL. The article must link the actual release and its code/data availability statement.
5. The public package contains scientific scripts, exact data, figures and tests. Keep private correction registers, drafts, development history and author-facing audits outside it. Obtain any required permissions and comply with the target journal's disclosure requirements for research software and other assistance.
