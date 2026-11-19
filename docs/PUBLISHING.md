# Publishing a versioned research artifact

The initial scientific snapshot is hosted at https://github.com/jusahil/ccz-conversion. The author-approved software license is MIT and the citation author is Sahil Imtiyaz. The initial public version is planned as v1.0.0.

1. Before public release, confirm the repository is cleared for public distribution, including institutional permissions where applicable. Review `LICENSE`, `CITATION.cff` and the scientific files.
2. On the final commit, confirm the GitHub Actions workflow succeeds. Inspect uploaded evidence and compare its JSON and data files to the reference outputs. If a platform-specific stdout difference occurs, diagnose it; do not overwrite reference results without review. For intentional file changes, use `python scripts/update_manifest.py --confirm-reviewed-change` before running the verification again.
3. When the author approves public release, change repository visibility to public. Sign in to Zenodo, connect the GitHub account, sync repositories and enable archival for this repository **before** publishing a GitHub release.
4. From GitHub Releases, create tag `v1.0.0` targeting the verified commit and publish the release with the title "Exact CCZ Conversion under Conserved Binary Charges: Research Software v1.0.0". Verify Zenodo archives the tagged snapshot and check its author, title, version and license metadata.
5. Cite the issued version-specific DOI and repository address in the article. Do not cite an uncreated DOI. If adding a DOI badge or metadata to GitHub later, regenerate the manifest and keep the archived tag unchanged.
6. The public package contains scientific scripts, exact data, figures and tests. Keep private correction registers, drafts, development history and author-facing audits outside it. Obtain any required permissions and comply with the target journal's disclosure requirements for research software and other assistance.
