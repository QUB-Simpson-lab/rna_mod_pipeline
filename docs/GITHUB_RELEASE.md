# Publish the source after the laptop trial

Target: <https://github.com/QUB-Simpson-lab/rna_mod_pipeline>.

The supported distribution is the full source checkout, installed with the
supplied helpers. A standalone app is not needed for this source release.

1. Finish the native Mac and Windows checks in `CROSS_PLATFORM_ACCEPTANCE.md`.
   Record the tested source checksum, Python version and any failures. Review
   the final contributor list in `CITATION.cff` and the MIT licence with the
   lab; MIT was selected for this preparation by Alexandru on 8 September 2026.
2. Clone the lab repository into a separate local folder using Git or GitHub
   Desktop. Preserve any existing lab files and history. The large research
   workspace currently has a different remote and is not the upload source.
3. Copy the **contents** of the tested source folder into that checkout.
   Include hidden `.github/`, `.gitignore` and `.gitattributes`. The root
   should contain `README.md`, `pyproject.toml`, `scripts/`, `src/`, `config/`,
   `tests/`, `examples/` and the install/run helpers. Do not nest everything
   under another `refactored_code/` folder.
4. Review the changes before committing. Do not add `.venv`, `rna_mod_demo`,
   experimental inputs, real results, caches, temporary builds or local
   resource/dataset profile JSON files. The synthetic fixture generator and
   resource checksum inventory belong in the source; the generated fixtures
   and large external resources do not.
5. Commit and push to the branch agreed with the lab. There is no need to
   force-push or rewrite the history. If the lab uses pull requests, follow
   that workflow; otherwise a normal main-branch push can use the same
   reviewed source tree.
6. Open the repository's **Actions** tab. Require the Windows, macOS, Linux
   and package checks to finish successfully. Merely including a workflow
   file does not mean these jobs have passed.
7. Record the final commit in the acceptance notes. Create a tagged release
   only after any failures are fixed and the tested version is the one being
   published. Keep future fixes as new commits/versions.

The prepared `LICENSE` uses MIT for software. External databases retain their
own terms, and no database redistribution permission is implied. The citation
currently contains Alexandru Zob's verified name; add further authors and
identifiers only when the lab agrees their attribution.

The old research repository can be backed up separately. Uploading this source
does not automatically back up its raw data or historical results.
