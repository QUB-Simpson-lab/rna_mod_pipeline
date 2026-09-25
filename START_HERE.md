# Test the pipeline on your laptop

Use this source folder for the Mac and Windows trial. Python 3.11 is
recommended; 3.10–3.12 are supported. An internet connection is needed for
the first installation. You do not need the large research databases for
the synthetic demonstration.

## 1. Put the code in its own workspace

Extract the supplied source ZIP so the layout is:

```text
RNA mod test/
└── rna_mod_pipeline/
    ├── START_HERE.md
    ├── install_macos.command
    ├── install_windows.bat
    ├── scripts/
    └── src/
```

Keep all files together, including hidden `.github`, `.gitignore` and
`.gitattributes` files. Do not copy a `.venv` from another computer.

On Windows, extract the code to a local drive, for example
`C:\RNA mod test\rna_mod_pipeline`. The helpers do not run directly from a
UNC network location such as `\\server\share`; they stop if they cannot open
their own source folder. Copy the source locally before installation.

## 2. Install, test and open the example

| Action | macOS | Windows |
|---|---|---|
| Install dependencies | Double-click `install_macos.command` | Double-click `install_windows.bat` |
| Run the automated tests and synthetic example | Double-click `test_demo_macos.command` | Double-click `test_demo_windows.bat` |
| Inspect the synthetic workspace in the GUI | Double-click `run_demo_macos.command` | Double-click `run_demo_windows.bat` |
| Open your normal research workspace | Double-click `run_macos.command` | Double-click `run_windows.bat` |

The example creates a sibling `rna_mod_demo/` directory with artificial
reference, annotation, binding tracks and one synthetic m6A dataset. It checks
analytically known overlap counts and odds ratios. These are software test
data, not real genes, RBP interactions or biological findings. No download of
hg38, GENCODE, ENCORI or POSTAR3 is needed for this example.

The test helper refuses to replace an existing demonstration directory.
For a repeat, keep/rename the previous demo or use the documented command
with a different `--workspace` path. Your research workspace is separate.

If macOS says a helper is not executable, open Terminal in the code folder
and run `chmod +x *.command` once. The installer checks standard Homebrew and
python.org locations as well as the command path for supported Python.

## 3. Record the laptop checks

Follow [CROSS_PLATFORM_ACCEPTANCE.md](docs/CROSS_PLATFORM_ACCEPTANCE.md).
Check that figures and tables open, saved datasets reload, a disposable run
can be cancelled, and results stay outside the source folder. The automated
example does not record a manual GUI check on your behalf.

Keep the test reports from both laptops. Statistical table values should
agree; PNG file hashes can differ because fonts and rendering differ.
Exact example commands and the expected-answer definitions are in
[examples/acceptance/README.md](examples/acceptance/README.md).

## 4. Use real data

Launch the normal GUI, choose **Resources…**, then add your own supported
bedMethyl as a dataset. Follow [RESOURCE_SETUP.md](docs/RESOURCE_SETUP.md) for
the external resources. The exact existing research snapshot can be checked
using [RESOURCE_SNAPSHOT.md](docs/RESOURCE_SNAPSHOT.md).

The synthetic profile must never be reused for an experimental bedMethyl
file. Select the matching real reference and annotation for real analyses.

## 5. Upload after testing

Follow [GITHUB_RELEASE.md](docs/GITHUB_RELEASE.md). The complete source folder
belongs at the root of the lab repository. Keep your virtual environment,
research data and generated results outside the upload. Standalone `.app`
and `.exe` bundles are a later step.
