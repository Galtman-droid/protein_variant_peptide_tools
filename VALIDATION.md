# Validation record

Release: 1.0.0. Checked on 2026-10-05.

## Executed checks

* Python 3.12.14, Linux, openpyxl 3.1.5.
* Quick environment checker: passed both command-line scripts using temporary artificial data, the Excel handoff, exact expected sequences and catalog, pool summary, and red amino-acid highlighting.
* 23 unittest cases passed, including all 380 possible single substitutions of the 20-residue artificial reference to the other 19 standard amino acids; numbering and chain validation; malformed inputs; overlapping repeated windows; pool selection; deduplication across groups; and CSV/Excel handoff equivalence.
* An additional compatibility check used the original builder to generate an Excel workbook from artificial sequences; the new pool reader accepted and independently validated both variants.
* A Python wheel built successfully with the supplied packaging metadata.
* The final ZIP was assembled from an explicit file allowlist. It contains no original uploaded scripts, real protein structures, research datasets, generated workbooks, private paths, or Git history. A text scan checked for known source-specific names, example sequences, mutation tokens, and personal path patterns.

## Documentation revision

README.md was expanded with a worked example, a first-run walkthrough, workbook interpretation, command options, and troubleshooting. Python code and dependency requirements are unchanged. The example tables were checked against generated CSV outputs.

The user supplied a successful quick-check transcript from Windows with Python 3.11.9 and openpyxl 3.1.5. This is a user-reported smoke-test result, not a full test-suite execution on Windows by the release preparer.

## Not executed here

* Windows and Python 3.10 runtime checks. GitHub Actions is configured to perform these after upload.
* EvoEF2 energy or structure calculations: those are outside this program's scope.
* Biological validation of any user input, peptide, mutation, or experimental design.

## Deliberate changes from the earlier scripts

Input files are explicit command-line arguments. A shared validator checks chain identity, wild-type residues, unique mutation positions, optional constant numbering offset, a single FASTA record, and full sequence agreement. Terminal semicolons are accepted, and the hardcoded six-mutation cap is replaced by an optional cap. Invalid input now stops processing instead of silently dropping rows. Ambiguous windows require explicit coordinates. All output directories must be new.

The default affected-core rule is preserved. An optional changed-windows mode narrows selection. The synthesis catalog deduplicates changed peptides globally by sequence, while window provenance stays in the pool table. Native reference windows remain a separate output; controls are not automatically added to a pool.

Input and output formats, limitations, error handling, and licensing status are documented in README.md.

The v3 distribution uses standard Markdown documentation only. The complete Python program, tests, configuration, and artificial examples are retained.
