# Protein variant and peptide tools

Validate protein substitutions, build mutant sequences, and generate traceable peptide pools from user-defined sequence windows. Export Excel reports with changed amino acids in red, machine-readable CSV tables, mutant FASTA records, and input checksums.

The repository contains only code and artificial examples. Supply your own reference sequence, chosen mutations, and window definitions.

## Start here: what problem does this solve?

Suppose you have a reference protein and a list of substitution combinations you want to study. You need to reconstruct the full protein sequence for each combination, then extract the short peptide sequences that cover the changed regions. Several variants may share the same peptide, so you also need to know which peptides are unique and which pools use them.

This project handles that bookkeeping in two stages. The mutant builder checks your mutation list against your reference and constructs each variant sequence. The pool generator then uses your explicitly defined peptide windows to extract the corresponding sequences from each variant. It produces both a detailed record of where every peptide came from and a deduplicated catalog of changed peptide sequences.

You choose the substitutions, window locations, window lengths, and groupings. The program applies those choices consistently and records the results. A “pool” here is a computational list associated with a variant; the software does not physically mix peptides or calculate reagent amounts.

### How this relates to EvoEF2

If you used EvoEF2 or another tool to select variants, bring the selected **mutation strings** into this workflow. You do not need to provide EvoEF2 energy tables, model structures, or install its executable. An energy-results spreadsheet is not a supported mutation-list input: put each chosen mutation combination on its own line in the variant file.

For example, two single mutations and their combined variant require three separate lines. A single line containing both mutations produces only the combined variant. The software does not expand that line into all subsets.

### What you need before starting

| File you supply | What you decide | Why it is needed |
| --- | --- | --- |
| Reference FASTA or plain-text sequence | The exact protein sequence and chain being analyzed | Establishes the wild-type residue at every position |
| Variant text file | The substitutions that belong together in each variant | Defines the complete sequence to construct for each record |
| Window CSV | The native peptide segments and how they are grouped | Defines which segments to extract and include in each pool |

The first two files are enough to build full-length mutant sequences. You need the window CSV only when generating peptide pools. All three files must describe the same reference and coordinate system.

## What the workflow does

| Step | Input | Result |
| --- | --- | --- |
| Build variants | One reference sequence and one mutation combination per line | Validated full-length mutant sequences |
| Generate pools | Builder CSV or Excel output, the same reference, and window definitions | Peptide-window rows and a summary for every variant |
| Deduplicate | Changed peptide rows | One catalog entry per distinct amino-acid sequence, linked back to its variants and windows |

This is a **sequence-processing workflow**. It accepts mutation strings using EvoEF2 notation. It does not run EvoEF2, model structures, calculate energies, discover mutations, predict epitopes, assess immunogenicity, or establish biological function. There is no automatic generation of all mutation combinations: each input line specifies one variant.

## Requirements and installation

Python 3.10 or newer and `openpyxl==3.1.5`. Neither R nor EvoEF2 is required to run the Python scripts. Internet access is needed to install missing dependencies; the workflow itself has no network calls.

Extract the repository, open a terminal in its top-level folder, and run:

```bash
python -m pip install -r requirements.txt
python check_environment.py
```

On Windows, use `py` in place of `python` if that is how Python is installed. An isolated virtual environment is recommended; run the installation and scripts with the same interpreter.

The checker runs both scripts on tiny artificial data in a temporary directory. It verifies exact output sequences, Excel handoff, mutation highlighting, pool summaries, and global peptide deduplication. It leaves user data untouched and removes its temporary files. A successful check tests this small path in your environment, not the correctness of your own scientific inputs.

Optional package installation gives the `variant-peptides` command:

```bash
python -m pip install -e .
variant-peptides --help
```

## Run the bundled example

Run the following commands from the repository root. Output directories must not exist yet; choose new directory names when repeating a run.

```bash
python mutant_builder.py --sequence examples/reference.fasta --variants examples/variants.txt --chain A --output results/example_mutants
```

```bash
python peptide_pool_generator.py --sequence examples/reference.fasta --mutants results/example_mutants/mutant_sequences.xlsx --windows examples/windows.csv --chain A --output results/example_pools
```

You may substitute `mutant_sequences.csv` for the Excel input. Both handoffs are checked against the reference and mutation strings.

Expected example results:

* Four variants, including one unchanged `WT` record.
* Four window rows under the default selection rule.
* Two distinct changed peptides in the synthesis catalog: `AADEF` and `FGHIA`.
* Two variants reported as `NO_AFFECTED_WINDOWS` in the pool summary.

The example reference is an alphabetical arrangement of the 20 standard amino acids. It is not a biological protein or a proposed experimental sequence.

## A worked example: from mutations to peptides

The bundled artificial reference is:

```text
ACDEFGHIKLMNPQRSTVWY
```

Positions are counted from 1. The second residue is C, the ninth is K, and the twentieth is Y. The example variant file defines:

| Variant ID | Mutation string | Meaning |
| --- | --- | --- |
| `VAR001` | `CA2A` | Change C to A at sequence position 2 on chain A |
| `VAR002` | `CA2A,KA9A` | Apply that change plus K to A at position 9 in the same protein |
| `VAR003` | `YA20F` | Change Y to F at position 20 |
| `VAR004` | `WT` | Retain the reference sequence |

The builder constructs the whole sequence, not just the changed segment. For example, `VAR002` becomes `AADEFGHIALMNPQRSTVWY`. The changed residues appear in red in the Excel sequence column. FASTA and CSV contain the same sequence without color formatting.

The example windows are:

| Group | Window | Reference positions | Native peptide |
| --- | --- | --- | --- |
| Group1 | 1 | 1–5 | `ACDEF` |
| Group1 | 2 | 5–9 | `FGHIK` |
| Group2 | 1 | 10–14 | `LMNPQ` |

Under the default affected-core rule, a change in either Group1 window selects **both** Group1 windows. The pool table therefore contains:

| Variant | Window | Native peptide | Extracted peptide | Status | Catalog ID |
| --- | --- | --- | --- | --- | --- |
| VAR001 | Group1:1 | `ACDEF` | `AADEF` | CHANGED | PEP0001 |
| VAR001 | Group1:2 | `FGHIK` | `FGHIK` | UNCHANGED | Blank |
| VAR002 | Group1:1 | `ACDEF` | `AADEF` | CHANGED | PEP0001 |
| VAR002 | Group1:2 | `FGHIK` | `FGHIA` | CHANGED | PEP0002 |

Group2 is not selected because neither variant changes positions 10–14. `VAR003` changes a position outside every supplied window, and `VAR004` has no changes. Both still appear in the summary, with `NO_AFFECTED_WINDOWS`.

There are four window rows but only three changed rows and **two distinct changed sequences**. `AADEF` is needed by both VAR001 and VAR002, so it appears once in the synthesis catalog and is linked to both variants. The unchanged `FGHIK` row remains visible for context but is not automatically included in that catalog.

This distinction is central to interpreting the output: the pool table describes selected windows; the synthesis catalog lists unique changed sequences. It is not a complete experimental ordering sheet, because native controls and peptide chemistry must be specified separately.

## Use your own data

Create an `inputs/` directory and keep your input files there. Run:

```bash
python mutant_builder.py --sequence inputs/reference.fasta --variants inputs/variants.txt --chain A --output results/my_mutants
```

```bash
python peptide_pool_generator.py --sequence inputs/reference.fasta --mutants results/my_mutants/mutant_sequences.csv --windows inputs/windows.csv --chain A --output results/my_pools
```

Paths are relative to your terminal's current directory, not the script location. Quote any path that contains spaces. The equivalent combined entry points are `python variant_peptide_tools.py mutants ...` and `python variant_peptide_tools.py pools ...`.

### Prepare your first real run

1. Put your exact reference sequence in `inputs/reference.fasta`.
2. Put the selected substitution combinations in `inputs/variants.txt`, one combination per line.
3. Run the mutant builder and open `mutant_sequences.xlsx`. Confirm that the mutation strings and full-length sequences match your intended variants.
4. Define the peptide segments in `inputs/windows.csv`. The software does not automatically tile the protein or choose window lengths.
5. Run the pool generator with the same reference, chain, and numbering offset used for the builder.
6. Open `peptide_pools.xlsx`. Start with **Pool summary**, inspect **Peptide pools**, then use **Synthesis catalog** to identify shared changed sequences.
7. Retain the inputs and both output directories together so each result can be traced to the files and settings that produced it.

### Reference sequence

Supply exactly one FASTA record or a plain-text protein sequence. Whitespace and lowercase letters are normalized. Only the 20 standard amino acids are accepted. Multiple FASTA records, empty sequences, stop symbols, ambiguous letters, and nonstandard residues are rejected. The sequence length is limited to 32,767 characters for Excel compatibility.

### Variant list

A plain-text file contains one variant per nonblank line:

```text
# Artificial examples matching examples/reference.fasta
CA2A;
CA2A,KA9A;
WT;
```

`CA2A` means change C at residue number 2 on chain A to A. Commas separate substitutions within a variant; a terminal semicolon is accepted but optional. Blank lines and comments starting with `#` are ignored. `WT` explicitly represents the unchanged reference.

Chain IDs are case-sensitive single letters or digits. Every mutation must match the selected `--chain`; this workflow processes one chain per run. The tool rejects wrong reference residues, out-of-range positions, duplicate positions within a variant, malformed tokens, and substitutions that do not change an amino acid. Repeated variant lines remain separate records with separate IDs. IDs are assigned in input-record order as `VAR001`, `VAR002`, and so on.

There is no six-mutation limit. If needed, set a deliberate cap:

```bash
python mutant_builder.py --sequence inputs/reference.fasta --variants inputs/variants.txt --max-mutations 6 --output results/capped_mutants
```

### Residue numbering

By default, token residue number 1 means the first amino acid of the supplied sequence. The program does not read PDB files or infer structure-to-sequence numbering.

For a contiguous numbering shift, use `--position-offset` on **both** commands:

```text
sequence_position = mutation_number - position_offset
```

For example, with `--position-offset 100`, token number 102 refers to sequence position 2. Windows always use one-based sequence coordinates, regardless of this option. Numbering gaps, insertion codes, indels, multi-chain variants, and general residue maps are unsupported. Map such mutation coordinates to your reference before using this program; a constant offset is not sufficient for arbitrary PDB numbering.

### Window definitions

Supply comma-separated CSV with these columns:

| Column | Required | Meaning |
| --- | --- | --- |
| `core_id` | Yes | User-defined group label, nonempty and not `NONE` |
| `window` | Yes | Positive window identifier, unique within its group |
| `window_len` | Yes | Number of amino acids in `native_window` |
| `native_window` | Yes | Unmodified peptide exactly matching the reference |
| `start` | No | One-based sequence start; use it to disambiguate repeated sequences |
| `core_seq` | No | Optional annotation retained in window reports; not used to select windows |

```csv
core_id,core_seq,window,window_len,native_window,start
Group1,CDE,1,5,ACDEF,1
Group1,CDE,2,5,FGHIK,5
Group2,MNP,1,5,LMNPQ,10
```

When `start` is omitted or blank, the native window must occur exactly once in the reference. Overlapping repeated occurrences count as repeats. With an explicit start, the window must match the reference at that exact position. End coordinates are inclusive.

`core_seq` is an annotation for the group, not a requirement that each window contain the entire core. If supplied, its amino-acid characters must be valid and nonblank annotations for a group must agree. Biological relevance of groups and windows is determined by the user.

Every window must validate. Invalid rows stop the run rather than silently producing incomplete pools.

## Pool-selection rules

Default: `--selection affected-core`.

If a mutation lies inside any window of a group, **all windows in that group** appear in that variant's pool table. This preserves group-level selection even when a neighboring window remains unchanged. Each row is marked `CHANGED` or `UNCHANGED`.

For only windows that directly contain a mutation, use:

```bash
python peptide_pool_generator.py --sequence inputs/reference.fasta --mutants results/my_mutants/mutant_sequences.csv --windows inputs/windows.csv --selection changed-windows --output results/changed_windows_only
```

A variant with no mutations in any defined window still appears in the pool summary as `NO_AFFECTED_WINDOWS`; it has no peptide-window rows. `WT` consequently does not create a pool automatically. Native windows are exported separately so users can plan reference controls explicitly.

### What “core,” “window,” and “pool” mean in this program

A **core** is a group label supplied in `core_id`. It lets you treat several windows as a unit. The software does not infer a binding core or predict a biological epitope from that label.

A **window** is one exact contiguous segment of the reference, with a start, length, and native sequence. Windows may overlap or differ in length. Their coordinates remain fixed when substitutions are applied, because substitutions do not change sequence length.

A **pool** is the set of selected window records for one variant. The default rule selects entire affected groups; `changed-windows` selects only windows that directly contain changes. In both modes, only changed sequences enter the deduplicated catalog.

Choose `affected-core` when your analysis is organized around groups of related windows. Choose `changed-windows` when you want the report restricted to segments directly changed by a substitution. The choice changes which window rows are retained, not the full-length mutant sequences.

## What each output means

### Variant output directory

* `mutant_sequences.xlsx`: validated variants, red changed residues, and run metadata.
* `mutant_sequences.csv`: the same records as plain data.
* `mutant_sequences.fasta`: full-length sequences with variant IDs and mutation strings.
* `manifest.json`: input filenames and SHA-256 checksums, software versions, and numbering settings.

### Pool output directory

* `peptide_pools.xlsx`: peptide rows, pool summary, native windows, synthesis catalog, and run metadata.
* `peptide_pools.csv`: one row per selected variant/window combination, including start/end, mutation profile, status, and catalog peptide ID.
* `pool_summary.csv`: one row for every input variant, including variants without affected windows.
* `native_windows.csv`: all supplied windows, their coordinates, and annotations.
* `synthesis_catalog.csv`: one row per distinct **changed** peptide sequence across all variants and windows.
* `manifest.json`: input checksums and run settings.

The catalog deduplicates globally by amino-acid sequence. It retains variant and window provenance; the same catalog ID can occur in many pool rows. Unchanged rows have no catalog ID. The catalog does not automatically include native controls, modification chemistry, quantities, purity requirements, or aliquoting instructions. Deduplication assumes chemically identical unmodified sequences; different terminal modifications or labels require separate ordering decisions.

Counts distinguish selected window rows, changed window rows, and unique changed sequences. A pool can have multiple window rows pointing to one physical peptide. Catalog IDs are deterministic for a given set of sequences but can change if that set changes. Variant IDs likewise depend on input order; retain each run's inputs and outputs together.

Excel cells containing user labels are stored as text. CSV files preserve raw labels; import them as text if using a spreadsheet application that automatically interprets formulas.

### Reading the Excel workbook

| Sheet | What to look for |
| --- | --- |
| Mutant sequences, in the builder workbook | Full-length sequence for each input variant; red letters mark substitutions |
| Peptide pools | Native and mutant peptide side by side, sequence coordinates, mutations within the window, and catalog ID |
| Pool summary | Every variant, affected groups, selected window count, changed window count, and distinct changed sequence count |
| Native windows | All original window definitions, including windows that no variant affected |
| Synthesis catalog | One row per distinct changed sequence, its shared peptide ID, and the variants/windows using it |
| Run metadata | Software versions, input checksums, chain, offset, and selection settings |

In the pool sheet, red letters identify the positions that differ between the native and mutant peptides. Red formatting is a visual aid; it does not indicate an energy score, a predicted effect, or a failed validation.

The `mutation_profile` column lists only mutations inside that particular window. The `mutation_string` column retains the complete variant definition. These columns can differ when a variant contains substitutions in distant parts of the protein.

### Common result patterns

| Result | Interpretation | Next step |
| --- | --- | --- |
| CHANGED | At least one residue in the window differs from its native sequence | Follow its peptide ID to the catalog |
| UNCHANGED | This window was selected through its group, but its own sequence did not change | Decide whether it is relevant as context or a separately planned control |
| NO_AFFECTED_WINDOWS | None of this variant's changes falls inside any supplied window | Check whether your window definitions cover the intended positions |
| Same peptide ID on several rows | The exact changed sequence is shared | Preserve each row's pool association; the catalog lists the sequence once |
| Empty catalog with a completed run | No selected peptide has a sequence change | Inspect the summary and coverage; completion alone does not imply a useful pool |

## Command options at a glance

| Option | Applies to | Meaning |
| --- | --- | --- |
| `--sequence` | Both scripts | Required reference FASTA or plain-text sequence path |
| `--variants` | Builder | Required text file of chosen mutation combinations |
| `--mutants` | Pool generator | Required builder CSV or compatible Excel workbook |
| `--windows` | Pool generator | Required window-definition CSV |
| `--output` | Both scripts | Required new directory for that stage's results |
| `--chain` | Both scripts | Expected token chain, default `A`; case-sensitive |
| `--position-offset` | Both scripts | Constant subtracted from token residue numbers, default `0` |
| `--max-mutations` | Builder | Optional positive cap per variant; no cap by default |
| `--selection` | Pool generator | `affected-core` by default, or `changed-windows` |

For the help text built into each script:

```bash
python mutant_builder.py --help
python peptide_pool_generator.py --help
```

## Troubleshooting

| Message or symptom | Likely explanation | What to do |
| --- | --- | --- |
| `python` is not recognized | Python is unavailable under that command | Try `py` on Windows, or select the installed interpreter |
| `No module named openpyxl` | Dependencies were installed into a different environment or not installed | Run `python -m pip install -r requirements.txt` with the interpreter used for the scripts |
| File not found | The relative path is being resolved from a different current folder | Run from the repository root or supply a quoted absolute path |
| Output directory already exists | The program is protecting an earlier run | Use a fresh output directory |
| Mutation chain differs from selected chain | Mutation tokens and `--chain` disagree | Select the intended chain consistently in both stages |
| Reference residue does not match | The mutation number, reference sequence, or offset is wrong | Reconcile the numbering with the actual reference; do not simply substitute a different wild-type letter |
| Native window occurs more than once | The same peptide has multiple possible locations, including overlapping repeats | Provide its intended `start` coordinate |
| Explicit window coordinates do not match | The native sequence or length does not match that segment | Check the one-based start and exact reference substring |
| Mutant sequence does not match substitutions | The builder table was altered or belongs to another reference/settings | Rebuild from the intended input files and use its unedited output |
| No affected windows | Mutations fall outside the definitions, or the variant is WT | Review window coverage and the numbering convention |

The quick environment check does not need your data. If it passes but your run fails, start with the reported input error and the corresponding file rather than reinstalling everything.

## Validation and error handling

All variants are checked before output creation. Any invalid variant stops the run. The pool generator independently rebuilds each mutant from the reference and mutation string and compares it against the supplied sequence; an `OK` flag alone is insufficient.

The pool generator also accepts the earlier builder's five-column Excel layout with headers at row 4 and data from row 5. Invalid rows and altered sequences are rejected rather than skipped. CSV generated by this repository is the simplest portable handoff. A different workbook layout requires conversion, not a guessed column match.

Exit status `0` means the operation completed, including cases where no windows were affected. Status `2` indicates invalid input, command-line usage, or a handled file/dependency error. Unexpected failures may return another nonzero status. If a disk or file-writing failure leaves a partial directory, inspect it and choose a fresh output directory for the next run.

Existing output directories are never overwritten. Keep the input files: checksums identify the exact bytes used but do not replace backups or establish biological validity.

## Tests and reproducibility

```bash
python check_environment.py
python -m unittest discover -s tests -v
```

Tests use only artificial sequences and temporary files. They check substitutions against independently constructed expected sequences, malformed inputs, chain and numbering rules, repeated windows, pool selection, global deduplication, Excel/CSV round trips, and overwrite protection.

The synthetic end-to-end checker also passed on a user-reported Windows environment with Python 3.11.9 and openpyxl 3.1.5. This confirms that smoke-test path, not an independent full Windows test-suite run.

The included GitHub Actions workflow is configured for Windows and Ubuntu with Python 3.10 and 3.12. Those hosted runs happen after upload; see `VALIDATION.md` for checks actually executed while preparing this release.

## Repository contents and sharing

The `.gitignore` excludes local input and output directories, common sequence/data files, workbooks, environments, and caches. Only the named synthetic example data are allowlisted. Ignore rules do not remove files already tracked by Git. Review `git status` and staged changes before committing your own work.

The original uploaded files and private research data are not part of this repository. Generic program logic is necessarily included because it implements the workflow. No project-specific sequences, selected mutations, targets, results, personal paths, or vendor information are included.

No license has been selected for this code. Add your chosen license before offering it under open-source reuse terms.

## Reference

EvoEF2 mutation-file notation is documented in the [official EvoEF2 repository](https://github.com/tommyhuangthu/EvoEF2#4-buildmutant). EvoEF2 is a separate project and is not bundled or invoked here. Its notation identifies the wild-type residue, chain, residue number, and substituted residue; each mutation combination can end in a semicolon.
