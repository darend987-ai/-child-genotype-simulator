# Child Genotype Simulator

A Streamlit application that generates **200 unique**, Mendelian-consistent simulated STR child profiles for one selected mother/father pair, while excluding profiles that exactly match uploaded real children.

## Input format: instrument matrix (Investigator 24plex / GlobalFiler)

Upload the laboratory export directly. The app expects the same general matrix form as the supplied `Parent Child Test.xlsx`:

- **Row 1** contains the locus names, beginning in column B.
- **Column A** contains sample IDs.
- Each later column contains the reported genotype for one locus.
- Control, ladder, negative, and positive-control rows can remain in the file; simply do not select them as parents or known children.

The application does **not** hard-code an Investigator 24plex locus list. It takes every nonblank locus name from row 1, in that order. This supports GlobalFiler or another compatible kit with a different locus panel and allele values.

### Selecting people in the app

1. Upload a matrix export.
2. Choose the sample ID for the mother and father.
3. Select one or more real children whose complete genotype profiles must be excluded.
4. Set the case/family ID, target count (default 200), and optional random seed.
5. Generate and download the results.

## Output format

The Excel and CSV exports use the same matrix layout as the upload:

- blank first header cell;
- sample ID in the first column;
- the original locus names and their original order across the remaining columns;
- one row per generated child, named `<case_ID>_SIM_001` through `<case_ID>_SIM_200`.

The export contains the generated children only, not controls or original samples. This avoids accidentally treating a simulated result as an original laboratory record.

## Genetic rules and exclusions

- Autosomal loci use one transmitted allele from each selected parent.
- Allele order is normalized internally, so `15, 16` and `16, 15` are identical profiles.
- A known child blocks a simulated child only when the **entire multilocus profile** matches.
- `Amel`/`Amelogenin` accepts `X` as `X,X` internally and exports `X` in the usual instrument style.
- `QS1` and `QS2` are PCR quality-control markers, not inherited DNA loci. The simulator always writes their successful-run single calls (`QS1 = 1`, `QS2 = 2`), does not require usable parent calls in those columns, and excludes them from whole-profile identity comparisons.
- Loci beginning with `DYS` and `Yindel` are treated as Y-linked: they are blank for an XX child and use the paternal Y genotype for an XY child. Where `Amel` is included, this is tied to the child’s simulated Amel result.
- The app refuses to produce a partial output when fewer than the requested number of unique permissible profiles exist.

## Run locally

```bash
python -m venv .venv
# macOS/Linux
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1
pip install -r requirements.txt
streamlit run app.py
```

## Deployment

Deploy the repository to Streamlit Community Cloud with `app.py` as the entry point.

## Scope and validation

This application is a simulation tool. It does not model mutation, linkage, allele frequencies, population substructure, dropout, stutter, mixtures, sample contamination, probabilistic genotyping, or kinship likelihood ratios. Validate any laboratory or forensic workflow under your laboratory’s approved procedures and quality system before operational use.
