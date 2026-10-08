# SKINK TIME limb-gene training pipeline

**Gene list → NCBI Datasets RefSeq CDS across Lepidosauria → longest eligible
CDS per species → MACSE → IQ-TREE.** Lepidosauria includes squamates and tuatara.

Step 01 processes `genes.txt` in a batch. Steps 02–04 prompt for one gene at a
time so Nick can inspect sequence selection and alignment before building trees.
`species.txt` is an optional reference panel; it does not restrict retrieval.

## Files

- `genes.txt`: one gene symbol per line; SHH is a starter example.
- `scripts/01_download_cds.sh`: downloads one Datasets package per listed gene.
- `scripts/02_extract_longest_refseq.sh`: selects one eligible CDS per species.
- `scripts/03_macse.sh`: aligns CDS and saves raw/cleaned alignments.
- `scripts/04_build_tree.sh`: fits a DNA model and infers a gene tree.
- `scripts/ncbi_cds.py`: package handling, selection and provenance.
- `scripts/common.sh`, `scripts/pipeline.py`: prompts and validation.
- `environment.yml`: software specification.
- `metadata/gene_tracking.tsv`: one reviewed row per gene.
- `genes/.gitkeep`: placeholder; generated gene directories are ignored.
- `tests/test_pipeline.py`: offline checks.

## Create the environment in Monsoon scratch

From the repository root on Monsoon:

```bash
module avail mambaforge
module load mambaforge
source "$(conda info --base)/etc/profile.d/conda.sh"

export SKINK_SCRATCH="/scratch/$USER"  # Adjust if your assigned path differs.
export CONDA_ENVS_PATH="$SKINK_SCRATCH/conda/envs"
export CONDA_PKGS_DIRS="$SKINK_SCRATCH/conda/pkgs"
mkdir -p "$CONDA_ENVS_PATH" "$CONDA_PKGS_DIRS"

export CONDARC="$SKINK_SCRATCH/conda/condarc"
cat > "$CONDARC" <<CONFIG
envs_dirs:
  - $SKINK_SCRATCH/conda/envs
pkgs_dirs:
  - $SKINK_SCRATCH/conda/pkgs
channel_priority: strict
CONFIG

conda env create --prefix "$CONDA_ENVS_PATH/skink-time-limb-genes" --file environment.yml
conda activate "$CONDA_ENVS_PATH/skink-time-limb-genes"
conda config --show envs_dirs pkgs_dirs
```

NAU also documents an `anaconda3` module; inspect available modules if the name
changes. Both environments and package caches go in scratch; the centrally
installed base conda stays in its existing location. Repeat module loading,
`source`, exports and activation in each new session/job.

**Already created an earlier environment?** Activate it and update:

```bash
conda env update --prefix "$CONDA_ENVS_PATH/skink-time-limb-genes" --file environment.yml
```

NCBI Datasets and Biopython are required. Python handles ZIP extraction;
a separate unzip package is unnecessary. Version constraints are not a complete
transitive dependency lock. After a successful Monsoon installation, save an
exact Linux package snapshot:

```bash
conda list --explicit > metadata/conda-linux-64-explicit.txt
# Later recreate the resolved environment:
# conda create --prefix "$CONDA_ENVS_PATH/skink-time-limb-genes-locked" \
#   --file metadata/conda-linux-64-explicit.txt
```

Review the snapshot for private channel credentials before committing it.
Back up scripts/specifications on GitHub so scratch environments can be rebuilt.

## Download the gene list

Edit `genes.txt` with your chosen gene symbols, one per line. Blank lines and
`#` comments are allowed. Symbols normalize to uppercase; invalid or duplicate
symbols fail before any download. The initial SHH is an example, not an assumed
final panel of limb genes.

```bash
bash scripts/01_download_cds.sh
# Or use a different list:
# bash scripts/01_download_cds.sh /path/to/my_genes.txt
```

**No email prompt, E-utilities client or NCBI login is used.** For each gene,
the downloader runs the equivalent of:

```bash
datasets download gene symbol SHH --taxon Lepidosauria \
  --include cds,product-report \
  --filename genes/SHH/downloads/SHH_lepidosauria.zip
```

Saves the ZIP, extracted package and download log in `genes/GENE/downloads/`.
The package includes `ncbi_dataset/data/cds.fna`, `data_report.jsonl`, and product
metadata. `download_manifest.json` records the command, date, CDS count, metadata
count and file hashes; completion requires a valid ZIP, parseable report and at
least one CDS record. It keeps all downloaded isoforms; selection happens later.
There is no manually imposed species cap or result-count truncation.

`genes/download_batch_status.tsv` records each gene's result. A failure or
zero-CDS package is recorded while later genes are still attempted. The batch
exits with an error if any gene failed. Matching complete downloads with verified
file hashes are skipped on rerun. Move incomplete/old-scope gene directories aside
before retrying; automated partial-download resume is not provided.

“All available” means the CDS returned by Datasets for matching gene symbols
within Lepidosauria, including predicted RefSeq transcripts. This is not a BLAST
search for every homolog. Missing annotations, different names or species without
available RefSeq CDS can limit coverage. Check `download.log` when results are
unexpected. Symbol queries can match synonyms and do not establish orthology.

## Process one gene at a time

```bash
bash scripts/02_extract_longest_refseq.sh
# Enter SHH, then inspect genes/SHH/candidate_cds.tsv and selection.tsv.
```

Step 02 parses the Datasets CDS headers and joins them to the gene report by
GeneID. Uses metadata/header organism names, groups subspecies under the species
binomial, and retains original names/accessions. Requires an NM_/XM_ RefSeq
coding-transcript accession, matching gene symbol or documented report synonym,
consistent species metadata and a protein-coding gene record where supplied.
Unresolved species, invalid/uninformative sequences, lengths not divisible by
three and internal stops are excluded and documented. If multiple eligible
GeneIDs remain for one species, that species is excluded for manual review.
Ambiguous bases, missing terminal stops, non-ATG starts and synonym matches are
flagged for review. FASTA CDS does not expose all GenBank feature quality flags;
these checks do not guarantee completeness or correct annotation.

Among eligible isoforms, longest CDS wins; accession/header breaks ties. There
is no curated-NM_ preference over a longer predicted XM_ transcript. Writes:

- `candidate_cds.tsv`: every downloaded CDS, species, taxid, transcript accession,
  GeneID, symbol, length, original header, eligibility and notes.
- `selection.tsv`: selected candidate for each included species.
- `extraction_summary.json`: CDS and species counts.
- `GENE_longest_cds.fasta`: one selected sequence per species, with Species_name labels.

Inspect lengths, isoforms, exclusions and coverage before alignment. The CDS are
already extracted by Datasets; introns are not passed to MACSE. Longest selection
and mechanical checks do not verify orthology or the preferred biological isoform.
Metadata/schema changes cause missing or mismatched records to be excluded rather
than silently guessed.

Inside a Monsoon compute allocation, run:

```bash
bash scripts/03_macse.sh
# Inspect alignment and logs before proceeding.
bash scripts/04_build_tree.sh
```

Use your lab's current Slurm allocation instructions; do not run alignment/tree
inference on a login node. Step 04 uses `SLURM_CPUS_PER_TASK`, or one thread if
unset. For a job, `printf 'SHH\n' | bash scripts/03_macse.sh` supplies the symbol;
likewise for step 04.

MACSE requires at least two sequences and saves raw nucleotide/amino-acid
alignments and logs under `alignment/`. Export masks internal/terminal stop
codons and internal frameshift codons with NNN, and remaining frameshift
characters with a gap. Inspect raw !/* characters and logs. Cleaned alignments
must have equal lengths divisible by three; masking is not biological validation.

IQ-TREE requires at least four species. Fits a DNA model with ModelFinder,
1,000 ultrafast bootstrap and SH-aLRT replicates, seed 12345, and results under
`iqtree/`. Inspect `GENE.treefile`. Trees are unrooted; root on Sphenodon_punctatus
in a viewer if present and appropriate. Gene trees may differ from species trees;
outputs are not automatically approved for HyPhy or codon-model analyses.

Stages refuse to overwrite outputs. **Old E-utilities/species-list/GenBank runs
are incompatible with these Datasets packages.** Move old `genes/GENE/` directories
outside the checkout before downloading that gene fresh. The unchanged steps
03–04 expect `GENE_longest_cds.fasta` and the existing cleaned-alignment layout.

## Track and publish

Maintain one reviewed row per gene in `metadata/gene_tracking.tsv`.
`species_requested` can be all_RefSeq_Lepidosauria; `species_downloaded` counts
named species represented in candidates; `n_species_final` counts reviewed
alignment taxa. Record alignment length, MACSE issues, relative tree path,
obvious tree problems and exclusions. Keep accession/provenance tables with
scratch outputs and summarize important findings in the tracked table.

Generated contents under `genes/`, except .gitkeep, are ignored. Do not commit
ZIP packages, FASTAs, alignments or IQ-TREE outputs, or force-add generated files.

To update an existing Monsoon checkout, transfer the latest source ZIP, extract
it into a separate staging folder, and copy its source files into the checkout.
Keep the checkout's `.git` directory and any customized `genes.txt` list.
After updating the active conda environment and reviewing files:

```bash
git add README.md environment.yml genes.txt species.txt scripts tests metadata .gitignore genes/.gitkeep
git diff --cached --stat
git commit -m "Use NCBI Datasets for RefSeq Lepidosauria gene-list downloads"
git push
```

For a new checkout extracted from ZIP, first run `git init -b main`, then commit
as above, set `git remote add origin https://github.com/YOUR-ACCOUNT/skink-time-limb-genes.git`
and push with `git push -u origin main`.

## Validation and documentation

Run `python -m unittest discover -s tests -v`. Offline checks cover mocked
Datasets commands/packages, source filtering, longest isoforms, tuatara grouping,
ambiguous GeneIDs, failed/empty downloads, gene-list validation, completed-download
skipping, package modification detection, alignment validation and stage
orchestration with mock alignment/tree tools. Live Datasets retrieval, the Monsoon
conda solve and actual MACSE/IQ-TREE inference still require verification there.

- [NCBI Datasets gene downloads](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/gene/datasets_download_gene_symbol/)
- [Datasets package and CDS header format](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/gene-package/)
- [Gene report schema](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/gene/)
- [Conda environment/cache locations](https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/custom-env-and-pkg-locations.html)
- [Monsoon software modules](https://in.nau.edu/arc/installing-software-packages/)
- [Monsoon scratch storage](https://in.nau.edu/arc/overview/file-management/)
- [MACSE documentation](https://www.agap-ge2pop.org/wp-content/uploads/macse/doc/doc_MACSE_v2.03.pdf)
- [IQ-TREE tutorial](https://www.iqtree.org/doc/Tutorial)
