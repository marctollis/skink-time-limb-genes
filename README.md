# SKINK TIME limb-gene training pipeline

Interactive training pipeline: **gene → NCBI annotated CDS → one CDS per species
→ MACSE alignment → IQ-TREE gene tree**.

Provide gene symbols in **`genes.txt`**, then step 01 retrieves all matching
**RefSeq** Nucleotide records across **Lepidosauria**, including squamates and
tuatara. RefSeq predictions are included. There is no ten-species limit.
`species.txt` is an optional reference panel and does not restrict downloads.
Steps 02–04 still prompt for one gene at a time.

“All” means all records returned by the annotated gene-symbol search, not every
homolog in every genome. Unannotated genes and annotations under other symbols
may be missed. This workflow downloads CDS-bearing source records; it does not
use BLAST to recover genes. Longest CDS alone does not establish orthology.

## Files

- `scripts/01_download_cds.sh`: batch RefSeq Lepidosauria retrieval from `genes.txt`.
- `scripts/02_extract_longest_refseq.sh`: annotated CDS extraction and selection.
- `scripts/03_macse.sh`: raw and cleaned codon-aware alignments.
- `scripts/04_build_tree.sh`: DNA model selection and gene tree inference.
- `scripts/ncbi_cds.py`: retrieval, annotation parsing and provenance tables.
- `scripts/pipeline.py`, `scripts/common.sh`: validation and shared prompts.
- `environment.yml`: software specification.
- `genes.txt`: one gene symbol per line; starts with SHH as an example.
- `species.txt`: optional reference panel, independent of retrieval.
- `metadata/gene_tracking.tsv`: one manually reviewed row per gene.
- `genes/.gitkeep`: placeholder; all generated gene directories are ignored.
- `tests/test_pipeline.py`: offline checks with synthetic annotations/mock tools.

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

Biopython is now required. The previous Datasets CLI/unzip packages may remain
installed but are no longer used. Version constraints are not a complete
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

## Download the gene list, then process each gene

Edit `genes.txt` with your chosen symbols, one per line. The supplied `SHH` is
an example, not an assumed final limb-gene panel. Blank lines and `#` comments
are allowed. Symbols normalize to uppercase; duplicates and unsafe symbols fail
before any network requests. Every gene keeps its own `genes/GENE/` directory.

Set your NCBI contact email, then run the list downloader:

```bash
export NCBI_EMAIL='your-real-email@example.org'
bash scripts/01_download_cds.sh
# Optional alternative list:
# bash scripts/01_download_cds.sh /path/to/my_genes.txt
```

Step 01 prompts for email if unset, but does not ask for a single gene. It runs
sequentially through the list and writes `genes/download_batch_status.tsv`.
A failed/zero-result gene is recorded while later genes are still attempted;
the batch exits with an error if any failed. Complete matching downloads are
skipped on rerun. Incomplete or older-scope directories must be moved aside
before retrying that gene. The list itself never chooses a preferred species.

An optional `NCBI_API_KEY` environment variable is supported. Neither the key
nor contact email is saved in run metadata. Biopython handles NCBI request pacing
and transient retries. Do not run many download processes simultaneously:
NCBI rate limits can apply across a shared IP address.

For each downloaded gene, run step 02 and enter its symbol:

```bash
bash scripts/02_extract_longest_refseq.sh
```

Inspect `genes/SHH/candidate_cds.tsv` and `genes/SHH/selection.tsv`, then run these
inside a Monsoon compute allocation using your lab's current Slurm instructions:

```bash
bash scripts/03_macse.sh
bash scripts/04_build_tree.sh
```

Do not run alignment/tree inference on a login node. Step 04 uses
`SLURM_CPUS_PER_TASK`, or one thread if unset. For jobs, supply the gene on stdin:
`printf 'SHH\n' | bash scripts/03_macse.sh` (likewise for step 04).

### 01: retrieve all matching Nucleotide records

Uses NCBI ESearch/EFetch with a query such as:

```text
("SHH"[Gene Name]) AND Lepidosauria[Organism] AND srcdb_refseq[PROP]
```

Retrieves all search results through server history in batches of 20, saving
GenBank-format records under `genes/GENE/downloads/`. `download_manifest.json`
records the query, translated query, UTC retrieval date, counts and batch files.
Records may be transcripts, individual genomic sequences or annotated
chromosomes. Large genomic records can make retrieval slow and disk-intensive.
There is no fixed result-count truncation. Some records will have no matching
usable CDS; download counts are **not species counts**.

Each batch is parsed and its record count verified. The manifest is marked
complete only when every batch succeeds. Incomplete downloads cannot proceed to
extraction. A zero-result gene is reported as failed and later genes are still attempted. If interrupted or failed, move
the incomplete gene directory aside and restart; automatic resume is not provided.

### 02: extract and select CDS from annotations

Reads matching CDS feature `/gene` or `/gene_synonym` qualifiers, using exact
case-insensitive symbol matches. A record-level search hit is insufficient:
unrelated CDS on the same chromosome are not selected. Parses exon joins and
reverse-strand locations using Biopython. Only named species within the target
lineage are eligible; subspecies are grouped under the species binomial, while
original organism names and source taxids remain in the tables.

Writes:

- `candidate_cds.tsv`: all matching CDS candidates, source/accessions, gene IDs,
  feature location, eligibility, selection and exclusion/review notes.
- `selection.tsv`: selected candidate for each included species.
- `extraction_summary.json`: retrieved records, candidates and selected species,
  including source counts and records without a matching CDS feature.
- `GENE_longest_cds.fasta`: one sequence per species, labeled `Species_name`.

Default selection excludes fuzzy/partial CDS annotations, noninitial reading
frames, pseudogenes, remote feature locations requiring other records,
translation exceptions, nonstandard genetic codes, invalid sequences, lengths
not divisible by three and internal stop codons. These exclusions remain visible
in the candidate table. Full gene sequences with introns are not aligned as CDS.
Ambiguous bases, missing terminal stops and non-ATG starts are flagged for review.
Excluded candidates may have a blank CDS length if extraction was not attempted.

Only RefSeq candidates are eligible. Among them, longest CDS wins;
accession/location breaks equal-length ties. RefSeq includes predicted XM_/XP_ records.
Different GeneIDs within a species are flagged for review. **Eligible/selected
means passed these mechanical checks, not verified orthology or completeness.**
Check isoforms, paralogs, duplicate assemblies and unusual lengths. A simple gene
name search can miss LOC-only annotations, alternative names, unannotated genomes
and records outside the searched index. This pipeline does not promise every
NCBI species or every deposited sequence for the gene.

### 03: MACSE

Requires at least two sequences. Saves original nucleotide/amino-acid alignments
and logs under `alignment/`. MACSE export masks internal/terminal stop codons and
internal frameshift codons with `NNN`, and remaining frameshift characters with
`-`. Inspect raw `!`/`*` characters and logs; masking is not biological validation.
The cleaned alignment must have equal sequence lengths divisible by three.

### 04: IQ-TREE

Requires at least four species. Uses a DNA model (`-st DNA`), ModelFinder (`-m
MFP`), 1,000 ultrafast bootstrap and SH-aLRT replicates, and seed 12345. Outputs
are under `iqtree/`; inspect `GENE.treefile`. Trees are unrooted; root using
`Sphenodon_punctatus` in a viewer only if present and appropriate. A gene tree
may differ from the species tree. Outputs are not automatically approved for
HyPhy or codon-model analysis.

Stages refuse to overwrite existing outputs. **Earlier species-list or GenBank-inclusive runs are incompatible
with this RefSeq Lepidosauria scope.** Move old `genes/GENE/` directories
outside the checkout if retaining them, then start fresh from step 01. Do not
mix old/new runs or change the search scope between retrieval and extraction.

## Tracking, updating and publishing

Maintain one reviewed row per gene in `metadata/gene_tracking.tsv`. For this
clade-wide version, `species_requested` can be `all_RefSeq_Lepidosauria`;
`species_downloaded` counts distinct named species in matching candidate rows,
and `n_species_final` counts taxa in the final reviewed alignment. Record alignment
length, MACSE issues, relative tree path, tree problems and exclusions in notes.
The generated provenance tables stay with the scratch results; summarize key
findings in the tracked table.

All per-gene content except `.gitkeep` is ignored. Do not commit downloaded NCBI
records, FASTAs, alignments or IQ-TREE results, or force-add generated files.
For a new checkout extracted from the source ZIP:

```bash
git init -b main
git add README.md environment.yml genes.txt species.txt scripts tests metadata .gitignore genes/.gitkeep
git diff --cached --stat
git commit -m "Batch download RefSeq Lepidosauria CDS by gene list"
git remote add origin https://github.com/YOUR-ACCOUNT/skink-time-limb-genes.git
git push -u origin main
```

If transferring these updated files into your **existing** Monsoon checkout,
keep that checkout's `.git` directory, replace only source files, then run the
`git add`, review, commit and push commands above. Skip `git init` and
`git remote add` when they are already configured. Do not extract under a nested
`skink-time-limb-genes/skink-time-limb-genes` directory by accident.

## Checks and documentation

```bash
python -m unittest discover -s tests -v
```

Offline checks cover annotated GenBank/RefSeq CDS, exon joins, reverse strands,
partial/pseudogene/paralog handling, gene-list validation and batch failure reporting,
retrieval pagination/incomplete downloads,
alignment validation, prompt/path behavior and stage orchestration with mock
alignment/tree tools. No live NCBI retrieval, Monsoon conda solve or actual
MACSE/IQ-TREE inference is claimed by local checks.

- [NCBI Nucleotide sources](https://www.ncbi.nlm.nih.gov/books/NBK44863/)
- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/)
- [Biopython feature extraction](https://biopython.org/docs/latest/Tutorial/chapter_seq_annot.html)
- [Conda environment/cache locations](https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/custom-env-and-pkg-locations.html)
- [Monsoon modules](https://in.nau.edu/arc/installing-software-packages/)
- [Monsoon scratch storage](https://in.nau.edu/arc/overview/file-management/)
- [MACSE documentation](https://www.agap-ge2pop.org/wp-content/uploads/macse/doc/doc_MACSE_v2.03.pdf)
- [IQ-TREE tutorial](https://www.iqtree.org/doc/Tutorial)
