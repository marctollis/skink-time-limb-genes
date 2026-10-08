# SKINK TIME limb-gene training pipeline

A small teaching workflow for Nick: **choose gene → choose species → download
CDS → select one sequence per species → align → build and inspect a gene tree**.
Run one gene at a time. This is a training pipeline, not an orthology or selection
analysis. Longest CDS selection does not establish orthology or identify the
biologically preferred isoform.

## Repository layout

```text
environment.yml
species.txt
scripts/
  01_download_cds.sh
  02_extract_longest_refseq.sh
  03_macse.sh
  04_build_tree.sh
  common.sh
  pipeline.py
metadata/gene_tracking.tsv
genes/.gitkeep
```

The earlier conversation chose a hand-curated species list over a clade-wide
download. The exact original list and script bodies were not available in the
retrieved chat. `species.txt` is an explicitly labeled starter panel; replace or
extend it with the agreed taxa, especially focal skinks. Presence in this list
does not guarantee an NCBI Gene annotation or CDS for a particular gene.

## Set up conda on Monsoon in scratch

Log in to Monsoon and obtain this repository (clone your eventual GitHub URL,
or transfer this folder). Keep the checkout in `/scratch/$USER/` too if desired.
NAU currently documents `mambaforge` and `anaconda3` modules; inspect available
modules if a name changes. Run the following **from the repository root**:

```bash
module avail mambaforge
module load mambaforge
source "$(conda info --base)/etc/profile.d/conda.sh"

# Adjust if your assigned scratch path differs.
export SKINK_SCRATCH="/scratch/$USER"
mkdir -p "$SKINK_SCRATCH/conda/envs" "$SKINK_SCRATCH/conda/pkgs"

# Set both locations explicitly; do not fall back to home or a shared cache.
export CONDA_ENVS_PATH="$SKINK_SCRATCH/conda/envs"
export CONDA_PKGS_DIRS="$SKINK_SCRATCH/conda/pkgs"

# Optional scratch-local configuration, without editing ~/.condarc.
export CONDARC="$SKINK_SCRATCH/conda/condarc"
cat > "$CONDARC" <<EOF
envs_dirs:
  - $SKINK_SCRATCH/conda/envs
pkgs_dirs:
  - $SKINK_SCRATCH/conda/pkgs
channel_priority: strict
EOF

conda env create --prefix "$CONDA_ENVS_PATH/skink-time-limb-genes" --file environment.yml
conda activate "$CONDA_ENVS_PATH/skink-time-limb-genes"
conda config --show envs_dirs pkgs_dirs
conda list
```

Repeat the module loading, `source`, exports and activation in each new session
or job. The environment and package cache are both in scratch on the same
filesystem. The centrally supplied base conda installation stays where NAU
installed it. Do not install a separate base distribution in home.

`environment.yml` constrains tool versions and channels; Python, Java and the
Datasets CLI allow patch/minor updates within the specified release series.
It is a reproducible environment specification, not a complete dependency lock.
For an exact resolved **Monsoon Linux** snapshot, after successful installation:

```bash
conda list --explicit > metadata/conda-linux-64-explicit.txt
# Recreate that exact resolved environment later:
# conda create --prefix "$CONDA_ENVS_PATH/skink-time-limb-genes-locked" \
#   --file metadata/conda-linux-64-explicit.txt
```

Commit that small snapshot after reviewing it for private channel credentials.
Scratch is working storage; keep the repository and environment specifications
backed up on GitHub so they can be rebuilt.

## Run one gene

First edit `species.txt`: one scientific binomial per line, with optional `#`
comments. Each script prompts for the same gene symbol (e.g. `SHH`). Symbols are
normalized to uppercase and all paths are relative to the script location.

```bash
bash scripts/01_download_cds.sh
bash scripts/02_extract_longest_refseq.sh
bash scripts/03_macse.sh
bash scripts/04_build_tree.sh
```

Run alignment and tree inference inside a Monsoon compute allocation, not on a
login node. Request CPUs/memory/time using your lab's current Slurm allocation
instructions. Step 04 uses `SLURM_CPUS_PER_TASK`, or one thread if unset; it
does not guess the available CPUs. For a noninteractive job, feed a gene on stdin:
`printf 'SHH\n' | bash scripts/03_macse.sh`. Repeat for step 04 in that job.

### 01: download

Uses `datasets download gene symbol GENE --taxon "Species name" --include cds`
once per species. Stores the ZIP, extracted package, NCBI metadata and log under
`genes/GENE/downloads/Species_name/`. Saves the actual species list and a status
table in that gene directory. Attempts all species and exits with an error if
any download/unzip fails. Review logs: network errors and unavailable genes must
not be treated as the same biological result. You may continue to step 02 with
the successful packages after reviewing and documenting failures.

### 02: select the longest RefSeq CDS

Reads downloaded `cds.fna` files. Retains CDS records with NM_/XM_ transcript or
NP_/XP_ protein accessions, including predicted RefSeqs. Where a header includes
`[gene=...]`, its symbol must match. Packages with multiple gene IDs are excluded
as ambiguous because symbol queries can match synonyms; resolve those manually.
Selects the longest eligible CDS per species, breaking ties by the original
header. It does not prefer curated NM_ over longer predicted XM_ records.
Non-ACGT IUPAC ambiguity bases are allowed; invalid sequences are excluded.
Lengths not divisible by three are flagged for review, not silently discarded.

Produces `GENE_longest_refseq.fasta` with unique `Species_name` labels and
`selection.tsv` with selected accession, length, original header, candidate
count and exclusions. The Python helper implements selection; SeqKit reports
FASTA statistics. Inspect unusual lengths, ambiguity, duplicate gene matches,
and missing taxa before proceeding. Header/schema changes may require adapting
the helper; zero eligible sequences causes a clear failure.

### 03: MACSE alignment

Requires at least two species. Saves raw nucleotide/amino-acid alignments and
logs under `alignment/`. Uses MACSE's `exportAlignment` to mask internal and
terminal stops with `NNN`, internal frameshift codons with `NNN`, and remaining
frameshift characters with `-`. The raw alignment remains available to inspect
`!` and `*` characters. This is a convenient training DNA alignment; masking is
not biological validation. Review substantial issues before inferring a tree.

### 04: IQ-TREE

Requires at least four species and a valid equal-length cleaned alignment. Fits
a **DNA** model with ModelFinder (`-m MFP`), 1,000 ultrafast bootstrap replicates,
1,000 SH-aLRT replicates and seed 12345. Results go under `iqtree/`; open
`GENE.treefile` in your preferred tree viewer. The tree is unrooted. Root on
`Sphenodon_punctatus` during inspection only if it survived selection and is an
appropriate outgroup. These are gene trees, which may differ from a species tree.
The cleaned alignment is not automatically approved for downstream HyPhy or
codon-model analyses.

All stages refuse to overwrite existing stage outputs. To restart a gene, move
its entire directory aside (outside this checkout if retaining it), then run
from step 01. Do not mix outputs from different species lists or runs.

## Tracking and GitHub

Maintain **one row per gene** in `metadata/gene_tracking.tsv`. Requested species
come from `species.tsv`; downloaded counts come from `download_status.tsv`
(successful packages may still contain no usable CDS). Final counts and alignment
length come from the reviewed alignment/IQ-TREE report. Record MACSE issues,
relative tree path, obvious tree problems and your interpretation in notes.
Keep unavailable genes, download failures and excluded sequences distinct.
Per-species provenance is in generated `selection.tsv`; summarize important
exclusions in the tracked notes.

The entire contents of `genes/`, except `.gitkeep`, are ignored. ZIPs, FASTAs,
alignments and common IQ-TREE outputs are also ignored elsewhere. Version-control
scripts, species lists, environment specifications and small tracking tables.
Keep bulk data in scratch; do not force-add generated files.

After review, create an empty GitHub repository called `skink-time-limb-genes`,
then publish from this folder (substitute your own URL):

```bash
git init -b main  # Only needed if you extracted the ZIP rather than using the initialized checkout.
git add README.md environment.yml species.txt scripts tests metadata .gitignore genes/.gitkeep
git diff --cached --stat
git commit -m "Add SKINK TIME limb-gene training pipeline"
git remote add origin https://github.com/YOUR-ACCOUNT/skink-time-limb-genes.git
git push -u origin main
```

## Documentation and validation limits

- [NCBI gene package downloads](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-gene-data-package/)
- [Conda environment and package locations](https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/custom-env-and-pkg-locations.html)
- [NAU software installation and modules](https://in.nau.edu/arc/installing-software-packages/)
- [NAU storage and scratch paths](https://in.nau.edu/arc/overview/file-management/)
- [MACSE documentation](https://www.agap-ge2pop.org/wp-content/uploads/macse/doc/doc_MACSE_v2.03.pdf)
- [IQ-TREE tutorial](https://www.iqtree.org/doc/Tutorial)

Local validation covers shell syntax, selection/alignment checks with synthetic
fixtures, prompt/path behavior, stage orchestration with mock tools, and Git
ignore rules. The conda solve, live NCBI downloads, MACSE and IQ-TREE must still
be exercised on Monsoon; no real biological outputs or exact dependency lock
are claimed by this initial repository.

Run the offline checks with `python -m unittest discover -s tests -v`.
