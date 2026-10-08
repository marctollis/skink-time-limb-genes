# SKINK TIME: limb-gene training pipeline

Use this workflow to download RefSeq coding sequences (CDS), select one sequence
per species, align sequences, and build a gene tree. The goal is to compare a
mix of limbed and limb-reduced squamates, including iguanians.

## Set up on Monsoon

From the repository folder:

```bash
module load mambaforge
source "$(conda info --base)/etc/profile.d/conda.sh"

export CONDA_ENVS_PATH="/scratch/$USER/conda/envs"
export CONDA_PKGS_DIRS="/scratch/$USER/conda/pkgs"
mkdir -p "$CONDA_ENVS_PATH" "$CONDA_PKGS_DIRS"

conda env create --prefix "$CONDA_ENVS_PATH/skink-time-limb-genes" \
  --file environment.yml
conda activate "$CONDA_ENVS_PATH/skink-time-limb-genes"
```

Both the environment and package cache are stored in scratch. In later sessions,
load the module, set the two paths, and activate the existing environment.
To update an existing environment after changes to `environment.yml`:

```bash
conda env update --prefix "$CONDA_ENVS_PATH/skink-time-limb-genes" \
  --file environment.yml
```

## 1. Download CDS

Edit `genes.txt` with one gene symbol per line. Start with `SHH` to test the
workflow, then add other genes.

```bash
bash scripts/01_download_cds.sh
```

This step searches NCBI Gene for matching genes across Squamata, then downloads
CDS packages using NCBI Datasets. It processes the entire gene list without a
gene-name prompt. Available species will differ among genes. Predicted RefSeq
transcripts are included.

Results are saved under `genes/<GENE>/`. Check `genes/download_batch_status.tsv`
for failures. Completed matching downloads are skipped when the script is rerun.
`species.txt` is a reference checklist; it does not limit downloads.

## 2. Select one CDS per species

```bash
bash scripts/02_extract_longest_refseq.sh
```

Enter a gene symbol when prompted. This step selects the longest eligible
RefSeq CDS per species and creates:

- `<GENE>_longest_cds.fasta`: selected sequences.
- `candidate_cds.tsv`: downloaded candidates and review notes.
- `selection.tsv`: selected accessions and sequence lengths.

Review species coverage, unusual lengths, and exclusions before alignment.
Longest CDS is a starting rule, not proof of the best isoform or correct orthology.

## 3. Align and clean

Run MACSE on a compute node. Supply the gene explicitly so the job does not
wait for an interactive prompt:

```bash
printf 'SHH\n' | srun --cpus-per-task=1 --mem=4G --time=01:00:00 \
  bash scripts/03_macse.sh
```

Replace `SHH` with the gene being processed. Files and logs are saved under
`genes/<GENE>/alignment/`. Review the raw nucleotide and amino-acid alignments.
The cleaned nucleotide alignment masks stop and frameshift codons; it does not
repair questionable sequences.

## 4. Build a gene tree

After reviewing the alignment:

```bash
printf 'SHH\n' | srun --cpus-per-task=2 --mem=4G --time=01:00:00 \
  bash scripts/04_build_tree.sh
```

IQ-TREE selects a DNA model and estimates branch support. The tree is saved as
`genes/<GENE>/iqtree/<GENE>.treefile`. Inspect its topology and branch lengths.
The tree is unrooted; a gene tree may differ from the species tree.

The resource requests above are starting values. Larger or more difficult genes
may require additional memory or time. Alignment and tree building must run on
compute nodes, not login nodes.

## Record and review results

Maintain one row per gene in `metadata/gene_tracking.tsv`. Record species counts,
alignment length, alignment issues, tree location, and review notes. Check that
each gene includes limbed and limb-reduced taxa and iguanians.

The SHH training run retained 42 species and produced a cleaned alignment of
1,911 bases. The selected *Paroedura picta* isoform had three MACSE frameshift
markers; an alternative isoform is available. Resolve this issue before
codon-based selection analyses. The training tree uses the cleaned alignment.

Keep scripts, gene lists, environment specifications, and tracking metadata in
GitHub. Downloaded packages, FASTAs, alignments, and tree outputs are ignored.

Existing extraction, alignment, and tree outputs are not overwritten. Move old
outputs aside before rerunning a stage with changed inputs. Automated Slurm
array alignment is planned; the current alignment step runs one gene at a time.
