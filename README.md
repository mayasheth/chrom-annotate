# chrom-annotate

Snakemake pipeline for annotating genomic elements with chromatin signals from ChIP-seq, ATAC-seq, and DNase-seq. Given a list of genomic elements (e.g., enhancers, CRISPR-tested regions) and a set of chromatin assays, the pipeline computes:

- **Peak overlap** — binary indicator of whether an element overlaps a ChIP-seq peak
- **RPM** — read count (reads per million) in the element's original coordinates
- **RPM (expanded region)** — RPM in coordinates expanded by a configurable number of base pairs
- **Fold-change signal** — average fold-change bigWig signal over the element

## Pipeline overview

```
Input elements (TSV)  +  Chromatin data (BAMs / peaks / bigWigs)
        |
        v
  [split by cell type]
        |
        +---> [peak overlap]  bedtools intersect → binary 0/1 per assay
        |
        +---> [read counting] run.neighborhoods.extended.py → RPM per assay
        |                     (original coordinates + expanded coordinates)
        +---> [fold-change]   run.neighborhoods.extended.py on fold-change bigWig
        |
        v
  [gather per cell type]  R script joins all signals onto element table
        |
        v
  [merge across cell types]  R script concatenates cell types
        |
        v
  Output TSV: original columns + chromatin annotation columns
```

## Pipeline modes

The pipeline has three entry points for different use cases.

### Main mode — `workflow/Snakefile`

Annotates arbitrary element lists (TSV or TSV.gz files) with chromatin signals. Suitable for any genomic region set: MPRA elements, finemo peaks, EnhancerLists, etc.

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config_JT.yml \
  --profile <slurm_profile> \
  --use-conda
```

### CRISPR mode — `workflow/Snakefile_CRISPR`

Annotates both ABC model EnhancerList predictions and CRISPR benchmark datasets simultaneously across multiple cell types. Also optionally annotates H3K27me3 peaks with H3K4me1 signal (for bivalent peak analysis).

```bash
snakemake --snakefile workflow/Snakefile_CRISPR \
  --configfile config/config_CRISPR.yml \
  --profile <slurm_profile> \
  --use-conda
```

## Configuration

Each run requires a YAML config file. See `config/config_JT.yml` for a complete example.

### Main mode config fields (`Snakefile`)

| Field | Description |
|---|---|
| `results_dir` | Output directory |
| `scratch_dir` | Temp directory for intermediate files and ENCODE downloads |
| `base_dir` | Absolute path to the repo root (used to resolve relative paths) |
| `chr_sizes` | Path to hg38 chromosome sizes file |
| `metadata_file` | ENCODE metadata TSV (required when using `experiment_accessions`) |
| `elements_to_annotate` | Dict mapping label → path to element TSV (plain or .gz) |
| `element_cell_types` | Dict mapping label → list of cell types to annotate |
| `chr_columns` | Dict mapping label → chromosome column name in element file |
| `start_columns` | Dict mapping label → start coordinate column name |
| `end_columns` | Dict mapping label → end coordinate column name |
| `cell_type_columns` | Dict mapping label → cell type column name, or `NONE` if all elements share one cell type |
| `peak_overlap_assays` | List of assays for binary peak overlap |
| `RPM_assays` | List of assays for RPM counting at element coordinates |
| `RPM_expanded_assays` | List of assays for RPM counting at expanded coordinates |
| `FC_assays` | List of assays for fold-change bigWig signal |
| `element_trim_size` | BP to trim from each side of element for baseline counts (0 = no trimming) |
| `element_ext_size` | BP to extend each side of element for expanded-region counts (default: 150 bp) |
| `peak_ext_size` | Dict mapping assay → BP to extend peaks before overlap; use 175 bp for broad histone marks (H3K27me3, H3K4me1, H3K27ac), 0 for TFs (e.g., CTCF) |

### Specifying chromatin data (main mode)

For each cell type listed in `element_cell_types`, add a config section keyed by cell type name. Data can be provided three ways (checked in order):

**Pre-downloaded files (preferred):**
```yaml
K562:
  processed_files:
    H3K27ac:
      reads: [/path/to/rep1.bam, /path/to/rep2.bam]
      peaks: /path/to/peaks.bed.gz  # only needed if assay is in peak_overlap_assays
```

**Auto-download from ENCODE (by experiment accession):**
```yaml
K562:
  experiment_accessions:
    GABPA: ENCSR000BLO   # pipeline parses metadata TSV to find BAM accessions, then downloads
```

**Fold-change bigWigs (for `FC_assays` only):**
```yaml
K562:
  fold_change_bws:
    H3K27ac: https://www.encodeproject.org/files/ENCFF.../@@download/...bigWig
```

### CRISPR mode config fields (`Snakefile_CRISPR`)

The CRISPR mode uses a different config schema. It annotates both ABC model EnhancerList outputs and CRISPR benchmark TSVs.

| Field | Description |
|---|---|
| `results_dir` | Output directory |
| `scratch_dir` | Temp directory |
| `chr_sizes` | Path to hg38 chromosome sizes file |
| `pred_cell_types` | List of cell types to annotate ABC predictions for |
| `peak_overlap_assays` | List of assays for binary peak overlap |
| `RPM_assays` | List of assays for RPM counting |
| `RPM_expanded_assays` | List of assays for RPM at expanded coordinates |
| `pred_trim_size` | BP to trim prediction element boundaries (to match CRISPR element sizes) |
| `pred_ext_size` | BP to extend prediction elements for expanded-region counts (default: 150 bp) |
| `peak_ext_size` | Dict mapping assay → BP to extend peaks before overlap; use 175 bp for broad histone marks (H3K27me3, H3K4me1, H3K27ac), 0 for TFs (e.g., CTCF) |
| `CRISPR_benchmark` | Dict mapping dataset label → path to CRISPR benchmark TSV.gz |
| `CRISPR_cell_types` | Dict mapping dataset label → list of cell types present in that dataset |
| `chr_columns` / `start_columns` / `end_columns` / `cell_type_columns` | Column name mappings, keyed by dataset label and `predictions` |
| `WTC11_cell_type` | Cell type name to use for WTC11 data lookup (allows aliasing if CRISPR data uses a different name) |

Per-cell-type data in CRISPR mode uses `E2G_results`, `bam`, and `bed` keys:
```yaml
K562:
  E2G_results: /path/to/ABC/results/K562_run   # directory containing Neighborhoods/EnhancerList.{bed,txt}
  bam:
    H3K27ac: [/path/rep1.bam, /path/rep2.bam]
    CTCF: [/path/rep1.bam]
  bed:
    H3K27ac: /path/to/peaks.bed.gz
    CTCF: /path/to/peaks.bed.gz
```

## Input format

Element files must be tab-separated with at minimum chromosome, start, and end columns. Column names are configured per-label via `chr_columns`, `start_columns`, `end_columns`. Files may be plain TSV or gzip-compressed.

Set `cell_type_columns` to `NONE` for a given label if all elements in that file belong to one cell type (no filtering applied). Otherwise specify the column name containing cell type values — the pipeline will subset to rows matching each cell type in `element_cell_types`.

## Output format

One output TSV per element label at `{results_dir}/{label}.chromatin_annotations.tsv`. Contains all original columns plus appended annotation columns:

| Column pattern | Description |
|---|---|
| `{assay}_peak_overlap` | 1 if element overlaps a peak, 0 otherwise |
| `{assay}.RPM` | Reads per million in element coordinates |
| `{assay}.RPM.expandedRegion` | RPM in element ± `element_ext_size` bp |
| `{assay}_fc.RPM` | Fold-change bigWig signal |

## Dependencies and environment

A pixi environment (`pixi.toml`) is provided with all required dependencies:

```bash
module load pixi
pixi install
```

Run the pipeline within the pixi environment (omit `--use-conda`, since all deps are bundled):

```bash
pixi run snakemake --snakefile workflow/Snakefile \
  --configfile config/config_JT.yml \
  --profile <slurm_profile>
```

**Included packages:**
- Python: `numpy`, `pandas`, `pysam`, `pyBigWig`, `pyranges`, `scipy`
- R: `dplyr`, `tidyr`, `data.table`, `stringr`
- CLI: `samtools` (1.23), `bedtools` (2.31), `csvtk`, `curl`
- Workflow: `snakemake` (>=7)

Read counting reuses scripts from the [ABC model pipeline](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction).

## Acknowledgments

The following components are reused/adapted from the [ABC model pipeline](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction):
- `workflow/scripts/neighborhoods.py`, `workflow/scripts/run.neighborhoods.extended.py`, and `workflow/scripts/tools.py` for quantifying chromatin signals

All reused code is licensed under the [original license](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/LICENSE), and proper attribution is maintained.
