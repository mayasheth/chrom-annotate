# chrom-annotate

Snakemake pipeline for annotating genomic elements with chromatin signals from ChIP-seq, ATAC-seq, and DNase-seq. Given a list of genomic elements (e.g., enhancers, CRISPR-tested regions) and a set of chromatin assays, the pipeline computes:

- **Peak overlap** — binary indicator of whether an element overlaps a ChIP-seq peak
- **RPM** — read count (reads per million) in the element's original coordinates
- **RPM (expanded region)** — RPM in coordinates expanded by a configurable number of base pairs
- **Fold-change signal** — average fold-change bigWig signal over the element (deprecatd)

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

Annotates arbitrary element lists (TSV or TSV.gz files) with chromatin signals. Suitable for any genomic region set specified by chromosome, start, and end location. 

```bash
pixi run snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yml \
  --profile profiles/slurm
```

### CRISPR mode (deprecated) — `workflow/Snakefile_CRISPR`

Annotates both ABC model EnhancerList predictions and CRISPR benchmark datasets simultaneously across multiple cell types. Also optionally annotates H3K27me3 peaks with H3K4me1 signal (for bivalent peak analysis).

```bash
pixi run snakemake --snakefile workflow/Snakefile_CRISPR \
  --configfile config/config_CRISPR.yml \
  --profile profiles/slurm
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

### Element files

Tab-separated (TSV or TSV.gz) with at minimum chromosome, start, and end columns. Column names are mapped via `chr_columns`, `start_columns`, `end_columns` in the config. All other columns are passed through unchanged to the output.

If elements span multiple cell types in one file, set `cell_type_columns` to the column name containing cell type labels. The pipeline will subset rows per cell type. Set to `NONE` if all elements belong to a single cell type.

### BAM and tagAlign files

- **BAM** (ChIP-seq, DNase-seq): must be sorted (coordinate order) and indexed (`.bai` file alongside). Pre-filtered BAMs are preferred: for single-end data, ENCODE recommends filtering with flags `-F 780` and `MAPQ ≥ 30`; paired-end data is used as-is.
- **tagAlign** (ATAC-seq): provide `.tagAlign.gz` files directly as the `reads` list. The pipeline auto-detects the format based on filename and uses the appropriate counting method. No BAM indexing is required.
- Multiple files per assay are supported (replicates are counted jointly and averaged).

```yaml
K562:
  processed_files:
    ATAC:
      reads: [/path/to/rep1.tagAlign.gz, /path/to/rep2.tagAlign.gz]
```

### Peak files

- BED format, at minimum 3 columns (chr, start, end); additional columns are ignored
- Plain BED or gzip-compressed (`.bed.gz`) are both accepted
- Peaks are extended by `peak_ext_size[assay]` bp before overlap; use larger values (e.g., 175 bp) for broad marks (H3K27ac, H3K27me3, H3K4me1) and 0 for sharp TF peaks (CTCF, EP300)

## Annotation types and output format

One output TSV per element label at `{results_dir}/{label}.chromatin_annotations.tsv`, containing all original columns plus the annotation columns described below. An assay can appear in multiple lists (e.g., H3K27ac in both `RPM_assays` and `RPM_expanded_assays`).

### Peak overlap — `peak_overlap_assays`

Binary indicator of whether the element overlaps a ChIP-seq or ATAC-seq peak.

- Output column: `{assay}_peak_overlap` (1 or 0)
- Requires: a peak BED file for the assay (the `peaks` key under `processed_files`)
- Peaks are extended by `peak_ext_size[assay]` bp (symmetric) before the overlap is tested
- Use for: detecting whether an element falls within a called peak region (e.g., CTCF binding, H3K27me3 domains)

### RPM at element coordinates — `RPM_assays`

Read count (reads per million mapped reads) in the element's own coordinates, after optional trimming.

- Output column: `{assay}.RPM`
- Requires: one or more BAM files for the assay
- `element_trim_size` shrinks the element by this many bp from each side before counting (set to 0 to use original coordinates)
- Use for: quantifying signal at the element itself (e.g., H3K27ac, DNase/ATAC accessibility)

### RPM at expanded coordinates — `RPM_expanded_assays`

Same as RPM but at coordinates expanded by `element_ext_size` bp on each side (default: 150 bp).

- Output column: `{assay}.RPM.expandedRegion`
- Requires: BAM files
- Captures signal in the neighborhood around the element, which is more robust for broad marks and for small elements
- Use for: H3K27ac signal used in chromatin categorization (the expanded window smooths over peak boundaries)

### Fold-change signal — `FC_assays`

Average fold-change over input signal from a pre-computed bigWig, at element coordinates.

- Output column: `{assay}_fc.RPM`
- Requires: a bigWig file (local path or URL) under `fold_change_bws` in the cell-type config
- The bigWig is downloaded at runtime if a URL is given
- Use for: normalized signal tracks from ENCODE (fold-change over control), useful when BAMs are not available or when you want a signal that is already input-normalized

## Chromatin categorization

`workflow/scripts/chromatin_categories.R` provides helper functions to assign each annotated element to one of five mutually exclusive chromatin categories. This categorization was developed for the [DC-TAP paper](https://github.com/EngreitzLab/DC_TAP_Paper) and is adapted here for general use.

### Categories

| Category | Criteria |
|---|---|
| High H3K27ac | H3K27ac.RPM.expandedRegion ≥ 90th percentile (among H3K27ac-peak elements), or overlaps H3K27ac peak and is above that threshold |
| H3K27ac | Overlaps H3K27ac peak, or H3K27ac.RPM.expandedRegion ≥ 50th percentile |
| CTCF element | Overlaps CTCF peak; no appreciable H3K27ac |
| H3K27me3 element | Overlaps H3K27me3 peak; no appreciable H3K27ac or CTCF |
| No H3K27ac | None of the above |

Thresholds are computed from the genome-wide distribution of candidate elements within each cell type, so they are data-driven and cell-type-specific.

### Recommended workflow

The categorization thresholds should be computed from a **genome-wide reference set** — typically all non-promoter candidate elements from the E2G universe — and then applied to whatever target set you want to categorize (e.g., CRISPR-tested pairs, E-G predictions). This ensures the quantile cutoffs reflect the full chromatin landscape rather than a biased subset.

**Step 1: Annotate the genome-wide reference set**

Run chrom-annotate on all E2G candidate elements (one TSV per cell type, or a combined multi-cell-type TSV) using the required assays listed below. This is typically the `EnhancerList` output from the ABC pipeline.

**Step 2: Annotate your target set**

Run chrom-annotate on the elements you want to categorize (e.g., CRISPR benchmark elements, E-G pairs). This can be the same run as step 1 if your target elements are a subset of the universe.

**Step 3: Compute thresholds from the reference set**

```r
library(dplyr)
library(data.table)
source("workflow/scripts/chromatin_categories.R")

# Load genome-wide annotations; add cell_type column if not present
enh_universe <- fread("results/my_run/EnhancerList.chromatin_annotations.tsv") %>%
    mutate(cell_type = "K562") %>%   # omit if cell_type column already exists
    filter(class != "promoter") %>%  # exclude promoters before computing thresholds
    replace(is.na(.), 0)

thresholds <- get_category_thresholds(enh_universe, quantiles = c(0.5, 0.9))
```

**Step 4: Apply categories to the target set**

```r
pairs <- fread("results/my_run/crispr_pairs.chromatin_annotations.tsv") %>%
    mutate(cell_type = "K562") %>%
    replace(is.na(.), 0)

pairs_cat <- categorize_elements(pairs, thresholds)
table(pairs_cat$element_category)
```

For multiple cell types, bind rows from all cell types before calling `get_category_thresholds` — thresholds are computed and applied per cell type via the `cell_type` column.

**Step 5: Plot**

```r
library(ggplot2)

ggplot(pairs_cat, aes(x = element_category, fill = element_category)) +
    geom_bar() +
    scale_fill_manual(values = CATEGORY_COLORS) +
    scale_x_discrete(limits = CATEGORY_ORDER) +
    coord_flip() +
    labs(x = NULL, y = "Number of elements") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"))
```

### Required chrom-annotate assays

The following assays must be run to use the categorization functions:

```yaml
peak_overlap_assays: [H3K27ac, CTCF, H3K27me3]
RPM_assays:          [H3K27ac, CTCF, DHS]   # use ATAC instead of DHS if providing tagAlign files
RPM_expanded_assays: [H3K27ac, H3K27me3]
```

`get_category_thresholds` accepts either `DHS.RPM` or `ATAC.RPM` as the accessibility column — whichever is present in the input data.

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
  --profile profiles/slurm
```

Before running, update `profiles/slurm/config.yaml` with your own SLURM settings — at minimum set `slurm_account` to your group account (e.g. run `sacctmgr show user $USER withassoc format=account -n` to find it) and adjust `slurm_partition`, `mem_mb`, and `runtime` defaults as needed.

**Included packages:**
- Python: `numpy`, `pandas`, `pysam`, `pyBigWig`, `pyranges`, `scipy`
- R: `dplyr`, `tidyr`, `data.table`, `stringr`, `ggplot2`
- CLI: `samtools` (1.23), `bedtools` (2.31), `csvtk`, `curl`
- Workflow: `snakemake` (>=9), `snakemake-executor-plugin-slurm`

Read counting reuses scripts from the [ABC model pipeline](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction).

## Acknowledgments

The following components are reused/adapted from the [ABC model pipeline](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction):
- `workflow/scripts/neighborhoods.py`, `workflow/scripts/run.neighborhoods.extended.py`, and `workflow/scripts/tools.py` for quantifying chromatin signals

All reused code is licensed under the [original license](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/LICENSE), and proper attribution is maintained.
