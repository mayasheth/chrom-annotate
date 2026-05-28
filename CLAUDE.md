# chrom-annotate — CLAUDE.md

Reference for AI-assisted development. Keep this up to date.

## Project overview

Snakemake pipeline that annotates genomic elements (enhancers, CRISPR-tested regions, etc.) with chromatin signals from ChIP-seq, ATAC-seq, and DNase-seq. Supports three operational modes via separate Snakefiles.

## Repository layout

```
chrom-annotate/
├── workflow/
│   ├── Snakefile                  # main mode: annotate arbitrary element lists
│   ├── Snakefile_CRISPR           # CRISPR mode: annotate ABC predictions + CRISPR benchmarks
│   ├── Snakefile_neighborhoods    # neighborhoods mode (separate)
│   ├── rules/
│   │   ├── annotate_elements.smk  # core annotation rules
│   │   ├── download.smk           # ENCODE BAM/peak download + preprocessing
│   │   └── utils.smk              # Python helper functions + metadata parsing
│   └── scripts/
│       ├── gather_annotations_for_one_cell_type.R  # merges all annotations per cell type
│       ├── gather_across_cell_types.R               # merges across cell types into final output
│       ├── gather_all_crispr_data.R                 # merges CRISPR benchmark annotations
│       ├── run.neighborhoods.extended.py            # BAM read counting (from ABC model)
│       ├── neighborhoods.py                         # neighborhoods library (from ABC model)
│       └── tools.py                                 # utilities (from ABC model)
├── config/
│   ├── config_JT.yml              # MPRA elements annotation (main Snakefile)
│   ├── config_p300.yml            # TF/p300 peaks annotation (main Snakefile)
│   ├── config_CRISPR.yml          # CRISPR benchmark annotation (Snakefile_CRISPR)
│   ├── ENCODE_K562_ChIP_metadata.tsv          # full ENCODE metadata dump
│   └── filtered_ENCODE_K562_ChIP_metadata.tsv # filtered version used by pipeline
├── resources/
│   ├── allGoodMPRA_EorP_hg38.tsv  # MPRA elements input file
│   └── metadata/
│       └── log.sh
└── results/
    └── 2025_1022_JT/              # completed run output
        └── MPRA_EorP.chromatin_annotations.tsv
```

## Pipeline modes

### Main mode (`Snakefile`)
Annotates user-defined element lists (TSV) with chromatin signals in one or more cell types.

**Invocation:**
```bash
snakemake --snakefile workflow/Snakefile --configfile config/config_JT.yml \
  --profile <slurm_profile> --use-conda
```

**Output:** `{results_dir}/{data_cat}.chromatin_annotations.tsv` — original element table + annotation columns appended.

### CRISPR mode (`Snakefile_CRISPR`)
Annotates both ABC EnhancerList predictions and CRISPR benchmark datasets in parallel.

**Output:** `{results_dir}/{cell_type}/EnhancerList.extended.tsv`, `{results_dir}/CRISPR_data/{data_cat}.extended.tsv`, optionally `annotated_H3K27me3_peaks.tsv.gz`.

## Config schema (main Snakefile)

```yaml
results_dir: <path>
metadata_file: <path to filtered ENCODE metadata TSV>
scratch_dir: <path>             # temp BAM/intermediate storage
base_dir: <project root>
chr_sizes: <path to .chrom.sizes>


# Annotation types (lists of assay names)
peak_overlap_assays: [EP300, CTCF, ...]   # binary peak overlap (0/1)
RPM_assays: [H3K27ac, DNase, ...]         # RPM at element coordinates
RPM_expanded_assays: [H3K27ac]            # RPM at element ± element_ext_size bp
FC_assays: []                             # fold-change bigWig signal

# Elements to annotate
elements_to_annotate:
  <label>: <path to TSV/TSV.gz>
element_cell_types:
  <label>: [<cell_type>, ...]

# Column name mapping (per label)
chr_columns / start_columns / end_columns / cell_type_columns: ...

# Resizing params
element_trim_size: 0   # bp to trim element boundaries for baseline counts
element_ext_size: 150  # bp to expand for expandedRegion counts
peak_ext_size:
  <assay>: <bp>        # slopbed applied to peaks before overlap

# Per-cell-type data (one key per cell type, matching element_cell_types values)
<CellType>:
  processed_files:           # pre-downloaded BAMs/peaks (preferred)
    <assay>:
      reads: [<bam>, ...]
      peaks: <bed.gz>
  experiment_accessions:     # fallback: download from ENCODE
    <assay>: <ENCSR...>
  fold_change_bws:
    <assay>: <URL or path>
```

## Output columns

For each annotated element set, the output TSV contains all original columns plus:
- `{assay}_peak_overlap` — 1/0 binary overlap with peak BED for `peak_overlap_assays`
- `{assay}.RPM` — reads per million in element coordinates for `RPM_assays`
- `{assay}.RPM.expandedRegion` — RPM in ±150 bp expanded region for `RPM_expanded_assays`
- `{assay}_fc.RPM` — fold-change bigWig signal for `FC_assays`
- `elementName` — `chr:start-end` identifier (intermediate, may be dropped in final output)

## Key dependencies

- **Environment:** pixi (`pixi.toml`). Run with `pixi run snakemake` — no `--use-conda` needed.
- **Read counting:** Adapted from ABC model's `run.neighborhoods.extended.py` — uses `samtools view` + bedtools to count BAM reads in regions, normalizes to RPM.
- **ENCODE metadata:** Pipeline parses ENCODE metadata TSV to map experiment accessions to BAM/peak file accessions. Metadata must include both fastq (for run type) and BAM rows for each experiment.

## Data sources

- ENCODE portal: K562, GM12878, HCT116, Jurkat, WTC11, H9, H1 ChIP-seq/DNase/ATAC
- Pre-downloaded BAMs live at `/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/<CellType>/`
- hg38 chrom sizes: `/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes`

## Important patterns

- `process_encode_metadata()` in `utils.smk` resolves assay → BAM/peak paths. Checks `processed_files` first, then `experiment_accessions` (downloads from ENCODE).
- For paired-end BAMs: sorts only. For single-end: filters (flags -780, MAPQ≥30) then sorts.
- Peak extension (`peak_ext_size`) is applied via `bedtools slop` before overlap, not after.
- `cell_type_col: NONE` means all elements belong to one cell type (no filtering by cell type column).
- The CRISPR Snakefile has a WTC11→WTC11 alias mechanism (`get_data_cell_type`) for when CRISPR cell type names differ from data cell type names.

## Environment

Uses pixi (`pixi.toml`). Run `module load pixi && pixi install` once, then invoke with `pixi run snakemake ...` (no `--use-conda`).
