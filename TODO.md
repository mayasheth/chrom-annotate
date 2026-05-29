# TODO

## Active tasks

- [x] **Document pipeline** — README has pipeline overview, all config fields, invocation examples (pixi run), input/output format, annotation types, chromatin categorization workflow
- [x] **Create pixi environment for pipeline and test that it works** — `pixi.toml` created with Sherlock system requirements; all Python/R/CLI deps verified; use `pixi run snakemake` without `--use-conda`
- [x] **Describe options/capacity of annotations and input/output file formats** — README documents all annotation types (peak overlap, RPM, RPM expanded, fold-change), column naming, BAM/tagAlign/peak input formats
- [ ] **Organize metadata and available genomic data files** — audit which cell types and assays are available in `/oak/.../Data/ENCODE/`; create a manifest/table of available BAMs, peaks, and cell types; clean up config metadata files
- [x] **Add documentation and helper scripts for chromatin categorization** — `workflow/scripts/chromatin_categories.R` with `get_category_thresholds`, `get_threshold_key`, `categorize_elements` + canonical colors/order; README section with required assays and usage example

## Chromatin categorization improvements

- [ ] Remove `CTCF.H3K27ac.ratio` — confirmed unused: computed in `get_category_thresholds` as `CTCF.H3K27ac.ratio.CTCF_peak` but never retrieved by `get_threshold_key` or referenced in `categorize_elements`; remove from `chromatin_categories.R`, README example, and any upstream `mutate()` calls
- [ ] General cleanup pass on `chromatin_categories.R` — verify all computed thresholds are actually used, check edge cases (e.g. cell types missing from thresholds table)
- [x] Ensure pixi environment supports categorization code — added `r-ggplot2` to `pixi.toml`; `rlang` available transitively via `dplyr`
- [x] Add ATAC-seq support as alternative to DHS/DNase for accessibility metric — `get_category_thresholds` now accepts `DHS.RPM` or `ATAC.RPM`; tagAlign files already handled by `neighborhoods.py` (auto-detected by filename); provide `.tagAlign.gz` under `reads` in config

## Code review fixes

- [x] `gather_across_cell_types.R`: always drop `elementName` before writing output
- [x] `annotate_elements.smk`: add guard when cell type column name not found in header
- [x] `annotate_elements.smk`: validate col3 > col2 after trim_size applied
- [x] `download.smk`: add `resources:` block to `download_peaks`
- [x] `download.smk`: fix `=` vs `==` in bash conditional
- [x] `utils.smk`: deduplicate `RPM_assays + RPM_expanded_assays`
- [x] `utils.smk`: add explicit check that each cell type has a config section
- [x] `utils.smk` `get_metadata_for_rpm_assay`: fall back to "alignments" if "unfiltered alignments" not found for SE
- [x] `utils.smk` / `Snakefile`: drop unused return values from `process_encode_metadata`
- [x] `neighborhoods.py` `count_bam`: use context manager for pysam file handle
- [x] `gather_annotations_for_one_cell_type.R`: remove debug `print(assays)` from `get_assays()`
- PE BAM filtering: intentionally not applied — SE and PE treated differently per ENCODE guidelines

## Completed

<!-- move items here when done -->
