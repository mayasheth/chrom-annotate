# TODO

## Active tasks

- [ ] **Document pipeline** — expand README and inline comments; add usage examples and a diagram of the workflow DAG; document all config fields
- [x] **Create pixi environment for pipeline and test that it works** — `pixi.toml` created with Sherlock system requirements; all Python/R/CLI deps verified; use `pixi run snakemake` without `--use-conda`
- [ ] **Describe options/capacity of annotations and input/output file formats** — document all annotation types (peak overlap, RPM, RPM expanded, fold-change), their column naming conventions, and accepted input formats (TSV/TSV.gz, BED/BED.gz); add to README
- [ ] **Organize metadata and available genomic data files** — audit which cell types and assays are available in `/oak/.../Data/ENCODE/`; create a manifest/table of available BAMs, peaks, and cell types; clean up config metadata files
- [ ] **Add documentation and helper scripts for chromatin categorization** — document the chromatin state categories used downstream (active enhancer, bivalent, Polycomb, etc.); add helper R/Python scripts to classify elements from annotation output columns

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
