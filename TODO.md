# TODO

## Active tasks

- [ ] **Document pipeline** — expand README and inline comments; add usage examples and a diagram of the workflow DAG; document all config fields
- [x] **Create pixi environment for pipeline and test that it works** — `pixi.toml` created with Sherlock system requirements; all Python/R/CLI deps verified; use `pixi run snakemake` without `--use-conda`
- [ ] **Describe options/capacity of annotations and input/output file formats** — document all annotation types (peak overlap, RPM, RPM expanded, fold-change), their column naming conventions, and accepted input formats (TSV/TSV.gz, BED/BED.gz); add to README
- [ ] **Organize metadata and available genomic data files** — audit which cell types and assays are available in `/oak/.../Data/ENCODE/`; create a manifest/table of available BAMs, peaks, and cell types; clean up config metadata files
- [ ] **Add documentation and helper scripts for chromatin categorization** — document the chromatin state categories used downstream (active enhancer, bivalent, Polycomb, etc.); add helper R/Python scripts to classify elements from annotation output columns

## Completed

<!-- move items here when done -->
