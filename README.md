# chrom-annotate
Annotate genomic elements with various chromatin signals (including direct download from the ENCODE portal, if required!)


# Acknowledgments
The following components of this pipeline are resued/adapted from the [ABC model pipeline](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction):
- `workflow/scripts/neighborhoods.py`, `workflow/scripts/run.neighborhoods.extended.py`, and `workflow/scripts/tools.py` for quantifying chromatin signals.
- `utilities.py` for signal normalization

All reused code is licensed under [original license](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/LICENSE), and proper attribution is maintained.
