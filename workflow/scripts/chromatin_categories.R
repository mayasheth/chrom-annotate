## Chromatin categorization functions for chrom-annotate output
##
## These functions assign each genomic element to one of five mutually exclusive
## chromatin categories based on peak overlap and quantile-normalized signal values.
## Thresholds are computed from the genome-wide distribution of candidate elements
## within each cell type, making the categorization data-driven and cell-type-specific.
##
## Categories (in priority order):
##   "High H3K27ac"    -- strong H3K27ac signal (above 90th percentile by default)
##   "H3K27ac"         -- moderate H3K27ac signal or overlaps H3K27ac peak
##   "CTCF element"    -- overlaps CTCF peak, no appreciable H3K27ac
##   "H3K27me3 element"-- overlaps H3K27me3 peak, no appreciable H3K27ac or CTCF
##   "No H3K27ac"      -- none of the above
##
## Source: adapted from EngreitzLab/DC_TAP_Paper (interpretation_analysis/)
## Required chrom-annotate assays:
##   peak_overlap_assays: [H3K27ac, CTCF, H3K27me3]
##   RPM_assays:          [H3K27ac, CTCF, DHS]  # DHS or ATAC accepted
##   RPM_expanded_assays: [H3K27ac, H3K27me3]

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})


## Compute per-cell-type quantile thresholds from a genome-wide element table.
##
## Args:
##   enh       data.frame with one row per element, must contain columns:
##               cell_type, CTCF_peak_overlap, CTCF.H3K27ac.ratio,
##               H3K27ac_peak_overlap, H3K27ac.RPM,
##               H3K27me3_peak_overlap, H3K27me3.RPM.expandedRegion,
##               CTCF.RPM, H3K27ac.RPM.expandedRegion, DHS.RPM or ATAC.RPM
##   quantiles numeric vector of quantile probabilities to compute (e.g. c(0.5, 0.9))
##
## Returns:
##   long-format data.frame: cell_type | feature | quantile | value
get_category_thresholds <- function(enh, quantiles) {
    enh_ctcf <- enh %>% filter(CTCF_peak_overlap == 1) %>%
        select(cell_type, CTCF.H3K27ac.ratio.CTCF_peak = CTCF.H3K27ac.ratio) %>%
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>%
        group_by(cell_type, feature) %>%
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    enh_h3k27ac <- enh %>% filter(H3K27ac_peak_overlap == 1) %>%
        select(cell_type, H3K27ac.RPM.H3K27ac_peak = H3K27ac.RPM) %>%
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>%
        group_by(cell_type, feature) %>%
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    enh_h3k27me3 <- enh %>% filter(H3K27me3_peak_overlap == 1) %>%
        select(cell_type, H3K27me3.RPM.expandedRegion.H3K27me3_peak = H3K27me3.RPM.expandedRegion) %>%
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>%
        group_by(cell_type, feature) %>%
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    acc_col <- if ("DHS.RPM" %in% names(enh)) "DHS.RPM" else "ATAC.RPM"
    enh_other <- enh %>%
        select(cell_type, CTCF.RPM, H3K27ac.RPM, H3K27ac.RPM.expandedRegion, !!rlang::sym(acc_col)) %>%
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>%
        group_by(cell_type, feature) %>%
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    rbind(enh_ctcf, enh_h3k27ac, enh_h3k27me3, enh_other)
}


## Look up threshold values for one feature and quantile across cell types.
##
## Args:
##   thresholds   output of get_category_thresholds()
##   feature_col  character, name of the feature (e.g. "H3K27ac.RPM.expandedRegion")
##   quantile_this numeric, quantile probability to retrieve (e.g. 0.9)
##
## Returns:
##   named numeric vector: names are cell types, values are threshold values
get_threshold_key <- function(thresholds, feature_col, quantile_this) {
    filt <- thresholds %>% filter(feature == feature_col, quantile == quantile_this)
    setNames(filt$value, filt$cell_type)
}


## Assign each element to a chromatin category.
##
## Categorization priority (first matching rule wins):
##   1. H3K27ac_peak_overlap == 1 AND H3K27ac.RPM.expandedRegion >= H3K27ac_q_high  -->  "High H3K27ac"
##   2. H3K27ac_peak_overlap == 1                                                    -->  "H3K27ac"
##   3. H3K27ac.RPM.expandedRegion >= H3K27ac_q_high                               -->  "High H3K27ac"
##   4. H3K27ac.RPM.expandedRegion >= H3K27ac_q_low                                -->  "H3K27ac"
##   5. CTCF_peak_overlap == 1                                                       -->  "CTCF element"
##   6. H3K27me3_peak_overlap == 1                                                   -->  "H3K27me3 element"
##   7. (default)                                                                    -->  "No H3K27ac"
##
## Thresholds (rules 1, 3, 4) are looked up per-cell-type from the thresholds table,
## using the H3K27ac.RPM.expandedRegion distribution among H3K27ac-peak elements.
##
## Args:
##   enh            data.frame with element annotations; must have a `cell_type` column
##                  and the columns listed in get_category_thresholds()
##   thresholds     output of get_category_thresholds()
##   H3K27ac_q_high quantile for "High H3K27ac" cutoff (default: 0.9)
##   H3K27ac_q_low  quantile for "H3K27ac" cutoff (default: 0.5)
##
## Returns:
##   enh with an added `element_category` character column
categorize_elements <- function(enh, thresholds, H3K27ac_q_high = 0.9, H3K27ac_q_low = 0.5) {
    key_high <- get_threshold_key(thresholds, "H3K27ac.RPM.expandedRegion", H3K27ac_q_high)
    key_low  <- get_threshold_key(thresholds, "H3K27ac.RPM.expandedRegion", H3K27ac_q_low)

    enh %>% mutate(element_category = case_when(
        H3K27ac_peak_overlap == 1 & H3K27ac.RPM.expandedRegion >= key_high[cell_type] ~ "High H3K27ac",
        H3K27ac_peak_overlap == 1                                                      ~ "H3K27ac",
        H3K27ac.RPM.expandedRegion >= key_high[cell_type]                             ~ "High H3K27ac",
        H3K27ac.RPM.expandedRegion >= key_low[cell_type]                              ~ "H3K27ac",
        CTCF_peak_overlap == 1                                                         ~ "CTCF element",
        H3K27me3_peak_overlap == 1                                                     ~ "H3K27me3 element",
        TRUE                                                                           ~ "No H3K27ac"
    ))
}


## Recommended color palette for plotting the five categories
CATEGORY_COLORS <- c(
    "High H3K27ac"     = "#c5373d",
    "H3K27ac"          = "#d9694a",
    "CTCF element"     = "#49bcbc",
    "H3K27me3 element" = "#429130",
    "No H3K27ac"       = "#c5cad7"
)

## Canonical category order (left to right in horizontal bar plots)
CATEGORY_ORDER <- c("High H3K27ac", "H3K27ac", "No H3K27ac", "CTCF element", "H3K27me3 element")
