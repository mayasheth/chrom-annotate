suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
})

## INPUTS
peak_overlap_assays <- snakemake@params$peak_overlap_assays %>% strsplit(" ") %>% unlist()
rpm_assays <- snakemake@params$rpm_assays %>% strsplit(" ") %>% unlist()
rpm_expanded_assays <- snakemake@params$rpm_expanded_assays %>% strsplit(" ") %>% unlist()
fc_assays <- snakemake@params$fc_assays %>% strsplit(" ") %>% unlist()

peak_overlap_files <- snakemake@input$peak_overlaps %>% strsplit(" ") %>% unlist() %>% setNames(peak_overlap_assays)
rpm_files <- snakemake@input$rpm_original %>% strsplit(" ") %>% unlist() %>% setNames(rpm_assays)
rpm_expanded_files <- snakemake@input$rpm_expanded %>% strsplit(" ") %>% unlist() %>% setNames(rpm_expanded_assays)
fc_files <- snakemake@input$fc_original %>% strsplit(" ") %>% unlist() %>% setNames(fc_assays)

chr_col <- snakemake@params$chr_col
start_col <- snakemake@params$start_col
end_col <- snakemake@params$end_col

out_file <- snakemake@output$enh_plus

## READ IN MAIN ENH LIST
enh_list <- fread(snakemake@input$enh_list, sep = "\t") %>% 
    mutate(elementName = paste0(!!sym(chr_col), ":", !!sym(start_col), "-", !!sym(end_col)))

## ADD PEAK OVERLAPS
for (assay in peak_overlap_assays) {
    this_overlap <- readLines(peak_overlap_files[assay])
    col_name <- paste0(assay, "_peak_overlap")
    if (col_name %in% colnames(enh_list)) {enh_list <- enh_list %>% select(-!!sym(col_name))}

    enh_list <- enh_list %>% mutate(!!sym(col_name) := ifelse(elementName %in% this_overlap, 1, 0))
}

## ADD ORIGINAL (OR RESIZED SMALLER) RPM COUNTS
for (assay in rpm_assays) {
    col_name <- paste0(assay, ".RPM")
    if (col_name %in% colnames(enh_list)) {enh_list <- enh_list %>% select(-!!sym(col_name))}

    this_rpm <- fread(rpm_files[assay], sep = "\t") %>%
        select(elementName = name, all_of(col_name))

    enh_list <- left_join(enh_list, this_rpm, by = "elementName")
}

## ADD ORIGINAL (OR RESIZED SMALLER) FOLD-CHANGE "COUNTS"
for (assay in fc_assays) {
    col_name <- paste0(assay, "_fc.RPM")
    if (col_name %in% colnames(enh_list)) {enh_list <- enh_list %>% select(-!!sym(col_name))}

    this_fc <- fread(fc_files[assay], sep = "\t") %>%
        select(elementName = name, all_of(col_name))

    enh_list <- left_join(enh_list, this_fc, by = "elementName")
}

## ADD EXPANDED RPM COUNTS
for (assay in rpm_expanded_assays) {
    input_col <- paste0(assay, ".RPM")
    output_col <- paste0(assay, ".RPM.expandedRegion")
    if (output_col %in% colnames(enh_list)) {enh_list <- enh_list %>% select(-!!sym(output_col))}

    this_rpm <- fread(rpm_expanded_files[assay], sep = "\t") %>%
        rename(!!sym(output_col) := input_col) %>%
        select(elementName = name, output_col)

    enh_list <- left_join(enh_list, this_rpm, by = "elementName")
}

fwrite(enh_list, out_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
