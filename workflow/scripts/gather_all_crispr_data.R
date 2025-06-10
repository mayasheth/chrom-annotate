suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
})

input_files <- snakemake@input$enh_plus %>% strsplit(" ") %>% unlist()
wtc11_ct <- snakemake@params$WTC11_ct
cell_type_col <- snakemake@params$cell_type_col

res <- lapply(input_files, fread, sep = "\t") %>% 
    rbindlist() %>% 
    as.data.frame() %>% 
    select(-elementName) %>% 
    mutate(!!sym(cell_type_col) := ifelse(!!sym(cell_type_col) == "WTC11", paste0(wtc11_ct, "_for_WTC11"), !!sym(cell_type_col)))

fwrite(res, snakemake@output$crispr_ext, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

