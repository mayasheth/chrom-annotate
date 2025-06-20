suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
})

input_files <- snakemake@input$enh_plus %>% strsplit(" ") %>% unlist()
cell_types <- snakemake@params$cell_types %>% strsplit(" ") %>% unlist()
output_file <- snakemake@output$elements_final

# only one cell type
if (snakemake@params$cell_type_col == "NONE") {
    res <- fread(input_files[1], sep = "\t")
} else {
    res <- lapply(input_files, fread, sep = "\t") %>% 
        rbindlist() %>% 
        as.data.frame() %>% 
        select(-elementName)
}

fwrite(res, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

