# by atom process symbols
library(biomaRt)
library(tidyverse)
GENES = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
GENES.ATTRIBUTES <- listAttributes(GENES)
GENES.FILTERS <- listFilters(GENES)

symbols <- c("ATG13", "ATG101", "ATG14", "DEPTOR", "VMP1", "RAB1C", "CTSL", "CTSV",
    "LAMTOR1", "LAMTOR2")
# RAB1C is a psudogene, need to be excluded.
symbols_entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
    filters = "hgnc_symbol", values = symbols, mart = GENES) %>% drop_na() %>% as_tibble() %>%
    mutate(entrez_id = as.character(entrezgene))

symbols_entrez
te <- expr %>% .[1, ] %>% unnest() %>% filter(entrez_id %in% symbols_entrez$entrez_id) %>%
    .[, c(1:3)] %>% left_join(symbols_entrez, by = "entrez_id")
te
change_symbol <- te$hgnc_symbol
names(change_symbol) <- te$symbol
change_symbol
tm <- te[, c(1:3)] %>% mutate(symbol = ifelse(symbol %in% names(change_symbol), plyr::revalue(symbol,
    replace = change_symbol), symbol))
change_name <- function(.x) {
    .x %>% mutate(symbol = ifelse(symbol %in% names(change_symbol), plyr::revalue(symbol,
        replace = change_symbol), symbol))
}
te <- expr %>% mutate(expr = map(expr, change_name))
te[1, ] %>% unnest() %>% filter(symbol %in% change_symbol) %>% .[, c(1:3)]
expr <- te
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr %>% readr::write_rds(file.path(tcga_path, "pancan_expr_20160513.rds.gz"), compress = "gz")

# save clinical data




















