library(googlesheets)
library(readr)
library(magrittr)
rppa_name_ss <- gs_title(x = "RPPA_Standard_Ab_List_Updated.xlsx")
rppa_list <- gs_read(ss = rppa_name_ss, ws = "Previous_Standard Ab List_218")
#rppa form MDAnderson
rppa_list %>% 
  dplyr::select(
    pre_protein = `Ab Name Reported on Dataset`,
    symbol = `Gene Name`,
    species = Species,
    status = `Validation Status*`
  ) %>% 
  dplyr::select(protein = pre_protein, symbol) %>% 
  dplyr::mutate(name = stringr::str_replace(protein, pattern = "-\\w-\\w$", "")) %>% 
  dplyr::mutate(name = stringr::str_to_upper(name)) %>% 
  dplyr::mutate(name = stringr::str_replace_all(name, pattern = "-|_", ""))-> mda


# tcpa portal
library(magrittr)
rppa_json_url <- "http://tcpaportal.org/tcpa/_design/basic/_show/annotation-antibody_list/annotation-antibody"

jsonlite::fromJSON(txt = rppa_json_url) %>% 
  dplyr::as_tibble() %>% 
  dplyr::select(`0`:`5`) %>% 
  dplyr::rename(
    protein = `0`,
    symbol = `1`,
    status = `2`,
    origin = `3`,
    source = `4`,
    catalog = `5`
  ) %>% 
  dplyr::mutate_all(.funs = dplyr::funs(stringr::str_trim)) %>% 
  dplyr::mutate(
    origin = stringr::str_sub(origin, end = 1),
    status = dplyr::recode(
      status,
      "Validated as ELISA" = "V",
      "Use with Caution" = "C",
      .default = "E"
    ),
    protein = stringr::str_replace(protein, pattern = "^x", "")
  ) %>% 
  dplyr::mutate(
    protein = stringr::str_c(protein, origin, status, sep = "-")
  ) %>% 
  dplyr::select(protein, symbol) %>% 
  dplyr::mutate(name = stringr::str_replace(protein, pattern = "-\\w-\\w$", "")) %>% 
  dplyr::mutate(name = stringr::str_to_upper(name)) %>% 
  dplyr::mutate(name = stringr::str_replace_all(name, pattern = "-|_", "")) -> tcpa_name

rds_url <- "https://github.com/chunjie-sam-liu/R_Leng_1/raw/master/autophagy_codes/01_expr/rppa_protein_name.rds"
rppa_name <- 
  tibble::tibble(protein = read_rds(path = url(rds_url))) %>% 
  dplyr::mutate(name = stringr::str_replace(protein, pattern = "-\\w-\\w$", "")) %>% 
  dplyr::mutate(name = stringr::str_to_upper(name)) %>% 
  dplyr::mutate(name = stringr::str_replace_all(name, pattern = "-|_", ""))

tibble::tribble(
  ~symbol, ~name,
  "EIF4EBP1", "4EBP1PT70",
  "PRKAA1", "AMPKPT172",
  "MAPK8", "JNKPT183PY185",
  "RPS6KB1", "P70S6K"
) -> sup

dplyr::bind_rows(tcpa_name, mda, sup) %>% 
  dplyr::select(symbol, name) %>% 
  dplyr::distinct(name, .keep_all = T) -> merged_name_symbol

rppa_name %>% 
  dplyr::left_join(merged_name_symbol, by = "name") %>% 
  dplyr::select(protein, symbol) %>% 
  dplyr::mutate(symbol = stringr::str_replace_all(symbol, " |,", ", ")) %>% print(n = Inf) %>% 
  readr::write_rds("rppa_name_symbol.rds.gz", compress = "gz")


