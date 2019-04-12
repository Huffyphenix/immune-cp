############################################
# get mutual exclusive mutation profile between all ICPs.
############################################

library(magrittr)
library(ggplot2)
out_path <- "/home/huff/project/immune_checkpoint/result_20171025"
snv_path <- file.path(out_path, "m_2_snv")


# LOAD DATA 
ICP_highGene_snv <- readr::read_rds(file.path(snv_path,"ICP_highGene_snv.rds.gz")) %>%
  dplyr::select(-highGene_SNV)

gene_list_snv.count_per.filter_hypermutation <-
  readr::read_rds(path = file.path(snv_path, "syn7824274_snv_gene_list_snv_count_per.filter_hypermutation.rds.gz"))

ICP_highGene_snv %>%
  dplyr::inner_join(gene_list_snv.count_per.filter_hypermutation,by="cancer_types") %>%
  dplyr::rename("sm_count"="res") -> ICP_SNV_combine_data

# GET EXCLUSIVE MUTATION PROFILE BETWEEN ICPS IN EACH CANCERS -----------------------
library(showtext)
font_add("Arial","ARIAL.TTF") # sudo进入container，将/home/huff/fonts/中的字体拷贝到/usr/share/fonts/,then do:fc-cache

library(export)

fn_draw_mutation_profile <- function(cancer_types,ICP_SNV,sm_count){
  if(nrow(sm_count)<1 | nrow(ICP_SNV)<1){
    print(paste(cancer_types, "NO data"))
  }else{
    ICP_SNV %>%
      dplyr::select(symbol) %>%
      dplyr::inner_join(sm_count,by = "symbol") %>%
      dplyr::arrange(desc(sm_count)) %>%
      dplyr::mutate(rank = 1:nrow(.)) %>%
      dplyr::filter(sm_count != 0) -> gene_rank
    
    # only TOP 10 mutated ICPs 
    gene_rank %>%
      head(10) -> gene_rank 
    
    per5 <- (ncol(ICP_SNV) - 1)*0.05
    sum_sm <- sum(gene_rank$sm_count)
    if(sum_sm >= per5){
      # fig_name <- paste(cancer_types,"SNV","ALL",sep = "_")
      fig_name <- paste(cancer_types,"SNV","TOP10",sep = "_")
      height <- round(nrow(gene_rank)*0.3)
      ICP_SNV %>%
        tidyr::gather(-symbol, key = "barcode", value="mut") %>%
        dplyr::mutate(mut = ifelse(is.na(mut), 0, mut)) %>%
        dplyr::mutate(mut = as.integer(mut)) -> plot_ready
      
      plot_ready %>%
        dplyr::inner_join(gene_rank, by = "symbol") %>%
        dplyr::select(barcode, mut, rank) -> rank_ready
      
      rank_ready %>%
        dplyr::filter(mut!=0) %>%
        dplyr::group_by(barcode) %>%
        dplyr::mutate(rank_n = min(rank)) %>%
        dplyr::arrange(rank_n) %>%
        dplyr::ungroup() %>%
        dplyr::select(barcode) %>%
        unique() -> sample_rank.mut
      
      rank_ready %>%
        dplyr::filter(!barcode %in% sample_rank.mut$barcode) %>%
        dplyr::select(barcode) %>%
        unique() -> sample_rank.nomut
      
      sample_rank.mut %>%
        rbind(sample_rank.nomut) %>%
        dplyr::pull(barcode) -> sample_rank
      
      plot_ready %>%
        dplyr::mutate(mut = ifelse(mut >= 2,2,mut)) %>%
        ggplot(aes(x = barcode, y = symbol, fill = as.factor(mut))) +
        geom_tile(color = "white",width = 0.9) +
        scale_x_discrete(limits = sample_rank) +
        scale_y_discrete(limits = gene_rank$symbol) +
        # facet_grid(symbol ~ .) +
        scale_fill_manual(
          name = "Mutation_type",
          limits = c("0", "1", "2"),
          label = c("NA", "1", ">=2"),
          values = c( "#D6D6D6", "#FF6A6A","#FF0000")) +
        labs(x = "", y = "", title = cancer_types) +
        theme(
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          
          text = element_text(size = 8),
          title = element_text(size = 10),
          # strip.text.y = element_text(angle = 0,hjust = 0,size = 8),
          # strip.text.x = element_text(size = 8,angle = 90,vjust = 0),
          # strip.background = element_blank(),
          
          legend.title = element_blank(),
          legend.position = "bottom",
          
          panel.background = element_blank()
          # panel.spacing.y  = unit(0, "lines")
        )  -> p;p
      ggsave(file.path(snv_path,"mutation_waterfall_plot",paste(fig_name,"pdf",sep = ".")),device = "pdf",width=6,height=4)
      ggsave(file.path(snv_path,"mutation_waterfall_plot",paste(fig_name,"png",sep = ".")),device = "png",width = 6,height = 4)
      
      # p <- p + theme(text = element_text(family = "Arial"))
      # graph2ppt(x = p,file = file.path(snv_path,"mutation_waterfall_plot",paste(fig_name,"pptx",sep = ".")),width = 6,height = 4)
      print(fig_name)
    }else{
      print(paste(cancer_types,"NO"))
    }
  }
  
}

ICP_SNV_combine_data %>%
  # head() %>%
  purrr::pmap(.f = fn_draw_mutation_profile)
