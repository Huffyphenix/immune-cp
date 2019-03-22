fn_survival <- function(data,title,color,sur_name,result_path,h,w,lx=0.8,ly=0.6){
  library(survival)
  library(survminer)
  fit <- survfit(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)
  diff <- survdiff(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)
  kmp <- 1 - pchisq(diff$chisq, df = length(levels(as.factor(data$group))) - 1)
  # legend <- data.frame(group=paste("C",sort(unique(data$group)),sep=""),n=fit$n)
  # legend %>%
  #   dplyr::mutate(
  #     label = purrr::map2(
  #       .x = group,
  #       .y = n,
  #       .f = function(.x,.y){
  #         latex2exp::TeX(glue::glue("<<.x>>, n = <<.y>>", .open = "<<", .close = ">>"))
  #       }
  #     )
  #   ) -> legend
  tibble::tibble(group = group, color = color) %>%
    dplyr::inner_join(data,by="group") %>%
    dplyr::mutate(group = paste(group,sep="")) %>%
    dplyr::select(group,color) %>%
    unique() -> color_paired
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        data = data,
                        surv.median.line = "hv",
                        title = paste(title,", p =", signif(kmp, 2)), # change it when doing diff data
                        xlab = "Survival in days",
                        ylab = 'Probability of survival',
                        # legend.title = "Methyla group:",
                        legend= c(lx,ly),
                        # ggtheme = theme_survminer(),
                        ggtheme = theme(
                          panel.border = element_blank(), panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", 
                                                                                       size = 0.5), 
                          panel.background = element_rect(fill = "white"),
                          legend.key = element_blank(),
                          legend.background = element_blank(),
                          legend.text = element_text(size = 8, colour = "black"),
                          axis.text = element_text(size = 12, colour = "black"),
                          legend.title = element_blank(),
                          axis.title = element_text(size = 12,color = "black")
                        )
   ) +
    scale_color_manual(
      values = color_paired$color,
      labels = color_paired$group
    )
  ggsave(filename =paste(sur_name,signif(kmp, 2),"png",sep="."), path = result_path,device = "png",height = h,width = w)
  ggsave(filename =paste(sur_name,signif(kmp, 2),"pdf",sep="."), path = result_path,device = "pdf",height = h,width = w)
  
}

fn_mutation_burden <- function(data,group,facet="~ cancer_types",value,color,xlab,comp_list,m_name,result_path,w=7,h=10){
  data %>%
    ggpubr::ggboxplot(x = group, y = value,fill="white",alpha = 0,width = 0.1,
                      color = group #add = "jitter",#, palette = "npg"
    ) +
    geom_violin(aes(fill = group),alpha = 0.5) +
    # geom_jitter(aes(color = group),alpha = 0.2,width = 0.1) +
    # geom_point(aes(x=as.numeric(group)+b,y=sm_count,color=group),alpha = 0.5) +
    scale_x_discrete(#breaks = c(1:10),
                     labels = unique(data$group)
                     # expand = c(0.2,0.2,0.2)
    ) +
    facet_wrap(as.formula(facet), strip.position = "bottom", scales = "free") +
    scale_color_manual(
      values = color
    )+
    # ylim(4,12) +
    ylab(xlab) +
    xlab("Group") +
    theme(legend.position = "none",
          # axis.title.x = element_blank(),
          strip.background = element_rect(fill = "white",colour = "white"),
          text = element_text(size = 10, colour = "black"),
          strip.text = element_text(size = 8)) +
    # ggpubr::stat_compare_means(label.y = 14,paired = TRUE) +
    ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif")
  
  ggsave(filename =paste(m_name,"png",sep="."), path = result_path,device = "png",width = w,height = h)
  ggsave(filename =paste(m_name,"pdf",sep="."), path = result_path,device = "pdf",width = w,height = h)
}

fn_mutation_burden_all <- function(data,group,value,color,xlab,comp_list,m_a_name,result_path,w=4,h=3){
  data %>%
    ggpubr::ggboxplot(x = group, y = value,fill="white",alpha = 0,width = 0.1,
                      color = group #add = "jitter",#, palette = "npg"
    ) +
    geom_violin(aes(fill = group),alpha = 0.5) +
    geom_jitter(aes(color = group),alpha = 0.2,width = 0.1,size=0.1) +
    geom_boxplot(fill="white",alpha = 0,width = 0.1) +
    scale_x_discrete(#breaks = c(1:10),
                     labels = unique(data$group)
                     # expand = c(0.2,0.2,0.2)
    ) +
    # facet_wrap(~ cancer_types, strip.position = "bottom", scales = "free") +
    scale_color_manual(
      values = color
    )+
    scale_fill_manual(
      values = color
    ) +
    # ylim(4,12) +
    ylab(xlab) +
    xlab("Group") +
    theme(legend.position = "none",
          # axis.title.x = element_blank(),
          strip.background = element_rect(fill = "white",colour = "white"),
          text = element_text(size = 12, colour = "black"),
          strip.text = element_text(size = 12))
    # ggpubr::stat_compare_means(label.y = 14,paired = TRUE) +
    # ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif")
  ggsave(filename =paste(m_a_name,"png",sep="."), path = result_path,device = "png",width = w,height = h)
  ggsave(filename =paste(m_a_name,"pdf",sep="."), path = result_path,device = "pdf",width = w,height = h)
}
