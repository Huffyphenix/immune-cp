fn_survival <- function(data,title,color){
  fit <- survfit(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)
  diff <- survdiff(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)
  kmp <- 1 - pchisq(diff$chisq, df = length(levels(as.factor(data$group))) - 1)
  legend <- data.frame(group=paste("C",unique(data$group),sep=""),n=fit$n)
  legend %>%
    dplyr::mutate(
      label = purrr::map2(
        .x = group,
        .y = n,
        .f = function(.x,.y){
          latex2exp::TeX(glue::glue("<<.x>>, n = <<.y>>", .open = "<<", .close = ">>"))
        }
      )
    ) -> legend
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        surv.median.line = "hv",
                        title = paste(title,", p =", signif(kmp, 2)), # change it when doing diff data
                        xlab = "Survival in days",
                        ylab = 'Probability of survival',
                        # legend.title = "Methyla group:",
                        legend= c(0.8,0.6),
                        # ggtheme = theme_survminer(),
                        ggtheme = theme(
                          panel.border = element_blank(), panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", 
                                                                                       size = 0.5), 
                          panel.background = element_rect(fill = "white"),
                          legend.key = element_blank(),
                          legend.background = element_blank(),
                          legend.text = element_text(size = 12, colour = "black"),
                          axis.text = element_text(size = 12, colour = "black"),
                          legend.title = element_blank(),
                          axis.title = element_text(size = 12,color = "black")
                        )
  ) +
    scale_color_manual(
      values = color,
      labels = legend$label
    ) 
}

fn_mutation_burden <- function(data,color,comp_list){
  data %>%
    dplyr::mutate(sm_count = log2(sm_count)) %>%
    ggpubr::ggboxplot(x = "group", y = "sm_count",
                      color = "group" #add = "jitter",#, palette = "npg"
    ) +
    # geom_point(aes(x=as.numeric(group)+b,y=sm_count,color=group),alpha = 0.5) +
    scale_x_discrete(breaks = c(1:10),
                     labels = paste("C",c(1:10),sep="")
                     # expand = c(0.2,0.2,0.2)
    ) +
    facet_wrap(~ cancer_types, strip.position = "bottom", scales = "free") +
    scale_color_manual(
      values = color
    )+
    # ylim(4,12) +
    ylab("log2(Mutation burden)") +
    xlab("Clusters") +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          strip.background = element_rect(fill = "white",colour = "white"),
          strip.text = element_text(size = 12)) +
    # ggpubr::stat_compare_means(label.y = 14,paired = TRUE) +
    ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif")
}

fn_mutation_burden_all <- function(data,color,comp_list){
  data %>%
    dplyr::mutate(sm_count = log2(sm_count)) %>%
    ggpubr::ggboxplot(x = "group", y = "sm_count",
                      color = "group" #add = "jitter",#, palette = "npg"
    ) +
    # geom_point(aes(x=as.numeric(group)+b,y=sm_count,color=group),alpha = 0.5) +
    scale_x_discrete(breaks = c(1:10),
                     labels = paste("C",c(1:10),sep="")
                     # expand = c(0.2,0.2,0.2)
    ) +
    # facet_wrap(~ cancer_types, strip.position = "bottom", scales = "free") +
    scale_color_manual(
      values = color
    )+
    # ylim(4,12) +
    ylab("log2(Mutation burden)") +
    xlab("Clusters") +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          strip.background = element_rect(fill = "white",colour = "white"),
          strip.text = element_text(size = 12)) +
    # ggpubr::stat_compare_means(label.y = 14,paired = TRUE) +
    ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif")
}