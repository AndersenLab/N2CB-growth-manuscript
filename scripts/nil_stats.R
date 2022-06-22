library(tidyverse)
library(cowplot)
library(ggsignif)


regress_stats <- function(phenodf, genodf, stages) {
  # regress out bleach effect
  regressed <- phenodf %>%
    dplyr::filter(stage %in% c(stages)) %>%
    dplyr::mutate(residual.meanTOF = residuals(aov(mean.TOF ~ bleach, data = .)),
                  residual.meanwidth = residuals(aov(mean.normEXT ~ bleach, data = .)))
  
  ### mean.TOF
  meanTOF.regressed <- regressed %>%
    summarize(comparison = row.names(TukeyHSD(aov(residual.meanTOF ~ strain))$strain),
              padj = TukeyHSD(aov(residual.meanTOF ~ strain))$strain[,4], .groups = "drop") 
  
  # ChrIV
  chrIV <- meanTOF.regressed %>%
    dplyr::filter(comparison %in% c("N2-ECA1064","N2-ECA597","CB4856-ECA575","CB4856-ECA599")) %>%
    tidyr::separate(comparison, into = c("Start", "End"), remove = FALSE) %>%
    dplyr::mutate(sig = dplyr::case_when(padj < 0.0001 ~ "****",
                                         padj < 0.001 ~ "***",
                                         padj < 0.01 ~ "**",
                                         padj < 0.05 ~ "*",
                                         TRUE ~ "ns"))
  # plot pheno
  chrIV.phenoplot <- regressed %>%
    dplyr::filter(strain %in% c("N2","CB4856","ECA1064","ECA597","ECA575","ECA599")) %>%
    ggplot(.) +
    aes(x = strain, y = residual.meanTOF, fill = strain) +
    geom_boxplot(outlier.color = NA, width = 0.6, alpha = 0.8, size = 0.2) +
    geom_jitter(width = 0.1, size = 0.2, alpha = 0.6) + 
    theme_cowplot() +
    scale_fill_manual(values = c("N2" = "#ffa500", "CB4856" = "#0000ff", "ECA575" = "gray","ECA1064" = "gray")) +
    facet_wrap(~stage, scales = "free") + 
    panel_border(color = "black", size = 0.8) +
    ggsignif::geom_signif(data = chrIV, manual = TRUE, inherit.aes = F,
                          aes(xmax = End, xmin = Start, y_position= c(82,73,82,73), annotations = sig)) +
    coord_flip() + facet_wrap(~"Phenotype") +
    scale_color_manual(values = c("N2" = "#ffa500", "CB4856" = "#0000ff")) + 
    guides(color = "none", size = "none", fill = "none") + ylab("Mean Animal Length") + xlab("") 
  
  # plot chrIV genotype
  chrIV.genoplot <- genodf %>%
    dplyr::ungroup() %>%
    dplyr::filter(chrom == "IV",
                  sample %in% c("N2","ECA1064","ECA597","CB4856","ECA575","ECA599"),
                  !end == 18000000) %>%
    dplyr::mutate(chrom = paste0("chr", "IV"),
                  sample = factor(sample, levels = rev(c("N2","ECA1064","ECA597","CB4856","ECA575","ECA599")))) %>%
    dplyr::distinct() %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name), size = 8) +
    facet_wrap(~chrom)+
    geom_vline(xintercept = 9.392639) +
    geom_vline(xintercept = c(6.211685,12.868784), linetype="dotted") +
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_cowplot() + panel_border(color = "black") +
    theme(axis.ticks.y = element_blank(),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    labs(x = "Genomic position (Mb)", y = "") 
  
  # plot genome background
  chrIV.backplot <- genodf %>%
    dplyr::ungroup() %>%
    dplyr::filter(sample %in% c("N2","ECA1064","ECA597","CB4856","ECA575","ECA599"),
                  chrom == ifelse("IV" == "III", "II", "III")) %>%
    dplyr::mutate(bp = end - start) %>%
    dplyr::group_by(sample, gt_name) %>%
    dplyr::summarize(max = sum(bp)) %>%
    dplyr::arrange(desc(max)) %>%
    dplyr::filter(max == first(max)) %>%
    dplyr::mutate(start = 0,
                  end = 15e6,
                  chr = "Genome")  %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sample = factor(sample, levels = rev(c("N2","ECA1064","ECA597","CB4856","ECA575","ECA599")))) %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name), size = 8)+
    facet_wrap(~chr)+
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_cowplot() + panel_border(color = "black") +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          strip.text.y = element_blank())+
    labs(x = "Genomic position (Mb)", y = "")
  
  allIV <- cowplot::plot_grid(chrIV.genoplot, chrIV.backplot, chrIV.phenoplot,
                              nrow = 1, ncol = 3, rel_widths = c(3, 1, 4.5), align = "h", axis = "b", labels = c("A", "", "B"))
  
  ### mean.normEXT
  meanwidth.regressed <- regressed %>%
    summarize(comparison = row.names(TukeyHSD(aov(residual.meanwidth ~ strain))$strain),
              padj = TukeyHSD(aov(residual.meanwidth ~ strain))$strain[,4], .groups = "drop")
  
  # ChrV
  chrV <- meanwidth.regressed %>%
    dplyr::filter(comparison %in% c("N2-ECA2006","N2-ECA232","CB4856-ECA1060","CB4856-ECA1058")) %>%
    tidyr::separate(comparison, into = c("Start", "End"), remove = FALSE) %>%
    dplyr::mutate(sig = dplyr::case_when(padj < 0.0001 ~ "****",
                                         padj < 0.001 ~ "***",
                                         padj < 0.01 ~ "**",
                                         padj < 0.05 ~ "*",
                                         TRUE ~ "ns"))
  chrV.phenoplot <- regressed %>%
    dplyr::filter(strain %in% c("N2","CB4856","ECA2006","ECA232","ECA1060","ECA1058")) %>%
    ggplot(.) +
    aes(x = strain, y = residual.meanwidth, fill = strain) +
    geom_boxplot(outlier.color = NA, width = 0.6, alpha = 0.8, size = 0.2) +
    geom_jitter(width = 0.1, size = 0.2, alpha = 0.6) + 
    theme_cowplot() +
    scale_fill_manual(values = c("N2" = "#ffa500", "CB4856" = "#0000ff", "ECA1060" = "gray","ECA2006" = "gray")) +
    #facet_wrap(~stage, scales = "free") + 
    panel_border(color = "black", size = 0.8) +
    ggsignif::geom_signif(data = chrV, manual = TRUE, inherit.aes = F,
                          aes(xmax = End, xmin = Start, y_position= c(8,7,8,7), annotations = sig)) +
    coord_flip() + facet_wrap(~"Phenotype") +
    scale_color_manual(values = c("N2" = "#ffa500", "CB4856" = "#0000ff")) + 
    guides(color = "none", size = "none", fill = "none") + ylab("Mean Animal Width") + xlab("") 
  
  # plot chrV genotype
  chrV.genoplot <- genodf %>%
    dplyr::ungroup() %>%
    dplyr::filter(chrom == "V",
                  sample %in% c("N2","CB4856","ECA2006","ECA232","ECA1060","ECA1058"),
                  !end == 18000000) %>%
    dplyr::mutate(chrom = paste0("chr", "V"),
                  sample = factor(sample, levels = rev(c("N2","ECA2006","ECA232","CB4856","ECA1060","ECA1058")))) %>%
    dplyr::distinct() %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name), size = 8) +
    facet_wrap(~chrom)+
    geom_vline(xintercept = 11.806498) +
    geom_vline(xintercept = c(5.371124,12.112105), linetype="dotted") +
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_cowplot() + panel_border(color = "black") +
    theme(axis.ticks.y = element_blank(),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    labs(x = "Genomic position (Mb)", y = "") 
  
  # plot genome background
  chrV.backplot <- genodf %>%
    dplyr::ungroup() %>%
    dplyr::filter(sample %in% c("N2","ECA2006","ECA232","CB4856","ECA1060","ECA1058"),
                  chrom == ifelse("V" == "III", "II", "III")) %>%
    dplyr::mutate(bp = end - start) %>%
    dplyr::group_by(sample, gt_name) %>%
    dplyr::summarize(max = sum(bp)) %>%
    dplyr::arrange(desc(max)) %>%
    dplyr::filter(max == first(max)) %>%
    dplyr::mutate(start = 0,
                  end = 15e6,
                  chr = "Genome")  %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sample = factor(sample, levels = rev(c("N2","ECA2006","ECA232","CB4856","ECA1060","ECA1058")))) %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name), size = 8)+
    facet_wrap(~chr)+
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_cowplot() + panel_border(color = "black") +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          strip.text.y = element_blank())+
    labs(x = "Genomic position (Mb)", y = "")
  
  allV <- cowplot::plot_grid(chrV.genoplot, chrV.backplot, chrV.phenoplot,
                             nrow = 1, ncol = 3, rel_widths = c(3, 1, 4.5), 
                             align = "h", axis = "b", labels = c("A", "", "B"))
  
  
  # ChrX
  chrX <- meanwidth.regressed %>%
    dplyr::filter(comparison %in% c("N2-ECA929","CB4856-ECA828")) %>%
    tidyr::separate(comparison, into = c("Start", "End"), remove = FALSE) %>%
    dplyr::mutate(sig = dplyr::case_when(padj < 0.0001 ~ "****",
                                         padj < 0.001 ~ "***",
                                         padj < 0.01 ~ "**",
                                         padj < 0.05 ~ "*",
                                         TRUE ~ "ns"))
  
  chrX.phenoplot <- regressed %>%
    dplyr::filter(strain %in% c("N2","ECA929","CB4856","ECA828")) %>%
    ggplot(.) +
    aes(x = strain, y = residual.meanwidth, fill = strain) +
    geom_boxplot(outlier.color = NA, width = 0.6, alpha = 0.8, size = 0.2) +
    geom_jitter(width = 0.1, size = 0.2, alpha = 0.6) + 
    theme_cowplot() +
    scale_fill_manual(values = c("N2" = "#ffa500", "CB4856" = "#0000ff", "ECA929" = "gray","ECA828" = "gray")) +
    #facet_wrap(~stage, scales = "free") + 
    panel_border(color = "black", size = 0.8) +
    ggsignif::geom_signif(data = chrX, manual = TRUE, inherit.aes = F,
                          aes(xmax = End, xmin = Start, y_position= c(6.5,6.5), annotations = sig)) +
    coord_flip() + facet_wrap(~"Phenotype") +
    scale_color_manual(values = c("N2" = "#ffa500", "CB4856" = "#0000ff")) + 
    guides(color = "none", size = "none", fill = "none") + ylab("Mean Animal Width") + xlab("") 
  
  # plot chrX genotype
  chrX.genoplot <- genodf %>%
    dplyr::ungroup() %>%
    dplyr::filter(chrom == "X",
                  sample %in% c("N2","ECA929","CB4856","ECA828"),
                  !end == 18000000) %>%
    dplyr::mutate(chrom = paste0("chr", "X"),
                  sample = factor(sample, levels = rev(c("N2","ECA929","CB4856","ECA828")))) %>%
    dplyr::distinct() %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name), size = 8) +
    facet_wrap(~chrom)+
    geom_vline(xintercept = 12.750794) +
    geom_vline(xintercept = c(12.565734,13.173080), linetype="dotted") +
    scale_color_manual(values = c("N2" = "#ffa500", "CB4856" = "#0000ff")) + 
    theme_cowplot() + panel_border(color = "black") +
    theme(axis.ticks.y = element_blank(),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    labs(x = "Genomic position (Mb)", y = "") 
  
  # plot genome background
  chrX.backplot <- genodf %>%
    dplyr::ungroup() %>%
    dplyr::filter(sample %in% c("N2","ECA929","CB4856","ECA828"),
                  chrom == ifelse("X" == "III", "II", "III")) %>%
    dplyr::mutate(bp = end - start) %>%
    dplyr::group_by(sample, gt_name) %>%
    dplyr::summarize(max = sum(bp)) %>%
    dplyr::arrange(desc(max)) %>%
    dplyr::filter(max == first(max)) %>%
    dplyr::mutate(start = 0,
                  end = 15e6,
                  chr = "Genome")  %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sample = factor(sample, levels = rev(c("N2","ECA929","CB4856","ECA828")))) %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name), size = 8)+
    facet_wrap(~chr) +
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_cowplot() + panel_border(color = "black") +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          strip.text.y = element_blank())+
    labs(x = "Genomic position (Mb)", y = "")
  
  allX <- cowplot::plot_grid(chrX.genoplot, chrX.backplot, chrX.phenoplot,
                             nrow = 1, ncol = 3, rel_widths = c(3, 1, 4.5), 
                             align = "h", axis = "b", labels = c("A", "", "B"))
  
  return(list(allIV, allV, allX))
}

# load("processed/nilgeno.Rda")
# L4s <- regress_stats(sumplate, nilgeno, "L4s")















