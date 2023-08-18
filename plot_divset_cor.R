library(tidyverse)


if(!'helperFuncsLoaded' %in% ls()){
  funcs_path <- "G:/Shared drives/Denby Lab Team Drive/Lab members/Harry/lettuce_data_shiny_app/plotting_funcs/helper_funcs.R"
  source(funcs_path)
}



divset_pheno <- db_query("SELECT sample_name, sqrt_lesion FROM divset_pheno") %>%
  tidyr::separate(sample_name, into=c('accession','biorep','fungi'),sep='_',remove = FALSE) %>%
  arrange(fungi, -sqrt_lesion)




plot_divset_cor <- function(GeneIDs = NULL,
                            At_orthologs = NULL,
                            GO_id = NULL,
                            protein_domain = NULL,
                            top_n_by_lesion_cor= NULL,
                            signif_cor_only = FALSE,
                            signif_cor_DEGs_filtered_only = FALSE,
                            fungi = c("Scl"),
                            single_panel = TRUE,
                            facet_scales = 'free_x',
                            facet_rows = NULL,
                            facet_title_size=11,
                            label_pos = c(0.05,0.9),
                            return_heatmap = TRUE){



  genes <- get_lettuce_genes_from_inputs(GeneIDs, At_orthologs,
                                         GO_id, protein_domain)

  if(length(genes)<1){return(ggplot())}





  ##if 8 or under genes, plot a scatterplot
  if(return_heatmap){
    require(pheatmap)

    genes_exp <- get_divset_exp(genes = genes,
                                fungi = fungi,
                                top_n_by_lesion_cor= top_n_by_lesion_cor,
                                signif_cor_only = signif_cor_only,
                                signif_cor_DEGs_filtered_only = signif_cor_DEGs_filtered_only,
                                return_long = FALSE)

    print(unique(genes_exp$GeneID))

    genes_exp <- genes_exp %>%
      merge(get_gene_names(genes_exp$GeneID),by='GeneID') %>%
      column_to_rownames('name')







    sample_annot_raw <- divset_pheno[divset_pheno$sample_name %in% colnames(genes_exp),c(1,4,5)]
    sample_annot <- sample_annot_raw %>%
      remove_rownames() %>%
      column_to_rownames('sample_name')

    if(length(fungi)>1){
      bot_scaled_lesion <- scale(sample_annot$sqrt_lesion[sample_annot$fungi=='Bot'])
      scl_scaled_lesion <- scale(sample_annot$sqrt_lesion[sample_annot$fungi=='Scl'])
      sample_annot$scaled_lesion <- c(bot_scaled_lesion,scl_scaled_lesion)
      sample_annot <- select(sample_annot, fungi, scaled_lesion) %>%
        arrange(-scaled_lesion)
    }else{
      sample_annot <- select(sample_annot, -fungi)
    }



    heatmap <- pheatmap(select(genes_exp,all_of(rownames(sample_annot))),
                        cluster_cols = F,
                        cluster_rows = T,
                        fontsize_row=facet_title_size*0.725,
                        color=colorRampPalette(c("navy", "white", "red"))(50),
                        scale = 'row',
                        annotation_col = sample_annot,
                        show_colnames = F,
                        border = NA)


    return(list(plot=heatmap,
                exp_data = genes_exp,
                lesion_data = divset_pheno))

  }
    # Merge the two data frames
    genes_exp <- get_divset_exp(genes = genes,
                                fungi = fungi,
                                top_n_by_lesion_cor= top_n_by_lesion_cor,
                                signif_cor_only = signif_cor_only,
                                signif_cor_DEGs_filtered_only = signif_cor_DEGs_filtered_only) %>%
      merge(get_gene_names(genes_exp$GeneID), by = 'GeneID')


    n_genes <- length(unique(genes_exp$GeneID))


    df2plot <- genes_exp %>%
      merge(divset_pheno, by=c('accession','biorep','fungi'),all.x=T) %>%
      mutate(pathogen = fungi %>% str_replace('Bot','B. cinerea') %>%
               str_replace('Scl','S. sclerotiorum'))

    if(is.null(facet_rows)){facet_rows <- length(fungi)}
    #make plot
    if(single_panel) {
      y_max <- round(max(df2plot$log2Exp))+(n_genes*0.14)
      print(head(df2plot))
      p <- ggplot(df2plot,
                  aes(x = sqrt_lesion,
                      y = log2Exp,
                      color = name)) +
        scale_y_continuous(limits = c(0,y_max),breaks = seq(0,y_max,2))+
        scale_color_discrete(name = "Gene")+
        facet_wrap( ~ pathogen, scale = facet_scales, nrow = 1)
    } else {
      if (is.null(facet_rows)) { facet_rows <- length(fungi) }
      p <- ggplot(df2plot, aes(x = sqrt_lesion, y = log2Exp)) +
        facet_wrap(pathogen ~ name, scale = facet_scales, nrow = facet_rows)
    }

    p <- p+
      labs(x = 'Sqrt Lesion Size (mm)', y = 'log2 expression') +
      geom_point(alpha=0.5) +
      geom_smooth(method = 'lm', formula = 'y ~ x') +
      ggpubr::stat_cor(aes(label = after_stat(r.label)),
                       method = 'spearman',
                       label.x.npc = label_pos[1],
                       label.y.npc = label_pos[2],
                       size = ifelse(!single_panel, 5.2, 5-(0.3*n_genes)),
                       show.legend = FALSE)+
      theme_bw()+
      theme(strip.text = ggtext::element_textbox( size = facet_title_size,
                                                  face='italic',
                                                  color = "black",
                                                  fill = "grey92",
                                                  box.color = "grey92",
                                                  halign = 0.5,
                                                  linetype = 1,
                                                  r = unit(5, "pt"),
                                                  width = unit(1, "npc"),
                                                  padding = margin(0.4, 0, 1, 0),
                                                  margin = margin(0.2,0.2,0.2,0.2)),
            strip.background = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size=11.5))


    print('plotted diversity set lesion correlation')

    return(list(plot=p,
                exp_data = genes_exp,
                lesion_data = divset_pheno))
}

# # GOid <- 'GO:0004672'
# plot_divset_cor(GO_id = 'GO:0004672', top_n_by_lesion_cor = 5,
#                 fungi = 'sclero', return_heatmap = FALSE,single_panel = FALSE)
