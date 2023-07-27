library(tidyverse)


funcs_path <- "G:/Shared drives/Denby Lab Team Drive/Lab members/Harry/lettuce_data_shiny_app/plotting_funcs/helper_funcs.R"
source(funcs_path)


divset_pheno <- db_query("SELECT sample_name, sqrt_lesion FROM divset_pheno") %>%
  tidyr::separate(sample_name, into=c('accession','biorep','fungi'),sep='_')



plot_divset_cor <- function(genes,
                        gene_names =NULL,
                        fungi = c("Scl"),
                        single_panel = FALSE,
                        facet_scales = 'free_x',
                        facet_rows = NULL,
                        label_size = 5-(0.3*length(genes)),
                        label_pos = c(0.05,1),
                        return_name_conversion = FALSE){

  ##get gene expression
  # Get the data for the given genes


  genes_exp <- get_divset_exp(genes = genes, fungi = fungi)


  #if no names provided, automatically generate them from annotations file
  if(is.null(gene_names)){
    gene_df <- db_query(paste0("SELECT GeneID,
  CASE
    WHEN (At_ShortName IS NULL OR At_ShortName LIKE 'AT%G%') THEN GeneID
    ELSE 'Ls' || At_ShortName
  END AS name
  FROM gene_annotations WHERE GeneID IN ",
                               paste0("('", paste(genes, collapse = "', '"), "')")))
  }else{ #use supplied names if provided
    gene_df <- data.frame(GeneID = genes, name = gene_names)
  }

  #if duplicate names, add letters at end to distinguish
  while(any(duplicated(gene_df$name))){
    dup_name <- gene_df$name [match(TRUE,duplicated(gene_df$name) )]
    dups <- filter(gene_df, name == dup_name)
    for( i in 1:nrow(dups)){
      row_i <- gene_df$GeneID == dups$GeneID[i]
      gene_df$name[row_i] <- paste0(gene_df$name[row_i], LETTERS[i])
    }
  }

  #print short name to geneID conversion if asked
  if (return_name_conversion){ print(gene_df)}





  # Merge the two data frames
  df2plot <- merge(genes_exp, gene_df, by = 'GeneID')%>%
    merge(divset_pheno, by=c('accession','biorep','fungi'),all.x=T) %>%
    mutate(pathogen = fungi %>% str_replace('Bot','*B. cinerea*') %>%
             str_replace('Scl','*S. sclerotiorum*'))

  if(is.null(facet_rows)){facet_rows <- length(fungi)}
  #make plot
  if(single_panel) {
    y_max <- round(max(df2plot$log2Exp))+(length(genes)*0.175)
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
    geom_point() +
    geom_smooth(method = 'lm', formula = 'y ~ x') +
    ggpubr::stat_cor(aes(label = after_stat(r.label)),
                     method = 'spearman',
                     label.x.npc = label_pos[1],
                     label.y.npc = label_pos[2],
                     size = label_size,
                     show.legend = FALSE)+
    theme_bw()+
    theme(strip.text = ggtext::element_markdown(size=11),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=11.5))

  return(p)
}


hubs <- db_print_head("grn_hubs")
genes<-  hubs$TF


plot_divset_cor(hubs$TF,single_panel = TRUE)
