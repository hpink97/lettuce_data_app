# Load the tidyverse library, which includes ggplot2 and dplyr among other packages
library(tidyverse)
library(RColorBrewer)
library(ggtext)
library(ComplexHeatmap)
library(circlize)


# Load helper functions if they are not already loaded
if(!'helperFuncsLoaded' %in% ls()){
  funcs_path <- "G:/Shared drives/Denby Lab Team Drive/Lab members/Harry/lettuce_data_shiny_app/plotting_funcs/helper_funcs.R"
  source(funcs_path)
}

z_score_clip <- function(x) {
  z <- (x - mean(x)) / sd(x)
  z[z > 1.999] <- 1.999
  z[z < -1.999] <- -1.999
  return(z)
}

# Function to plot time-series expression data for specific genes
plot_timeseries_expr <- function(GeneIDs = NULL,        # Vector of gene IDs
                                 At_orthologs = NULL,    # Vector of Arabidopsis orthologs
                                 GO_id = NULL,           # GO term ID
                                 protein_domain = NULL,  # Protein domain
                                 max_n = 1e6,            # Max number of genes to plot
                                 fungi = c("Bot", "Scl"),# Fungi names
                                 include_mock = TRUE,    # Whether to include mock treatment
                                 overlap_DEGs_only = FALSE, # Plot only overlapping DEGs
                                 up_down_only = 'both',  # 'up', 'down', or 'both' - for expression direction
                                 plot_type = 'line',
                                 strip.text.x_size = 9,  # Font size for facet labels
                                 x_axis_interval = 6,    # Interval for x-axis breaks
                                 single_panel = FALSE,
                                 facet_scales = 'free_y', # Scale type for facets
                                 facet_nrow = 2,         # Number of rows for facets
                                 return_name_conversion = FALSE){ # Return name conversion if needed


  ##check inputs are valid
  if(!(plot_type %in% c('heatmap','line'))){ return (ggplot())}
  if(!(up_down_only %in% c('both','up','down'))){ return (ggplot())}



  # Get gene IDs based on input parameters
  input_genes <- get_lettuce_genes_from_inputs(GeneIDs, At_orthologs, GO_id, protein_domain)

  # Return an empty plot if no genes are found
  if(length(input_genes) < 1){ return(ggplot())}

  fetch_mock_data = ((include_mock &(plot_type =='heatmap' |!single_panel)) | length(input_genes)>max_n)

  # Get time-series expression data for the selected genes
  TS_exp <- get_timeseries_exp(genes =input_genes,
                               fungi = fungi,
                               ##include mock data if it will be shown on plot
                               ##OR is required for log2FC filtering of genes
                               include_mock =  fetch_mock_data,
                               ts_DEGs_only = overlap_DEGs_only,
                               up_down_only = up_down_only,
                               include_raw_points = plot_type=='line',
                               include_conf_intervals = plot_type=='line'

  )

  # Get the number of unique genes in the dataset
  n_genes <- length(unique(TS_exp$GeneID))

  # If the number of genes exceeds max_n, filter to top max_n genes based on absolute log2 fold change
  if(n_genes > max_n){
    log2FC <- TS_exp %>%
      select(GeneID, Fungi, Treatment, mean_expr,hpi) %>%
      pivot_wider(names_from = 'Treatment',values_from = 'mean_expr') %>%
      mutate(log2FC = Infected-Mock) %>%
      group_by(GeneID) %>%
      summarise(mean_Abslog2FC = mean(abs(log2FC),na.rm = T),
                max_Abslog2FC = max(abs(log2FC),na.rm = T)) %>%
      arrange(-max_Abslog2FC)
    TS_exp <- filter(TS_exp, GeneID %in% log2FC$GeneID[1:max_n])
  }

  # Merge expression data with gene names and modify Fungi names
  TS_exp_w_names <- merge(TS_exp, get_gene_names(TS_exp$GeneID), by='GeneID',all.x=TRUE)

  # Set facet rows based on the number of genes
  facet_nrow = ifelse(facet_nrow >= n_genes, n_genes, facet_nrow)

  # Generate the base plot according to whether mock treatment is included

  if(plot_type=='heatmap'){
    exp_long<- TS_exp_w_names %>%
      mutate(hpi = str_pad(hpi, width = 2, pad = "0"),
             Treatment = factor(Treatment, levels=c('Mock','Infected'))) %>%
      arrange(Fungi,Treatment, hpi)

      exp_wide <- pivot_wider(exp_long,
                              names_from = c('Fungi','Treatment','hpi'),values_from = 'mean_expr') %>%
        column_to_rownames('name')

      exp_matrix <- exp_wide %>%
        dplyr::select(-GeneID) %>%
        as.matrix()



      z_score_scaled_data <- t(apply(exp_matrix, 1, z_score_clip))


      sample_annots <- data.frame(id = colnames(exp_matrix)) %>%
        separate(id, into = c('Fungi','Treatment','hpi'),sep='_',remove = FALSE) %>%
        mutate(hpi = as.numeric(hpi),
               Fungi = Fungi %>%
                 str_replace('Bot','Bc')%>%
                 str_replace('Scl','Ss'),
               Condition = paste0(Fungi, '-',Treatment)) %>%
        dplyr::select(-Fungi, -Treatment) %>%
        column_to_rownames('id')

      hpi_palette <- colorRamp2(seq(9, max(sample_annots$hpi), by = 3),
                                colorRampPalette(c("grey75", "black"))(length(seq(9,
                                                                                   max(sample_annots$hpi),
                                                                                   by = 3))))

      if(any(grepl('Mock', sample_annots$Condition))){
        cond_levels <- c()
        if(any(allowed_bot %in% tolower(fungi))){
          cond_levels <- c(cond_levels, c('Bc-Mock','Bc-Infected'))
        }
        if(any(allowed_sclero %in% tolower(fungi))){
          cond_levels <- c(cond_levels, c('Ss-Mock','Ss-Infected'))
        }

        sample_annots$Condition <- factor(sample_annots$Condition, levels = cond_levels)
      }

      n_cond <- length(unique(sample_annots$Condition))
      cond_palette <-brewer.pal(6, "Set1")[1:n_cond]
      names(cond_palette) <- unique(sample_annots$Condition)





      col_annot <- HeatmapAnnotation(df = sample_annots,
                                     col = list(hpi = hpi_palette,
                                                Condition=cond_palette),
                                     name = c("HPI", "Condition"))


      heatmap = Heatmap(z_score_scaled_data,
                        name = "log2 expression",
                        cluster_rows = TRUE,
                        cluster_columns = FALSE,
                        show_column_names = FALSE,
                        show_row_names = TRUE,
                        top_annotation = col_annot,
                        row_names_gp = gpar(fontsize=strip.text.x_size*0.8),
                        column_split = sample_annots$Condition,
                        heatmap_legend_param = list(direction = "horizontal"))

      x = list(data=exp_wide,
               plot = draw(heatmap, heatmap_legend_side = "top") )

      return(x)



  }else if(single_panel){
    if(n_genes <10){
      gene_colours <-  setNames(brewer.pal(n_genes, "Set1"), unique(TS_exp_w_names$name))
    }else{
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

      gene_colours <- setNames(sample(col_vector, n_genes), unique(TS_exp_w_names$name))

    }


    base_plot <- filter(TS_exp_w_names, Treatment =='Infected')%>%
      mutate(Fungi = Fungi %>% str_replace('Bot','B. cinerea')
             %>% str_replace('Scl','S. sclerotiorum'))%>%
      ggplot(data = ., aes(x = hpi, y = mean_expr)) +
      geom_line(aes(group = name, colour = name), linewidth = 1.325) +
      geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = name), alpha = 0.2) +
      geom_point(aes(x = hpi, y = rep1, colour = name), na.rm = TRUE, alpha = 0.6) +
      geom_point(aes(x = hpi, y = rep2, colour = name), na.rm = TRUE, alpha = 0.6) +
      geom_point(aes(x = hpi, y = rep3, colour = name), na.rm = TRUE, alpha = 0.6) +
      scale_color_manual(values = gene_colours) +
      scale_fill_manual(values = gene_colours)+
      facet_wrap(~Fungi)

    legend_fontface <- 'plain'

  }else if(!include_mock & length(fungi)>1){
    base_plot <- filter(TS_exp_w_names, Treatment =='Infected')%>%
      mutate(Fungi = Fungi %>% str_replace('Bot','B. cinerea')
             %>% str_replace('Scl','S. sclerotiorum'))%>%
      ggplot(data=., aes(x= hpi, y= mean_expr))+
      geom_line(aes(group = Fungi, colour = Fungi))+
      geom_ribbon(aes(ymin=LCL, ymax=UCL, fill = Fungi),alpha = 0.2) +
      geom_point(aes(x=hpi, y= rep1, colour = Fungi), na.rm=TRUE,alpha=0.6)+
      geom_point(aes(x=hpi, y= rep2,colour = Fungi), na.rm=TRUE,alpha=0.6)+
      geom_point(aes(x=hpi, y= rep3, colour = Fungi), na.rm=TRUE,alpha=0.6)+
      facet_wrap(~name, nrow=facet_nrow, scales = facet_scales)

    legend_fontface <- 'italic'
  }else{
    base_plot <- TS_exp_w_names%>%
      mutate(Fungi = Fungi %>% str_replace('Bot','B. cinerea')
             %>% str_replace('Scl','S. sclerotiorum')) %>%
      ggplot(data=., aes(x= hpi, y= mean_expr))+
      geom_line(aes(group = Treatment, colour = Treatment))+
      geom_ribbon(aes(ymin=LCL, ymax=UCL, fill = Treatment),alpha = 0.2) +
      geom_point(aes(x=hpi, y= rep1, colour = Treatment), na.rm=TRUE,alpha=0.6)+
      geom_point(aes(x=hpi, y= rep2,colour = Treatment), na.rm=TRUE,alpha=0.6)+
      geom_point(aes(x=hpi, y= rep3, colour = Treatment), na.rm=TRUE,alpha=0.6)+
      facet_wrap(name~Fungi, nrow=facet_nrow, scales = facet_scales)

    legend_fontface <- 'plain'
  }


  p <- base_plot+
    theme_bw()+
    labs(x = "Hours post Inoculation (HPI)",
         y= "Log2 Expression")+
    scale_x_continuous(breaks = seq(9,48,x_axis_interval))+
    theme(
      strip.text.x = ggtext::element_textbox( size = strip.text.x_size,
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
      legend.position="bottom",
      plot.title = element_blank(),
      legend.title = element_text(size =14, face = "bold"),
      legend.text = element_text(size = 13,
                                 face = legend_fontface))


  x = list(data=TS_exp_w_names,
           plot = p )


  return(x)



}


# plot <- plot_timeseries_expr(GO_id = 'response to jasmonic acid',
#                      #fungi = c('Bc'),
#                      max_n = 40,facet_nrow = 2,
#                      overlap_DEGs_only = TRUE,
#                      plot_type='heatmap',strip.text.x_size = 9)
#
#
# plot
