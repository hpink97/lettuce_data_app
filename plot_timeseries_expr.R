library(tidyverse)


if(!'helperFuncsLoaded' %in% ls()){
  funcs_path <- "G:/Shared drives/Denby Lab Team Drive/Lab members/Harry/lettuce_data_shiny_app/plotting_funcs/helper_funcs.R"
  source(funcs_path)
}

plot_timeseries_expr <- function(GeneIDs = NULL,
                                 At_orthologs = NULL,
                                 GO_id = NULL,
                                 protein_domain = NULL,
                    max_n = 1e6,
                    fungi = c("Bot", "Scl"),
                    include_mock = TRUE,
                    overlap_DEGs_only = FALSE,
                    up_down_only = 'both',
                    strip.text.x_size = 9,
                    x_axis_interval = 6,
                    facet_scales = 'free_y',
                    facet_nrow = 2,
                    return_name_conversion = FALSE){


  input_genes <- get_lettuce_genes_from_inputs(GeneIDs, At_orthologs, GO_id, protein_domain)
  if(length(input_genes) <1){ return(ggplot())}

  TS_exp <- get_timeseries_exp(genes =input_genes,
                     fungi = fungi,
                     include_raw_points = TRUE,
                     include_mock =  (include_mock | length(input_genes)>max_n),
                     ts_DEGs_only = overlap_DEGs_only,
                     up_down_only = up_down_only,
                     include_conf_intervals = TRUE
                     )
  n_genes <- length(unique(TS_exp$GeneID))


  if(n_genes > max_n){
    log2FC <- TS_exp %>%
      select(GeneID, Fungi, Treatment, mean_expr,hpi) %>%
      pivot_wider(names_from = 'Treatment',values_from = 'mean_expr') %>%
      mutate(log2FC = Infected-Mock) %>%
      group_by(GeneID) %>%
      summarise(mean_Abslog2FC = mean(abs(log2FC),na.rm = T), max_Abslog2FC = max(abs(log2FC),na.rm = T)) %>%
      arrange(-max_Abslog2FC)

    TS_exp <- filter(TS_exp, GeneID %in% log2FC$GeneID[1:max_n])




  }



  TS_exp_w_names <- merge(TS_exp, get_gene_names(TS_exp$GeneID), by='GeneID',all.x=TRUE) %>%
    mutate(Fungi = Fungi %>% str_replace('Bot','B. cinerea') %>% str_replace('Scl','S. sclerotiorum'))


  facet_nrow = ifelse(facet_nrow >=n_genes, n_genes, facet_nrow  )

  if(!include_mock & length(fungi)>1){
    base_plot <- filter(TS_exp_w_names, Treatment =='Infected')%>%
      ggplot(data=., aes(x= hpi, y= mean_expr))+
      geom_line(aes(group = Fungi, colour = Fungi))+
      geom_ribbon(aes(ymin=LCL, ymax=UCL, fill = Fungi),alpha = 0.2) +
      geom_point(aes(x=hpi, y= rep1, colour = Fungi), na.rm=TRUE,alpha=0.6)+
      geom_point(aes(x=hpi, y= rep2,colour = Fungi), na.rm=TRUE,alpha=0.6)+
      geom_point(aes(x=hpi, y= rep3, colour = Fungi), na.rm=TRUE,alpha=0.6)+
      facet_wrap(~name, nrow=facet_nrow, scales = facet_scales)

    legend_font_family <- 'italic'
  }else{
    base_plot <- ggplot(data=TS_exp_w_names, aes(x= hpi, y= mean_expr))+
      geom_line(aes(group = Treatment, colour = Treatment))+
      geom_ribbon(aes(ymin=LCL, ymax=UCL, fill = Treatment),alpha = 0.2) +
      geom_point(aes(x=hpi, y= rep1, colour = Treatment), na.rm=TRUE,alpha=0.6)+
      geom_point(aes(x=hpi, y= rep2,colour = Treatment), na.rm=TRUE,alpha=0.6)+
      geom_point(aes(x=hpi, y= rep3, colour = Treatment), na.rm=TRUE,alpha=0.6)+
      facet_wrap(name~Fungi, nrow=facet_nrow, scales = facet_scales)

    legend_font_family <- 'plain'
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
                                 face = legend_font_family))


  return(p)



}


plot_timeseries_expr(GO_id = 'jasmonic acid mediated signaling pathway',
                     fungi = c('Bc','Ss'),
                     max_n = 8,facet_nrow = 4,
                     include_mock = FALSE)
