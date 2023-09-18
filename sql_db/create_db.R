require(RSQLite)
require(DBI)
require(tidyverse)

setwd("G:/Shared drives/Denby Lab Team Drive/Lab members/Harry/R functions/Data")
rm(list=ls())
dir()

lettuce_divset <- read_rds("Lettuce_DivSet.rds") %>%
  mutate(sample_name = paste(sample, pathogen, sep='_'))

divset_pheno <- distinct(lettuce_divset, sample_name, SampleID, sqrt_lesion)

divset_expr <- dplyr::select(
  lettuce_divset,
  sample_name, GeneID, log_Exp
) %>%
  pivot_wider(
    names_from = 'sample_name', values_from = 'log_Exp'
)

divset_expr_Bot <- divset_expr[, grepl('GeneID|_Bot', names(divset_expr))]
divset_expr_Scl <- divset_expr[, grepl('GeneID|_Scl', names(divset_expr))]

head(lettuce_divset)

timeseries <- read_rds("Let_BotScl_TS.rds") %>%
  dplyr::filter(grepl(pattern = '^Lsat',x = GeneID)) %>%
  dplyr::select(GeneID, Fungi, Treatment,
                hpi=Time_Point,
                rep1_log2, rep2_log2,
                rep3_log2)


# bc_timeseries <- dplyr::filter(timeseries,
#                                Fungi =='Bot')
#
# ss_timeseries <- dplyr::filter(timeseries,
#                                Fungi =='Scl')

scl_divset_corr.DEGs_filtered <- readxl::read_excel("G:/Shared drives/Denby Lab Team Drive/Paper writing/Lettuce QTL paper/Supplemental Datasets/Supp_Dataset5_lesion_corr.xlsx",sheet = 3, skip = 1)

scl_divset_corr.signif <- readxl::read_excel("G:/Shared drives/Denby Lab Team Drive/Paper writing/Lettuce QTL paper/Supplemental Datasets/Supp_Dataset5_lesion_corr.xlsx",sheet = 1, skip = 1) %>%
  mutate(confounding_DE = ifelse(!(GeneID %in% scl_divset_corr.DEGs_filtered$GeneID) , TRUE, FALSE)) %>%
  select(GeneID, sclero_lesion_cor = Lesion_cor, sclero_lesion_cor_padj = p_adj, confounding_DE)




grn.path <- "G:/Shared drives/Denby Lab Team Drive/Paper writing/Lettuce Botrytis time series paper/Analysis/Host/Botrytis-Scelrotinia-combo/OutPredict_GRNs"

hubs <- read_csv(file.path(grn.path,"TF_outdegrees.csv"))

ww_modules <- readxl::read_excel("G:/Shared drives/Denby Lab Team Drive/Paper writing/Lettuce Botrytis time series paper/Submitted version/Supp Dataset 5 Modules.xlsx",sheet = 1,skip = 1) %>%
  dplyr::select(GeneID, wigwam_module= Module)

edges <- read_csv(file.path(grn.path,"high_confidence_edges.csv")) %>%
  merge(x=., y=ww_modules, by.x='TF',by.y='GeneID') %>%
  merge(x=., y=ww_modules, by.x='Target',by.y='GeneID',suffixes=c('_TF','_Target'))

#divset_samples <- colnames(divset_expr)[grepl('_Scl',colnames(divset_expr))]

divset_expressed <- unique(divset_expr$GeneID)

calculate_divset_correlation <- function(TF, target, divset_samples) {


  if(any( !c(TF, target) %in% divset_expressed)){
    return(NA)
  }




  TF_expr <- divset_expr[divset_expr$GeneID==TF, divset_samples] %>% unlist()
  target_expr <- divset_expr[divset_expr$GeneID==target, divset_samples] %>% unlist()

  both_not_NA <- !is.na(TF_expr) & !is.na(target_expr)

  cor_value <- cor(TF_expr[both_not_NA], target_expr[both_not_NA], method = "pearson") # Use complete.obs to handle any NAs

  return(cor_value)
}

edges$Ss_DivSet_coexp_cor <- mapply(calculate_divset_correlation,
                                        edges$TF, edges$Target,
                                        MoreArgs = list(divset_samples = colnames(divset_expr)[grepl('_Scl',colnames(divset_expr))]))

edges$Bc_DivSet_coexp_cor <- mapply(calculate_divset_correlation,
                                        edges$TF, edges$Target,
                                        MoreArgs = list(divset_samples = colnames(divset_expr)[grepl('_Bot',colnames(divset_expr))]))


prot_domain <- read_csv("G:/Shared drives/Denby Lab Team Drive/Paper writing/Lettuce Botrytis time series paper/Analysis/Host/Data/protein_domain_data/protein_domains.csv") %>%
  filter(database %in% c('Pfam','Panther')) %>%
  dplyr::select(GeneID=gene_id, id, description = desc) %>%
  #mutate(description = tolower(description)) %>%
  distinct() %>%
  drop_na()

# low_freq_domains <- filter(count(prot_domain, id),n<4)$id
# low_freq_domain_genes <- filter(prot_domain, id %in% low_freq_domains )$GeneID
# multiple_domains <- filter(prot_domain, !id %in%low_freq_domains )
#
# prot_domain_filtered <- prot_domain %>%
#   filter(!(GeneID %in% multiple_domains$GeneID & id %in% low_freq_domains))

count(prot_domain, id, description, sort=T) %>% head(100)

gene_annots <- read_csv("../../Data Analysis/Gene_Annots/Lettuce_Gene_Annots/Complete_Lettuce_Gene_Annots.csv") %>%
  mutate(AGI = str_sub(TAIR10.ID,1,9))

degs_tofde <- gene_annots %>%
  dplyr::filter(!is.na(Bot_Inf_De_Dir) & Scl_Inf_De_Dir==Bot_Inf_De_Dir) %>%
  dplyr::select(GeneID,Bot_Inf_De_Dir, Bot_Inf_TOFDE,Scl_Inf_De_Dir,Scl_Inf_TOFDE )

annots <- dplyr::select(gene_annots, GeneID, AGI,TAIR10.fescription,At_ShortName,Protein_Domains)


con <- dbConnect(RSQLite::SQLite(), '../../lettuce_data_shiny_app/sql_db/lettuce_data.sqlite')


dbWriteTable(con, "divset_pheno", divset_pheno,overwrite=TRUE)
dbWriteTable(con, "bot_divset_expr", divset_expr_Bot,overwrite=TRUE)
dbWriteTable(con, "sclero_divset_expr", divset_expr_Scl,overwrite=TRUE)
dbWriteTable(con, "timeseries", timeseries,overwrite=TRUE)
dbWriteTable(con, "sclero_divset_corr", scl_divset_corr.signif,overwrite=TRUE)
dbWriteTable(con, "grn_edges", edges,overwrite=TRUE)
dbWriteTable(con, "grn_hubs", hubs,overwrite=TRUE)
dbWriteTable(con, 'timeseries_overlap_degs',degs_tofde,overwrite=TRUE)
dbWriteTable(con, 'gene_annotations',annots,overwrite=TRUE)
dbWriteTable(con, 'domain_annotations',prot_domain,overwrite=TRUE)

con


# erf1 <- 'Lsat_1_v5_gn_3_121961'
#
# query = paste0("SELECT * FROM sclero_divset_expr ",
#                "WHERE GeneID ='",erf1,"'")
# dbGetQuery(con, query)



dbDisconnect(con)

