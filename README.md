

# Lettuce Data Explorer

Welcome to the Lettuce Data Explorer, a Shiny application for exploring lettuce transcriptomic datasets after infection with necrotrophic pathogens (*Botrytis cinerea* and *Sclerotinia sclerotiorum*) and gene regulatory network inferred from them!

This application is hosted on shinyapps.io;
## 🔗 **[Visit the Lettuce Data Explorer on Shinyapps.io](https://hpink97.shinyapps.io/Lettuce-Data/)** 🔗


## Overview
The Lettuce Data Explorer simplifies the exploration of lettuce transcriptomic datasets. It empowers users to swiftly identify lettuce genes using different selection criteria, such as Arabidopsis symbols/IDs, GO terms, protein domains, and more.

### Features:

1. **Data Selection:** Select from multiple datasets, such as:
   - **Time-series expression**: Dive into dynamic expression profiles.
   - **Gene Regulatory Network Analysis**: Examine predicted transcriptional regulation.
   - **Lesion Size Correlation**: Analyze gene expression's correlation with lesion size.
   
2. **Gene Selection:** Choose genes based on:
   - **Lettuce GeneID**: Input using Lettuce gene ID(s).
   - **Orthologues of Arabidopsis Genes**: Specify an Arabidopsis gene ID or symbol.
   - **Genes with GO-term**: Enter a GO term to filter genes.
   - **Genes with Protein Domain**: Discover genes based on their associated protein domains.
   
3. **Customization Options**:
   - Dataset-specific gene selection criteria.
   - Plot customization capabilities.
   
4. **Results Generation**: After making selections, click 'Generate Results' to view plots and tables.

5.  **Save results**: all generated plots, their underlying data, and GRN tables can be downloaded directly from the Shiny app

## Lettuce gene annotations
Gene annotations in the app, such as Arabidopsis orthologues and protein domains, are based on the work of [Reyes-Chin-Wo et al, 2017](https://doi.org/10.1038/ncomms14953). The GO-term annotations derive from the GO-term annotation of the closest Arabidopsis orthologue of each lettuce gene.


## Citations
For an in-depth understanding of the datasets, please consider reading the following publications:

- [Pink, H., Talbot, A., Graceson, A. et al. Identification of genetic loci in lettuce mediating quantitative resistance to fungal pathogens. Theor Appl Genet, 135, 2481–2500 (2022).](https://doi.org/10.1007/s00122-022-04129-5)
  
- [Pink, H., Talbot, A., Carter, R.,. et al. Identification of Lactuca sativa transcription factors impacting resistance to Botrytis cinerea through predictive network inference. bioRxiv (2023).](https://doi.org/10.1101/2023.07.19.549542)

