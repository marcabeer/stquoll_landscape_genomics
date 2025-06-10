# Spotted-tailed quoll landscape genomics
This repository contains scripts and data underlying population genomics analyses for a RADseq dataset of spotted-tailed quolls (Dasyurus maculatus) across Tasmania.

The two major subdirectories contain the following:
- data
  - genetic_data (3431 SNPs from 345 individual quolls in various formats)
  - metadata (data giving geographic coordinates and other sampling metadata)
- analyses
  - env_data
    - 1km_res (1km resolution environmental data used in ResistanceGA)
    	- all files are inputs and intermediates produced by the Rscript except for those in the subdirectory /resistanceGA_rasters
    - 100m_res (100m resolution environmental data used in tests for IBE and GEAs)
    - tasveg_reclassification_key.csv (key that relates TASVEG broad landcover categories to the reclassified numeric codes found in categorical TASVEG raster used in analyses)
  - pop_struct
    - genetic_clusters (scripts to determine the number of genetic clusters in the dataset)
      - dapc (discriminant analysis of principal components)
      - faststructure
    - ibe (scripts for analyses related to isolation-by-environment)
      - gdm (generalized dissimilarity modeling)
      - prda (partial redundancy analysis of population genomic structure)
    - ibr (scripts for analyses related to isolation-by-resistance)
      - eems (estimated effective migration surfaces)
      - resistanceGA (optimized, expert opinion-free parameterization of resistance surfaces)
  - gea_tests (scripts for genetic-environment association tests of divergent selection)
    - lfmm (latent factor mixed models)
    - prda (partial redundancy analysis of individual genotypes)
    - overlap (determine which SNPs overlap as detections by both LFMM and pRDA)
    - candidate_genes (determine genes in close proximity to SNPs with GEAs)



