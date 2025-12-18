#Drosophila fitness trade-offs
Repository of code for experiments assessing presence of fitness trade-offs using large, experimental populations of Drosophila melanogaster.

###Study Overview
The study consisted of two experiments:
1) An indoor experimental evolution study using four outbred, replicate popultions of the DGRP (http://dgrp2.gnets.ncsu.edu/data.html) whereby replicate popultions underwent sustained repredocution selection for nine, non-overlapping generations and a single bout of stress-tolerance, truncations selection. 
2) A paired indoor/outdoor mesocosm used to test whether SNPs identified via reproduction selection in an environmentally-controlled, lab setting displayed signals of selection and trade-offs in an outdoor mesocosm exposed to natural enviornmental fluctuations. Outdoor mesocosm data for this experiment were previously described and published by Bitter et al. 2024 (10.1038/s41586-024-07834-x).

###Repository description
This repository contains code for the analysis of allele frequency data from each of the two experiments described above.

The R scripts glm.FullExpansion.x2.R and glm.Truncation1.x2.R run a generlized linear model on allele frequency data from the expansion and truncation samples, respectively.

The notebook in IndoorCage_ExpansionGLM_Manhattans.html contains R code that takes the output of glm.FullExpansion.x2.R , identifies significant sites, and generates a manhattan plot from the output regression data.

The notebook in IndoorCage_Fst_MDS.html contains R code for all Fst-based analysis (e.g., trends of Fst throughout the experiment and Fst-MDS analysis)

The notebook FullExpansion_Truncation1_SiteComparison.html contains R code for assessing the dynamics of expansion-favored alleles throughout truncation.

The notebook IndoorCage_InversionAnalysis.html contains R code that visualizes the frequency of inversions throughout expansion and truncation.

The notebook Indoor_Outdoor_Mesocosm2021_Analysis.html contains R code that analyzes the paired indoor-outdoor mesocosm experiment conducted in 2021.

The directory SLiM_simulations contains the config files for the SLiM simulations run for the revised manuscript


###Data availability
Sequencing data from the DGRP lines used in the population expansion/truncation experiment (Experiment 1) are publicly available at http://dgrp2.gnets.ncsu.edu/data.html.  Founder line sequences for the paired indoor/outdoor mesocosm study (Experiment 2) are available at NCBI Accession PRJNA722305. Raw sequencing reads from pooled samples collected during Experiments 1 and 2 are available at NCBI Accessions PRJNA1390176 and PRJNA1031645 , respectively¬. All raw allele frequency data necessary to replicate the results presented here are publicly available on Dryad: https://doi.org/10.5061/dryad.hx3ffbgt1.<img width="468" height="85" alt="image" src="https://github.com/user-attachments/assets/c706599e-8d2e-42f6-b996-c43caae4523e" />

