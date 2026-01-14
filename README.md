#Drosophila fitness trade-offs
Repository of code for experiments assessing presence of fitness trade-offs using large, experimental populations of Drosophila melanogaster.

###Study Overview
The study consisted of two experiments:
1) An indoor experimental evolution study using four outbred, replicate popultions of the DGRP (http://dgrp2.gnets.ncsu.edu/data.html) whereby replicate popultions underwent sustained repredocution selection for nine, non-overlapping generations and a single bout of stress-tolerance, truncations selection. 
2) A paired indoor/outdoor mesocosm used to test whether SNPs identified via reproduction selection in an environmentally-controlled, lab setting displayed signals of selection and trade-offs in an outdoor mesocosm exposed to natural enviornmental fluctuations. Outdoor mesocosm data for this experiment were previously described and published by Bitter et al. 2024 (10.1038/s41586-024-07834-x).

###Repository description
This repository contains code for the analysis of allele frequency data from each of the two experiments described above.
It is important to note that the code notebooks are uploaded as .ipynb, .html, and .Rmd files. The code was orginally written using R and ipython notebooks and for reproducibility we have uploaded these notebooks as .html and .Rmd files.

For tidyness, functions called throughout the notebooks are from the R scripts indoor.cage.functions.R and general_cage_functions.R. I acknowledge co-author Sharon Greenblum for writing many of these functions.  

The notebook in Expansion_truncation_FstAnalyses.* contains R code for all Fst-based analysis (e.g., trends of Fst throughout the experiment and Fst-MDS analysis)

The notebook Expansion_Truncation_RawAlleleFrequencyAnalysis.* contains R code for identifying SNPs with evidence of linked selection during Experiment 1 and their relative behavior across expansion and truncation.

The directory SLiM_simulations contains the config files for the SLiM simulations exploring the likely proceses driving patterns of antagonistic pleiotropy observed during Ezpxeriment 1.


The notebook Indoor_Outdoor_Mesocosm2021_Analysis.html contains R code that analyzes the paired indoor-outdoor mesocosm experiment (conducted in 2021). This notebook draws upon a set of SNPs identified and reported previously by Bitter et al. 2021 (doi: 10.1038/s41586-024-07834-x).



###Data availability
Sequencing data from the DGRP lines used in the population expansion/truncation experiment (Experiment 1) are publicly available at NCBI accession PRJNA36679. Founder line sequences for the paired indoor/outdoor mesocosm study (Experiment 2) are available at NCBI Accession PRJNA722305. Raw sequencing reads from pooled samples collected during Experiments 1 and 2 are available at NCBI Accessions PRJNA1390176 and PRJNA1031645 , respectively. All raw allele frequency data necessary to replicate the results presented here are publicly available on Dryad: https://doi.org/10.5061/dryad.hx3ffbgt1

