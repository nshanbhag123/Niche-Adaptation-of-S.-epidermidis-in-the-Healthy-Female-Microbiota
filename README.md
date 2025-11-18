# Niche Adaptation of Staphylococcus epidermidis in the Healthy Female Microbiota

This is code is related to the paper: "Strain-Level Genomic Analysis Reveals Cross-Site Colonization and Niche Adaptation of Staphylococcus epidermidis in the Healthy Female Microbiota" conducted in the work of Sandra Jablonska, Niru Shanbhag, Alex Kula, and Catherine Putonti

## Abstract

Staphylococcus epidermidis is a ubiquitous skin commensal that also colonizes the nasal, oral, and urinary tracts. To investigate strain-level dynamics and niche specialization, we performed whole genome sequencing on 114 S. epidermidis isolates collected from the skin, nasal cavity, oral cavity, and urine of healthy female participants. Pairwise average nucleotide identity (ANI) revealed multiple instances of identical strains (>99.99% ANI) across different body sites within the same individual, suggesting a generalist colonization strategy. Accessory genome composition did not significantly differ by site, and intact prophage sequences were often shared across anatomical sites and individuals, highlighting pervasive gene flow mediated by mobile elements. Nevertheless, multinomial logistic regression using gene-cluster presence/absence identified site-associated functions—including amidases linked to skin isolates and LysR-type regulators and AraJ transporters enriched in nasal isolates—indicating subtle but detectable niche-adaptive signatures. Together, these findings show that S. epidermidis combines a generalist colonization strategy with targeted functional adaptations, supported by a highly fluid accessory genome. This study advances our understanding of S. epidermidis population structure in healthy hosts and provides a genomic framework for distinguishing commensal adaptation from pathogenic potential.

## Code Details
This repo includes two scripts to generate 1) multinomial logistic regression machine learning model to identify gene clusters associated with specific niches in S. epidermidis 2) PcoA plot generated using a COG presence/absence table of the S. epi accessory genome

### Packages for Machine Learning (Python)
* pandas
* sklearn
* numpy
* seaborn
* matplotlib

### Data for Machine Learning (Python)
* private_gc_matrix.csv - gene cluster presence/absence matrix for S. epidermidis samples collected from Loyola University Chicago
* private_acc2iso.csv - maps accession IDs from LUC samples to their isolation source
* public_gc_matrix.csv - gene cluster presence/absence matrix for S. epidermidis samples from publicly available sources
* public_acc2iso.csv - maps accession IDs from publicly available S. epi samples to their isolation source

### Packages for PcoA (R)
* vegan
* ggplot

### Data for PcoA (R)
* private_gc_matrix.csv - gene cluster presence/absence matrix for S. epidermidis samples collected from Loyola University Chicago

## Authors

* Sandra Jablonska (R code) - inquiries can be made at: sandrajablo@gmail.com
* Niru Shanbhag (Python code) - inquiries can be made at: shanbhagn123@gmail.com

## Acknowledgments
Thank you to SJ, who spent 2 years generating all this lovely data :) 
