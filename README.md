# Strain-Level Genomic Analysis of Staphylococcus epidermidis Across Multiple Body Sites in Healthy Females
This is code is related to the paper: "Strain-Level Genomic Analysis Reveals Cross-Site Colonization and Niche Adaptation of Staphylococcus epidermidis in the Healthy Female Microbiota" conducted in the work of Sandra Jablonska, Niru Shanbhag, Alex Kula, and Catherine Putonti

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
