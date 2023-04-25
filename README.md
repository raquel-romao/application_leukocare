# Coding MWE 

For a minimal working example I thought it would make sense to present two distinct projects I have done in the past. One of them a project in Urban Data Science in Python (version 3.8.4) and other a Case Study about differential expression of a Cancer Dataset in R (version 4.1.1).

## Project in _Urban_ Data Science

As one of the coding MWE, I decided to present a previous project that I did for one of the courses I completed in the Delft University of Technology, Introduction to _Urban_ Data Science, where I chose to step out of my confort zone and deal with non-biological data.

The project aimed to formulate and test an hypothesis about an Urban area. I hypothesized that the neighbourhoods of The Hague, The Netherlands, with **higher average incomes** are more prone to **emit $CO_2$**. Later, I also analysed if the bike lane length per neighbourhood could be related to it, hypotesizing that a bigger amount of **bike lanes** available could translate in a minor use of transportation by motorized vehicles and therefore a minor $CO_2$ emission, being 2016 the base year for my analysis.

The data used for this project was extracted from [Den Haag in Cijfers](https://denhaag.incijfers.nl/jive) as well as made available by the docent team, in particular a shapefile of The Hague's neighbourhoods (under the urbanDS_project/data folder in this repository).




## Case Study: Differential Expression of Cancer Dataset

As an additional MWE, I decided to include a case study that I solved in the past where, contrarily to the last MWE, I make use of R and biological data.

For this case study, I had to suppose that I was a scientist investigating which genes and pathways are differentially expressed in cancer. The gene expression of 3 affected patients (disease 1 - 3) and 3 controls (controls 1 - 3) of 1000 genes were measured. The resulting count matrix is stored in the file “expression_counts.txt” (under the CaseStudy folder in this repository). 

I had to call differential expression using the library DESeq2 between two sets of RNA-seq samples (control and disease) and then, obtain a list of differential expressed genes at False Discovery Rate (i.e. padj) after 0.05, and define which ones of them are up and which ones are down regulated (part 1 of the script).

Aditionally, I had to use the list of Tumour Suppressor Genes (TSGs) and oncogenes from the [COSMIC database](https://cancer.sanger.ac.uk/census) and find out if TSGs and/or oncogenes over-repressented among over- or under-expressed genes, while assessing statistical significance of the enrichment, if any (part 2 of the script).
