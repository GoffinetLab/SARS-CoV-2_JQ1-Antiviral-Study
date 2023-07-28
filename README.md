# SARS-CoV-2_JQ1-Antiviral-Study

## ABSTRACT

Inhibitors of bromodomain and extra-terminal proteins (iBETs), including JQ-1, have been suggested as potential prophylactics against SARS-CoV-2 infection. However, molecular mechanisms underlying JQ-1-induced antiviral activity and its susceptibility to viral antagonism remain incompletely understood. iBET pre-treatment of cells reduced their susceptibility to infection by SARS-CoV-2 variants and SARS-CoV, but not MERS-CoV, with IC50 concentrations of 0.2-1.3 µM. The antiviral effect manifested itself by reduced reporter expression of recombinant viruses, reduced abundance of viral RNA and lowered infectious titer in the culture supernatant. While we confirmed JQ-1-mediated downregulation of angiotensin-converting enzyme 2 (ACE2) expression, our multi-omics analysis addressing the chromatin accessibility, transcriptome and proteome of infected and uninfected cells in the absence and presence of JQ-1 uncovered induction of an antiviral NRF-2-mediated cytoprotective response as an additional antiviral component of JQ-1 treatment. Serial passaging of SARS-CoV-2 in the presence of JQ-1 resulted in predominance of ORF6-deficient variants, suggesting a minimised need for SARS-CoV-2 ORF6-mediated repression of IFN signaling. The induction and maintenance of the JQ-1-mediated antiviral state was prevented by active SARS-CoV-2 infection. We propose that JQ-1 exerts pleiotropic effects that collectively induce an antiviral state that is ultimately nullified by an established SARS-CoV-2 infection, raising questions about the clinical suitability of iBETs in the context of COVID-19.

## OVERVIEW

This is the GitHub repository for the manuscript "Pharmacological inhibition of bromodomain and extra-terminal proteins induces NRF-2-mediated inhibition of SARS-CoV-2 replication and is subject to viral antagonism" by Mhlekude _et al._, 2023. The code utilised for the analysis of the RNAseq and ATACseq data is deposited here. 

### RNAseq

Unix/BASH and R scripts are deposited in a directory structure that reflects the different compartments of the analysis. Code used for the initial alignment of raw reads to the reference genome and counting of gene features is deposited in the folder "0_ReferenceAlignment_FeatureCount". Preprocessing of the resulting count matrix is included in the folder "1_Preprocessing" and further downstream analyses are included in the folder "2_FurtherAnalyses".

### ATACseq

See [ATAC-seq analysis README file](ATACseq)

## REFERENCE GENOME

The analysis for the RNAseq reads was performed by aligning the reads to a combined human/SARS-CoV-2 genome. To create this reference, the SARS-CoV-2 genome (GenBank entry MN908947) was appended to the human genome assembly h38 (ENSEMBL v102), i.e. both the FASTA and GTF files were merged. The SARS-CoV-2 FASTA and GTF files are included in the folder "0_ReferenceAlignment_FeatureCount". 

For further information, please refer to the manuscript itself.
