## ATAC seq analysis

The ATAC seq data was analysed using the seasnap pipeline
[https://github.com/bihealth/seasnap-pipeline](https://github.com/bihealth/seasnap-pipeline), 
branch ATACseq. 

In brief:

 - adapters removed with cutadapt
 - reads were mapped using STAR
 - peaks detected with Macs2
 - peak counts were collected with R package DiffBind
 - differential expression with DESeq2
 - motif enrichment using the MEME suite
 - gene set enrichment using the R packages tmod (Zyla et al. 2019) and clusterProfiler
 - the pipeline config files are `mapping_config.yaml` and `DE_config.yaml`

The full pipeline output and all objects and log files created by the
pipeline are archived in the SODAR system (Nieminen et al. 2023) and available upon request
from january.weiner _at_ bih-charite.de.

## Visualization / figures

The R markdown file "figures.rmd" contains the code necessary to replicate
the exact figures from the manuscript. To run the code, you will need the
original objects generated by the pipeline (see above) as well as the
Rseasnap R package, available from 
[https://github.com/bihealth/Rseasnap](https://github.com/bihealth/Rseasnap).

## Literature

Nieminen  M, Stolpe  O, Kuhring  M, Weiner  J, Pett  P, Beule  D, Holtgrewe  M
(2023). “SODAR: managing multiomics study data and metadata.”
_Gigascience_, *12*, giad052. doi:10.1093/gigascience/giad052
 [☞ Link](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giad052/7232111)

Zyla J, Marczyk M, Domaszewska T, Kaufmann SH, Polanska J, Weiner 3rd J.
Gene set enrichment for reproducible science: comparison of CERNO and eight
other algorithms.
[Bioinformatics](https://academic.oup.com/bioinformatics/article-abstract/35/24/5146/5511403).
2019 Dec 15;35(24):5146-54.
