# Evaluation of TopDownCrawl (TDC)
Here we describe a detailed workflow for reproducing the tests used to evaluate the TopDownCrawl (TDC) method available at topdowncrawl.usc.edu. The method is used for the alignment of sequences with quantitative binding metrics from DNA binding experiments such as SELEX-seq.

## Datasets Used for Evaluation
We evaluated TDC using 10 previously published SELEX-seq datasets released with the following publications. 

*Abe, N., Dror, I., Yang, L., Slattery, M., Zhou, T., Bussemaker, H. J., Rohs, R., & Mann, R. S. (2015). Deconvolving the recognition of DNA shape from sequence. Cell, 161(2), 307-318.*

*Dantas Machado, A. C., Cooper, B. H., Lei, X., Di Felice, R., Chen, L., & Rohs, R. (2020). Landscape of DNA binding signatures of myocyte enhancer factor-2B reveals a unique interplay of base and shape readout. Nucleic acids research, 48(15), 8529-8544.*

*Zhang, L., Martini, G. D., Rube, H. T., Kribelbauer, J. F., Rastogi, C., FitzPatrick, V. D., Houtman, J. C., Bussemaker, H. J., & Pufall, M. A. (2018). SelexGLM differentiates androgen and glucocorticoid receptor DNA-binding preference over an extended binding site. Genome research, 28(1), 111-121.*

<details><summary style="font-size:14px">SRA Details
| SRA | Renamed File |
| --- | --- |
| SRR5340724 | AR_R0.fastq.gz |
| SRR5340729 | AR_R4.fastq.gz |
| SRR5340724 | GR_R0.fastq.gz |
| SRR5340720 | GR_R4.fastq.gz |
| SRR7450249 | MEF2B_R0.fastq.gz |
| SRR7450251 | MEF2B_R2.fastq.gz |
| SRR1765757 | AbdA-Exd_R0.fastq.gz |
| SRR1765754 | AbdA-Exd_R3.fastq.gz |
| SRR1765757 | Dfd-Exd_R0.fastq.gz |
| SRR1765752 | Dfd-Exd_R3.fastq.gz |
| SRR1765757 | Lab-Exd_R0.fastq.gz |
| SRR1765751 | Lab-Exd_R3.fastq.gz |
| SRR1765759 | PbFl-Exd_R0.fastq.gz |
| SRR1765746 | PbFl-Exd_R3.fastq.gz |
| SRR1765756 | Scr-Exd_R0.fastq.gz |
| SRR1765733 | Scr-Exd_R3.fastq.gz |
| SRR1765756 | UbxIa-Exd_R0.fastq.gz |
| SRR1765750 | UbxIa-Exd_R3.fastq.gz |
| SRR1765757 | UbxIVa-Exd_R0.fastq.gz |
| SRR1765753 | UbxIVa-Exd_R3.fastq.gz |
</summary></details>