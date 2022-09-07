# Evaluation of TopDownCrawl (TDC)
Here we describe a detailed workflow for reproducing the tests used to evaluate the TopDownCrawl (TDC) method available at topdowncrawl.usc.edu. The method is used for the alignment of sequences with quantitative binding metrics from DNA binding experiments such as SELEX-seq.

Many popular tools for processing SELEX-seq data are designed to create a position weight matrix (PWM) to represent the binding preferences of a given DNA-binding protein, which assumes all positions contribute independently to binding. Although k-mer level enrichment is a highly reproducible measurement of binding affinity, it is difficult to know which region of the binding site is covered by each k-mer. For this purpose, we developed TopDownCrawl for aligning k-mers using their associated binding metrics. We evaluate the quality of alignments based on the performance of a multiple linear regression (MLR) model which is inherently sensitive to the alignment of the input sequences. To compare with PWM-based methods, we align the same set of k-mers to the generated PWMs and then perform MLR as was done for the TDC-based alignment. 

## Datasets Used for Evaluation
We evaluated TDC using 10 previously published SELEX-seq datasets released with the following publications. All reads were trimmed to only include the variable region.

*Abe, N., Dror, I., Yang, L., Slattery, M., Zhou, T., Bussemaker, H. J., Rohs, R., & Mann, R. S. (2015). Deconvolving the recognition of DNA shape from sequence. Cell, 161(2), 307-318.*

*Dantas Machado, A. C., Cooper, B. H., Lei, X., Di Felice, R., Chen, L., & Rohs, R. (2020). Landscape of DNA binding signatures of myocyte enhancer factor-2B reveals a unique interplay of base and shape readout. Nucleic acids research, 48(15), 8529-8544.*

*Zhang, L., Martini, G. D., Rube, H. T., Kribelbauer, J. F., Rastogi, C., FitzPatrick, V. D., Houtman, J. C., Bussemaker, H. J., & Pufall, M. A. (2018). SelexGLM differentiates androgen and glucocorticoid receptor DNA-binding preference over an extended binding site. Genome research, 28(1), 111-121.*

<details><summary style="font-size:14px">SRA Details</summary>

| SRA | Renamed File |
| --- | --- |
| SRR5340724 | AR_R0.fastq.gz |
| SRR5340730 | AR_R3.fastq.gz |
| SRR5340729 | AR_R4.fastq.gz |
| SRR5340724 | GR_R0.fastq.gz |
| SRR5340721 | GR_R3.fastq.gz |
| SRR5340720 | GR_R4.fastq.gz |
| SRR7450249 | MEF2B_R0.fastq.gz |
| SRR7450250 | MEF2B_R1.fastq.gz |
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
</details>

## Calculation of Relative Enrichment
The relative enrichments for 10-mers were calculated using a modified R script provided for the SELEX package available on bioconductor. Additionally, the script depends on the R package, 'stringr'.

https://bioconductor.org/packages/release/bioc/html/SELEX.html

```
./calculateEnrichment.R <R0 input> <R# input> <round #> <k length>
./calculateEnrichment.R AbdA-Exd_R0.fastq.gz AbdA-Exd_R3.fastq.gz 3 10
# Output: AbdA-Exd_R3_k10.tsv
```

## Alignment of 10-mers using TDC
The table of 10-mers with their assciated relative enrichments can be directly aligned using the TDC algorithm provided through our webserver at topdowncrawl.usc.edu or using the our python package available through pip.

```
pip install TopDownCrawl
TopDownCrawl AbdA-Exd_R3_k10.tsv
# Output: AbdA-Exd_R3_k10_aligned.tsv
```

## Alignment of 10-mers using PWM-based Methods

### BEESEM
BEESEM is available as a python package through a github page maintained by the original creator. In the original publication, the authors use the sequenceing results from the previous round of selection as a prior, but this data was not published for the Hox-Exd datasets. Instead, we used the initial library as the prior. Although this may affect the scaling of the output motif, I do not expect it to have an impact on the predicted ranking of differing binding sites, and should therefore not impact algnment. **For the AR and GR datasets, R3 sequencing was used as a prior for R4. For the MEF2B dataset, R1 was used as the prior for R2.** 

https://github.com/sx-ruan/BEESEM

*Ruan S, Swamidass SJ, Stormo GD. BEESEM: estimation of binding energy models using HT-SELEX data. Bioinformatics (2017)*

```
./countUnique.py <input fastq.gz>
./countUnique.py AbdA-Exd_R0.fastq.gz
./countUnique.py AbdA-Exd_R3.fastq.gz

beesem.py -s <most enriched 10-mer as seed> <output name> <processed prior> <processed R# input>
beesem.py -s ATGATTTATG beesem_k10 AbdA-Exd_R0.tsv AbdA-Exd_R3.tsv
# Output: ./beesem_k10_rep=1_phs=10/results/pfm_r0_p9_w10.txt

./alignWithBEESEM.py <input sequences> <beesem generated PFM>
./alignWithBEESEM.py AbdA-Exd_R3_k10.tsv ./beesem_k10_rep=1_phs=10/results/pfm_r0_p9_w10.txt
# Output: ./AbdA-Exd_R3_k10_BEESEM.tsv
``` 

### SelexGLM
SelexGLM is available as an R package through a github page maintained by the original lab. The included runSELEXGLM.R script is based on the script published for the original AR / GR paper by Zhang et al. For consistency accross datasets, symmetry was not enforced. 

https://github.com/BussemakerLab/SelexGLM/

<details><summary style="font-size:14px">SELEX Run Details (Used for SelexGLM configuration)</summary>

```
# Hox-Exd
varLen = 16
leftFixed = "CCGACGATCTGG"
rightFixed = "CCAGCTGTCGTAT"

# AR / GR
varLen = 23
leftFixed = "GTTCAGAGTTCTACAGTCCGACGATC"
rightFixed = "TGGAATTCTCGGGTGCCAAGG"

# MEF2B
varLen = 16
leftFixed = "GAGTTCTACAGTCCGACGATCCGC"
rightFixed = "CCTGGAATTCTCGGGTGCCA"
```
</details>

```
./runSelexGLM.R <R0 fastq.gz> <R# fastq.gz> <R#> <most enriched 10-mer as seed> <variable region length> <left adapter> <right adapter>
./runSelexGLM.R AbdA-Exd_R0.fastq.gz AbdA-Exd_R3.fastq.gz 3 ATGATTTATG 16 CCGACGATCTGG CCAGCTGTCGTAT
# Output: SelexGLM_k10_R3.pwm

./alignWithSelexGLM.py <input sequences> <SelexGLM generated PWM>
./alignWithSelexGLM.py AbdA-Exd_R3_k10.tsv SelexGLM_k10_R3.pwm
# Output: ./AbdA-Exd_R3_k10_SelexGLM.tsv

```

### MEME
MEME is a very popular and powerful tool for the alignment of DNA sequences, but inherently weights all sequences equally. This restricts the ability to use quantitative binding measurements in the alignment process. For this reason, we must restrict the table of 10-mers to only 10-mers which might create a useful binding motif. In this case, we arbitratily chose 10-mers with a log enrichment 2 standard deviations above the mean. This set of 10-mers is padded with k-1 ambiguous N bases and used to generate a single motif, which is then used for alignment of all 10-mers. Padding helps generate larger and more representative motifs. MEME can be accessed through a publicly available webserver or downloaded as part of the meme-suite package available on conda. 

https://meme-suite.org/meme/tools/meme
https://anaconda.org/bioconda/meme

```
./Zfilt.py AbdA-Exd_R3_k10.tsv
./tsvToFa.py AbdA-Exd_R3_k10_Z2.tsv
meme -revcomp -csites <# greater than the input length> -searchsize <# chars greater than input length> -dna <filtered input fa>
meme -revcomp -csites 10000 -searchsize 1000000 -dna AbdA-Exd_R3_k10_Z2.fa 
meme-get-motif -all meme_out/meme.txt > meme_out/motif.txt
# Output: ./meme_out/motif.txt

./alignWithMEME.py <input sequences> <MEME generated PFM>
./alignWithMEME.py AbdA-Exd_R3_k10.tsv meme_out/motif.txt
# Output: ./AbdA-Exd_R3_k10_MEME.tsv
```

## Multiple Linear Regression (MLR)
Multiple linear regression is performed with elastic net regularization using the sci-kit learn package in python. Models are trained and evaluated using 5-fold cross validation and the median R² is reported. Models are trained to predict the log enrichment of the aligned input using 1-hot encoded sequences, in which case gaps are encoded as an empty vector. Additionally, minor groove width (MGW) and electrostatic potential (EP) features are added to the input using feature tables derived from the DNAshapeR package available on bioconductor. A random seed of 0 is used arbitratily for consistence accross datasets. Since low-affinity binding sites are less likely to have a meaningful alignment, we chose to evaluate only sequences with a log enrichment at least 2 standard deviations above the mean. 

```
./Zfilt.py <input alignment>
./Zfilt.py AbdA-Exd_R3_k10_aligned.tsv
./Zfilt.py AbdA-Exd_R3_k10_BEESEM.tsv
./Zfilt.py AbdA-Exd_R3_k10_SelexGLM.tsv
./Zfilt.py AbdA-Exd_R3_k10_MEME.tsv

./MLR.py <filtered alignment>
./MLR.py AbdA-Exd_R3_k10_aligned_Z2.tsv
./MLR.py AbdA-Exd_R3_k10_BEESEM_Z2.tsv
./MLR.py AbdA-Exd_R3_k10_SelexGLM_Z2.tsv
./MLR.py AbdA-Exd_R3_k10_MEME_Z2.tsv
```