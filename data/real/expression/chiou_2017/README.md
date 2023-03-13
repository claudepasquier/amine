# DATASETS

This directory contains the results of differential gene expression analyses performed on the datasets [GSE90625](https://www.omicsdi.org/dataset/geo/GSE90625) and [GSE90824](https://www.omicsdi.org/dataset/geo/GSE90824) associated with the article:
> Chiou SH, Risca VI, Wang GX, Yang D et al. BLIMP1 Induces Transient Metastatic Heterogeneity in Pancreatic Cancer. Cancer Discovery 2017 Oct;7(10):1184-1199. PMID: [28790031](https://www.ncbi.nlm.nih.gov/pubmed/28790031).

The analysis has been performed with the DESeq2 R package using default configuration.

## 1 - differential expression analysis of GSE90625

The file **GSE90625_all_counts_InVitro.txt.gz**, downloaded from [OmicsDI website(https://www.omicsdi.org/) contains the read counts of 30747 transcripts in 8 samples.
* shGFP_nor_1: Normoxia control KD rep1 
* shGFP_nor_2: Normoxia control KD rep2
* shGFP_hyp_1: Hypoxia control KD rep1
* shGFP_hyp_2: Hypoxia control KD rep2
* shBlimp1_nor_1: Normoxia Blimp1 KD rep1
* shBlimp1_nor_2: Normoxia Blimp1 KD rep2
* shBlimp1_hyp_1: Hypoxia Blimp1 KD rep1
* shBlimp1_hyp_2: Hypoxia Blimp1 KD rep2	

### 1.1 - shGFP (hypoxia vs normoxia): 
  
#### Input onditions:
|      Condition 0    |    Condition 1      |
| ------------------- | ------------------- |
|      shGFP_nor_1    |    shGFP_hyp_1      |
|      shGFP_nor_2    |    shGFP_hyp_2      |

#### Generated file:
* *shGFP_normoxia_vs_hypoxia.csv

### 1.2 - ShBlimp1 vs control (normoxia)
  
#### Input onditions:
|      Condition 0    |    Condition 1      |
| ------------------- | ------------------- |
|      shGFP_nor_1    |    shBlimp1_nor_1   |
|      shGFP_nor_2    |    shBlimp1_nor_2   |

#### Generated file:
* shGFP_vs_shBlimp1_normoxia.csv

### 1.2 - ShBlimp1 vs control (hypoxia)
  
#### Input onditions:
|      Condition 0    |    Condition 1      |
| ------------------- | ------------------- |
|      shGFP_hyp_1    |    shBlimp1_hyp_1   |
|      shGFP_hyp_2    |    shBlimp1_hyp_2   |

#### Generated file:
* shGFP_vs_shBlimp1_hypoxia.csv

## 2 - differential expression analysis of GSE90824

The file **GSE90824_All_ExVivo_counts.txt.gz**, downloaded from [OmicsDI website(https://www.omicsdi.org/) contains the read counts of 31468 transcripts in 12 samples.
* X0784.G.T: purified GFP/Hmga2-positive Tomato-positive PDAC cancer cells from mouse 0784
* X2317.G.T: purified GFP/Hmga2-positive Tomato-positive PDAC cancer cells from mouse 2317
* X2518.G.T: purified GFP/Hmga2-positive Tomato-positive PDAC cancer cells from mouse 2518
* X2542.G.T: purified GFP/Hmga2-positive Tomato-positive PDAC cancer cells from mouse 2542
* X2689.G.T: purified EpCAM-negative GFP/Hmga2-positive Tomato-positive PDAC cancer cells from mouse 2689
* X2691.G.T: purified GFP/Hmga2-positive Tomato-positive PDAC cancer cells from mouse 2691
* X0784.G.T..1: purified GFP/Hmga2-negative Tomato-positive PDAC cancer cells from mouse 0784
* X2317.G.T..1: purified GFP/Hmga2-negative Tomato-positive PDAC cancer cells from mouse 2317
* X2518.G.T..1: purified GFP/Hmga2-negative Tomato-positive PDAC cancer cells from mouse 2518
* X2542.G.T..1: purified GFP/Hmga2-negative Tomato-positive PDAC cancer cells from mouse 2542
* X2689.G.T..1: purified EpCAM-negative GFP/Hmga2-negative Tomato-positive PDAC cancer cells from mouse 2689
* X2691.G.T..1: purified GFP/Hmga2-negative Tomato-positive PDAC cancer cells from mouse 2691
 
### 1.2 - Hmga2_positive_vs_negative
  
#### Input onditions:
|      Codition 0     |      Condition 1    |
| ------------------- | ------------------- |
|      X0784.G.T      |    X0784.G.T..1     |
|      X2317.G.T      |    X2317.G.T..1     |
|      X2518.G.T      |    X2518.G.T..1     |
|      X2542.G.T      |    X2542.G.T..1     |
|      X2689.G.T      |    X2689.G.T..1     |
|      X2691.G.T      |    X2691.G.T..1     |

#### Generated file:
* Hmga2_positive_vs_negative.csv
