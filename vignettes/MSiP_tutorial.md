---
title: "Mass Spectrometry interaction Prediction (MSiP)"
author: | 
 | Matineh Rahmatbakhsh
 | matinerb.94@gmail.com
email: "matinerb.94@gmail.com"
package: "MSiP"
output: html_vignette
vignette: >
    %\VignetteEngine{knitr::knitr}
    %\VignetteIndexEntry{MSiP tutorial}
    %\usepackage[UTF-8]{inputenc}
---




#Introduction
The MSiP is a computational approach to predict protein-protein
interactions (PPIs) from large-scale affinity purification mass spectrometry (AP-MS) data. This approach includes both spoke and matrix models for interpreting AP-MS data in a network context. The "spoke" model considers only bait-prey interactions, whereas the "matrix" model assumes that each of the identified proteins (baits and prey) in a given AP-MS experiment interacts with each of the others. The spoke model has a high false-negative rate, whereas the matrix model has a high false-positive rate. Thus, although both statistical models have merits, a combination of both models has shown to increase the performance of machine learning classifiers in terms of their capabilities in discrimination between true and false positive interactions [Drew et al., 2017](https://www.embopress.org/doi/full/10.15252/msb.20167490/). 

###Load the package

```r
library(MSiP)
```

###Sample Data Description:
A demo AP-MS proteomics dataset is provided in this package to guide the users about data structure.

```r
data("SampleDatInput")
head(SampleDatInput)
```

```
##    Experiment.id Replicate  Bait   Prey counts
## 1:            14    BN2198  NSP5 Q07065     17
## 2:            23    BN2214 ORF7A O75947      3
## 3:            16    BN2173  NSP7 Q8TCU6      2
## 4:             7    BN2177 NSP11 P30153      2
## 5:            22    BN2186  ORF6 Q8NBM4      1
## 6:             4    BN2190     N P51571      1
```

###Scoring based on "spoke-model":
Comparative Proteomic Analysis Software Suite (CompPASS) is a robust statistical scoring scheme for assigning scores to bait-prey interactions [Sowa et al., 2009](https://pubmed.ncbi.nlm.nih.gov/19615732/). The output from CompPASS scoring includes Z-score, S-score, D-score, WD-score and other features. This function was optimized from the [source code](https://github.com/zqzneptune/SMAD).  


















