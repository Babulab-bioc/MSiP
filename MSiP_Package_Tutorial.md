---
    title: "MSiP tutorial"
author: "Matineh Rahmatbakhsh"
date: "June 11, 2021"
output: github_document
vignette: >
    %\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Title of your vignette}
%\usepackage[UTF-8]{inputenc}
---
    
    
    #Introduction
    The MSiP is a computational approach to predict protein-protein
interactions (PPIs) from large-scale affinity purification mass spectrometry (AP-MS) data. This approach includes both spoke and matrix models for interpreting AP-MS data in a network context. The "spoke" model considers only bait-prey interactions, whereas the "matrix" model assumes that each of the identified proteins (baits and prey) in a given AP-MS experiment interacts with each of the others. The "spoke" model has a high false-negative rate, whereas a high positive rate has resulted from the matrix model. Thus, although both statistical models have merits, a combination of both spoke and matrix models has shown to increase the performance of machine learning classifiers in terms of their capabilities in discrimination between true and false positive interactions [Drew et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5488662). 

###Load the package

```r
library(MSiP)
```

```
## Loading required package: magrittr
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
Comparative Proteomic Analysis Software Suite (CompPASS) is a robust statistical scoring scheme for assigning confidence scores to bait-prey interactions [Sowa et al., 2009](https://www.cell.com/cell/fulltext/S0092-8674(09)00503-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867409005030%3Fshowall%3Dtrue). The output from CompPASS scoring includes Z-score, S-score, D-score, WD-score and other features. This function was optimized from the [source code](https://github.com/zqzneptune/SMAD).  


```r
datScoring <- 
    cPASS(SampleDatInput)
head(datScoring)
```

```
##   Bait   Prey AvePSM TotalPSM Ratio_totalPSM      Ratio   S_score  Z_score
## 1    E A5YKK6      4        4      1.0000000 0.03846154 10.198039 4.899136
## 2    E O15144      5       17      0.2941176 0.11538462  6.582806 2.344863
## 3    E O95299      1        1      1.0000000 0.03846154  5.099020 4.899136
## 4    E P04899      7        7      1.0000000 0.03846154 13.490738 4.899136
## 5    E P07355      1        1      1.0000000 0.03846154  5.099020 4.899136
## 6    E P07954      4        4      1.0000000 0.03846154 10.198039 4.899136
##     D_score  WD_score Entropy
## 1 10.198039 0.6014181       0
## 2  6.582806 0.2893455       0
## 3  5.099020 0.3007091       0
## 4 13.490738 0.7956014       0
## 5  5.099020 0.3007091       0
## 6 10.198039 0.6014181       0
```

###Scoring based on "matrix-model":
The Dice coefficient was first applied by [Zhang et al., 2008](https://academic.oup.com/bioinformatics/article/24/7/979/296061) to score interaction between all identified proteins (baits and preys) in a given AP-MS expriment.

```r
datScoring <- 
    diceCoefficient(SampleDatInput)
head(datScoring)
```

```
##             BPI      Dice
## 1 A4D1P6~A5YKK6 0.0000000
## 2      A4D1P6~E 0.0000000
## 3 A4D1P6~E9PAV3 0.0000000
## 4      A4D1P6~M 0.0000000
## 5      A4D1P6~N 0.0000000
## 6   A4D1P6~NSP1 0.6666667
```

Alternatively, Jaccard, Simpson, and Overlap scores can be used to score the interaction between all the identified proteins in a given AP-MS experiment. 

```r
#Jaccard coefficient
datScoring <- 
    jaccardCoefficient(SampleDatInput)
head(datScoring)
```

```
##             BPI Jaccard
## 1 A4D1P6~A5YKK6     0.0
## 2      A4D1P6~E     0.0
## 3 A4D1P6~E9PAV3     0.0
## 4      A4D1P6~M     0.0
## 5      A4D1P6~N     0.0
## 6   A4D1P6~NSP1     0.5
```

```r
#Simpson coefficient
datScoring <- 
    simpsonCoefficient(SampleDatInput)
head(datScoring)
```

```
##             BPI Simpson
## 1 A4D1P6~A5YKK6       0
## 2      A4D1P6~E       0
## 3 A4D1P6~E9PAV3       0
## 4      A4D1P6~M       0
## 5      A4D1P6~N       0
## 6   A4D1P6~NSP1       1
```

```r
#Overlap score
datScoring <- 
    simpsonCoefficient(SampleDatInput)
head(datScoring)
```

```
##             BPI Simpson
## 1 A4D1P6~A5YKK6       0
## 2      A4D1P6~E       0
## 3 A4D1P6~E9PAV3       0
## 4      A4D1P6~M       0
## 5      A4D1P6~N       0
## 6   A4D1P6~NSP1       1
```

Finally, a weighted matrix model [Drew et al., 2017](https://www.embopress.org/doi/full/10.15252/msb.20167490) can also be employed to score interactions between identified proteins in a given AP-MS experiment. The output of the weighted matrix model includes the number of experiments for which the pair of proteins is co-purified (i.e., k) and $-1$*log(P-value) of the hypergeometric test (i.e., logHG) given the experimental overlap value, each protein's total number of observed experiments, and the total number of experiments.


```r
datScoring <- 
Weighted.matrixModel(SampleDatInput)
```

```
## Joining, by = "UniprotID"
## Joining, by = "UniprotID"
```

```r
head(datScoring)
```

```
##             BPI k    logHG
## 1 A4D1P6~Q12931 2 4.762174
## 2 O00160~P53396 2 3.690590
## 3 O00160~Q9UBT2 2 4.762174
## 4 O00231~P28331 2 3.558201
## 5 O15144~Q9P2D7 2 3.690590
## 6 O43776~P12004 2 4.762174
```

###Assign a confidence score to each instances using classifiers:
The labeled feature matrix can be used as input for Support Vector Machine (SVM) or random forest (RF) classifiers. The classifier then assigns each bait-prey pair a confidence score, indicating the level of support for that pair of proteins to interact. Hyperparameter optimization can also be performed to select a set of parameters that maximizes the model's performance. The RF and the SVM functions provided in this package also computes the areas under the precision-recall (PR) and ROC curve to evalute the performance of the classifer. 

####Import the demo data:

```r
data("testdfClassifier")
head(testdfClassifier)
```

```
##             BPI      Dice   Jaccard   Overlap   Simpson  k     logHG target
## 1 O00232~P35998 0.3478261 0.2105263 0.2105263 1.0000000  4 1.5102503      1
## 2 P18077~P42766 0.7272727 0.5714286 0.5538462 0.9230769 12 3.0401840      1
## 3 Q08722~Q96N66 0.6153846 0.4444444 0.4444444 1.0000000  8 3.9266218      0
## 4 O43684~P42224 0.6363636 0.4666667 0.4083333 0.7000000  7 3.0103311      0
## 5 Q13838~Q14498 0.5882353 0.4166667 0.3571429 0.7142857  5 3.1524199      1
## 6 P22314~P61224 0.3333333 0.2000000 0.1142857 0.4000000  2 0.9471408      0
```

####Run the RF classifier:

```
## named numeric(0)
```

![plot of chunk rfTrain output figure](figure/rfTrain output figure-1.png)

####Output from RF classifier:

```r
#positve score corresponds to  the level of support for the pair of proteins to be true positive
#negaitve score corresponds to the level of support for the pair of proteins to be true negative
head(predidcted_RF)
```

```
##    fulldat[, 1]   Positive  Negative
## 1 O00232~P35998 0.67585256 0.3241474
## 2 P18077~P42766 0.73349942 0.2665006
## 3 Q08722~Q96N66 0.15404286 0.8459571
## 4 O43684~P42224 0.00950000 0.9905000
## 5 Q13838~Q14498 0.61530159 0.3846984
## 6 P22314~P61224 0.07906212 0.9209379
```


####Run the SVM classifier:

```r
#only generate the ROC curve
predidcted_SVM <- 
    svmTrain(testdfClassifier,impute = FALSE,p = 0.3,parameterTuning = TRUE,
        cost = seq(from = 2, to = 10, by = 2),
        gamma = seq(from = 0.01, to = 0.10, by = 0.02),
        kernel = "radial",ncross = 10,
        pr.plot = FALSE, roc.plot = TRUE
    ) 
```

```
## named numeric(0)
```

```
## cost range:
```

```
## [1]  2  4  6  8 10
```

```
## gamma range:
```

```
## [1] 0.01 0.03 0.05 0.07 0.09
```

```
## Setting levels: control = Positive, case = Negative
```

```
## Setting direction: controls < cases
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

####Output from SVM classifier:

```r
#positve score corresponds to  the level of support for the pair of proteins to be true positive
#negaitve score corresponds to the level of support for the pair of proteins to be true negative
head(predidcted_SVM)
```

```
##                   Positive            Negative           
## 1 "O00232~P35998" "0.444645689445843" "0.555354310554157"
## 2 "P18077~P42766" "0.493806183074298" "0.506193816925702"
## 3 "Q08722~Q96N66" "0.406379056468524" "0.593620943531476"
## 4 "O43684~P42224" "0.406951170135447" "0.593048829864553"
## 5 "Q13838~Q14498" "0.446726482567619" "0.553273517432381"
## 6 "P22314~P61224" "0.378309573314996" "0.621690426685004"
```






