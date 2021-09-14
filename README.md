<!-- badges: start -->
[![CRAN Version](https://www.r-pkg.org/badges/version/MSiP)](https://cran.r-project.org/package=MSiP)
[![Downloads from the RStudio CRAN mirror](https://cranlogs.r-pkg.org/badges/MSiP)](https://cranlogs.r-pkg.org/badges/MSiP)
<!-- badges: end -->

# Mass Spectrometry interaction Prediction (MSiP)
The MSiP is a computational approach to predict protein-protein interactions from largescale affinity purification mass spectrometry (AP-MS) data. This approach includes both spoke and matrix models for interpreting AP-MS data in a network context. The 'spoke' model considers only bait-prey interactions, whereas the 'matrix' model assumes that each of the identified proteins (baits and prey) in a given AP-MS experiment interacts with each of the others. The spoke model has a high false-negative rate, whereas the matrix model has a high false-positive rate. Although, both statistical models have merits, a combination of both models has shown to increase the performance of machine learning classifiers in terms of their capabilities in discrimination between true and false positive interactions [Drew et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5488662).  

### Installation from cran:
```{r}
install.packages('MSiP')
library(MSiP)
```
### Installation from github:
```{r}
install_github("mrbakhsh/MSiP")
library(MSiP)
```

### Sample Data Description:
A demo AP-MS proteomics dataset is provided in this package to guide the users about data structure.
```{r}
data("SampleDatInput")
head(SampleDatInput)
```
### Scoring based on "spoke-model":
Comparative Proteomic Analysis Software Suite (CompPASS) is a robust statistical scoring scheme for assigning confidence scores to bait-prey interactions [Sowa et al., 2009](https://www.cell.com/cell/fulltext/S0092-8674(09)00503-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867409005030%3Fshowall%3Dtrue). The output from CompPASS scoring includes Z-score, S-score, D-score, WD-score and other features.  

```{r}
datScoring <- 
    cPASS(SampleDatInput)
```
### Scoring based on "matrix-model":
The Dice coefficient was first applied by [Zhang et al., 2008](https://academic.oup.com/bioinformatics/article/24/7/979/296061) to score interaction between all identified proteins (baits and preys) in a given AP-MS expriment.
```{r}
datScoring <- 
    diceCoefficient(SampleDatInput)
```

Alternatively, Jaccard, Simpson, and Overlap scores can be used to score the interaction between all the identified proteins in a given AP-MS experiment. 
```{r}
#Jaccard coefficient
datScoring <- 
    jaccardCoefficient(SampleDatInput)

#Simpson coefficient
datScoring <- 
    simpsonCoefficient(SampleDatInput)

#Overlap score
datScoring <- 
    overlapCoefficient(SampleDatInput)
```

Finally, a weighted matrix model [Drew et al., 2017](https://www.embopress.org/doi/full/10.15252/msb.20167490) can also be employed to score interactions between identified proteins in a given AP-MS experiment. The output of the weighted matrix model includes the number of experiments for which the pair of proteins is co-purified (i.e., k) and $-1$*log(P-value) of the hypergeometric test (i.e., logHG) given the experimental overlap value, each protein's total number of observed experiments, and the total number of experiments.

```{r}
datScoring <- 
Weighted.matrixModel(SampleDatInput)
```

### Assign a confidence score to each instances using classifiers:
The labeled feature matrix can be used as input for Support Vector Machine (SVM) or Random Forest (RF) classifiers. The classifier then assigns each bait-prey pair a confidence score, indicating the level of support for that pair of proteins to interact. Hyperparameter optimization can also be performed to select a set of parameters that maximizes the model's performance. The RF and the SVM functions provided in this package also computes the areas under the precision-recall (PR) and ROC curve to evalute the performance of the classifier. 

#### Import the demo data:
```{r}
data("testdfClassifier")
head(testdfClassifier)

```

#### Run the RF classifier:
```{r rfTrain output figure, echo=FALSE, fig.height=4, fig.width=5, message=FALSE, warning=FALSE, paged.print=FALSE}
#only generate the pr.curve
predidcted_RF <- 
    rfTrain(testdfClassifier,impute = FALSE, p = 0.3, parameterTuning = FALSE,
        mtry  = seq(from = 1, to = 5, by = 1),
        min_node_size = seq(from = 1, to = 5, by = 1),
        splitrule =c("gini"),metric = "Accuracy",
        resampling.method = "repeatedcv",iter = 5,repeats = 5,
        pr.plot = TRUE, roc.plot = FALSE
    )
```

#### Run the SVM classifier:
```{r}
#only generate the ROC curve
predidcted_SVM <- 
    svmTrain(testdfClassifier,impute = FALSE,p = 0.3,parameterTuning = TRUE,
        cost = seq(from = 2, to = 10, by = 2),
        gamma = seq(from = 0.01, to = 0.10, by = 0.02),
        kernel = "radial",ncross = 10,
        pr.plot = FALSE, roc.plot = TRUE
    ) 
```




