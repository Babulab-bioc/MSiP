## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
#Jaccard coefficient
datScoring <- 
    jaccardCoefficient(SampleDatInput)
head(datScoring)

#Simpson coefficient
datScoring <- 
    simpsonCoefficient(SampleDatInput)
head(datScoring)

#Overlap score
datScoring <- 
    simpsonCoefficient(SampleDatInput)
head(datScoring)

## -----------------------------------------------------------------------------
datScoring <- 
Weighted.matrixModel(SampleDatInput)
head(datScoring)

## -----------------------------------------------------------------------------
data("testdfClassifier")
testdfClassifier <- #smaller data size 
dplyr::sample_n(testdfClassifier, 1000)
head(testdfClassifier)


## -----------------------------------------------------------------------------
#SVM
predidcted_SV <- 
svmTrain(testdfClassifier,impute = FALSE,p = 0.3, parameterTuning = FALSE,
cost = seq(from = 2, to = 5, by = 1),
gamma = seq(from = 0.01, to = 0.10, by = 0.02),
kernel = "radial",ncross = 5,
pr.plot = TRUE, roc.plot = TRUE
)
head(predidcted_SV)

## -----------------------------------------------------------------------------
#RF
predidcted_RF <- 
rfTrain(testdfClassifier,impute = FALSE, p = 0.3, parameterTuning = FALSE,
mtry  = seq(from = 1, to = 5, by = 1),
min_node_size = seq(from = 1, to = 5, by = 1),
splitrule =c("gini"),metric = "Accuracy",
resampling.method = "repeatedcv",iter = 5,repeats = 5,
pr.plot = TRUE, roc.plot = TRUE
)
head(predidcted_RF)

