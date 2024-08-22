### Load Packages ----
library(tidyverse)
library(yardstick)
library(skimr)
library(caret)

### Load Data ----
result_data <- read_csv("data/result.csv")

### Analysis ----
calc_metrics <- function(dat) {
  accuracy <- c()
  sensitivity <- c()
  specificity <- c()
  ppv <- c()
  npv <- c()
  mcc <- c()
  for (i in 5:11) {
    x <- confusionMatrix(factor(as.double(unlist(dat[i]))),
                         factor(dat$label),
                         positive = "1")
    y <- mcc(dat,
             truth = factor(dat$label),
             estimate = factor(unlist(dat[i])))
    accuracy <- c(accuracy, x[[3]][[1]])
    sensitivity <- c(sensitivity, x[[4]][[1]])
    specificity <- c(specificity, x[[4]][[2]])
    ppv <- c(ppv, x[[4]][[3]])
    npv <- c(npv, x[[4]][[4]])
    mcc <- c(mcc, y[[3]])
    }
  return(tibble(Tool = c("PolyPhen-2 HumDiv", "PolyPhen-2 HumVar", 
                         "MutationTaster2021", "Align GVGD",
                         "REVEL", "CADD", "FATHMM"),
                Accuracy = accuracy,
                Sensitivity = sensitivity,
                Specificity = specificity,
                PPV = ppv,
                NPV = npv,
                MCC = mcc))
}

metrics <- calc_metrics(result_data)

# save result
save(metrics, file = "data/metrics.rds")
