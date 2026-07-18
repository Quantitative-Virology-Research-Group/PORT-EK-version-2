## Construction of logistic regression-based models
To examine the prediction power in the classification of the origin of enriched *k*-mers (DNA versus RNA *k*-mers), 1,000 random enriched *k*-mers from each individual pool of enriched *k*-mers were bootstrapped with replacement and divided into a training set (80% of the dataset) and a testing set (20% of the dataset) for logistic regression running on R. The logistic regression model was fitted using the function glm() with the argument family specified as “na.exclude” and ‘‘binomial’’ in the R package ‘‘stats’’ (https://www.r-project.org/). Classifiers associated with the five mentioned features (i.e., “*k*-mer average count”, “*k*-mer RMSE”, “*k*-mer weight”, “subtype *k*-mer count”, and “isolate *k*-mer count”) were independently constructed; the “response” (i.e., the origin of enriched *k*-mers), coupled with a corresponding *k*-mer feature, was thus included in an object of class “formula”, independently. Receiver operating characteristic (ROC) and the area under the curve (AUC) were calculated using the functions multiclass.roc() and auc() in the R package **pROC**, respectively. Given that there is no great variety of the input data for modeling, we have chosen ROC-AUC for evaluating a model's overall ranking ability across all thresholds. The whole procedure was repeated 1,000 times for statistical robustness.

### R code
```
ogistic.Regression.Kmer.average.count <- function(df) {
  # 1. feature selection
  df.kmer.type <- df %>% dplyr::select(kmer.type, kmer.average.count, RMSE, kmer.weight) 
  df.kmer.type$kmer.type <- as.factor(df.kmer.type$kmer.type)
  df.kmer.type$kmer.averaage.count <- as.numeric(df.kmer.type$kmer.average.count)
  df.kmer.type$RMSE <- as.numeric(df.kmer.type$RMSE)
  df.kmer.type$kmer.weight <- as.numeric(df.kmer.type$kmer.weight)
  
  df.kmer.type.b <- sample_n(df.kmer.type, 1000, replace = T)
  
  # 2. create training and test samples
  Sample <- sample(c(TRUE, FALSE), nrow(df.kmer.type.b), replace = TRUE, prob = c(0.8,0.2))
  train <- df.kmer.type.b[Sample,]
  test <- df.kmer.type.b[!Sample,]
  
  # 3. fit the LR model
  model <- glm(kmer.type ~ kmer.average.count, na.action = na.exclude, family = "binomial", data = train)
  
  # 4. prediction
  prediction <- predict(model, test, type = "response")
  
  # 5. create roc curve
  roc <- multiclass.roc(test$kmer.type, prediction)
  
  auc <- auc(roc)
  
  return(auc)
}

Logistic.Regression.Kmer.RMSE <- function(df) {
  # 1. feature selection
  df.kmer.type <- df %>% dplyr::select(kmer.type, kmer.average.count, RMSE, kmer.weight) 
  df.kmer.type$kmer.type <- as.factor(df.kmer.type$kmer.type)
  df.kmer.type$kmer.averaage.count <- as.numeric(df.kmer.type$kmer.average.count)
  df.kmer.type$RMSE <- as.numeric(df.kmer.type$RMSE)
  df.kmer.type$kmer.weight <- as.numeric(df.kmer.type$kmer.weight)
  
  df.kmer.type.b <- sample_n(df.kmer.type, 1000, replace = T)
  
  # 2. create training and test samples
  Sample <- sample(c(TRUE, FALSE), nrow(df.kmer.type.b), replace = TRUE, prob = c(0.8,0.2))
  train <- df.kmer.type.b[Sample,]
  test <- df.kmer.type.b[!Sample,]
  
  # 3. fit the LR model
  model <- glm(kmer.type ~ RMSE, na.action = na.exclude, family = "binomial", data = train)
  
  # 4. prediction
  prediction <- predict(model, test, type = "response")
  
  # 5. create roc curve
  roc <- multiclass.roc(test$kmer.type, prediction)
  
  auc <- auc(roc)
  
  return(auc)
}

Logistic.Regression.Kmer.weight <- function(df) {
  # 1. feature selection
  df.kmer.type <- df %>% dplyr::select(kmer.type, kmer.average.count, RMSE, kmer.weight) 
  df.kmer.type$kmer.type <- as.factor(df.kmer.type$kmer.type)
  df.kmer.type$kmer.averaage.count <- as.numeric(df.kmer.type$kmer.average.count)
  df.kmer.type$RMSE <- as.numeric(df.kmer.type$RMSE)
  df.kmer.type$kmer.weight <- as.numeric(df.kmer.type$kmer.weight)
  
  df.kmer.type.b <- sample_n(df.kmer.type, 1000, replace = T)
  
  # 2. create training and test samples
  Sample <- sample(c(TRUE, FALSE), nrow(df.kmer.type.b), replace = TRUE, prob = c(0.8,0.2))
  train <- df.kmer.type.b[Sample,]
  test <- df.kmer.type.b[!Sample,]
  
  # 3. fit the LR model
  model <- glm(kmer.type ~ kmer.weight, na.action = na.exclude, family = "binomial", data = train)
  
  # 4. prediction
  prediction <- predict(model, test, type = "response")
  
  # 5. create roc curve
  roc <- multiclass.roc(test$kmer.type, prediction)
  
  auc <- auc(roc)
  
  return(auc)
}

Logistic.Regression.Kmer.count <- function(df) {
  # 1. feature selection
  df.kmer.type <- df %>% dplyr::select(kmer.type, kmer.average.count, RMSE, kmer.weight) 
  df.kmer.type$kmer.type <- as.factor(df.kmer.type$kmer.type)
  df.kmer.type$kmer.count <- as.numeric(df.kmer.type$kmer.count)
  
  df.kmer.type.b <- sample_n(df.kmer.type, 1000, replace = T)
  
  # 2. create training and test samples
  Sample <- sample(c(TRUE, FALSE), nrow(df.kmer.type.b), replace = TRUE, prob = c(0.8,0.2))
  train <- df.kmer.type.b[Sample,]
  test <- df.kmer.type.b[!Sample,]
  
  # 3. fit the LR model
  model <- glm(kmer.type ~ kmer.count, na.action = na.exclude, family = "binomial", data = train)
  
  # 4. prediction
  prediction <- predict(model, test, type = "response")
  
  # 5. create roc curve
  roc <- multiclass.roc(test$kmer.type, prediction)
  
  auc <- auc(roc)
  
  return(auc)
}

Logistic.Regression.subtype.kmer.count <- function(df.dna, df.rna) {
  # 1. feature selection
  df <- dplyr::bind_rows(df.dna, df.rna)
  
  df.kmer.type <- df %>% dplyr::select(type, kmer.count) 
  df.kmer.type$type <- as.factor(df.kmer.type$type)
  df.kmer.type$kmer.count <- as.numeric(df.kmer.type$kmer.count)
  
  df.kmer.type.b <- sample_n(df.kmer.type, 1000, replace = T)
  
  # 2. create training and test samples
  Sample <- sample(c(TRUE, FALSE), nrow(df.kmer.type.b), replace = TRUE, prob = c(0.8,0.2))
  train <- df.kmer.type.b[Sample,]
  test <- df.kmer.type.b[!Sample,]
  
  # 3. fit the LR model
  model <- glm(type ~ kmer.count, na.action = na.exclude, family = "binomial", data = train)
  
  # 4. prediction
  prediction <- predict(model, test, type = "response")
  
  # 5. create roc curve
  roc <- multiclass.roc(test$type, prediction)
  
  auc <- auc(roc)
  
  return(auc)
}
```
