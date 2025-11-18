Logistic.Regression.Kmer.average.count <- function(df) {
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
  df.kmer.type$kmer.kmer.count <- as.numeric(df.kmer.type$kmer.count)
  
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
