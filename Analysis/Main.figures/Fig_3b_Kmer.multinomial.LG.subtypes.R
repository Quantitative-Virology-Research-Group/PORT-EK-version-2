Logistic.Regression.subtype.kmer.average.count <- function(df) {
  # 1. feature selection
  df.subtype <- df %>% dplyr::select(group, kmer.average.count, RMSE, kmer.weight) 
  df.subtype$group <- as.factor(df.subtype$group)
  df.subtype$kmer.average.count <- as.numeric(df.subtype$kmer.average.count)
  df.subtype$RMSE <- as.numeric(df.subtype$RMSE)
  df.subtype$kmer.weight <- as.numeric(df.subtype$kmer.weight)
  
  df.subtype.b <- sample_n(df.subtype, 1000, replace = T)
  
  # 2. create training and test samples
  Sample <- sample(c(TRUE, FALSE), nrow(df.subtype.b), replace = TRUE, prob = c(0.8,0.2))
  train <- df.subtype.b[Sample,]
  test <- df.subtype.b[!Sample,]
  
  # 3. fit the LR model
  model <- multinom(group ~ kmer.average.count, data = train)
  
  # 4. prediction
  prediction <- as.data.frame(predict(model, test, type = "probs"))
  
  value.A <- mean(prediction$A)
  value.B <- mean(prediction$B)
  value.C <- mean(prediction$C)
  value.D <- mean(prediction$D)
  value.R <- mean(prediction$R)
  
  df.out <- data.frame(value = c(value.A, value.B, value.C, value.D, value.R), subtype = c("A", "B", "C", "D", "R"))
  
  return(df.out)
}

Logistic.Regression.subtype.kmer.RMSE <- function(df) {
  # 1. feature selection
  df.subtype <- df %>% dplyr::select(group, kmer.average.count, RMSE, kmer.weight) 
  df.subtype$group <- as.factor(df.subtype$group)
  df.subtype$kmer.average.count <- as.numeric(df.subtype$kmer.average.count)
  df.subtype$RMSE <- as.numeric(df.subtype$RMSE)
  df.subtype$kmer.weight <- as.numeric(df.subtype$kmer.weight)
  
  df.subtype.b <- sample_n(df.subtype, 1000, replace = T)
  
  # 2. create training and test samples
  Sample <- sample(c(TRUE, FALSE), nrow(df.subtype.b), replace = TRUE, prob = c(0.8,0.2))
  train <- df.subtype.b[Sample,]
  test <- df.subtype.b[!Sample,]
  
  # 3. fit the LR model
  model <- multinom(group ~ RMSE, data = train)
  
  # 4. prediction
  prediction <- as.data.frame(predict(model, test, type = "probs"))
  
  value.A <- mean(prediction$A)
  value.B <- mean(prediction$B)
  value.C <- mean(prediction$C)
  value.D <- mean(prediction$D)
  value.R <- mean(prediction$R)
  
  df.out <- data.frame(value = c(value.A, value.B, value.C, value.D, value.R), subtype = c("A", "B", "C", "D", "R"))
  
  return(df.out)
}

Logistic.Regression.subtype.kmer.weight <- function(df) {
  # 1. feature selection
  df.subtype <- df %>% dplyr::select(group, kmer.average.count, RMSE, kmer.weight) 
  df.subtype$group <- as.factor(df.subtype$group)
  df.subtype$kmer.average.count <- as.numeric(df.subtype$kmer.average.count)
  df.subtype$RMSE <- as.numeric(df.subtype$RMSE)
  df.subtype$kmer.weight <- as.numeric(df.subtype$kmer.weight)
  
  df.subtype.b <- sample_n(df.subtype, 1000, replace = T)
  
  # 2. create training and test samples
  Sample <- sample(c(TRUE, FALSE), nrow(df.subtype.b), replace = TRUE, prob = c(0.8,0.2))
  train <- df.subtype.b[Sample,]
  test <- df.subtype.b[!Sample,]
  
  # 3. fit the LR model
  model <- multinom(group ~ kmer.weight, data = train)
  
  # 4. prediction
  prediction <- as.data.frame(predict(model, test, type = "probs"))
  
  value.A <- mean(prediction$A)
  value.B <- mean(prediction$B)
  value.C <- mean(prediction$C)
  value.D <- mean(prediction$D)
  value.R <- mean(prediction$R)
  
  df.out <- data.frame(value = c(value.A, value.B, value.C, value.D, value.R), subtype = c("A", "B", "C", "D", "R"))
  
  return(df.out)
}

Logistic.Regression.subtype.kmer.count <- function(df) {
  # 1. feature selection
  df.subtype <- df %>% dplyr::select(group, kmer.count) 
  df.subtype$group <- as.factor(df.subtype$group)
  df.subtype$kmer.count <- as.numeric(df.subtype$kmer.count)
  
  df.subtype.b <- sample_n(df.subtype, 1000, replace = T)
  
  # 2. create training and test samples
  Sample <- sample(c(TRUE, FALSE), nrow(df.subtype.b), replace = TRUE, prob = c(0.8,0.2))
  train <- df.subtype.b[Sample,]
  test <- df.subtype.b[!Sample,]
  
  # 3. fit the LR model
  model <- multinom(group ~ kmer.count, data = train)
  
  # 4. prediction
  prediction <- as.data.frame(predict(model, test, type = "probs"))
  
  value.A <- mean(prediction$A)
  value.B <- mean(prediction$B)
  value.C <- mean(prediction$C)
  value.D <- mean(prediction$D)
  value.R <- mean(prediction$R)
  
  df.out <- data.frame(value = c(value.A, value.B, value.C, value.D, value.R), subtype = c("A", "B", "C", "D", "R"))
  
  return(df.out)
}

Logistic.Regression.sub.subtype.kmer.count <- function(df) {
  # 1. feature selection
  df.subtype <- df %>% dplyr::select(type, kmer.count, subtype) 
  df.subtype$type <- as.factor(df.subtype$type)
  df.subtype$kmer.count <- as.numeric(df.subtype$kmer.count)
  df.subtype$subtype <- as.factor(df.subtype$subtype)
  
  df.subtype.b <- sample_n(df.subtype, 1000, replace = T)
  
  # 2. create training and test samples
  Sample <- sample(c(TRUE, FALSE), nrow(df.subtype.b), replace = TRUE, prob = c(0.8,0.2))
  train <- df.subtype.b[Sample,]
  test <- df.subtype.b[!Sample,]
  
  # 3. fit the LR model
  model <- multinom(subtype ~ kmer.count, data = train)
  
  # 4. prediction
  prediction <- as.data.frame(predict(model, test, type = "probs"))
   
   value.A <- mean(prediction$A)
   value.B <- mean(prediction$B)
   value.C <- mean(prediction$C)
   value.D <- mean(prediction$D)
   value.R <- mean(prediction$R)
   
   df.out <- data.frame(value = c(value.A, value.B, value.C, value.D, value.R), subtype = c("A", "B", "C", "D", "R"))
   
   return(df.out)
}

combine.dataframe.list <- function(df) {
  df.split <- lapply(df, as.data.frame.list)
  df.merge <- do.call(what = "rbind", df.split)
  df.merge[is.na(df.merge)] <- 0
  
  return(df.merge)
}
