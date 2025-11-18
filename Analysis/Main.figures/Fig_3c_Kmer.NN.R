Neural.Network.subtype.kmer.average.count <- function(df) {
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
  model <- neuralnet::neuralnet(group ~ kmer.average.count, data = train)
  
  # 4. prediction
  prediction <- predict(model, test)
  
  # 5. calculate accuracy
  ## To check the accuracy, we have to first convert actual categorical values  into numerical ones and compare them with predicted values. As a result, we will receive a list of boolean values. 
  
  ##We can use the `sum` function to find the number of `TRUE` values and divide it  by the total number of samples to get the accuracy. 
  
  check = as.numeric(test$group) == max.col(prediction)
  accuracy = as.numeric(sum(check)/nrow(test))*100
  
  return(accuracy)
}

Neural.Network.subtype.kmer.RMSE <- function(df) {
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
  model <- neuralnet::neuralnet(group ~ RMSE, data = train)
  
  # 4. prediction
  prediction <- predict(model, test)
  
  # 5. calculate accuracy
  ## To check the accuracy, we have to first convert actual categorical values  into numerical ones and compare them with predicted values. As a result, we will receive a list of boolean values. 
  
  ##We can use the `sum` function to find the number of `TRUE` values and divide it  by the total number of samples to get the accuracy. 
  
  check = as.numeric(test$group) == max.col(prediction)
  accuracy = as.numeric(sum(check)/nrow(test))*100
  
  return(accuracy)
}

Neural.Network.subtype.kmer.weight <- function(df) {
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
  model <- neuralnet::neuralnet(group ~ kmer.weight, data = train)
  
  # 4. prediction
  prediction <- predict(model, test)
  
  # 5. calculate accuracy
  ## To check the accuracy, we have to first convert actual categorical values  into numerical ones and compare them with predicted values. As a result, we will receive a list of boolean values. 
  
  ##We can use the `sum` function to find the number of `TRUE` values and divide it  by the total number of samples to get the accuracy. 
  
  check = as.numeric(test$group) == max.col(prediction)
  accuracy = as.numeric(sum(check)/nrow(test))*100
  
  return(accuracy)
}

Neural.Network.subtype.kmer.count <- function(df) {
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
  model <- neuralnet::neuralnet(group ~ kmer.count, data = train)
  
  # 4. prediction
  prediction <- predict(model, test)
  
  # 5. calculate accuracy
  ## To check the accuracy, we have to first convert actual categorical values  into numerical ones and compare them with predicted values. As a result, we will receive a list of boolean values. 
  
  ##We can use the `sum` function to find the number of `TRUE` values and divide it  by the total number of samples to get the accuracy. 
  
  check = as.numeric(test$group) == max.col(prediction)
  accuracy = as.numeric(sum(check)/nrow(test))*100
  
  return(accuracy)
}

Neural.Network.sub.subtype.kmer.count <- function(df) {
  # 1. feature selection
  df.subtype <- df %>% dplyr::select(subtype, kmer.count) 
  df.subtype$subtype <- as.factor(df.subtype$subtype)
  df.subtype$kmer.count <- as.numeric(df.subtype$kmer.count)

  df.subtype.b <- sample_n(df.subtype, 1000, replace = T)
  
  # 2. create training and test samples
  Sample <- sample(c(TRUE, FALSE), nrow(df.subtype.b), replace = TRUE, prob = c(0.8,0.2))
  train <- df.subtype.b[Sample,]
  test <- df.subtype.b[!Sample,]
  
  # 3. fit the LR model
  model <- neuralnet::neuralnet(subtype ~ kmer.count, data = train)
  
  # 4. prediction
  prediction <- predict(model, test)
  
  # 5. calculate accuracy
  check = as.numeric(test$group) == max.col(prediction)
  accuracy = as.numeric(sum(check)/nrow(test))*100
  
  return(accuracy)
}
