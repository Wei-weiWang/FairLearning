
##Glm law school dataset with vairbles used from GLM paper(The one we use)
lawschool = read.csv('lawschool(glm).csv')
# we choose the covariates fair glm uses
lawschool <- lawschool[ , c('age', 'fam_inc', 'lsat', 'ugpa', 'gender', 'fulltime', 'zgpa', 'race1')]
# Following GLM paper, we remove NAs
lawschool = na.omit(lawschool)


#gender
unique(lawschool$gender) # "female" "male"   "" We delete two rows where gender==""
lawschool <- lawschool[-which(lawschool$gender == ""), ]
unique(lawschool$gender)
gender_num <- as.integer(lawschool$gender == "female")+1  # 2 if female, 1 is male
lawschool$gender  <- factor(gender_num)
#race
lawschool <- lawschool[-which(lawschool$race1 == ""), ]
lawschool$race1 <- factor(lawschool$race1, levels = c("white", "black", "hisp","asian", "other"), labels = 1:5)
#fulltime
lawschool$fulltime <- factor(lawschool$fulltime)


set.seed(50)
# Pre-allocate result vectors (30 iterations)
n_iter <- 30
QR_Pinball_T <- numeric(n_iter)
QR_KS_T      <- numeric(n_iter)
QR_Pinball_L <- numeric(n_iter)
QR_KS_L      <- numeric(n_iter)

QRI_Pinball_T <- numeric(n_iter)
QRI_KS_T      <- numeric(n_iter)
QRI_Pinball_L <- numeric(n_iter)
QRI_KS_L      <- numeric(n_iter)

QRPA_Pinball_T <- numeric(n_iter)
QRPA_KS_T      <- numeric(n_iter)
QRPA_Pinball_L <- numeric(n_iter)
QRPA_KS_L      <- numeric(n_iter)

qu <- 0.75
sample_size <- floor(0.8 * nrow(lawschool))  # 80% training, 20% test

for (k in 1:n_iter) {
  cat("Iteration:", k, "\n")
  
  # Create train/test split
  train_indices <- sample(seq_len(nrow(lawschool)), size = sample_size)
  write.table(lawschool[train_indices,], "lawtrain75.csv", sep = ",", 
              row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE)
  write.table(lawschool[-train_indices,], "lawtest75.csv", sep = ",", 
              row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE)
  
  train_data <- lawschool[train_indices, ]
  test_data  <- lawschool[-train_indices, ] 
  
  # Separate covariates and response and split by group (by race)
  idx_cov = 6
  idx_y = 7
  idx_race = 8
  
  
  #--------------------------
  # Plain Quantile Regression
  #--------------------------
  quantile_model <- rq(lawschool[train_indices, idx_y] ~ ., 
                       data = lawschool[train_indices, c(1:idx_cov, idx_race)], 
                       tau = qu) 
  
  QR_Pinball_L[k] <- quantile_loss(lawschool[train_indices, idx_y], quantile_model$fitted.values, qu) 
  QR_KS_L[k] <- f_ks(quantile_model$fitted.values, lawschool[train_indices, idx_race]) 
  
  predicted_y <- predict(quantile_model, newdata = lawschool[-train_indices, c(1:idx_cov, idx_race)])
  QR_Pinball_T[k] <- quantile_loss(lawschool[-train_indices, idx_y], predicted_y, qu) 
  QR_KS_T[k] <- f_ks(predicted_y, lawschool[-train_indices, idx_race])
  
  
  
  # Fair regression methods
  split_by_race_train <- split(train_data, train_data[,idx_race])
  
  split_by_race_test <- split(test_data, test_data[,idx_race])
  
  # For train
  Ytrain = c()
  Strain = c()
  
  LYpred = c()
  rankYLF = c()
  
  # For test  
  Ytest = c()
  Stest = c()
  
  rankYTF = c()
  
  for (i in 1:length(unique(train_data[,idx_race])) ) {
    train_cov_split = split_by_race_train[[i]][1:idx_cov]
    Ytrain_split = split_by_race_train[[i]][idx_y]
    Ytrain = c(Ytrain, Ytrain_split[[1]])
    Strain = c(Strain, rep(i, nrow(train_cov_split)))
    
    
    test_cov_split = split_by_race_test[[i]][1:idx_cov]
    Ytest_split = split_by_race_test[[i]][idx_y]
    Ytest = c(Ytest, Ytest_split[[1]])
    Stest = c(Stest, rep(i, nrow(test_cov_split)))
    
    model = rq(Ytrain_split[[1]] ~ ., data = train_cov_split, tau = qu)
    LYpred_split = predict(model, train_cov_split)
    TYpred_split = predict(model, test_cov_split)

    
    #print( quantile_loss( split_by_race_train[[i]][idx_y][[1]], LYpred_split, 0.75)  )
    #print( quantile_loss( split_by_race_test[[i]][idx_y][[1]], TYpred_split, 0.75)  )
    
    LYpred = c(LYpred, LYpred_split)
    
    rankYLF_split = rank(LYpred_split, ties.method = "max")/length(LYpred_split) 
    rankYLF = c(rankYLF, rankYLF_split) 
    
    rankYTF_split = sapply(TYpred_split, function(y) sum(LYpred_split <= y)/length(LYpred_split))
    rankYTF = c(rankYTF, rankYTF_split)
    
    
  }
  names(rankYLF) <- NULL
  
  #--------------------------
  # Ispline Method
  #--------------------------
  
  
  # Fit the Ispline model on training data
  ispline_fit <- fit_ispline_model_Quantile(rankYLF, Ytrain, tau = qu, 
                                   knots = c(0.25, 0.5, 0.75), degree = 2)  
  train_pred <- ispline_predict(rankYLF, ispline_fit$coeff, 
                                knots = c(0.25, 0.5, 0.75), degree = 2)  
  QRI_Pinball_L[k] <- quantile_loss(Ytrain, train_pred, qu) 
  QRI_KS_L[k] <- f_ks(train_pred, Strain) 
  
  
  test_pred <- ispline_predict(rankYTF, ispline_fit$coeff, 
                               knots = c(0.25, 0.5, 0.75), degree = 2)
  QRI_Pinball_T[k] <- quantile_loss(Ytest, test_pred, qu)
  QRI_KS_T[k] <- f_ks(test_pred, Stest) 
  
  #--------------------------
  # PAVA Method
  #--------------------------
  pava_fit <- fit_pava_model_Quantile(rankYLF, Ytrain, tau = qu, solverinput = weighted.fractile ) 
  gpava_pred_train <- gpava(rankYLF, Ytrain, solver = weighted.fractile, 
                            ties = "secondary", p = qu)$x 
  QRPA_Pinball_L[k] <- quantile_loss(Ytrain, gpava_pred_train, qu) 
  QRPA_KS_L[k] <- f_ks(gpava_pred_train, Strain)
  
  # For test predictions, use our wrapper pava_predict:
  test_pred_pava <- pava_predict(pava_fit, rankYTF) 
  QRPA_Pinball_T[k] <- quantile_loss(Ytest, test_pred_pava, qu) 
  QRPA_KS_T[k] <- f_ks(test_pred_pava, Stest) 
}

#--------------------------
# Summarize Results
#--------------------------
cat("Means:\n")
cat("QR_Pinball_T:", mean(QR_Pinball_T), "\n")
cat("QR_KS_T:", mean(QR_KS_T), "\n")
cat("QR_Pinball_L:", mean(QR_Pinball_L), "\n")
cat("QR_KS_L:", mean(QR_KS_L), "\n")
cat("QRI_Pinball_T:", mean(QRI_Pinball_T), "\n")
cat("QRI_KS_T:", mean(QRI_KS_T), "\n")
cat("QRI_Pinball_L:", mean(QRI_Pinball_L), "\n")
cat("QRI_KS_L:", mean(QRI_KS_L), "\n")
cat("QRPA_Pinball_T:", mean(QRPA_Pinball_T), "\n")
cat("QRPA_KS_T:", mean(QRPA_KS_T), "\n")
cat("QRPA_Pinball_L:", mean(QRPA_Pinball_L), "\n")
cat("QRPA_KS_L:", mean(QRPA_KS_L), "\n")

cat("\nStandard Errors:\n")
cat("QR_Pinball_T:", sd(QR_Pinball_T)/sqrt(n_iter), "\n")
cat("QR_KS_T:", sd(QR_KS_T)/sqrt(n_iter), "\n")
cat("QR_Pinball_L:", sd(QR_Pinball_L)/sqrt(n_iter), "\n")
cat("QR_KS_L:", sd(QR_KS_L)/sqrt(n_iter), "\n")
cat("QRI_Pinball_T:", sd(QRI_Pinball_T)/sqrt(n_iter), "\n")
cat("QRI_KS_T:", sd(QRI_KS_T)/sqrt(n_iter), "\n")
cat("QRI_Pinball_L:", sd(QRI_Pinball_L)/sqrt(n_iter), "\n")
cat("QRI_KS_L:", sd(QRI_KS_L)/sqrt(n_iter), "\n")
cat("QRPA_Pinball_T:", sd(QRPA_Pinball_T)/sqrt(n_iter), "\n")
cat("QRPA_KS_T:", sd(QRPA_KS_T)/sqrt(n_iter), "\n")
cat("QRPA_Pinball_L:", sd(QRPA_Pinball_L)/sqrt(n_iter), "\n")
cat("QRPA_KS_L:", sd(QRPA_KS_L)/sqrt(n_iter), "\n")
