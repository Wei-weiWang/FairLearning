
# Law School Dataset
library(fairml)

#https://www.jstor.org/stable/pdf/40040209.pdf cite 1998 . Linda F. Wightman, LSAC National Longitudinal Bar Passage Study
#They used the dataset to analyze.
#Description from 1998. Linda: This study presents national longitudinal bar passage data gathered from the class that started law school in fall 1991.
#Description from Law School Admission Council data: Survey among students attending law school in the U.S. in 1991. Thus, we believe they are the same data.
lawschool = as.data.frame(law.school.admissions)

# Check NA
any(is.na(lawschool)) 

# Vairbles used
# conformal fair paper uses covariates as 'race', 'cluster', 'lsat', 'zfygpa', 'zgpa', 'fulltime', 'fam_inc', 'age', 'gender'(sensitive feature)
# Fair GLM uses 'race1': Cat, 'gender': Cat, 'age': Num, 'fam_inc': Num, 'fulltime': Cat, 'zgpa': Num, 'ugpa': Num, 'lsat': Num
# We use 'race1': Cat, 'gender': Cat(Sensitive), 'age': Num, 'fam_inc': Num, 'fulltime': Cat, 'ugpa': Num(response), 'lsat': Num
# Here, we follow the variable logic of conformal fair paper but use the dataset from R.
lawschool = lawschool[,c(-2, -3, -11)]
# Change positions of columns
lawschool <- lawschool[ , c('age', 'fam_inc', 'lsat', 'race1', 'cluster', 'fulltime', 'ugpa', 'gender')]

gender_num <- as.integer(lawschool$gender == "female")+1  # 2 if female, 1 is male
lawschool$gender  <- factor(gender_num)

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
  
  train_data <- lawschool[train_indices, ]
  test_data  <- lawschool[-train_indices, ] 
  
  # Separate covariates and response and split by group (by race)
  idx_cov = 6
  idx_y = 7
  idx_gender = 8
  
  
  #--------------------------
  # Plain Quantile Regression
  #--------------------------
  quantile_model <- rq(lawschool[train_indices, idx_y] ~ ., 
                       data = lawschool[train_indices, c(1:idx_cov, idx_gender)], 
                       tau = qu) 
  
  QR_Pinball_L[k] <- quantile_loss(lawschool[train_indices, idx_y], quantile_model$fitted.values, qu) 
  QR_KS_L[k] <- f_ks(quantile_model$fitted.values, lawschool[train_indices, idx_gender]) 
  
  predicted_y <- predict(quantile_model, newdata = lawschool[-train_indices, c(1:idx_cov, idx_gender)])
  QR_Pinball_T[k] <- quantile_loss(lawschool[-train_indices, idx_y], predicted_y, qu) 
  QR_KS_T[k] <- f_ks(predicted_y, lawschool[-train_indices, idx_gender])
  
  
  
  # Fair regression methods
  split_by_race_train <- split(train_data, train_data$gender)
  
  split_by_race_test <- split(test_data, test_data$gender)
  
  # For train
  Ytrain = c()
  Strain = c()
  
  LYpred = c()
  rankYLF = c()
  
  # For test  
  Ytest = c()
  Stest = c()
  
  rankYTF = c()
  
  for (i in 1:length(unique(train_data$gender)) ) {
    train_cov_split = split_by_race_train[[i]][1:idx_cov]
    Ytrain_split = split_by_race_train[[i]][idx_y]
    Ytrain = c(Ytrain, Ytrain_split[[1]])
    Strain = c(Strain, rep(i, nrow(train_cov_split)))
    
    
    test_cov_split = split_by_race_test[[i]][1:idx_cov]
    Ytest_split = split_by_race_test[[i]][idx_y]
    Ytest = c(Ytest, Ytest_split[[1]])
    Stest = c(Stest, rep(i, nrow(test_cov_split)))
    
    LYpred_split = predict(rq(Ytrain_split[[1]] ~ ., data = train_cov_split, tau = qu), train_cov_split)
    TYpred_split = predict(rq(Ytrain_split[[1]] ~ ., data = train_cov_split, tau = qu), test_cov_split)

    
    #print( quantile_loss( split_by_race_train[[i]][idx_y][[1]], LYpred_split, 0.75)  )
    #print( quantile_loss( split_by_race_test[[i]][idx_y][[1]], TYpred_split, 0.75)  )
    
    LYpred = c(LYpred, LYpred_split)
    
    rankYLF_split = rank(LYpred_split, ties.method = "max")/length(LYpred_split) 
    rankYLF = c(rankYLF, rankYLF_split) 
    
    rankYTF_split = sapply(TYpred_split, function(y) sum(LYpred_split <= y)/length(LYpred_split))
    rankYTF = c(rankYTF, rankYTF_split)
    
    
  }
  
  
  #--------------------------
  # Ispline Method
  #--------------------------
  
  
  # Fit the Ispline model on training data
  ispline_fit <- fit_ispline_model(rankYLF, Ytrain, tau = qu, 
                                   knots = c(0.2, 0.5, 0.75), degree = 2)  
  train_pred <- ispline_predict(rankYLF, ispline_fit$coeff, 
                                knots = c(0.2, 0.5, 0.75), degree = 2)  
  QRI_Pinball_L[k] <- quantile_loss(Ytrain, train_pred, qu) 
  QRI_KS_L[k] <- f_ks(train_pred, Strain) 
  
  
  test_pred <- ispline_predict(rankYTF, ispline_fit$coeff, 
                               knots = c(0.2, 0.5, 0.75), degree = 2)
  QRI_Pinball_T[k] <- quantile_loss(Ytest, test_pred, qu)
  QRI_KS_T[k] <- f_ks(test_pred, Stest) 
  
  #--------------------------
  # PAVA Method
  #--------------------------
  pava_fit <- fit_pava_model(rankYLF, Ytrain, tau = qu, solverinput = weighted.fractile ) 
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

