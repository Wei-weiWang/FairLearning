hrs = as.data.frame(health.retirement)#https://search.r-project.org/CRAN/refmans/fairml/html/health.retirement.html
hrs <- na.omit(hrs) 

# Categorical variables
nummar = ifelse(hrs$marriage == "Married/Partner", 1, 2)
hrs$marriage = factor(nummar)

numgender = ifelse(hrs$gender == 'Female', 1, 2)
hrs$gender = factor(numgender)

hrs$race.ethnicity =  factor(hrs$race.ethnicity, levels = c('NHW','NHB','Hispanic','Other'), labels = 1:4)
hrs$race =  factor(hrs$race, levels = c('Black','White','Other'), labels = 1:3)

# Change column positions
cols <- names(hrs)
# find their positions
iA <- which(cols == "score")
iB <- which(cols == "race")
iC <- which(cols == "race.ethnicity")
iD <- which(cols == "marriage")
# swap them
cols[c(iA, iB, iC, iD)] <- cols[c(iB, iA, iD, iC)]

# reindex the df
hrs <- hrs[ , cols]

set.seed(50)
# Pre-allocate result vectors (30 iterations)
n_iter <- 30
PR_PL_T <- numeric(n_iter)
PR_KS_T <- numeric(n_iter)
PR_PL_L <- numeric(n_iter)
PR_KS_L <- numeric(n_iter)

PRI_PL_T <- numeric(n_iter)
PRI_KS_T <- numeric(n_iter)
PRI_PL_L <- numeric(n_iter)
PRI_KS_L <- numeric(n_iter)

PRPA_PL_T <- numeric(n_iter)
PRPA_KS_T <- numeric(n_iter)
PRPA_PL_L <- numeric(n_iter)
PRPA_KS_L <- numeric(n_iter)

sample_size <- floor(0.8 * nrow(hrs))  # 80% training, 20% test

for (k in 1:n_iter) {
  cat("Iteration:", k, "\n")

  train_indices <- sample(seq_len(nrow(hrs)), size = sample_size)
  train_data <- hrs[train_indices, ]
  test_data  <- hrs[-train_indices, ] 
  
  idx_cov = 25
  idx_y = 26
  idx_sensi = 27
  
  #--------------------------
  # Plain Poisson Regression
  #--------------------------
  
  poisson_model <- glm(hrs$score[train_indices] ~ ., data =  hrs[train_indices, c(-idx_y)], family = poisson())
  
  PR_PL_L[k] = poissonloss(hrs$score[train_indices], poisson_model$fitted.values  )
  PR_KS_L[k] = f_ks(poisson_model$fitted.values, hrs[train_indices, idx_sensi])
  
  predicted_y <- predict(poisson_model, newdata =  hrs[-train_indices, c(-idx_y)], type = 'response')
  PR_PL_T[k] = poissonloss(hrs$score[-train_indices], predicted_y )
  PR_KS_T[k] =  f_ks(predicted_y, hrs[-train_indices, idx_sensi])
  
  # Fair regression methods
  split_by_gender_train <- split(train_data, train_data[,idx_sensi])
  split_by_gender_test <- split(test_data, test_data[,idx_sensi])
  
  # For train
  Ytrain = c()
  Strain = c()
  
  LYpred = c()
  rankYLF = c()
  
  # For test  
  Ytest = c()
  Stest = c()
  
  rankYTF = c()
  
  for (i in 1:length(unique(train_data[,idx_sensi])) ) {
    train_cov_split = split_by_gender_train[[i]][1:idx_cov]
    Ytrain_split = split_by_gender_train[[i]][idx_y]
    Ytrain = c(Ytrain, Ytrain_split[[1]])
    Strain = c(Strain, rep(i, nrow(train_cov_split)))
    
    
    test_cov_split = split_by_gender_test[[i]][1:idx_cov]
    Ytest_split = split_by_gender_test[[i]][idx_y]
    Ytest = c(Ytest, Ytest_split[[1]])
    Stest = c(Stest, rep(i, nrow(test_cov_split)))
    
   
    model = glm(Ytrain_split[[1]] ~ ., data = train_cov_split,  family = poisson())
    LYpred_split = predict(model, train_cov_split)
    TYpred_split = predict(model, test_cov_split)
    
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
  ispline_fit <- fit_ispline_model_Poisson(rankYLF, Ytrain, 
                                   knots = c(0.2, 0.5, 0.75), degree = 2)  
  train_pred <- ispline_predict(rankYLF, ispline_fit$coeff, 
                                knots = c(0.2, 0.5, 0.75), degree = 2)  
  PRI_PL_L[k] <- poissonloss(Ytrain, exp(train_pred)) 
  PRI_KS_L[k] <- f_ks(train_pred, Strain) 

                           
  test_pred <- ispline_predict(rankYTF, ispline_fit$coeff, 
                               knots = c(0.2, 0.5, 0.75), degree = 2)
  PRI_PL_T[k] <- poissonloss(Ytest, exp(test_pred))
  PRI_KS_T[k] <- f_ks(test_pred, Stest) 
  
  #--------------------------
  # PAVA Method
  #--------------------------
  pava_fit <- fit_pava_model_Poisson(rankYLF, Ytrain,  solverinput = poisson_solver ) 
  gpava_pred_train <- pava_predict(pava_fit, rankYLF) 
  PRPA_PL_L[k] <- poissonloss(Ytrain, gpava_pred_train) 
  PRPA_KS_L[k] <- f_ks(gpava_pred_train, Strain)
  
  # For test predictions, use our wrapper pava_predict:
  test_pred_pava <- pava_predict(pava_fit, rankYTF) 
  PRPA_PL_T[k] <- poissonloss(Ytest, test_pred_pava) 
  PRPA_KS_T[k] <- f_ks(test_pred_pava, Stest) 
}

#--------------------------
# Summarize Results
#--------------------------
cat("Means:\n")
cat("PR_PL_T:", mean(PR_PL_T), "\n")
cat("PR_KS_T:", mean(PR_KS_T), "\n")
cat("PR_PL_L:", mean(PR_PL_L), "\n")
cat("PR_KS_L:", mean(PR_KS_L), "\n")
cat("PRI_PL_T:", mean(PRI_PL_T), "\n")
cat("PRI_KS_T:", mean(PRI_KS_T), "\n")
cat("PRI_PL_L:", mean(PRI_PL_L), "\n")
cat("PRI_KS_L:", mean(PRI_KS_L), "\n")
cat("ORPA_PL_T:", mean(PRPA_PL_T), "\n")
cat("PRPA_KS_T:", mean(PRPA_KS_T), "\n")
cat("PRPA_PL_L:", mean(PRPA_PL_L), "\n")
cat("PRPA_KS_L:", mean(PRPA_KS_L), "\n")

cat("\nStandard Errors:\n")
cat("PR_PL_T:", sd(PR_PL_T)/sqrt(n_iter), "\n")
cat("PR_KS_T:", sd(PR_KS_T)/sqrt(n_iter), "\n")
cat("PR_PL_L:", sd(PR_PL_L)/sqrt(n_iter), "\n")
cat("PR_KS_L:", sd(PR_KS_L)/sqrt(n_iter), "\n")
cat("PRI_PL_T:", sd(PRI_PL_T)/sqrt(n_iter), "\n")
cat("PRI_KS_T:", sd(PRI_KS_T)/sqrt(n_iter), "\n")
cat("PRI_PL_L:", sd(PRI_PL_L)/sqrt(n_iter), "\n")
cat("PRI_KS_L:", sd(PRI_KS_L)/sqrt(n_iter), "\n")
cat("PRPA_PL_T:", sd(PRPA_PL_T)/sqrt(n_iter), "\n")
cat("PRPA_KS_T:", sd(PRPA_KS_T)/sqrt(n_iter), "\n")
cat("PRPA_PL_L:", sd(PRPA_PL_L)/sqrt(n_iter), "\n")
cat("PRPA_KS_L:", sd(PRPA_KS_L)/sqrt(n_iter), "\n")
