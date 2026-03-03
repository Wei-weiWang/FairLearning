
#------------------------------------
# Main Loop 
#------------------------------------
set.seed(45)
n_iter <- 200
# Pre-allocate result vectors
RR_HL_T <- numeric(n_iter)
RR_KS_T <- numeric(n_iter)
RR_HL_L <- numeric(n_iter)
RR_KS_L <- numeric(n_iter)

RRI_HL_T <- numeric(n_iter)
RRI_KS_T <- numeric(n_iter)
RRI_HL_L <- numeric(n_iter)
RRI_KS_L <- numeric(n_iter)

RRPA_HL_T <- numeric(n_iter)
RRPA_KS_T <- numeric(n_iter)
RRPA_HL_L <- numeric(n_iter)
RRPA_KS_L <- numeric(n_iter)

#The parameter represents the level of bias. The larger the parameter, the higher the bias level.
a = 1
#The parameter represents the cutoff in the Huber loss.
delta <- 9.67*1.345

for (k in 1:n_iter) {
  cat("Iteration:", k, "\n")
  
  # Train data generation
  # Sample size and sensitive attribute S 
  n = 500
  S = rbinom(n, size = 1, prob = 0.5)
  # Number of samples with S=0 and S=1
  n_train1 <- sum(S==0)
  n_train2 <- sum(S==1)
  
  # Generate features for group S=0 
  mu = rep(0, 8)
  Sigma <- toeplitz(   0.5^(0:(8-1))  )
  X1 = rmvnorm(n_train1, mu, Sigma)
  beta1 = c(3, 1.5, 0, 0, 2, 0, 0, 0) 
  
  # Generate outcome for group S=0
  Y1_clean <- 1 + X1%*%beta1 
  
  # Generate and standardize noise  
  z1  <- rbinom(n_train1, 1, 0.9)                  
  esd1 <- ifelse(z1 == 1, 1, 15)
  V1  <- rnorm(n_train1, mean = 0, sd = esd1)
  sdV <- sqrt(0.9 * 1^2 + 0.1 * 15^2)  
  eps1 <- V1 / sdV
  # Final Y for group S=0
  Y1  <- Y1_clean + 9.67*eps1
  
  # Generate features for group S=1
  X2 = rmvnorm(n_train2, mu, Sigma)
  beta2 = c(3, 1.5, 0, 0, 2, 0, 0, 0)
  
  # Generate outcome for group S=1
  Y2_clean <- 1 + X2%*%beta2 + a
  
  z2  <- rbinom(n_train2, 1, 0.9)                  
  esd2 <- ifelse(z2 == 1, 1, 15)
  V2  <- rnorm(n_train2, mean = 0, sd = esd2)
  eps2 <- V2 / sdV
  Y2       <- Y2_clean + 9.67*eps2
  
  # Combine into training dataset
  X <- rbind(X1, X2)      
  Y   <- c(Y1, Y2)       
  S   <- c(rep(0, n_train1),
           rep(1, n_train2)) 
  
  # Training data frame
  train_df   <- data.frame(X = X, Y = Y, S = as.factor(S))
  
  # Test data generation
  n_test1 <- 5000
  n_test2 <- 5000
  
  # Generate features for group S=0 in test
  tX1 = rmvnorm(n_test1, mu, Sigma)
  
  # Generate outcomes for group S=0 in test
  tY1_clean <-  1 + tX1%*%beta1 
  tz1  <- rbinom(n_test1, 1, 0.9)                  
  tesd1 <- ifelse(tz1 == 1, 1, 15)
  tV1  <- rnorm(n_test1, mean = 0, sd = tesd1)
  teps1 <- tV1 / sdV
  tY1  <- tY1_clean + 9.67*teps1
  
  # Generate features for group S=1 in test
  tX2 = rmvnorm(n_test2, mu, Sigma)
  
  # Generate outcomes for group S=1 in test
  tY2_clean <-  1 + tX2%*%beta2 + a
  tz2  <- rbinom(n_test2, 1, 0.9)                  
  tesd2 <- ifelse(tz2 == 1, 1, 15)
  tV2  <- rnorm(n_test2, mean = 0, sd = tesd2)
  teps2 <- tV2 / sdV
  tY2  <- tY2_clean + 9.67*teps2
  
  # Combine into test dataset
  tX <- rbind(tX1, tX2)     
  tY   <- c(tY1,  tY2)      
  tS   <- c(rep(0, n_test1),
            rep(1, n_test2))  
  
  # Test data frame
  test_df   <- data.frame(X = tX, Y = tY, S = as.factor(tS))
  
  # Column indices: covariates, response, and group indicator 
  idx_cov = 8
  idx_y = 9
  idx_sensi = 10
  
  # Plain Robust Regression
  # Fit robust linear regression by minimizing Huber loss with cutoff delta
  Rfit  <- huber_fit_delta(Y, cbind(X, S), delta = delta)
  
  # Predictions on training and test
  LYpredictions = Rfit$coef[1]  + cbind(X, S)%*%Rfit$coef[-1] 
  TYpredictions = Rfit$coef[1]  + cbind(tX, tS)%*%Rfit$coef[-1] 
  
  # Store performance (Huber loss) and fairness (KS statistic) on train/test
  RR_HL_L[k] =  mean(huber_loss(LYpredictions-Y, delta))
  RR_KS_L[k] =  f_ks(LYpredictions,  S)
  
  RR_HL_T[k] =  mean(huber_loss(TYpredictions-tY, delta))
  RR_KS_T[k] =  f_ks(TYpredictions,  tS)
  
  # Attach fitted values for later group-wise ranking
  train_fit = cbind(train_df, yfit = LYpredictions)
  test_fit = cbind(test_df, yfit = TYpredictions)
  # Column index of yfit
  idx_yfit = ncol(train_fit)
  
  # Split data by sensitive group
  split_by_sensi_train <- split(train_fit, train_fit[, idx_sensi]) 
  split_by_sensi_test <- split(test_fit, test_fit[, idx_sensi])
  
  # Estimate F*_S(f*(X,S)) on training and test dataset
  # Train pooled vectors
  Ytrain = c()
  Strain = c()
  LYpred = c()
  rankYLF = c()
  
  # Test pooled vectors
  Ytest = c()
  Stest = c()
  rankYTF = c()
  racenum = c()
  
  for (i in 1:length(unique(train_fit[, idx_sensi])) ) {
    LYpred_split = split_by_sensi_train[[i]][idx_yfit][[1]] 
    Ytrain_split = split_by_sensi_train[[i]][idx_y]
    racenum[i] = nrow(Ytrain_split)
    Ytrain = c(Ytrain, Ytrain_split[[1]])
    Strain = c(Strain, rep(i, nrow(Ytrain_split)))
    
    TYpred_split = split_by_sensi_test[[i]][idx_yfit][[1]]
    Ytest_split = split_by_sensi_test[[i]][idx_y]
    Ytest = c(Ytest, Ytest_split[[1]])
    Stest = c(Stest, rep(i, nrow(Ytest_split)))
    
    LYpred = c(LYpred, LYpred_split)
    
    rankYLF_split = rank(LYpred_split, ties.method = "max")/length(LYpred_split)
    rankYLF = c(rankYLF, rankYLF_split)
    
    rankYTF_split = sapply(TYpred_split, function(y) sum(LYpred_split <= y)/length(LYpred_split))
    rankYTF = c(rankYTF, rankYTF_split)
  }
  
  #--------------------------
  # Ispline Method
  #--------------------------
  # make knots
  n_interior_vec <- 0:10
  knot_list <- lapply(n_interior_vec, function(m) make_knots(m))
  degree_vec <- 1:10
  names(knot_list) <- paste0("m", n_interior_vec)  
  
  # K fold cross validation
  K=5
  ntrain <- length(Ytrain)
  fold_id <- sample(rep(1:K, length.out = ntrain))
  
  # Storage for results
  cv_results <- data.frame(
    knots_id = character(),
    degree   = integer(),
    mean_huber = numeric(),
    sd_huber   = numeric(),
    mean_ks      = numeric(),
    sd_ks        = numeric()
  )
  
  # Main grid loop
  for(kn_id in names(knot_list)){
    knots_vec <- knot_list[[kn_id]]
    for(deg in degree_vec){
      
      huber_fold <- numeric(K)
      ks_fold      <- numeric(K)
      
      for(kk in 1:K){
        # Get training and evaluation data for cross validation
        idx_tr <- fold_id != kk
        idx_va <- !idx_tr
        
        R_tr <- rankYLF[idx_tr]
        Y_tr <- Ytrain[idx_tr]
        S_tr <- Strain[idx_tr]   
        R_va <- rankYLF[idx_va]
        Y_va <- Ytrain[idx_va]
        S_va <- Strain[idx_va]
        
        # Fit
        fit_k <- fit_ispline_model_Robust_Convex(
          R_tr, Y_tr,
          knots = knots_vec,
          degree = deg,
          delta = delta
        )
        
        # Predict
        va_pred <- ispline_predict(
          R_va, fit_k$coeff,
          knots = knots_vec,
          degree = deg
        )
        
        # Metrics
        huber_fold[kk] <- mean(huber_loss(Y_va- va_pred, delta) )
        ks_fold[kk]      <- f_ks(va_pred, S_va)
      }
      
      cv_results <- rbind(cv_results,
                          data.frame(
                            knots_id = kn_id,
                            degree = deg,
                            mean_huber = mean(huber_fold),
                            sd_huber   = sd(huber_fold),
                            mean_ks      = mean(ks_fold),
                            sd_ks        = sd(ks_fold)
                          ))
    }
  }
  
  # Get the best knot and degree values according to KS and pinball loss
  ks_threshold = sort( cv_results$mean_ks)[floor(length(cv_results$mean_ks)*0.1)]
  cv_results = cv_results[which(cv_results$mean_ks<=ks_threshold ), ]
  cv_results <- cv_results[order(cv_results$mean_huber), ]
  best_row <- cv_results[1, ]
  best_knots <- knot_list[[ best_row$knots_id ]]
  best_degree <- best_row$degree
  print(best_knots)
  print(best_degree)
  
  # Fit the Ispline model on training data
  ispline_fit <- fit_ispline_model_Robust_Convex(rankYLF, Ytrain, 
                                                 knots = best_knots, degree = best_degree, delta = delta)  
  train_pred <- ispline_predict(rankYLF, ispline_fit$coeff, 
                                knots = best_knots, degree = best_degree)  
  RRI_HL_L[k] <- mean(huber_loss(Ytrain- train_pred, delta) )
  RRI_KS_L[k] <- f_ks(train_pred, Strain) 
  
  # For test predictions, use our wrapper ispline_predict
  test_pred <- ispline_predict(rankYTF, ispline_fit$coeff, 
                               knots = best_knots, degree = best_degree)
  RRI_HL_T[k] <- mean(huber_loss(Ytest-test_pred,delta))
  RRI_KS_T[k] <- f_ks(test_pred, Stest) 
  
  #--------------------------
  # PAVA Method
  #--------------------------
  # Fit the PAVA model on training data
  pava_fit <- fit_pava_model_Robust(rankYLF, Ytrain, solverinput = huber_solver ) 
  gpava_pred_train <- pava_predict(pava_fit, rankYLF) 
  RRPA_HL_L[k] <- mean(huber_loss(Ytrain-gpava_pred_train, delta) ) 
  RRPA_KS_L[k] <- f_ks(gpava_pred_train, Strain)
  
  # For test predictions, use our wrapper pava_predict:
  test_pred_pava <- pava_predict(pava_fit, rankYTF) 
  RRPA_HL_T[k] <- mean(huber_loss(Ytest-test_pred_pava, delta)) 
  RRPA_KS_T[k] <- f_ks(test_pred_pava, Stest) 
}

#--------------------------
# Summarize Results
#--------------------------
cat("Means:\n")
cat("RR_HL_T:", mean(RR_HL_T), "\n")
cat("RR_KS_T:", mean(RR_KS_T), "\n")
cat("RR_HL_L:", mean(RR_HL_L), "\n")
cat("RR_KS_L:", mean(RR_KS_L), "\n")
cat("RRI_HL_T:", mean(RRI_HL_T), "\n")
cat("RRI_KS_T:", mean(RRI_KS_T), "\n")
cat("RRI_HL_L:", mean(RRI_HL_L), "\n")
cat("RRI_KS_L:", mean(RRI_KS_L), "\n")
cat("RRPA_HL_T:", mean(RRPA_HL_T), "\n")
cat("RRPA_KS_T:", mean(RRPA_KS_T), "\n")
cat("RRPA_HL_L:", mean(RRPA_HL_L), "\n")
cat("RRPA_KS_L:", mean(RRPA_KS_L), "\n")

cat("\nStandard Errors:\n")
cat("RR_HL_T:", sd(RR_HL_T)/sqrt(n_iter), "\n")
cat("RR_KS_T:", sd(RR_KS_T)/sqrt(n_iter), "\n")
cat("RR_HL_L:", sd(RR_HL_L)/sqrt(n_iter), "\n")
cat("RR_KS_L:", sd(RR_KS_L)/sqrt(n_iter), "\n")
cat("RRI_HL_T:", sd(RRI_HL_T)/sqrt(n_iter), "\n")
cat("RRI_KS_T:", sd(RRI_KS_T)/sqrt(n_iter), "\n")
cat("RRI_HL_L:", sd(RRI_HL_L)/sqrt(n_iter), "\n")
cat("RRI_KS_L:", sd(RRI_KS_L)/sqrt(n_iter), "\n")
cat("RRPA_HL_T:", sd(RRPA_HL_T)/sqrt(n_iter), "\n")
cat("RRPA_KS_T:", sd(RRPA_KS_T)/sqrt(n_iter), "\n")
cat("RRPA_HL_L:", sd(RRPA_HL_L)/sqrt(n_iter), "\n")
cat("RRPA_KS_L:", sd(RRPA_KS_L)/sqrt(n_iter), "\n")

