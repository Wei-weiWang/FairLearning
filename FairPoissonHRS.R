library(quantreg)
library(splines2)  
library(dfoptim)
library(Iso)   
library(isotone)
library(rqPen) 

set.seed(50)
n_iter <- 100
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


read_appended_csv <- function(path) {
  lines <- readLines(path)
  header <- lines[1]
  lines_clean <- c(header, lines[lines != header])  # drop repeated headers
  read.csv(text = paste(lines_clean, collapse = "\n"))
}


# read the data files
train_all <- read_appended_csv("100hrstrain75.csv")
test_all  <- read_appended_csv("100hrstest75.csv")


n_train_each <- nrow(train_all) / n_iter
n_test_each  <- nrow(test_all)  / n_iter

# split into a list of 100 data frames
train_list <- split(train_all, rep(seq_len(n_iter), each = n_train_each))
test_list  <- split(test_all,  rep(seq_len(n_iter), each = n_test_each))


for (k in 1:n_iter) {
  cat("Iteration:", k, "\n")
  
  
  train_data <- train_list[[k]]
  test_data  <- test_list[[k]]
  
  
  idx_cov = 25
  idx_y = 26
  idx_sensi = 27
  
  train_data$marriage <- factor(train_data$marriage)
  test_data$marriage <- factor(test_data$marriage)
  

  # Plain Poisson Regression
  
  poisson_model <- glm(train_data$score ~ ., data =  train_data[, c(-idx_y)], family = poisson())
  
  
  PR_PL_L[k] = poissonloss(train_data$score, poisson_model$fitted.values  )
  PR_KS_L[k] = f_ks(poisson_model$fitted.values, train_data[, idx_sensi])
  
  
  predicted_y <- predict(poisson_model, newdata =  test_data[, c(-idx_y)], type = 'response')
  PR_PL_T[k] = poissonloss(test_data$score, predicted_y )
  PR_KS_T[k] =  f_ks(predicted_y, test_data[, idx_sensi])
  
  
  train_fit = cbind(train_data, yfit = poisson_model$fitted.values)
  test_fit = cbind(test_data, yfit = predicted_y)
  idx_yfit = ncol(train_fit)
  
  split_by_sensi_train <- split(train_fit, train_fit[, idx_sensi]) 
  
  split_by_sensi_test <- split(test_fit, test_fit[, idx_sensi])
  
  # For train
  Ytrain = c()
  Strain = c()
  
  LYpred = c()
  rankYLF = c()
  
  # For test  
  Ytest = c()
  Stest = c()
  
  rankYTF = c()
  racenum = c()
  
  for (i in 1:length(unique(train_fit[, idx_sensi])) ) {
    LYpred_split = split_by_sensi_train[[i]][idx_yfit][[1]] #note: get columns
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
  
 
    
  
  names(rankYLF) <- NULL
  

  # Ispline Method
  # make knots
  n_interior_vec <- 0:10
  
  # Build interior knots at equally spaced quantile positions of X
  make_knots <- function(n_interior){
    if(n_interior == 0) return(numeric(0))
    probs <- seq(0, 1, length.out = n_interior + 2)[-c(1, n_interior + 2)]
    return(probs)
  }
  
  knot_list <- lapply(n_interior_vec, function(m) make_knots(m))
  degree_vec <- 1:10
  names(knot_list) <- paste0("m", n_interior_vec)  
  
  # K fold cross validation
  K <- 5
  ntrain <- length(Ytrain)
  fold_id <- sample(rep(1:K, length.out = ntrain))
  
  # Storage for results
  cv_results <- data.frame(
    knots_id = character(),
    degree   = integer(),
    mean_poisson = numeric(),
    sd_poisson   = numeric(),
    mean_ks      = numeric(),
    sd_ks        = numeric()
  )
  
  # Main grid loop
  for(kn_id in names(knot_list)){
    #print(kn_id)
    knots_vec <- knot_list[[kn_id]]
    #print(knots_vec)
    for(deg in degree_vec){
      
      poisson_fold <- numeric(K)
      ks_fold      <- numeric(K)
      
      for(kk in 1:K){
        idx_tr <- fold_id != kk
        idx_va <- !idx_tr
        
        R_tr <- rankYLF[idx_tr]
        Y_tr <- Ytrain[idx_tr]
        S_tr <- Strain[idx_tr]   
        R_va <- rankYLF[idx_va]
        Y_va <- Ytrain[idx_va]
        S_va <- Strain[idx_va]
        
        # Fit
        fit_k <- fit_ispline_model_Poisson(
          R_tr, Y_tr,
          knots = knots_vec,
          degree = deg
        )
        
        # Predict
        va_pred <- ispline_predict(
          R_va, fit_k$coeff,
          knots = knots_vec,
          degree = deg
        )
        
        poisson_fold[kk] <- poissonloss(Y_va, exp(va_pred))
        ks_fold[kk]      <- f_ks(va_pred, S_va)
      }
      
      cv_results <- rbind(cv_results,
                          data.frame(
                            knots_id = kn_id,
                            degree = deg,
                            mean_poisson = mean(poisson_fold),
                            sd_poisson   = sd(poisson_fold),
                            mean_ks      = mean(ks_fold),
                            sd_ks        = sd(ks_fold)
                          ))
    }
  }
  
  
  
  ks_threshold = sort( cv_results$mean_ks)[floor(length(cv_results$mean_ks)*0.1)]
  cv_results = cv_results[which(cv_results$mean_ks<=ks_threshold ), ]
  
  
  cv_results <- cv_results[order(cv_results$mean_poisson), ]
  best_row <- cv_results[1, ]
  
  
  best_knots <- knot_list[[ best_row$knots_id ]]
  best_degree <- best_row$degree
  print(best_knots)
  print(best_degree)
  
  
  
  # Fit the Ispline model on training data
  ispline_fit <- fit_ispline_model_Poisson(rankYLF, Ytrain, 
                                           knots = best_knots, degree = best_degree)  
  train_pred <- ispline_predict(rankYLF, ispline_fit$coeff, 
                                knots = best_knots, degree = best_degree)  
  PRI_PL_L[k] <- poissonloss(Ytrain, exp(train_pred)) 
  PRI_KS_L[k] <- f_ks(train_pred, Strain) 
  
  
  test_pred <- ispline_predict(rankYTF, ispline_fit$coeff, 
                               knots = best_knots, degree = best_degree)
  PRI_PL_T[k] <- poissonloss(Ytest, exp(test_pred))
  PRI_KS_T[k] <- f_ks(test_pred, Stest) 
  

  # PAVA Method
  pava_fit <- fit_pava_model_Poisson(rankYLF, Ytrain,  solverinput = poisson_solver ) 
  gpava_pred_train <- pava_predict(pava_fit, rankYLF) 
  PRPA_PL_L[k] <- poissonloss(Ytrain, gpava_pred_train) 
  PRPA_KS_L[k] <- f_ks(gpava_pred_train, Strain)
  
  # Test predictions:
  test_pred_pava <- pava_predict(pava_fit, rankYTF) 
  PRPA_PL_T[k] <- poissonloss(Ytest, test_pred_pava) 
  PRPA_KS_T[k] <- f_ks(test_pred_pava, Stest) 
}


# Results of each iteration

cat("PR_PL_T:", PR_PL_T, "\n")
cat("PR_KS_T:", PR_KS_T, "\n")
cat("PR_PL_L:", PR_PL_L, "\n")
cat("PR_KS_L:", PR_KS_L, "\n")
cat("PRI_PL_T:", PRI_PL_T, "\n")
cat("PRI_KS_T:", PRI_KS_T, "\n")
cat("PRI_PL_L:", PRI_PL_L, "\n")
cat("PRI_KS_L:", PRI_KS_L, "\n")
cat("PRPA_PL_T:", PRPA_PL_T, "\n")
cat("PRPA_KS_T:", PRPA_KS_T, "\n")
cat("PRPA_PL_L:", PRPA_PL_L, "\n")
cat("PRPA_KS_L:", PRPA_KS_L, "\n")
