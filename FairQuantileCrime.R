library(quantreg)
library(splines2)  
library(dfoptim)
library(Iso)   
library(isotone)
library(rqPen) 

set.seed(50)
n_iter <- 100
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

qu <- 0.25


#sample_size <- floor(0.8 * nrow(community))  # 80% training, 20% test


# Get data from community dataset for 100 repetitions

#for (k in 1:100) {

  # Create train/test split
#  train_indices <- sample(seq_len(nrow(community)), size = sample_size)
#  write.table(community[train_indices,], "100crimetrain75.csv", sep = ",", 
#              row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE)
#  write.table(community[-train_indices,], "100crimetest75.csv", sep = ",", 
#              row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE)
#}

read_appended_csv <- function(path) {
  lines <- readLines(path)
  header <- lines[1]
  lines_clean <- c(header, lines[lines != header])  # drop repeated headers
  read.csv(text = paste(lines_clean, collapse = "\n"))
}

# read the data files
train_all <- read_appended_csv("100crimetrain75.csv")
test_all  <- read_appended_csv("100crimetest75.csv")


n_train_each <- nrow(train_all) / n_iter
n_test_each  <- nrow(test_all)  / n_iter

# split into a list of 100 data frames
train_list <- split(train_all, rep(seq_len(n_iter), each = n_train_each))
test_list  <- split(test_all,  rep(seq_len(n_iter), each = n_test_each))




for (k in 1:n_iter) {
  cat("Iteration:", k, "\n")
  
  
  train_data <- train_list[[k]]
  test_data  <- test_list[[k]]
  

  idx_cov = 96
  idx_y = 97
  idx_race = 98
  
  train_data$race <- factor(train_data$race)
  test_data$race <- factor(test_data$race)


  # Plain Quantile Regression

  quantile_model <- rq(train_data[, idx_y] ~ ., 
                       data = train_data[, c(1:idx_cov, idx_race)], 
                       tau = qu) 
  
  QR_Pinball_L[k] <- quantile_loss(train_data[, idx_y], quantile_model$fitted.values, qu) 
  QR_KS_L[k] <- f_ks(quantile_model$fitted.values, train_data[, idx_race]) 
  
  predicted_y <- predict(quantile_model, newdata = test_data[, c(1:idx_cov, idx_race)])
  QR_Pinball_T[k] <- quantile_loss(test_data[, idx_y], predicted_y, qu) 
  QR_KS_T[k] <- f_ks(predicted_y, test_data[, idx_race])
  
  
 
  
  train_fit = cbind(train_data, yfit = quantile_model$fitted.values)
  test_fit = cbind(test_data, yfit = predicted_y)
  idx_yfit = ncol(train_fit)
  
  split_by_race_train <- split(train_fit, train_fit$race) 
  
  split_by_race_test <- split(test_fit, test_fit$race)
  
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
  
  for (i in 1:length(unique(train_data$race)) ) {
    LYpred_split = split_by_race_train[[i]][idx_yfit][[1]] 
    Ytrain_split = split_by_race_train[[i]][idx_y]
    racenum[i] = nrow(Ytrain_split)
    Ytrain = c(Ytrain, Ytrain_split[[1]])
    Strain = c(Strain, rep(i, nrow(Ytrain_split)))
    
    
    TYpred_split = split_by_race_test[[i]][idx_yfit][[1]]
    Ytest_split = split_by_race_test[[i]][idx_y]
    Ytest = c(Ytest, Ytest_split[[1]])
    Stest = c(Stest, rep(i, nrow(Ytest_split)))
    
    
    
    LYpred = c(LYpred, LYpred_split)
    
    rankYLF_split = rank(LYpred_split, ties.method = "max")/length(LYpred_split) 
    rankYLF = c(rankYLF, rankYLF_split) 
    
    rankYTF_split = sapply(TYpred_split, function(y) sum(LYpred_split <= y)/length(LYpred_split))
    rankYTF = c(rankYTF, rankYTF_split)
  }
   
  

  # Ispline Method
  
  
  # make knots
  n_interior_vec <- 0:10
  
  # build interior knots at equally spaced quantile positions of X
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
    mean_pinball = numeric(),
    sd_pinball   = numeric(),
    mean_ks      = numeric(),
    sd_ks        = numeric()
  )
  
  # Main grid loop
  for(kn_id in names(knot_list)){
    knots_vec <- knot_list[[kn_id]]
    for(deg in degree_vec){
      
      pinball_fold <- numeric(K)
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
        fit_k <- fit_ispline_model_Quantile(
          R_tr, Y_tr, tau = qu,
          knots = knots_vec,
          degree = deg
        )
        
        # Predict
        va_pred <- ispline_predict(
          R_va, fit_k$coeff,
          knots = knots_vec,
          degree = deg
        )
        
        pinball_fold[kk] <- quantile_loss(Y_va, va_pred, qu)
        ks_fold[kk]      <- f_ks(va_pred, S_va)
      }
      
      cv_results <- rbind(cv_results,
                          data.frame(
                            knots_id = kn_id,
                            degree = deg,
                            mean_pinball = mean(pinball_fold),
                            sd_pinball   = sd(pinball_fold),
                            mean_ks      = mean(ks_fold),
                            sd_ks        = sd(ks_fold)
                          ))
    }
  }
  
  ks_threshold = sort( cv_results$mean_ks)[floor(length(cv_results$mean_ks)*0.1)]
  cv_results = cv_results[which(cv_results$mean_ks<=ks_threshold ), ]
  
  
  cv_results <- cv_results[order(cv_results$mean_pinball), ]
  best_row <- cv_results[1, ]

  
  best_knots <- knot_list[[ best_row$knots_id ]]
  best_degree <- best_row$degree
  print(best_knots)
  print(best_degree)
  
 
    
  
  # Fit the Ispline model on training data
  ispline_fit <- fit_ispline_model_Quantile(rankYLF, Ytrain, tau = qu, 
                                            knots = best_knots, degree = best_degree)  
  train_pred <- ispline_predict(rankYLF, ispline_fit$coeff, 
                                knots = best_knots, degree = best_degree)  
  QRI_Pinball_L[k] <- quantile_loss(Ytrain, train_pred, qu) 
  QRI_KS_L[k] <- f_ks(train_pred, Strain) 
  
  
  test_pred <- ispline_predict(rankYTF, ispline_fit$coeff, 
                               knots = best_knots, degree = best_degree)
  QRI_Pinball_T[k] <- quantile_loss(Ytest, test_pred, qu)
  QRI_KS_T[k] <- f_ks(test_pred, Stest) 
  

  # PAVA Method
  pava_fit <- fit_pava_model_Quantile(rankYLF, Ytrain, tau = qu, solverinput = weighted.fractile )  
  gpava_pred_train <- pava_predict(pava_fit, rankYLF) 
  QRPA_Pinball_L[k] <- quantile_loss(Ytrain, gpava_pred_train, qu) 
  QRPA_KS_L[k] <- f_ks(gpava_pred_train, Strain)
  
  # Test predictions:
  test_pred_pava <- pava_predict(pava_fit, rankYTF) 
  QRPA_Pinball_T[k] <- quantile_loss(Ytest, test_pred_pava, qu) 
  QRPA_KS_T[k] <- f_ks(test_pred_pava, Stest) 
}


# Summarize Results

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



