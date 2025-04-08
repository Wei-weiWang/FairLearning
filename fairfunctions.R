# Solver for gpava
poisson_solver = function(Y, w){ 
  
  return(   sum(Y*w)/sum(w)     )  
  
}

poissonloss = function( Y, expfhat  ){ 
  
  ans <- ifelse(Y == 0, expfhat, expfhat - Y* log(expfhat )  )
  return( mean (ans) ) 
  
}

quantile_loss <- function(y, y_hat, tau) {
  loss <- numeric(length(y))
  for (i in seq_along(y)) {
    if (y[i] < y_hat[i]) {
      loss[i] <- (tau - 1) * (y[i] - y_hat[i])
    } else {
      loss[i] <- tau * (y[i] - y_hat[i])
    }
  }
  return(mean(loss))
}

f_ks <- function(Y, S) {
  # Generate the sequence of thresholds
  tt <- seq(min(Y), max(Y), length.out = 2000)
  
  # Get unique classes
  classes <- unique(S)
  num_classes <- length(classes)
  
  # Initialize the maximum KS statistic
  max_fai <- 0
  
  # Loop through all pairs of classes
  for (i in 1:(num_classes - 1)) {
    for (j in (i + 1):num_classes) {
      class_i <- classes[i]
      class_j <- classes[j]
      
      # Get the counts for each class
      nn_i <- sum(S == class_i)
      nn_j <- sum(S == class_j)
      
      # Compute KS statistic for this pair of classes
      fai <- 0
      for (t in tt) {
        tmp1 <- sum(Y[S == class_i] <= t) / nn_i
        tmp2 <- sum(Y[S == class_j] <= t) / nn_j
        fai <- max(fai, abs(tmp1 - tmp2))
      }
      
      # Update the maximum KS statistic
      max_fai <- max(max_fai, fai)
    }
  }
  
  return(max_fai)
}


# Note that the training and test data should have been rearranged according to sensitive variable, respectively, before using this Ifair and Pfair.
# Here, we require that the input vectors Strain and Stest are sorted in increasing order (e.g., c(0, 0, 1, 1, 1)).
Ifairpoi = function(fLYpredictions, fTYpredictions, fStrain, fStest, fYtrain){
  
  YLF0 = fLYpredictions[which(fStrain==0)]
  rankYLF0 = rank(YLF0, ties.method = "max")/length(YLF0)
  YLF1 = fLYpredictions[which(fStrain==1)]
  rankYLF1 = rank(YLF1, ties.method = "max")/length(YLF1)
  
  YLT01 = fYtrain#YLT list
  rankYLF01 = c(rankYLF0, rankYLF1)# Rank of each element in YLT list(decided by YLF)
  
  isknots = c(0.2,0.5,0.75)
  degree = 2
  
  
  is_basis <- iSpline(rankYLF01, knots = isknots, degree = degree)
  
  poisson_loss_tie <- function(par) {
    
    
    column_of_ones <- rep(1, nrow(is_basis))
    
    # Add the column of 1s to the matrix
    is_basis_wiin <- cbind(is_basis, column_of_ones)
    
    # Calculate the predicted y values
    y_pred = is_basis_wiin%*%par
    
    # Calculate Pinball loss
    Ploss = poissonloss( YLT01, exp(y_pred))
    
    return(Ploss)
  }
  
  initial_par <- rep(1, ncol(is_basis)+1 )
  
  optim_result = hjkb(initial_par, poisson_loss_tie, lower = c(rep(0, ncol(is_basis)) ,-Inf) )
  
  # Predict on the train data using the optimal coefficients
  column_of_ones <- rep(1, nrow(is_basis))
  train_pred <- cbind(is_basis, column_of_ones) %*% optim_result$par
  
  
  # Predict on the test data
  YTF0 = fTYpredictions[which(fStest==0)]
  nYLF0 = length(YLF0)
  rankYTF0 <- sapply(YTF0, function(x) sum(YLF0 <= x)) / nYLF0
  
  YTF1 = fTYpredictions[which(fStest==1)]
  nYLF1 = length(YLF1)
  rankYTF1 <- sapply(YTF1, function(x) sum(YLF1 <= x)) / nYLF1
  
  
  rankYTF01 = c(rankYTF0, rankYTF1)
  
  i_spline_basis_test <- iSpline(rankYTF01, knots = isknots, degree = degree)
  column_of_ones <- rep(1, nrow(i_spline_basis_test))
  test_pred <- cbind(i_spline_basis_test, column_of_ones) %*% optim_result$par
  
  return(list(fairtrain = exp(train_pred), fairtest = exp(test_pred)))
  
  
  
}
Pfair = function(fLYpredictions, fTYpredictions, fStrain, fStest, fYtrain, method = poisson_solver, quantile = NA){
  
  YLF0 = fLYpredictions[which(fStrain==0)]
  rankYLF0 = rank(YLF0, ties.method = "max")/length(YLF0)
  YLF1 = fLYpredictions[which(fStrain==1)]
  rankYLF1 = rank(YLF1, ties.method = "max")/length(YLF1)
  
  YLT01 = fYtrain#YLT list
  rankYLF01 = c(rankYLF0, rankYLF1)# Rank of each element in YLT list(decided by YLF)
  
  
  # predictions of train and test data by pava
  gpavapre = gpava(rankYLF01, YLT01, solver = method, ties = "secondary", p =quantile)
  gpavapre = gpavapre$x
  
  
  # Predict on the test data
  YTF0 = fTYpredictions[which(fStest==0)]
  nYLF0 = length(YLF0)
  rankYTF0 <- sapply(YTF0, function(x) sum(YLF0 <= x)) / nYLF0
  
  YTF1 = fTYpredictions[which(fStest==1)]
  nYLF1 = length(YLF1)
  rankYTF1 <- sapply(YTF1, function(x) sum(YLF1 <= x)) / nYLF1
  rankYTF01 = c(rankYTF0, rankYTF1)
  
  
  
  # For test data to use PAVA
  gpavapreuni = c()
  rankYLF01uni = c()
  
  for (i in 1:length(rankYLF01)) {
    if(rankYLF01[i] %in% rankYLF01uni == TRUE){
      next
      
    }else{
      
      rankYLF01uni = c(rankYLF01uni, rankYLF01[i])
      gpavapreuni = c(gpavapreuni, gpavapre[i])
      
    }
    
  }
  orderrankYLF01uni = order(rankYLF01uni)
  rankYLF01uniordered = rankYLF01uni[orderrankYLF01uni]
  gpavapreuniordered = gpavapreuni[orderrankYLF01uni]
  
  
  knots = c()
  
  for (i in 2:length(rankYLF01uni)) {
    if(gpavapreuniordered[i]==gpavapreuniordered[i-1]){
      
      
    }else{
      
      knots = c(knots, i)
    }
    
  }
  
  if(gpavapreuniordered[1]!=0){
    y0 <- c(gpavapreuniordered[1], gpavapreuniordered[knots])
  }else{
    
    y0 <- c(gpavapreuniordered[knots[1]], gpavapreuniordered[knots])
    
  }
  
  
  pavatrainedYhatfun <- stepfun(knots, y0)
  
  
  
  YT0byPAVA = c()
  for (i in 1:length(YTF0)) {
    YT0byPAVA[i] = pavatrainedYhatfun( sum( rankYTF0[i]>= rankYLF01uni) ) }
  
  YT1byPAVA = c()
  for (i in 1:length(YTF1)) {
    YT1byPAVA[i] = pavatrainedYhatfun( sum( rankYTF1[i]>= rankYLF01uni) ) }
  
  
  # Combine the results
  YtesthatPAVA = c(YT0byPAVA, YT1byPAVA)
  
  
  return(list(fairtrain = gpavapre, fairtest = YtesthatPAVA ))
}
Ifairqua = function(fLYpredictions, fTYpredictions, fStrain, fStest, fYtrain){
  
  YLF0 = fLYpredictions[which(fStrain==0)]
  rankYLF0 = rank(YLF0, ties.method = "max")/length(YLF0)
  YLF1 = fLYpredictions[which(fStrain==1)]
  rankYLF1 = rank(YLF1, ties.method = "max")/length(YLF1)
  
  YLT01 = fYtrain#YLT list
  rankYLF01 = c(rankYLF0, rankYLF1)# Rank of each element in YLT list(decided by YLF)
  
  isknots = c(0.2,0.5,0.75)
  degree = 2
  
  
  is_basis <- iSpline(rankYLF01, knots = isknots, degree = degree)
  
  quantile_loss_tie <- function(par) {
    
    
    column_of_ones <- rep(1, nrow(is_basis))
    
    # Add the column of 1s to the matrix
    is_basis_wiin <- cbind(is_basis, column_of_ones)
    
    # Calculate the predicted y values
    y_pred = is_basis_wiin%*%par
    
    # Calculate Pinball loss
    Ploss = quantile_loss( YLT01, y_pred, qu)*length(y_pred)
    
    return(Ploss)
  }
  
  initial_par <- rep(1, ncol(is_basis)+1 )
  
  optim_result = hjkb(initial_par, quantile_loss_tie, lower = c(rep(0, ncol(is_basis)) ,-Inf) )
  
  # Predict on the train data using the optimal coefficients
  column_of_ones <- rep(1, nrow(is_basis))
  train_pred <- cbind(is_basis, column_of_ones) %*% optim_result$par
  
  
  # Predict on the test data
  YTF0 = fTYpredictions[which(fStest==0)]
  nYLF0 = length(YLF0)
  rankYTF0 <- sapply(YTF0, function(x) sum(YLF0 <= x)) / nYLF0
  
  YTF1 = fTYpredictions[which(fStest==1)]
  nYLF1 = length(YLF1)
  rankYTF1 <- sapply(YTF1, function(x) sum(YLF1 <= x)) / nYLF1
  
  
  rankYTF01 = c(rankYTF0, rankYTF1)
  
  i_spline_basis_test <- iSpline(rankYTF01, knots = isknots, degree = degree)
  column_of_ones <- rep(1, nrow(i_spline_basis_test))
  test_pred <- cbind(i_spline_basis_test, column_of_ones) %*% optim_result$par
  
  return(list(fairtrain = train_pred, fairtest = test_pred))
  
  
  
}
