
quantile_loss <- function(y, y_hat, tau) {
  loss <- ifelse(y < y_hat, (tau - 1) * (y - y_hat), tau * (y - y_hat))
  mean(loss)
}
# Test: 
# 1. This function is aligned with definition in https://www.lokad.com/pinball-loss-function-definition/
# 2.quantile_loss(c(1,2,3,4), c(0.3,1,2,4),0.9) > quantile_loss(c(1,2,3,4), c(5,6,7,8),0.9)

# KS statistic
f_ks <- function(Y, S) {
  ts <- sort(unique(Y))
  
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
      for (t in ts) {
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
# Test: 
# 1.definition of each line
# 2.examples: f_ks(1:5, c(1,1,1,2,2)) = 1 f_ks(1:5, c(1,2,1,2,1)) = 0.3333333
# f_ks(1:10, c(1,2,3,3,2,1,1,2,3,2)) = 0.5 f_ks(1:10, c(1,1,1,2,2,2,3,3,3)) = 1

#------------------------------------
# ISPLINE METHOD FUNCTIONS
#------------------------------------

# https://www.fon.hum.uva.nl/praat/manual/spline.html Mspline is decided by t. t is decided by order and k-sigh_i
# From R package: if both knots and Boundary.knots are supplied, the basis parameters do not depend on x
fit_ispline_model <- function(rank_train, Y_train, tau = 0.75, 
                              knots = c(0.2, 0.5, 0.75),  degree = 2) {
  is_basis <- iSpline(rank_train, knots = knots, Boundary.knots = c(0,1),  degree = degree) # Generate iSpline 
  
  quantile_loss_objective <- function(par) {
    y_pred <- cbind(is_basis, 1) %*% par  # add intercept to the last column
    quantile_loss(Y_train, y_pred, tau) # Calculate the loss
  }
  initial_par <- rep(1, ncol(is_basis) + 1) # Initialize parameter
  optim_result <- hjkb(initial_par, quantile_loss_objective, 
                       lower = c(rep(0, ncol(is_basis)), -Inf))# Update the parameter to get the best one
  list(coeff = optim_result$par)
}

# Test: is_basis is plotted, which looks good.

# Given a (possibly new) vector of ranks, predict responses using the fitted Ispline model.
ispline_predict <- function(rank_new, coeff, 
                            knots = c(0.2, 0.5, 0.75), degree = 2) {
  new_basis <- iSpline(rank_new, knots = knots,  Boundary.knots = c(0,1), degree = degree)
  as.vector(cbind(new_basis, 1) %*% coeff) # To do the prediction
}


#------------------------------------
# PAVA METHOD FUNCTIONS
#------------------------------------

fit_pava_model <- function(rank_train, Y_train, tau = 0.75, solverinput = weighted.fractile) {
  gpava_result <- gpava(rank_train, Y_train, 
                        solver = solverinput, ties = "secondary", p = tau)
  fitted_values <- gpava_result$x
  
  # Remove duplicates based on training ranks
  unique_idx <- !duplicated(rank_train)
  rank_unique <- rank_train[unique_idx]
  gpava_unique <- fitted_values[unique_idx]
  
  # Order the unique values
  ord <- order(rank_unique)
  rank_unique <- rank_unique[ord]
  gpava_unique <- gpava_unique[ord]
  
  # Determine knot positions: indices where the fitted values change
  knots_idx <- which(diff(gpava_unique) != 0) + 1
  y0 <- c(gpava_unique[1], gpava_unique[knots_idx])
  # Build a step function; here we use the index positions as knots.
  model_fun <- stepfun(knots_idx, y0)
  list(rank_unique = rank_unique, model = model_fun)
}

# Predict new values from the fitted PAVA model.
pava_predict <- function(pava_model, rank_new) {
  sapply(rank_new, function(x) {
    idx <- sum(x >= pava_model$rank_unique)
    pava_model$model(idx)
  })
}
