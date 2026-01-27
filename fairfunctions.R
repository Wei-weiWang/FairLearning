library(quantreg)
library(splines2)  
library(dfoptim)
library(Iso)   
library(isotone)
library(rqPen) 


MSE<- function(YT, YF) {
  err <- sum((YT - YF)^2)/length(YT)
  return(err)
}

quantile_loss <- function(y, y_hat, tau) {
  loss <- ifelse(y < y_hat, (tau - 1) * (y - y_hat), tau * (y - y_hat))
  mean(loss)
}


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
      

      nn_i <- sum(S == class_i)
      nn_j <- sum(S == class_j)
      

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


poisson_solver = function(Y, w){ 
  
  return(   sum(Y*w)/sum(w)     )  
  
}

poissonloss = function( Y, expfhat  ){ 
  
  ans = expfhat - Y* log(expfhat)
  return( mean (ans) ) 
  
} 

L_obj <- function(m, y, w, M) sum(w * huber_loss(y - m, M))

huber_solver <- function(y, w = NULL, M = delta) {
  if (is.null(w)) w <- rep(1, length(y))
  rng <- range(y) + c(-10*M, 10*M)       # wide bracket
  res <- optimize(L_obj, interval = rng, y = y, w = w, M = M)
  res$minimum
}

huber_loss <- function(r, delta) {ifelse(abs(r) <= delta,
                                        0.5 * r^2,
                                        delta * abs(r) - 0.5 * delta^2)}

huber_psi <- function(r, delta) pmax(pmin(r, delta), -delta)  # derivative

huber_fit_delta <- function(y, X, delta, intercept = TRUE) {
  X <- as.matrix(X)
  if (intercept && (ncol(X) == 0 || any(X[,1] != 1))) X <- cbind(1, X)
  p <- ncol(X)
  
  obj <- function(b) sum(huber_loss(y - drop(X %*% b), delta))
  grad <- function(b) {
    r <- y - drop(X %*% b)
    - crossprod(X, huber_psi(r, delta))
  }
  
  b0 <- rep(0, p)                   
  sol <- optim(b0, obj, grad, method = "BFGS", control = list(reltol = 1e-10))
  list(coef = sol$par, convergence = sol$convergence)
}


# I-spline METHOD FUNCTIONS

fit_ispline_model_Quantile <- function(rank_train, Y_train, tau = 0.75, 
                                       knots = c(0.25, 0.5, 0.75),  degree = 2) {
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

fit_ispline_model_Poisson <- function(rank_train, Y_train,
                                      knots = c(0.25, 0.5, 0.75),  degree = 2) {
  
  is_basis <- iSpline(rank_train, knots = knots, Boundary.knots = c(0,1),  degree = degree) # Generate iSpline 
  
  poisson_loss_tie <- function(par) {
    column_of_ones <- rep(1, nrow(is_basis))
    y_pred = cbind(is_basis, column_of_ones)%*%par
    Ploss = poissonloss(Y_train, exp(y_pred))
    
    return(Ploss)
  }
  initial_par <- rep(1, ncol(is_basis) + 1) # Initialize parameter
  optim_result <- hjkb(initial_par, poisson_loss_tie, 
                       lower = c(rep(0, ncol(is_basis)), -Inf))# Update the parameter to get the best one
  list(coeff = optim_result$par)
}

fit_ispline_model_Linear <- function(rank_train, Y_train,
                                     knots = c(0.25, 0.5, 0.75),  degree = 2) {
  
  is_basis <- iSpline(rank_train, knots = knots, Boundary.knots = c(0,1),  degree = degree) # Generate iSpline 
  
  linear_loss_tie <- function(par) {
    column_of_ones <- rep(1, nrow(is_basis))
    y_pred = cbind(is_basis, column_of_ones)%*%par
    return( MSE(Y_train, y_pred) )
    
    
  }
  initial_par <- rep(1, ncol(is_basis) + 1) # Initialize parameter
  optim_result <- hjkb(initial_par, linear_loss_tie, 
                       lower = c(rep(0, ncol(is_basis)), -Inf))# Update the parameter to get the best one
  list(coeff = optim_result$par)
}



fit_ispline_model_Robust <- function(rank_train, Y_train,
                                     knots = c(0.25, 0.5, 0.75),  degree = 2, delta) {
  
  is_basis <- iSpline(rank_train, knots = knots, Boundary.knots = c(0,1),  degree = degree) # Generate iSpline 
  
  huber_loss_tie <- function(par) {
    column_of_ones <- rep(1, nrow(is_basis))
    y_pred = cbind(is_basis, column_of_ones)%*%par
    return( sum(huber_loss(Y_train-y_pred, delta)) )
    
    
  }
  initial_par <- rep(1, ncol(is_basis) + 1) # Initialize parameter
  optim_result <- hjkb(initial_par, huber_loss_tie, 
                       lower = c(rep(0, ncol(is_basis)), -Inf))# Update the parameter to get the best one
  list(coeff = optim_result$par)
}



# Given a vector of ranks, predict responses using the fitted Ispline model.
ispline_predict <- function(rank_new, coeff, 
                            knots = c(0.25, 0.5, 0.75), degree = 2) {
  new_basis <- iSpline(rank_new, knots = knots,  Boundary.knots = c(0,1), degree = degree)
  as.vector(cbind(new_basis, 1) %*% coeff) 
}


# PAVA METHOD FUNCTIONS

fit_pava_model_Quantile <- function(rank_train, Y_train, tau = NA, solverinput = weighted.fractile) {
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


fit_pava_model_Poisson <- function(rank_train, Y_train, solverinput = poisson_solver) {
  gpava_result <- gpava(rank_train, Y_train, 
                        solver = solverinput, ties = "secondary")
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
  if(gpava_unique[1] == 0){
    gpava_unique[1] = unique(gpava_unique)[2]
    
  }
  y0 <- c(gpava_unique[1], gpava_unique[knots_idx])
  # Build a step function; here we use the index positions as knots.
  model_fun <- stepfun(knots_idx, y0)
  list(rank_unique = rank_unique, model = model_fun)
}



fit_pava_model_Linear <- function(rank_train, Y_train, solverinput = weighted.mean) {
  gpava_result <- gpava(rank_train, Y_train, 
                        solver = solverinput, ties = "secondary")
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

fit_pava_model_Robust <- function(rank_train, Y_train, solverinput = huber_solver) {
  gpava_result <- gpava(rank_train, Y_train, 
                        solver = solverinput, ties = "secondary")
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
# following gpava prediction method
pava_predict <- function(pava_model, rank_new) {
  sapply(rank_new, function(x) {
    idx <- sum(x >= pava_model$rank_unique)
    pava_model$model(idx)
  })
}

