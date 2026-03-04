# FairLearning

R implementation for the paper **Fair Regression under Demographic Parity: A Unified Framework**.

We propose a unified framework for fair regression formulated as risk minimization subject to a demographic parity constraint. Unlike many existing approaches that are tied to specific loss functions or require challenging non-convex optimization, our framework applies to a wide range of tasks, including squared-loss regression, binary classification with cross-entropy, quantile regression with pinball loss, and robust regression with Huber loss.
Our key result is a characterization of the fair risk minimizer, which leads to a computationally efficient procedure for solving general fair regression problems. We establish asymptotic consistency and convergence rates under mild assumptions, and we demonstrate the method’s versatility through both theoretical discussions and numerical experiments, showing that it can effectively minimize risk while satisfying fairness constraints across diverse regression settings.

---

## Method (two-step)

**Step 1.** Fit a base predictor: $\hat f(x,s)$.

**Step 2.** Compute the rank transform $\hat U = \hat F_s(\hat f(x,s))$, then learn a **shared monotone calibration** function $\hat Q$ on $[0,1]$. The final fair predictor is
\[
\hat f_{\text{fair}}(x,s)=\hat Q\left(\hat F_s(\hat f(x,s))\right).
\]

This enforces demographic parity by construction while optimizing the chosen loss.

We provide two solvers for Step 2:
- **PAVA isotonic** (produces a piecewise-constant monotone $\hat Q$)
- **I-spline** (produces a smooth monotone $\hat Q$)
  
---

## Descriptions of main R functions 

We learn a non-decreasing function $\hat Q$ that maps $\hat U$ to predictions, using either a PAVA (step function) solver or a smooth I-spline solver.

### Solver A: PAVA (piecewise-constant monotone Q)
- `fit_pava_model_Quantile(rank_train, Y_train, tau, solverinput=weighted.fractile)`  
  Fit monotone $\hat Q$ under pinball loss (quantile regression) via generalized PAVA.
- `fit_pava_model_Poisson(rank_train, Y_train, solverinput=poisson_solver)`  
  Fit monotone $\hat Q$ under Poisson loss via generalized PAVA.
- `fit_pava_model_Robust(rank_train, Y_train, solverinput=huber_solver)`  
  Fit monotone $\hat Q$ under Huber loss via generalized PAVA.
- `pava_predict(pava_model, rank_new)`  
  Predict $\hat Q$(rank_new) from a fitted PAVA model.

### Solver B: I-spline (smooth monotone Q)
- `fit_ispline_model_Quantile_LP(rank_train, Y_train, tau, knots, degree)`  
   Fit $\hat Q$ under pinball loss with I-spline basis and nonnegative spline coefficients.
- `fit_ispline_model_Poisson_Convex(rank_train, Y_train, knots, degree, solver)`  
   Fit $\hat Q$ under Poisson loss with I-spline basis and nonnegative spline coefficients.
- `fit_ispline_model_Robust_Convex(rank_train, Y_train, knots, degree, delta, solver)`  
   Fit $\hat Q$ under Huber loss with I-spline basis and nonnegative spline coefficients.
- `ispline_predict(rank_new, coeff, knots, degree)`  
  Predict $\hat Q$ (rank_new) given fitted spline coefficients.

---

## Files

- `fairfunctions.R` — core functions for the fair processing pipeline
- `FairQuantileCrime.R` — fair quantile regression demo (with CRIME dataset)
- `FairPoissonHRS.R` — fair Poisson regression demo (with HRS dataset)
- `FairRobustSimulation.R` — fair Huber robust regression simulation
- `data/` —including CRIME and HRS datasets we use

