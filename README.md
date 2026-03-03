# FairLearning

Code for our paper **“Fair Regression under Demographic Parity: A Unified Framework”**.  
We study fair prediction under **demographic parity** for **general regression tasks/losses** (e.g., squared, pinball/quantile, Poisson NLL, Huber).

## Method (two-step post-processing)

1) Fit a base predictor: `\hat f(x,s)`  
2) Build a **shared monotone calibration** `\hat Q` using `\hat U = \hat F_s(\hat f(x,s))`, then output  
`f_fair(x,s) = \hat Q( \hat F_s(\hat f(x,s)) )`  
This enforces demographic parity by construction while optimizing the target loss.

## Files

- `fairfunctions.R` — core functions for the fair processing pipeline
- `FairQuantileCrime.R` — fair quantile regression demo (with CRIME dataset)
- `FairPoissonHRS.R` — fair Poisson regression demo (with HRS dataset)
- `FairRobustSimulation.R` — fair Huber robust regression simulation
- `data/` —including CRIME and HRS datasets we use
