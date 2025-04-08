library('quantreg')
library('splines2')
library('dfoptim')
library('isotone')


hrs = read.csv('HRS_ADL_IADL.csv')
hrs <- hrs[, -1]# delete ID column
hrs <- na.omit(hrs)
sample_size <- floor(0.8 * nrow(hrs))

numgender = ifelse(hrs$gender == "Female", 0, 1)
hrs$gender = numgender

numrace = ifelse(hrs$race == "Black", 0, 1)
hrs$race = numrace

nummar = ifelse(hrs$marriage == "Married/Partner", 0, 1)
hrs$marriage = nummar



# We scale the covariates according to the original dataset
hrs$AGE    = scale(hrs$AGE   )
hrs$educa = scale(hrs$educa )
hrs$networth = scale(hrs$networth)
hrs$cognition_catnew = scale(hrs$cognition_catnew )
hrs$bmi   = scale(hrs$bmi)
hrs$hlthrte = scale(hrs$hlthrte)
hrs$bloodp = scale(hrs$bloodp)
hrs$diabetes = scale(hrs$diabetes)
hrs$cancer = scale(hrs$cancer)
hrs$lung = scale(hrs$lung)
hrs$heart = scale(hrs$heart)
hrs$stroke = scale(hrs$stroke)
hrs$pchiat = scale(hrs$pchiat)
hrs$arthrit = scale(hrs$arthrit)
hrs$fall  = scale(hrs$fall )
hrs$pain = scale(hrs$pain)
hrs$A1c_adj = scale(hrs$A1c_adj)
hrs$CRP_adj = scale(hrs$CRP_adj)
hrs$CYSC_adj = scale(hrs$CYSC_adj)
hrs$HDL_adj = scale(hrs$HDL_adj)
hrs$TC_adj = scale(hrs$TC_adj)


train_indices <- sample(seq_len(nrow(hrs)), size = sample_size)
  
# Create training set
train_data <- hrs[train_indices, ]
traindataframe <- hrs[train_indices, c(-22)]
gscoretrain = hrs$score[train_indices ]
  
traindataframe0 = train_data[train_data$race==0, c(-25,-22) ]
gscore0train = train_data[train_data$race==0, 22]
  
traindataframe1 = train_data[train_data$race==1, c(-25,-22) ]
gscore1train = train_data[train_data$race==1, 22]
  
# rearrange the data
traindataframe = rbind(traindataframe0, traindataframe1)
gscoretrain = c(gscore0train, gscore1train)
  
Strain = c(rep(0, length(traindataframe0[,1]) ), rep(1,  length(traindataframe1[,1]) ) )
Ytrain = gscoretrain
  
  
# Create test set using the remaining rows
test_data <- hrs[-train_indices, ]
  
testdataframe <- hrs[-train_indices, c(-25,-22) ]
gscoretest = hrs[-train_indices, 22]
  
  
testdataframe0 = test_data[test_data$race==0, c(-25,-22) ]
gscore0test = test_data[test_data$race==0, 22]
  
testdataframe1 = test_data[test_data$race==1, c(-25,-22)]
gscore1test = test_data[test_data$race==1, 22]
  
  
testdataframe <- rbind(testdataframe0, testdataframe1)
gscoretest = c(gscore0test, gscore1test)
  
Stest= c(rep(0, length(testdataframe0[,1]) ), rep(1,  length(testdataframe1[,1]) ) )
Ytest = gscoretest
  
  
quantile_model0 =  glm(gscore0train ~ ., data = traindataframe0, family = poisson())
quantile_model1 =  glm(gscore1train ~ ., data = traindataframe1, family = poisson())
  
  
LYpredictions0 = predict(quantile_model0,  traindataframe0, type = "response")
LYpredictions1 = predict(quantile_model1,  traindataframe1, type = "response")
LYpredictions = c(LYpredictions0, LYpredictions1)
  
  
TYpredictions0 <- predict(quantile_model0,  newdata = testdataframe0, type = "response")
TYpredictions1 <- predict(quantile_model1,  newdata = testdataframe1, type = "response")
TYpredictions = c(TYpredictions0, TYpredictions1)


  
# Note that the training and test data should have been rearranged according to sensitive varible, respectively, before using this function.
# Here, we require that the input vectors Strain and Stest are sorted in increasing order (e.g., c(0, 0, 1, 1, 1)).
HRSI = Ifairpoi(LYpredictions, TYpredictions, Strain, Stest, Ytrain)
poissonloss(Ytrain,HRSI$fairtrain)
f_ks(HRSI$fairtrain, Strain )


# Note that the training and test data should have been rearranged according to sensitive varible, respectively, before using this function.
# Here, we require that the input vectors Strain and Stest are sorted in increasing order (e.g., c(0, 0, 1, 1, 1)).
HRSP = Pfair(LYpredictions, TYpredictions, Strain, Stest, Ytrain, method = poisson_solver)
poissonloss(Ytrain, HRSP$fairtrain)
f_ks(HRSP$fairtrain, Strain)
  


  
  

