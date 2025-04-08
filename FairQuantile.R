#library('quantreg')
#library('splines2')
#library('dfoptim')
#library('isotone')

community = read.csv('communityours.csv')


qu=0.75

sample_size <- floor(0.8 * nrow(community)) # around 80% train 20% test

train_indices <- sample(seq_len(nrow(community)), size = sample_size)# To get some training data



# Create training set
train_data <- community[train_indices, ]
traindata_corvariates <- community[train_indices, 1:96 ]


traindata_corvariates_0 = train_data[train_data$race==0, 1:96 ]
ViolentCrimesPerPoptrain_0 = train_data[train_data$race==0, 97 ]

traindata_corvariates_1 = train_data[train_data$race==1, 1:96 ]
ViolentCrimesPerPoptrain_1 = train_data[train_data$race==1, 97 ]


# rearrange the data
traindata_corvariates = rbind(traindata_corvariates_0, traindata_corvariates_1)
ViolentCrimesPerPoptrain = c(ViolentCrimesPerPoptrain_0, ViolentCrimesPerPoptrain_1 )

Strain = c(rep(0, length(traindata_corvariates_0[,1]) ), rep(1,  length(traindata_corvariates_1[,1]) ) )
Ytrain = c(ViolentCrimesPerPoptrain_0, ViolentCrimesPerPoptrain_1)


# Create test set using the remaining rows
test_data <- community[-train_indices, ]
testdata_covariates <- community[-train_indices, 1:96 ]


testdata_corvariates_0 = test_data[test_data$race==0, 1:96 ]
ViolentCrimesPerPoptest_0 = test_data[test_data$race==0, 97 ]

testdata_corvariates_1 = test_data[test_data$race==1, 1:96 ]
ViolentCrimesPerPoptest_1 = test_data[test_data$race==1, 97 ]


testdata_covariates <- rbind( testdata_corvariates_0, testdata_corvariates_1  )
ViolentCrimesPerPoptest = c(ViolentCrimesPerPoptest_0, ViolentCrimesPerPoptest_1)

Stest= c(rep(0, length(testdata_corvariates_0[,1]) ), rep(1,  length(testdata_corvariates_1[,1]) ) )
Ytest = c( ViolentCrimesPerPoptest_0, ViolentCrimesPerPoptest_1)



# To estimate f* for s = 0 and s = 1 respectively
quantile_model0 =  rq(ViolentCrimesPerPoptrain_0 ~ ., data = traindata_corvariates_0, tau = qu)
quantile_model1 =  rq(ViolentCrimesPerPoptrain_1 ~ ., data = traindata_corvariates_1, tau = qu)

LYpredictions0 = predict(quantile_model0,  traindata_corvariates_0)
LYpredictions1 = predict(quantile_model1,  traindata_corvariates_1)
LYpredictions = c(LYpredictions0, LYpredictions1)

TYpredictions0 <- predict(quantile_model0,  newdata = testdata_corvariates_0)
TYpredictions1 <- predict(quantile_model1,  newdata = testdata_corvariates_1)
TYpredictions = c(TYpredictions0, TYpredictions1)





# Note that the training and test data should have been rearranged according to sensitive variable, respectively, before using this function.
# Here, we require that the input vectors Strain and Stest are sorted in increasing order (e.g., c(0, 0, 1, 1, 1)).
CommunityI = Ifairqua(LYpredictions, TYpredictions, Strain, Stest, Ytrain)
quantile_loss(Ytest, CommunityI$fairtest, qu)
f_ks(CommunityI$fairtrain, Strain)


# Note that the training and test data should have been rearranged according to sensitive variable, respectively, before using this function.
# Here, we require that the input vectors Strain and Stest are sorted in increasing order (e.g., c(0, 0, 1, 1, 1)).
CommunityP = Pfair(LYpredictions, TYpredictions, Strain, Stest, Ytrain, method = weighted.fractile, quantile = qu)
quantile_loss(Ytest, CommunityP$fairtest, qu)
f_ks(CommunityP$fairtest, Stest)
