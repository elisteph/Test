setwd ("/Users/sgyaase/Documents/Machine_learning")
###just testing
library(tidyverse) # dataset & viz tools
library(magrittr) # pipe operator
library(MASS) # statistical tools, such as lda
library(parallel) # parallel computing
library(ggplot2)
library(caret) # machine learning toolkit
library(ROCR) # ROC and AUC
library(Hmisc)
library(rms) # calibration slope
library(rmda) # decision curve
library(OmicsPLS)  # data integration toolkit
library(Matrix)
library(glmnet)
library(kernlab)
library(boot)
library(glmnet)  # for ridge regression
library(psych)   # for function tr() to compute trace of a matrix
library(pROC)
#install.packages("tree")
library(tree)



############################# DAY 1 ######################################
load("/Users/sgyaase/Downloads/DownSyndrome_day1.RData")

ls() # ls shows all objects in your workspace
## str gives an overview of all kinds of objects
cat("\nGlycomics data:\n")
str(glycomics_train)
#####Introduction
##### Before Analysis make a description of your data
boxplot(glycomics_train) #Boxplot
table(ClinicalVars_Day1$RPF)
summary(glycomics_train) #To check on the summary of the data
############# MACHINE LEARNING METHODS VRS. LOGISTIC##############
outc <- factor(ClinicalVars_Day1$RPF)

cat("----------- Logistic Regression\n----------------------\n")
fit_glm <- train(x = glycomics_train, outc, method="glmnet", family="binomial", 
                  tuneGrid = expand.grid(alpha = 0, lambda = 0))
fit_glm

cat("----------- Support Vector Machine\n----------------------\n")
(fit_svm <- train(x = glycomics_train, outc, method="svmLinear"))
confusionMatrix(outc, predict(fit_svm, glycomics_train)) ###This gives the prediction based on the train data

cat("----------- Neural Network\n----------------------\n")
(fit_nnet <- train(x = glycomics_train, outc, method="nnet", trace=F, preProcess = "range"))
confusionMatrix(outc, predict(fit_nnet, glycomics_train))
confusionMatrix(outc, predict(fit_glm,  n))

cat("----------- k-Nearest Neighbors\n----------------------\n")
(fit_knn <- train(x = glycomics_train, outc, method="knn"))
confusionMatrix(outc, predict(fit_knn, glycomics_train))

####### DAY 3:PROBABILITIES,PREDICTION AND DECISION ANALYSIS#######
load("/Users/sgyaase/Downloads/DownSyndrome_day1.RData")
#### Modelling Model: RPF ~ P1 + .. + P10
set.seed(7385872)
## Declare the outcome as a factor
outc <- factor(ClinicalVars_Day1$RPF)
outc_num <- as.numeric(outc) - 1

## Draw indices for test data
tst_i <- createDataPartition(outc, p = .25, list = FALSE)

## Fit glmnet and calculate predicted probabilities
fit_glm <- train(x = glycomics_train[-tst_i,], outc[-tst_i], method="glmnet", 
                 family="binomial", 
                 tuneGrid = expand.grid(alpha = 1, lambda=0))
# confusionMatrix(outc[tst_i], predict(fit_glm, glycomics_train[tst_i,]))
glm_predprob <- predict(fit_glm, glycomics_train[tst_i,], "prob")[,2]

## Convert it into a prediction object and calculate TPR/FPR 
glm.pred = prediction(glm_predprob, outc_num[tst_i])
glm.perf = performance(glm.pred,"tpr","fpr")
glm.perf

##compute area under ROC curve using TPR and FPR
glm.auc <- performance(glm.pred,"auc")
glm.auc <- unlist(slot(glm.auc, "y.values"))
plot(glm.perf, main="ROC Curve for logistic model on the test data",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")

## Calibration plot through the function val.prob from the rms package
v <- val.prob(p=glm_predprob,y=outc_num[tst_i])
c(auc=glm.auc,cal.slope=v["Slope"])

#testdata <- as_tibble(glycomics_train) %>% bind_cols(outc = outc)
#testdata$predicted.hd <-predict(heart.fit,newdata=testdata,type="response")
glm_dcdata <- tibble(y = outc_num[tst_i], p = glm_predprob)
glm_dc <- decision_curve(y~p, data = glm_dcdata, fitted.risk = TRUE, 
                         thresholds = seq(0, 1, by = .05), bootstraps = 50)
plot_decision_curve(glm_dc, legend.position = "topright")

## attr(calibrate(lrm(y ~ p, data = data.frame(p=glm_predprob,y=outc_num[tst_i]), x=T, y=T)), "predicted") 
## as the new predicted probabilities in val.prob.
####################### DAY 4 ################################
##glycomics_train <- scale(glycomics_train, scale = FALSE)
methylation_train <- scale(methylation_train, scale = FALSE)

boxplot(glycomics_train)

boxplot(methylation_train[,1:100])
print(ClinicalVars_train)
table(ClinicalVars_train$group_ds)

###### CROSS VALIDATION OF JOINT AND SPECIFIC COMPONENTS ######
out_cv <- crossval_o2m_adjR2(methylation_train, glycomics_train, 
                             a = 1:9, ax = 0:9, ay = 0:9, nr_folds = 5, nr_cores = 5)
print(out_cv)

#### Scree Plots
# par(mfrow=c(1,3))
# plot(svd(crossprod(methylation_train,glycomics_train),0,0)$d^2 %>% (function(e) e/sum(e)), main='Joint Scree plot')
# plot(svd(tcrossprod(methylation_train),0,0)$d %>% (function(e) e/sum(e)), main="Methylation Scree plot") 
# plot(svd(crossprod(glycomics_train),0,0)$d %>% (function(e) e/sum(e)), main="Glycomics Scree plot")
# par(mfrow=c(1,1))
# ## -> 3 6 4

r <- 2; rx <- 1; ry <- 0

crossval_sparsity(methylation_train, glycomics_train, r,rx,ry, 10, 
                  groupx = CpG_groups, keepx_seq = 1:10*40, keepy_seq = 10)

##### GO2PLS CROSS-VALIDATION ######
fit <- o2m(methylation_train, glycomics_train, r, rx, ry, sparse = TRUE, 
           groupx = CpG_groups, keepx = c(240, 120))
summary(fit)
### INSPECT LOADINGS #####
plot(fit,loading_name = "Yj",i=1,j=2,label = "col")
plot(fit,loading_name = "gr_Xj",i=1,j=2,label = "col", alpha = (apply(fit$W_gr,1,prod) > 0)) 

#### CHECKING IF JOINT COMPONENTS ARE ASSOCIATED #####
plot(data.frame(Tt=fit$Tt,U=fit$U), col = as.numeric(as.factor(ClinicalVars_train$group_ds)))
data.frame(Group = ClinicalVars_train$group_ds, JPC = scores(fit, "Xjoint")) %>% 
  pivot_longer(-Group,
               names_to = "Comp", 
               values_to = "LoadingValue") %>% 
  ggplot(aes(x=Comp, y=LoadingValue, col=Group)) + 
  geom_boxplot() + 
  theme_bw()

###### LOGISTIC REGRESSION
glm(as.factor(ClinicalVars_train$group_ds) ~ ., 
    data = data.frame(JPC=scores(fit, "Xjoint")), 
    family = "binomial") %>% 
  summary()

#### USING THE JOINT COMPONENTS TO TRAIN ML MODEL ####
(fit_caret <- train(data.frame(scores(fit, "Xjoint")), ClinicalVars_train$group_ds, 
                    method="knn"))
confusionMatrix(as.factor(ClinicalVars_train$group_ds), 
                predict(fit_caret, data.frame(scores(fit, "Xjoint"))), 
                positive = "yes")



################################## DAY 2 ###################################
#load("/Users/sgyaase/Downloads/Data_day2")

###### Generating a Data containing 8predictors for multivariate normal distribution for 100
#### subject, setting correlation to 0.5
corX <- matrix(0.5, 8, 8); diag(corX) =1;
nsize=100; mu=c(rep(0,8)); covX = corX;
X = mvrnorm(nsize, mu, covX)
mean(X); cov(X);
# Generate outcome variable y with effect size = (3,1.5,0,0,2,0,0,0), and error with sd(error) = 3.
Betas = c(3,1.5,0,0,2,0,0,0)
error <- rnorm(100, 0, 3)
y = X%*%Betas + error
hist(y)

dat = data.frame(cbind(X,y))
names(dat) <- c(paste("X",1:8,sep=""), "y")
write.csv (dat, "/Users/sgyaase/Documents/Machine_learning/dat1.csv", row.names=F)

########## Running Linear Regressions and checking for Assumptions#######
# Run the linear regression analysis using the generated/simulated data
dat <- read.csv("/Users/sgyaase/Documents/Machine_learning/dat1.csv")
lm.fit =lm(y~.,data=dat)
summary(lm.fit) ##Variables found to be significant at 0.05 are x1,x2 and x5
# Check if error is normally distributed.
err = y-predict(lm.fit)
plot(err)
qqnorm(err, pch = 1, frame = FALSE)
qqline(err, col = "steelblue", lwd = 2)
######## Comparing the correlation of X and the heat map for X
###the correlation is given in off the diagonal elements. The heatmap shows strong 
### color towards red, which indicates strong correlation.
cor(X)
heatmap(X)
###### Removing all variable lists and the generating with 10000 records
rm(list=ls())
corX <- matrix(0.5, 8, 8); diag(corX) =1;
nsize=10000; mu=c(rep(0,8)); covX = corX;
X = mvrnorm(nsize, mu, covX)
# Generate outcome variable y with effect size = (3,1.5,0,0,2,0,0,0), and error with sd(error) = 3.
Betas = c(3,1.5,0,0,2,0,0,0)
error <- rnorm(nsize, 0, 3)
y = X%*%Betas + error
hist(y)

dat = data.frame(cbind(X,y))
names(dat) <- c(paste("X",1:8,sep=""), "y")
write.csv (dat, "/Users/sgyaase/Documents/Machine_learning/dat2.csv", row.names=F)

#2b. Run the linear regression analysis using the generated/simulated data
dat <- read.csv("/Users/sgyaase/Documents/Machine_learning/dat2.csv")
lm.fit =lm(y~.,data=dat)
summary(lm.fit)
# Check if error is normally distributed.
err = y-predict(lm.fit)
plot(err)
qqnorm(err, pch = 1, frame = FALSE)
qqline(err, col = "steelblue", lwd = 2)


###### Validation Set Approach#####
#### Split the data into a training (2/3 of the whole data) and validation (1/3) sets.
dat <- read.csv("/Users/sgyaase/Documents/Machine_learning/dat1.csv")
train=sample(1:nsize, 0.66*nsize) #2/3=0.66
lm.fit.1 =lm(y~.,data=dat, subset=train)
summary(lm.fit.1)
mean((dat$y-predict(lm.fit.1,dat))[-train]^2)

###### Compare LOOCV and k-fold CV
# You can also use glm function. 
glm.fit=glm(y~.,data=dat)
coef(glm.fit)
lm.fit=lm(y~.,data=dat)
coef(lm.fit)

# Now run LOOCV and compute cv.error --> For a large sample size, this can take a long time.
# install.packages("boot")
### The input parameter responsible in the glm is the parameter alpha
glm.fit=glm(y~.,data=dat)
cv.err=cv.glm(dat,glm.fit)
cv.err$delta

# k-Fold Cross-Validation: 
set.seed(20200623)
cv.error.5=rep(0,5)
for (i in 1:5){
  glm.fit=glm(y~.,data=dat)
  cv.error.5[i]=cv.glm(dat,glm.fit,K=10)$delta[1]
}
cv.error.5

# Writing a simple function: boot.fn(inp1, inp2). What would be the output?
boot.fn=function(data,index)
  return(coef(lm(y~.,data=dat,subset=index)))
boot.fn(dat,1:10)

# Next we use the boot() function to compute the SE of 1000 bootstrap estimates.
boot(dat,boot.fn,1000)
# Compare with the asymptotic standard error 
summary(lm(y~.,data=dat))$coef
#–> Resampling 10 times is not enough to get robust results. If you increase the 
# number of resampling, the coefficients estimates become similar
#–> With 1000 resampling, the coefficents estimates of two methods are in a similar 
#scale. Focusing on standard error estimates, asymptotic error is much smaller. 
#In this case, where many predictors are correlated, the bootstrap estimates are to be preferred.

################ RIDGE REGRESSION##############
rm(list=ls())

set.seed(20240618)    # seed for reproducibility
data("mtcars")
# Center y, X will be standardized in the modelling function
y <- mtcars %>% dplyr::select(mpg) %>% scale(center = TRUE, scale = FALSE) %>% as.matrix()
X <- mtcars %>% dplyr::select(-mpg) %>% as.matrix()

# Perform 10-fold cross-validation to select lambda 
lambdas_to_try <- 10^seq(-3, 5, length.out = 100)

# Setting alpha = 0 implements ridge regression
ridge_cv <- cv.glmnet(X, y, alpha = 0, lambda = lambdas_to_try,
                      standardize = TRUE, nfolds = 10)
# Plot cross-validation results
plot(ridge_cv)

# Best cross-validated lambda
lambda_cv <- ridge_cv$lambda.min
# Fit final model, get its sum of squared residuals and multiple R-squared
model_cv <- glmnet(X, y, alpha = 0, lambda = lambda_cv, standardize = TRUE)
# Lecture slide 26 (yhat)
y_hat_cv <- predict(model_cv, X)
RSS_cv <- t(y - y_hat_cv) %*% (y - y_hat_cv)
rsq_ridge_cv <- cor(y, y_hat_cv)^2

# See how increasing lambda shrinks the coefficients --------------------------
# Each line shows coefficients for one variables, for different lambdas.
# The higher the lambda, the more the coefficients are shrinked towards zero.
res <- glmnet(X, y, alpha = 0, lambda = lambdas_to_try, standardize = FALSE)
plot(res, xvar = "lambda")
legend("bottomright", lwd = 1, col = 1:6, legend = colnames(X), cex = .7)

################ LASSO REGRESSION ################
# Perform 10-fold cross-validation to select lambda ---------------------------
lambdas_to_try <- 10^seq(-3, 5, length.out = 100)
# Setting alpha = 1 implements lasso regression
lasso_cv <- cv.glmnet(X, y, alpha = 1, lambda = lambdas_to_try,
                      standardize = TRUE, nfolds = 10)
# Plot cross-validation results
plot(lasso_cv)
# Best cross-validated lambda
lambda_cv <- lasso_cv$lambda.min
# Fit final model, get its sum of squared residuals and multiple R-squared
model_cv <- glmnet(X, y, alpha = 1, lambda = lambda_cv, standardize = TRUE)
y_hat_cv <- predict(model_cv, X)
RSS_cv <- t(y - y_hat_cv) %*% (y - y_hat_cv)
rsq_lasso_cv <- cor(y, y_hat_cv)^2


# See how increasing lambda shrinks the coefficients --------------------------
# Each line shows coefficients for one variables, for different lambdas.
# The higher the lambda, the more the coefficients are shrinked towards zero.
res2 <- glmnet(X, y, alpha = 1, lambda = lambdas_to_try, standardize = FALSE)
plot(res2, xvar = "lambda")
legend("bottomright", lwd = 1, col = 1:6, legend = colnames(X), cex = .7)

### Comparing the two models
rsq <- cbind("R-squared" = c(rsq_ridge_cv, rsq_lasso_cv))
rownames(rsq) <- c("ridge cross-validated",  "lasso cross-validated")
print(rsq)

############## ELASTIC NET###############
# Set training control
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats = 5,
                              search = "random",
                              verboseIter = TRUE)
# Train the model
elastic_net_model <- train(mpg ~ .,
                           data = cbind(y, X),
                           method = "glmnet",
                           preProcess = c("center", "scale"),
                           tuneLength = 25,
                           trControl = train_control)

# Check multiple R-squared
y_hat_enet <- predict(elastic_net_model, X)
rsq_enet <- cor(y, y_hat_enet)^2

# Compare all three methods
rsq <- cbind("R-squared" = c(rsq_ridge_cv, rsq_lasso_cv, rsq_enet))
rownames(rsq) <- c("ridge cross-validated",  "lasso cross-validated", "elastic net cross-validated")
print(rsq)

################## LOGISTIC REGRESSION(HEART DISEASES) ###############
Heart <- read.csv("/Users/sgyaase/Documents/Machine_learning/Data_day2/heart_comp.csv")
dim(Heart)
library(Hmisc)
Hmisc::describe(Heart) ### Performs like the describe function in STATA

heart.fit=glm(AHD~.,data=Heart, family="binomial")
summary(heart.fit)

###### LOOCV (Leave-one-out), 10-fold cross-validation.Two numbers in the delta
###First:standard k-fold CV estimate and 2nd is bias corrected version
cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)
LOOCV.err <- cv.glm(Heart,heart.fit, cost, K = nrow(Heart))$delta
LOOCV.err
cv.10.err <- cv.glm(Heart,heart.fit, cost, K = 10)$delta
cv.10.err

###### ROC and AUC ##############
#make this example reproducible
set.seed(1)
#use 70% of dataset as training set and 30% as test set
df <- Heart
sample <- sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.7,0.3))

df_train  <- df[sample, ]
df_test   <- df[!sample, ]

# fit logistic model
df_model <- glm(AHD~.,  family="binomial", data=df_train)

# predicted data on the test set
test_prediction <- predict(df_model, df_test, type="response")

# create roc curve
test_roc_object <- roc( df_test$AHD, test_prediction, plot=T, print.auc=T)
## Setting levels: control = 0, case = 1
## Setting direction: controls < cases

# predicted data on the training set
train_prediction <- predict(df_model, df_train, type="response")

# create roc curve
train_roc_object <- roc( df_train$AHD, train_prediction, plot=T, print.auc=T)
## Setting levels: control = 0, case = 1
## Setting direction: controls < cases

########################### TREE - BASED METHODS #####################
# Read the "heart_all.csv" data, and describe the data
Heart <- read.csv("heart_all.csv", sep=",", header=T, as.is = FALSE)
Hmisc::describe(Heart)

### Use the tree() function to fit a classification tree in order to predict Yes (here AHD = 1) in using all variables but AHD.
# The unpruned tree
tree.Heart=tree(AHD~.,Heart) 
# If you want to know more about the function tree(),
?tree 
summary(tree.Heart)

##### Estimating the Test error by spliting the Observations into training set and a test set to evaluate
#### performance on the test data
###### Display the tree structure ########
plot(tree.Heart)
text(tree.Heart,pretty=0, cex=0.7)
set.seed(20230620)
train=sample(1:nrow(Heart), trunc(nrow(Heart)*0.7)) 
Heart.test=Heart[-train,]
AHD.test=Heart.test$AHD
tree.Heart=tree(AHD~.,Heart,subset=train)
tree.pred=predict(tree.Heart,Heart.test,type="class")
table(tree.pred,AHD.test)
accuracy_value <- (37+38)/91    #The values can be different according to the random seeds, or splitting percentage.
accuracy_value

############################## PRUNNING THE TREE ##############################
##  The function cv.tree() performs cross-validation in order to determine the 
## optimal level of tree complexity; cost complexity pruning is used in order 
## to select a sequence of trees for consideration.
set.seed(20230620)
cv.Heart=cv.tree(tree.Heart,FUN=prune.misclass)
names(cv.Heart)
cv.Heart

#### dev corresponds to the cross-validation error rate in this case
par(mfrow=c(1,2))
plot(cv.Heart$size,cv.Heart$dev,type="b")
plot(cv.Heart$k,cv.Heart$dev,type="b")

### Now apply the prune.misclass() function in order to prune the tree to obtain the nine-node tree.
set.seed(20230620)
prune.Heart=prune.misclass(tree.Heart,best=9)
plot(prune.Heart)
text(prune.Heart,pretty=0)

tree.pred=predict(prune.Heart,Heart.test,type="class")
tree.pred
table(tree.pred,AHD.test)
# Compute accuracy of predicted class
accuracy_value <- (37+38)/91    #The values can be different according to the random seeds, or splitting percentage.
accuracy_value

######################### FITTING REGRESSION TREES ######################
## We fit a regression tree to the Boston data set. First, create a training set, and fit the tree to the training data.
# Hmisc::describe(Boston)
library(MASS)
set.seed(1)
train = sample(1:nrow(Boston), nrow(Boston)/2)
tree.boston=tree(medv~.,Boston,subset=train)
summary(tree.boston) ### Variables used in constructing the tree are:rm, lstat,crim and age
plot(tree.boston) ## True. When (rm>=6.543 and lstat<14.405), mdev=27.73 (compared to 21.38).
text(tree.boston,pretty=0)

### PRUNNING THE TREE AND MSE ####
cv.boston=cv.tree(tree.boston)
plot(cv.boston$size,cv.boston$dev,type='b')

prune.boston=prune.tree(tree.boston,best=5)
plot(prune.boston)
text(prune.boston,pretty=0)

yhat=predict(tree.boston,newdata=Boston[-train,])
boston.test=Boston[-train,"medv"]
plot(yhat,boston.test)
abline(0,1)

# Compute MSE of the test set
mean((yhat-boston.test)^2)
## CONCLUSION: –> The test set MSE associated with the unpruned tree is 35.29, 
   ## and with the pruned tree 35.90. So, in this case, use the unpruned tree.

#################  BAGGING (AIR QUALITY DATA) #########################
set.seed(20230620)

# Reading data
bagging_data=data.table(airquality)
# Explore the data (i.e. make a plot)
ggplot(bagging_data,aes(Wind,Ozone))+geom_point()+ggtitle("Ozone vs wind speed")
## Warning: Removed 37 rows containing missing values or values outside the scale range
## (`geom_point()`).

data_test=na.omit(bagging_data[,.(Ozone,Wind)])

# Split the data into a training and test set
train_index=sample.int(nrow(data_test),size=round(nrow(data_test)*0.66),replace = F)
data_test[train_index,train:=TRUE][-train_index,train:=FALSE]

# Single regression tree (without bagging)
no_bag_model=rpart(Ozone~Wind,data_test[train_index],control=rpart.control(minsplit=6))
rpart.plot(no_bag_model)
# Training of the bagged model
n_model=50
bagged_models=list()
for (i in 1:n_model)
{
  new_sample=sample(train_index,size=length(train_index),replace=T)
  bagged_models=c(bagged_models,list(rpart(Ozone~Wind,data_test[new_sample],control=rpart.control(minsplit=6))))
}

# Getting estimate from the bagged model
bagged_result=NULL
i=0
for (from_bag_model in bagged_models)
{
  if (is.null(bagged_result))
    bagged_result=predict(from_bag_model,bagging_data)
  else
    bagged_result=(i*bagged_result+predict(from_bag_model,bagging_data))/(i+1)
  i=i+1
}

####### VISUALIZATION ############
### Let’s summarize the results of single tree and bagged model in one plot. 
## Using ggplot2 we build the plots step by step: 1) scatterplot of ‘Wind’ and ‘Ozone’; 2) 
## The results from each of 50 bags; 3) bagged model (green line); 4) add the results of single tree (green).
gg=ggplot(data_test,aes(Wind,Ozone))+geom_point(aes(color=train))
gg

# plots of the results from each of the bags.
for (tree_model in bagged_models[1:50])
{
  prediction=predict(tree_model,bagging_data)
  data_plot=data.table(Wind=bagging_data$Wind,Ozone=prediction)
  gg=gg+geom_line(data=data_plot[order(Wind)],aes(x=Wind,y=Ozone),alpha=0.2)
}
gg

# plot of the bagged model
data_bagged=data.table(Wind=bagging_data$Wind,Ozone=bagged_result)
gg=gg+geom_line(data=data_bagged[order(Wind)],aes(x=Wind,y=Ozone),color='green')
gg

# plot of the results form a single regression tree
data_no_bag=data.table(Wind=bagging_data$Wind,Ozone=result_no_bag)
gg=gg+geom_line(data=data_no_bag[order(Wind)],aes(x=Wind,y=Ozone),color='red')
gg

####### RANDOM FOREST #############
## The dataset is taken from UCI website (https://archive.ics.uci.edu/ml/machine-learning-databases/car/), 
## and prepared in a readable format in R: car.rds. The data contains 7 variables – six explanatory 
## (Buying Price, Maintenance, NumDoors, NumPersons, BootSpace, Safety) and 
## one response variable (Condition). All the variables are categorical in nature and have 3-4 factor levels in each.

### PREPARING THE DATA ######
data1 <- readRDS("car.rds")
head(data1)
str(data1)
summary(data1)

# Split into Train and Validation sets
# Training Set : Validation Set = 70 : 30 (random)
set.seed(20230620)
train <- sample(nrow(data1), 0.7*nrow(data1), replace = FALSE)
TrainSet <- data1[train,]
ValidSet <- data1[-train,]

#### FITTING RF MODEL ####
# Random Forest model with default parameters
model1 <- randomForest(Condition ~ ., data = TrainSet, importance = TRUE)
model1

# For example 
model2 <- randomForest(Condition ~ ., data = TrainSet, ntree = 500, mtry = 6, importance = TRUE)
model2

# Predicting on train set
predTrain <- predict(model2, TrainSet, type = "class")
# Checking classification accuracy
mean(predTrain == TrainSet$Condition)  
table(predTrain, TrainSet$Condition)  

# Predicting on Validation set
predValid <- predict(model2, ValidSet, type = "class")
# Checking classification accuracy
mean(predValid == ValidSet$Condition)                    
table(predValid,ValidSet$Condition)


# To check important variables
importance(model2)        
varImpPlot(model2) 

##### MISSING DATA INPUTATION USING MISSFOREST ######
## For this exercise we will use “some-data-1.Rda”. 1) Explore the data. Are there any missing values? 2)
## Do the linear regression with BMI as outcome and Cholesterol as a covariate. 3) 
## (Randomly) Generate 20 missing values in the covariate Cholesterol values. 
## Repeat the linear regression. Compare the results of the original data and the 
## data with missing values. 4) Generate missing values in the covariates Smoking, Education, and Age. Save the data as “dat.mis.Rda”
library(VIM) # optional
library(mice)
library(missForest)
load("some-data-1.Rda")

# Study the dataset
head(dat)
str(dat)

#Check the data for missing values.
sapply(dat, function(x) sum(is.na(x)))

# fit a linear regression
fit <- lm(BMI~Cholesterol, data=dat)
summary(fit)
dat.mis <- dat
set.seed(27032019)
# missing in each variable.
dat.mis[sample(1:nrow(dat.mis), 20), "Cholesterol"] <- NA
dat.mis[sample(1:nrow(dat.mis), 10), "Smoking"] <- NA
dat.mis[sample(1:nrow(dat.mis), 20), "Education"] <- NA
dat.mis[sample(1:nrow(dat.mis), 50), "Age"] <- NA

save(dat.mis,file="dat.mis.Rda")

######  MULTIPLE IMPUTATION USING THE MICE PACKAGE ####
## One popular package for multiple imputation (MI) is “mice”. 
## MICE creates several Imputed/Complete data sets, runs each analysis on each 
## imputed dataset, and combines the estimates by taking uncertainty into account 
## using the Rubin’s rule. See https://stefvanbuuren.name/fimd/sec-nutshell.html

# Load the data
load("dat.mis.Rda")
# Study patterns of missing
sapply(dat.mis, function(x) sum(is.na(x)))
md.pattern(dat.mis)

# library(VIM): optional
aggr(dat.mis, col=c('navyblue','yellow'),
     numbers=TRUE, sortVars=TRUE,
     labels=names(dat.mis), cex.axis=.7,
     gap=3, ylab=c("Missing data","Pattern"))

###############################################################################
set.seed(28032019)
# Creates m=5 imputed datasets.
dat.imp <- mice(dat.mis, m=5, maxit = 50, method = 'pmm', seed = 500)
summary(dat.imp)

# Build predictive model, runs each analysis on each imputed dataset
fit.imp <- with(data = dat.imp, exp = lm(BMI~Cholesterol)) 

#combine the results of all 5 models
combine <- pool(fit.imp)
summary(combine)

# Delete missing values
dat.comp <- dat.mis[complete.cases(dat.mis),]
fit.del <- lm(BMI~Cholesterol, data=dat.comp)
summary(fit.del)

###### IMPUTATION USING MSISSFOREST ######
load("some-data-1.Rda")

#seed 10% missing values
dat.mis2 <- prodNA(dat, noNA = 0.3)
summary(dat.mis2)

#impute missing values, using all parameters as default values
dat.imp2 <- missForest(dat.mis2)
#check imputed values
dat.imp2$ximp
#check imputation error
dat.imp2$OOBerror

#comparing actual data accuracy
dat.err <- mixError(dat.imp2$ximp, dat.mis2, dat)
dat.err

# Do the linear regression 
fit.rf <- lm(BMI~Cholesterol, data=dat.imp2$ximp)
summary(fit.rf)

# Compare the linear regression results using different imputation methods.
#1. Data without any missing
summary(fit)
#2. Data by removing nay missing
summary(fit.del)
#3. Imputeda data with missForest
summary(fit.rf)
#4. Optional: using mice, generated 5 datasets, run the analysis for each of 5 datasets, and combine the results.
summary(combine)

######### ELASTIC NET ###
library(caret)
library(glmnet)

# using microbiome data
bacteria <- readRDS("Bacteria.RDS")
#str(bacteria)

# Splitting of the data into a training and testing set.
set.seed(2019)
inTraining <- createDataPartition(bacteria$Immune, p = .70, list = FALSE)
training <- bacteria[ inTraining,]
testing  <- bacteria[-inTraining,]

# Check if the split is correct in both datasets regarding the outcome 'immune'
prop.table(table(training$Immune))
prop.table(table(testing$Immune))

# set-up the CV strategy: 5-fold CV
cv_5 = trainControl(method = "cv", number = 5)

# use train() with method="glm", which will fit elastic net
bac_elnet = train(Immune ~ ., data = training,
                  method = "glmnet",
                  trControl = cv_5
)
bac_elnet

##### FINE - TUNING #######
# More controls
myControl <- trainControl(
  method = "repeatedcv", number = 10, repeats = 10,
  summaryFunction = twoClassSummary,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE,
  savePredictions = TRUE,
)

# Fit the model
set.seed(20220621)

model1 <- train(Immune ~ ., training, method = "glmnet", trControl = myControl)
## Warning in train.default(x, y, weights = w, ...): The metric "Accuracy" was not
## in the result set. ROC will be used instead.
model1
## For visualizing the results of grid search, and computing predicted values and cofusion matrix
ggplot(model1)

model1Class = predict (model1, newdata = testing)
model1Probs <- predict(model1, newdata = testing, type = "prob")

# To evaluate prediction performance
confusionMatrix( as.factor(model1Class), as.factor(testing$Immune))

########### USING RANDOM FOREST ###########
# Now we run random forest 
# for random forest, we need 'ranger' package.
library(ranger)
## 
## Attaching package: 'ranger'
## The following object is masked from 'package:randomForest':
## 
##     importance

set.seed(201906)

model_rf <- train(
  Immune ~ ., training,
  metric = "ROC",
  method = "ranger",
  trControl = myControl
)

model_rf

# Make a list containing models to be compared

modelList <- list(
  glmnet1=model1,
  rf = model_rf
)
bwplot(resamples(modelList),
       metric = "ROC")

####### DATA PREPARATION & XGboost MODEL ######
library(xgboost)
## Attaching package: 'xgboost'
## The following object is masked from 'package:dplyr':
##     slice
library(caTools)
library(dplyr)
library(caret)

head(iris)

set.seed(42)
sample_split <- sample.split(Y = iris$Species, SplitRatio = 0.7)
train_set <- data.table(subset(x = iris, sample_split == TRUE))
test_set <- data.table(subset(x = iris, sample_split == FALSE))

y_train <- as.integer(train_set$Species) - 1
y_test <- as.integer(test_set$Species) - 1
X_train <- train_set[,-"Species"]
X_test <- test_set[,-"Species"]

### We also have to specify the parameters for the XGBoost model. You can learn 
## more about all the available parameters: see https://xgboost.readthedocs.io/en/latest/parameter.html
xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
xgb_params <- list(
  booster = "gbtree",
  eta = 0.01,
  max_depth = 8,
  gamma = 4,
  subsample = 0.75,
  colsample_bytree = 1,
  objective = "multi:softprob",
  eval_metric = "mlogloss",
  num_class = length(levels(iris$Species))
)

xgb_model <- xgb.train(
  params = xgb_params,
  data = xgb_train,
  nrounds = 5000,
  verbose = 1
)
xgb_model

#### PREDICTION & EVALUATION
## You can use the predict() function to make predictions with the XGBoost model.
xgb_preds <- predict(xgb_model, as.matrix(X_test), reshape = TRUE)
xgb_preds <- as.data.frame(xgb_preds)
colnames(xgb_preds) <- levels(iris$Species)
xgb_preds

xgb_preds$PredictedClass <- apply(xgb_preds, 1, function(y) colnames(xgb_preds)[which.max(y)])
xgb_preds$ActualClass <- levels(iris$Species)[y_test + 1]
xgb_preds

accuracy <- sum(xgb_preds$PredictedClass == xgb_preds$ActualClass) / nrow(xgb_preds)
accuracy

confusionMatrix(factor(xgb_preds$ActualClass), factor(xgb_preds$PredictedClass))





