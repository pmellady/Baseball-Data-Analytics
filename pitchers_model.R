# Import MASS for the polr function
library(MASS)
library(rpart)
library(randomForest)
library(rpart.plot)
library(gt)
library(ggplot2)

NCAA_tots<-read.csv("pitchers_data/wcl_pitchers.csv", header=TRUE)
NCAA_tots$Outcome<-factor(NCAA_tots$Outcome, levels=c("None", "Minors", "Majors"))


# Fitting and Prediction #######################################################
################################################################################
# Parametric Model #############################################################
# Split the data
train<-NCAA_tots[NCAA_tots$Age>20,]
test<-NCAA_tots[NCAA_tots$Age<=20,]

# Fit the model and predict on the training set
mod_age<-polr(Outcome~(Age+BF+BB+ER+H+IP+SO+HR)^3+(W+L)^2, data=train, Hess=TRUE)
temp<-predict(mod_age, train)

# Calculate the training error
train_err<-1-mean(temp==train$Outcome)

# How did we missclassify?
## How did we missclassify "None" players
n_n<-sum(train[temp=="None",]$Outcome=="None")/nrow(train[train$Outcome=="None",]) 
n_mn<-sum(train[temp=="Minors",]$Outcome=="None")/nrow(train[train$Outcome=="None",])
n_mj<-sum(train[temp=="Majors",]$Outcome=="None")/nrow(train[train$Outcome=="None",])

## How did we missclassify "Minors" players
mn_n<-sum(train[temp=="None",]$Outcome=="Minors")/nrow(train[train$Outcome=="Minors",])
mn_mn<-sum(train[temp=="Minors",]$Outcome=="Minors")/nrow(train[train$Outcome=="Minors",])
mn_mj<-sum(train[temp=="Majors",]$Outcome=="Minors")/nrow(train[train$Outcome=="Minors",])

## How did we missclassify "Majors" players
mj_n<-sum(train[temp=="None",]$Outcome=="Majors")/nrow(train[train$Outcome=="Majors",])
mj_mn<-sum(train[temp=="Minors",]$Outcome=="Majors")/nrow(train[train$Outcome=="Majors",])
mj_mj<-sum(train[temp=="Majors",]$Outcome=="Majors")/nrow(train[train$Outcome=="Majors",])

per<-function(x){
  paste0(100*round(x,3),"%")
}

errors<-data.frame(matrix(c(per(n_n), per(n_mn), per(n_mj), 
                            per(mn_n), per(mn_mn), per(mn_mj), 
                            per(mj_n), per(mj_mn), per(mj_mj)),nrow=3, byrow = TRUE))
errors<-cbind(c("Obs. None", "Obs. Minors", "Obs. Majors"), errors)
colnames(errors)<-c(" ", "Pred. None", "Pred. Minors", "Pred. Majors")
gt(errors)|>tab_header("Accuracy of PO Model for Pitchers")

# Where do we expect the newer players to end up?
preds<-predict(mod_age, test)

## Find the placement players who haven't been drafted yet
table<-cbind(test[,c(2,1)],preds)
noms<-c("Name", "Current Placement", "Projected Placement")
colnames(table)<-noms
fin<-cbind(train[,c(2,1)], temp)
colnames(fin)<-noms
table<-rbind(table,fin)

write.csv(table,"pitcher_PO_predictions.csv", row.names = FALSE)


# Non-Parametric Method ########################################################
## Random Forest
forest<-randomForest(Outcome~Age+BF+BB+ER+H+IP+SO+HR+W+L+ERA+WHIP+BB9+H9+HR9+SO9, data=train, ntree=10)
temp<-predict(forest, train)

# Calculate the training error
train_err<-1-mean(temp==train$Outcome)

# How did we missclassify?
## How did we missclassify "None" players
n_n<-sum(train[temp=="None",]$Outcome=="None")/nrow(train[train$Outcome=="None",]) 
n_mn<-sum(train[temp=="Minors",]$Outcome=="None")/nrow(train[train$Outcome=="None",])
n_mj<-sum(train[temp=="Majors",]$Outcome=="None")/nrow(train[train$Outcome=="None",])

## How did we missclassify "Minors" players
mn_n<-sum(train[temp=="None",]$Outcome=="Minors")/nrow(train[train$Outcome=="Minors",])
mn_mn<-sum(train[temp=="Minors",]$Outcome=="Minors")/nrow(train[train$Outcome=="Minors",])
mn_mj<-sum(train[temp=="Majors",]$Outcome=="Minors")/nrow(train[train$Outcome=="Minors",])

## How did we missclassify "Majors" players
mj_n<-sum(train[temp=="None",]$Outcome=="Majors")/nrow(train[train$Outcome=="Majors",])
mj_mn<-sum(train[temp=="Minors",]$Outcome=="Majors")/nrow(train[train$Outcome=="Majors",])
mj_mj<-sum(train[temp=="Majors",]$Outcome=="Majors")/nrow(train[train$Outcome=="Majors",])

per<-function(x){
  paste0(100*round(x,3),"%")
}

errors<-data.frame(matrix(c(per(n_n), per(n_mn), per(n_mj), 
                            per(mn_n), per(mn_mn), per(mn_mj), 
                            per(mj_n), per(mj_mn), per(mj_mj)),nrow=3, byrow = TRUE))
errors<-cbind(c("Obs. None", "Obs. Minors", "Obs. Majors"), errors)
colnames(errors)<-c(" ", "Pred. None", "Pred. Minors", "Pred. Majors")
gt(errors)|>tab_header("Accuracy of PO Model for Pitchers")

# Where do we expect the newer players to end up?
preds<-predict(forest, test)

## Find the placement players who haven't been drafted yet
table<-cbind(test[,c(2,1)],preds)
noms<-c("Name", "Current Placement", "Projected Placement")
colnames(table)<-noms
fin<-cbind(train[,c(2,1)], temp)
colnames(fin)<-noms
table<-rbind(table,fin)

write.csv(table,"predictions/pitcher_tree_predictions.csv", row.names = FALSE)

set.seed(922)
gt(table[sample(c(TRUE, FALSE), nrow(table), replace=TRUE, prob=c(.01, .999)),]) |>
  tab_header(title = "Players Actual vs Predicted Placements" )



