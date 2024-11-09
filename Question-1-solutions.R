library(glmnet)
library(egg)
library(ggplot2)
# Predicting strike given the batter did not swing #############################
## Read and clean the data
pitch<-read.csv("pitch_data.csv", header=TRUE)
noms<-colnames(pitch)
pitch[,3:17] <- lapply(pitch[,3:17] , factor)
pitch$spin_rate<-as.numeric(pitch$spin_rate)
pitch<-pitch[is.na(pitch$spin_rate)==FALSE,]

## Transform the data to add indicators for the events in question
strike<-pitch$is_strike
swing<-pitch$is_swing
A<-rep(0,nrow(pitch))
B<-rep(0,nrow(pitch))
pitch<-cbind(A, B, pitch)
colnames(pitch)<-c("A", "B", noms)

pitch[pitch$is_strike==1 & pitch$is_swing==0, 1]<-1
pitch[pitch$is_strike==0 & pitch$is_swing==0, 2]<-1
pitch<-pitch[,c(1,2,5:25)]

## Make the design matrices for the two models we need to fit
xfactors1<-model.matrix(A~.-1,pitch[, c(1, 3:17)])
X1<-as.matrix(data.frame(xfactors1, as.matrix(pitch[, 18:23]), 
                         as.matrix(pitch[, 18]^2), 
                         as.matrix(pitch[, 19]^2), 
                         as.matrix(pitch[, 18]*pitch[, 19])))

xfactors2<-model.matrix(B~.-1,pitch[, 2:17])
X2<-as.matrix(data.frame(xfactors2, as.matrix(pitch[, 18:23]), 
                         as.matrix(pitch[, 18]^2), 
                         as.matrix(pitch[, 19]^2), 
                         as.matrix(pitch[, 18]*pitch[, 19])))

## Extract the response vectors and fit the regularized models
y1<-as.vector(pitch$A)
mod1<-glmnet(X1, y1, family="binomial", alpha=1, standardize = TRUE, intercept=TRUE) 

y2<-as.vector(pitch$B)
mod2<-glmnet(X2, y2, family="binomial", alpha=1, standardize = TRUE, intercept=TRUE)

## Find and append the predicted probabilities to the data frame
p_10<-predict(mod1, X1, s=mod1$lambda[length(mod1$lambda)],type="response")
p_00<-predict(mod2, X2, s=mod1$lambda[length(mod1$lambda)], type="response")

pitch<-cbind(pitch, p_10/(p_10+p_00))
pitch$A<-strike
pitch$B<-swing
colnames(pitch)<-c("is_strike", "is_swing", 
                   colnames(pitch)[-c(1,2,length(colnames(pitch)))], "prob")

## Calculate the brier score and brier skill score
brier<-sum((pitch[pitch$is_swing==0,]$is_strike-pitch[pitch$is_swing==0,]$prob)^2)/
  nrow(pitch[pitch$is_swing==0,])
b_ref<-.25
bss<-1-brier/b_ref

## Calculate the error of the model
err<-sum(pitch[pitch$is_swing==0,1]!=round(pitch[pitch$is_swing==0,24]))/
  nrow(pitch[pitch$is_swing==0,])

# Projecting 2025 strikeout rate ###############################################
## Import the data
X<-data.frame(read.csv("pitcher_strikeout_data.csv", header=TRUE))

## Transform the data
noms<-colnames(X)
ids<-unique(X$pitcher_id)
years<-unique(X$year)
mat<-matrix(rep(0,8),nrow=1)
for(id in ids){
  for(yr in years){
    temp<-c(yr, id, mean(X[X$pitcher_id==id & X$year==yr,]$age), 
            mean(X[X$pitcher_id==id & X$year==yr,]$stuff),
            sum(X[X$pitcher_id==id & X$year==yr,]$is_strikeout), 
            sum(X[X$pitcher_id==id & X$year==yr,]$is_walk),
            length(X[X$pitcher_id==id & X$year==yr,]$age)-
              sum(X[X$pitcher_id==id & X$year==yr,]$is_strikeout)-
              sum(X[X$pitcher_id==id & X$year==yr,]$is_walk),
            length(X[X$pitcher_id==id & X$year==yr,]$age))
    mat<-rbind(mat, temp)
  }
}
mat<-data.frame(mat[-1,])
colnames(mat)<-c(noms, "other", "tot")
mat$pitcher_id<-as.factor(mat$pitcher_id)

## Make design matrix
Xr<-model.matrix(is_strikeout~.-1, mat[,c(2,5)])
Xr<-as.matrix(data.frame(Xr,as.matrix(mat[,c(1,3,4)])))

## Extract response vectors and fit the model
y<-matrix(c(as.vector(mat$is_strikeout), 
            as.vector(mat$is_walk), 
            as.vector(mat$other)), ncol=3)
mod<-glmnet(Xr, y, family="multinomial", 
            penalty.factor = c(rep(0.0001,200),0,0,0))

## Predict average stuff in 2025
pred_stuff<-c()
for(id in ids){
  x<-cbind(rep(1, length.out=length(X[X$pitcher_id==id,]$stuff)), 
           X[X$pitcher_id==id,]$year)
  beta<-solve(t(x)%*%x)%*%t(x)%*%(X[X$pitcher_id==id, ]$stuff)
  pred_stuff<-c(pred_stuff, t(beta)%*%c(1, 2025))
}

pred_rates<-matrix(rep(0,4), nrow=1)
for(i in 1:length(ids)){
  pitcher<-c(ids[i], 2025, max(mat[mat$pitcher_id==ids[i],]$age)+1, pred_stuff[i])
  pred_rates<-rbind(pred_rates, pitcher)
}

## Make design matrix for prediction
pred_rates<-pred_rates[-1,]
pred_rates<-data.frame(pred_rates)
colnames(pred_rates)<-c("pitcher_id", "year", "age", "stuff")
pred_rates$pitcher_id<-as.factor(pred_rates$pitcher_id)
Xp<-model.matrix(year~pitcher_id-1, pred_rates)
Xp<-as.matrix(data.frame(Xp,as.matrix(pred_rates[,c(2,3,4)])))

## Predict strikeout rates for 2025
rates<-predict(mod, Xp, s=mod$lambda[length(mod$lambda)], type="response")

## Visualize the predictions
### Extract each pitchers strikeout rate and store it
temp<-merge(rates, pred_rates, by=0)
colnames(temp)<-c("temp","K%", "BB%", "Other%", 
                  "pitcher_id", "year", "age", "stuff")
mat2<-matrix(c(0,0,0,0), nrow=1)
for(id in ids){
  mat2<-rbind(mat2,c(as.numeric(id), temp[temp$pitcher_id==id, ]$`K%`, 
                     temp[temp$pitcher_id==id, ]$`BB%`, 
                     temp[temp$pitcher_id==id, ]$`Other%`))
}
mat2<-data.frame(mat2[-1,])
colnames(mat2)<-c("pitcher_id", "K%", "BB%", "Other%")

### Select four random pitchers to visualize and collect their data
pits<-sample(ids, 4)
probs<-matrix(rep(0,4), nrow=1)
men<-matrix(rep(0,8), nrow=1)
for(id in pits){
  probs<-rbind(probs, as.matrix(mat2[mat2$pitcher_id==id,]))
  men<-rbind(men, as.matrix(mat[mat$pitcher_id==id,]))
}
probs<-data.frame(probs[-1,])
men<-men[-1,]
colnames(probs)<-c("pitcher_id", "K%", "BB%", "Other%")
colnames(men)<-c(noms, "other", "tot")

guys<-matrix(rep(0,ncol(men)*nrow(men)), nrow=nrow(men))
for(i in 1:ncol(men)){
  guys[,i]<-as.numeric(men[,i])
}

guys<-data.frame(cbind(guys[,1:4], round(guys[,5:7]/guys[,8],3)))
colnames(guys)<-c(noms, "other")

k_plot1<-c(guys[guys$pitcher_id==pits[1],]$is_strikeout, 
           probs[probs$pitcher_id==pits[1],]$'K%')
bb_plot1<-c(guys[guys$pitcher_id==pits[1],]$is_walk, 
            probs[probs$pitcher_id==pits[1],]$'BB%')

### Generate the plots for the four pitchers
plot1 <- ggplot() +
  geom_line(aes(2020:2025, k_plot1))+
  geom_line(aes(2020:2025, bb_plot1), color="orange")+
  ggtitle(paste("ID no. ", probs[probs$pitcher_id==pits[1],],"K and BB Rates by Year"))+
  xlab("Year")+ylab("Rate")

k_plot2<-c(guys[guys$pitcher_id==pits[2],]$is_strikeout, 
           probs[probs$pitcher_id==pits[2],]$'K%')
bb_plot2<-c(guys[guys$pitcher_id==pits[2],]$is_walk, 
            probs[probs$pitcher_id==pits[2],]$'BB%')

plot2 <- ggplot() +
  geom_line(aes(2020:2025, k_plot2))+
  geom_line(aes(2020:2025, bb_plot2), color="orange")+
  ggtitle(paste("ID no. ", probs[probs$pitcher_id==pits[2],],"K and BB Rates by Year"))+
  xlab("Year")+ylab("Rate")

k_plot3<-c(guys[guys$pitcher_id==pits[3],]$is_strikeout, 
           probs[probs$pitcher_id==pits[3],]$'K%')
bb_plot3<-c(guys[guys$pitcher_id==pits[3],]$is_walk, 
            probs[probs$pitcher_id==pits[3],]$'BB%')

plot3 <- ggplot() +
  geom_line(aes(2020:2025, k_plot3))+
  geom_line(aes(2020:2025, bb_plot3), color="orange")+
  ggtitle(paste("ID no. ", probs[probs$pitcher_id==pits[3],],"K and BB Rates by Year"))+
  xlab("Year")+ylab("Rate")

k_plot4<-c(guys[guys$pitcher_id==pits[4],]$is_strikeout, 
           probs[probs$pitcher_id==pits[4],]$'K%')
bb_plot4<-c(guys[guys$pitcher_id==pits[4],]$is_walk, 
            probs[probs$pitcher_id==pits[4],]$'BB%')

plot4 <- ggplot() +
  geom_line(aes(2020:2025, k_plot4))+
  geom_line(aes(2020:2025, bb_plot4), color="orange")+
  ggtitle(paste("ID no. ", probs[probs$pitcher_id==pits[4],],"K and BB Rates by Year"))+
  xlab("Year")+ylab("Rate")

### Plot the results
ggarrange(plot1, plot2, plot3, plot4)

## Test the model for 2024 vs career average strikeout rate
### Make design matrix for the model that predicts in 2024
mat1<-mat[mat$year!=2024,]
Xt<-model.matrix(is_strikeout~.-1, mat1[,c(2,5)])
Xt<-as.matrix(data.frame(Xt,as.matrix(mat1[,c(1,3,4)])))

### Extract the response vectors and fit the model
yt<-matrix(c(as.vector(mat1$is_strikeout),
             as.vector(mat1$is_walk), 
             as.vector(mat1$other)), ncol=3)
mod1<-glmnet(Xt, yt, family="multinomial", 
             penalty.factor = c(rep(0.0001,200),0,0,0))

### Predict stuff for 2024
pred_stuff1<-c()
for(id in ids){
  x<-cbind(rep(1, length.out=length(X[X$pitcher_id==id & X$year!=2024,]$stuff)),
           X[X$pitcher_id==id & X$year!=2024,]$year)
  beta<-solve(t(x)%*%x)%*%t(x)%*%(X[X$pitcher_id==id & X$year!=2024, ]$stuff)
  pred_stuff1<-c(pred_stuff1, t(beta)%*%c(1, 2024))
}

### Make design matrix to predict the rates for 2024
pred_rates1<-matrix(rep(0,4), nrow=1)
for(i in 1:length(ids)){
  pitcher<-c(ids[i], 2024, max(mat1[mat1$pitcher_id==ids[i],]$age)+1, pred_stuff1[i])
  pred_rates1<-rbind(pred_rates1, pitcher)
}
pred_rates1<-pred_rates1[-1,]
pred_rates1<-data.frame(pred_rates1)
colnames(pred_rates1)<-c("pitcher_id", "year", "age", "stuff")
pred_rates1$pitcher_id<-as.factor(pred_rates1$pitcher_id)
Xpt<-model.matrix(year~pitcher_id-1, pred_rates1)
Xpt<-as.matrix(data.frame(Xpt,as.matrix(pred_rates1[,c(2,3,4)])))

### Predict the 2024 rates and make a data frame
rates1<-predict(mod1, Xpt, s=mod1$lambda[length(mod1$lambda)], type="response")
rates1<-data.frame(rates1)
colnames(rates1)<-c("K%", "BB%", "Other%")

### Find each pitchers career average rates
k_ref<-c()
for(id in ids){
  k_ref<-c(k_ref, sum(mat1[mat1$pitcher_id==id,]$is_strikeout)/
             sum(mat1[mat1$pitcher_id==id,]$tot))
}

### Find each pitchers 2024 strikeout rates
k_2024<-mat[mat$year==2024,]$is_strikeout/mat[mat$year==2024,]$tot

### Find the ratio of MSEs for the two estimates
mean((k_2024-rates1[,1])^2)/mean((k_ref-k_2024)^2)








