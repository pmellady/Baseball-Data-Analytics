# Import packages and data for analysis
library(glmnet)
data<-read.csv("data.csv", header=TRUE)

# Partition the pitch types into the three desired categories
fb<-c("FC", "FF", "SI")
bb<-c("CU", "KC", "CS", "SL", "KN","ST", "SV", "FA")
os<-c("CH", "SC", "EP", "FS", "FO")

# Perform some data transformation by creating new indicators
WINNING<-rep(0, nrow(data)) # Indicator of batting team winning
COUNT<-rep(0, nrow(data))   # Indicator of hitter's count
RISP<-rep(0, nrow(data))    # Indicator if RISP
SHIFT<-rep(0, nrow(data))   # Indicator of a shift
FB<-rep(0, nrow(data))      # Indicator of a fastball
BB<-rep(0, nrow(data))      # Indicator of a breaking ball
OS<-rep(0, nrow(data))      # Indicator of an offspeed


# Find the value of the indicators
for(i in 1:nrow(data)){
  if(data[i,1] %in% fb){
    FB[i]<-1
  }
  if(data[i,1] %in% bb){
    BB[i]<-1
  }
  if(data[i,1] %in% os){
    OS[i]<-1
  }
  if(data[i,18]>data[i,19] & data[i,19]<2){
    COUNT[i]<-1
  }
  if(is.na(data[i,21])==FALSE | is.na(data[i,22])==FALSE){
    RISP[i]<-1
  }
  if(data[i,23]!="Standard" | data[i,24]!="Standard"){
    SHIFT[i]<-1
  }
  if(data[i,49]>data[i,50]){
    WINNING[i]<-1
  }
}


# Reformat the dataframe
data<-data[,c(1,4,6,7,9)] # Select the necessary columns from the data frame
data<-cbind(FB, BB, OS, data, COUNT, RISP, SHIFT, WINNING) # Append new columns to data frame
data<-data[data$PITCH_TYPE!="",-4] # Remove the PITCH_TYPE column as we now have the three categories we wish to predict
noms<-colnames(data) # Save the column names to reformat the dataframe later

# Find the total count values for each variable
batters<-unique(data$BATTER_ID) # Find how many distinct batters are in the data set
years<-unique(data$GAME_YEAR) # Find the distinct years in the data set
df<-matrix(rep(0, 11), nrow=1) # Make a matrix to store the aggregated data
# Loop through each year for each hitter to populated the aggregated data matrix
for(batter in batters){
  for(year in years){
    # Store the aggregated data for each hitter within each year in a temporary vector
    temp<-c(sum(data[data$BATTER_ID==batter & data$GAME_YEAR==year,1]), sum(data[data$BATTER_ID==batter & data$GAME_YEAR==year,2]),
            sum(data[data$BATTER_ID==batter & data$GAME_YEAR==year,3]), batter, 
            sum(data[data$BATTER_ID==batter & data$GAME_YEAR==year,]$BAT_SIDE=="R"), 
            sum(data[data$BATTER_ID==batter & data$GAME_YEAR==year,]$THROW_SIDE=="R"), year,
            sum(data[data$BATTER_ID==batter & data$GAME_YEAR==year,]$COUNT), sum(data[data$BATTER_ID==batter & data$GAME_YEAR==year,]$RISP),
            sum(data[data$BATTER_ID==batter & data$GAME_YEAR==year,]$SHIFT), sum(data[data$BATTER_ID==batter & data$GAME_YEAR==year,]$WINNING))
    # Append the vector to the aggregated data matrix
    df<-rbind(df, temp)
  }
}

# Convert the aggregated data matrix to a data frame
df<-data.frame(df[-1,]) # Drop the empty first row and format as a data frame
colnames(df)<-noms # Attach names to the columns
df[,4]<-as.factor(df[,4]) # Make BATTER_ID a factor


# Make design matrix with the newly formatted data
X<-model.matrix(FB~.-1, df[,c(1,4)]) # Makes the design matrix for the BATTER_ID factor
X<-as.matrix(data.frame(X,as.matrix(df[,-c(1:4)]))) # Appends the corresponding covariates to the design matrix
y<-matrix(c(as.vector(df$FB), as.vector(df$BB), as.vector(df$OS)), ncol=3) # Make a matrix of responses for the three pitch types
# Fit the multinomial model to the data. The value of penalty.factor drops the penalties from the covariates
mod<-glmnet(X, y, family="multinomial", penalty.factor = c(.001, rep(0,ncol(X)-1)), standardize=FALSE, intercept=TRUE)

# Define a function to find the median so we can predict future season situations
med <- function(x) {
  quantile(x,0.5)
}


# Make a covariate matrix for the 2024 season based on historical trends
df_pred<-matrix(rep(0, 8), nrow=1) # Make a matrix to store the predicted aggregated data
for(i in 1:length(batters)){
  # Store the predicted aggregated data in a temporary vector
  temp<-c(batters[i], med(df[df$BATTER_ID==batters[i],]$BAT_SIDE), med(df[df$BATTER_ID==batters[i],]$THROW_SIDE),
          2024, med(df[df$BATTER_ID==batters[i],]$COUNT), med(df[df$BATTER_ID==batters[i],]$RISP),
          med(df[df$BATTER_ID==batters[i],]$SHIFT), med(df[df$BATTER_ID==batters[i],]$WINNING))
  df_pred<-rbind(df_pred, temp) # Append the predicted aggregated data to the 2024 covariate data matrix
}
# Format the prediction data matrix as a dataframe
df_pred<-data.frame(df_pred[-1,]) # Drop the empty first row
df_pred<-cbind(sample(c(1,0),nrow(df_pred), replace=TRUE),df_pred) # Add a useless column to the df_pred dataframe to make the design matrix later
colnames(df_pred)<-noms[-(1:2)] # Attach names to the df_pred matrix to so the design matrix is easier to make
df_pred$BATTER_ID<-as.factor(df_pred$BATTER_ID) # Make BATTER_ID a factor

# Make design matrix for the prediction dataset
X_pred<-model.matrix(OS~.-1, df_pred[,c(1,2)]) # Make the design matrix for the BATTER_ID factor
X_pred<-as.matrix(data.frame(X_pred,as.matrix(df_pred[,-c(1,2)]))) # Append corresponding covariates to the design matrix

# Predict the values for the 2024 season
preds<-predict(mod, X_pred, s=mod$lambda[length(mod$lambda)], type="response")
colnames(preds)<-c("p_FB", "p_BB", "p_OS") # Add column names for ease of use later

# Some data processing to write to csv
temp<-merge(preds, df_pred, by=0) # Merge predictions with the old data frame on BATTER_ID
mat<-matrix(c(0,0,0,0), nrow=1) # Make a predicted data matrix
batters<-unique(temp$BATTER_ID)
# Match predicted probabilities to BATTER_ID
for(batter in batters){
  p_fb<-temp[temp$BATTER_ID==batter,]$p_FB 
  p_bb<-temp[temp$BATTER_ID==batter,]$p_BB
  p_os<-temp[temp$BATTER_ID==batter,]$p_OS
  mat<-rbind(mat,c(as.numeric(batter), p_fb, p_bb, p_os))
}
mat<-data.frame(mat[-1,])
colnames(mat)<-c("BATTER_ID", colnames(preds))

# Match the predicted probabilities and BATTER_ID to the player names
predictions<-read.csv("predictions.csv", header=TRUE)
for(batter in batters){
  predictions[predictions$BATTER_ID==batter,4]<-round(mat[mat$BATTER_ID==batter,2],3)
  predictions[predictions$BATTER_ID==batter,5]<-round(mat[mat$BATTER_ID==batter,3],3)
  predictions[predictions$BATTER_ID==batter,6]<-round(mat[mat$BATTER_ID==batter,4],3)
}

# Write the output to a csv
write.csv(predictions, "predictions.csv", row.names=FALSE)



