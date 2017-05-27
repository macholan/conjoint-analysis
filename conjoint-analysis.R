############################################################
# 1. Load the necessary packages for the Conjoint analysis #
############################################################

library(bayesm)
library(dummies)
library(pROC)

##################################
# 2. Load the survey data into R #
##################################

load("stc-cbc-respondents-v3.RData")
resp <- resp.data.v3
tasks <- read.csv("stc-dc-task-cbc-v3.csv")
extras3 <- read.csv("extra-scenarios-v3.csv")
extras2 <- read.csv("extra-scenarios-v3-2.csv")
extras1 <- read.csv("extra-scenarios.csv")
load("efCode.RData")

### Table of actual choices
table.actual <- apply(resp[4:39], 2, function(x){tabulate(na.omit(x))})
table.actual
write.csv(table.actual,"table_actual.csv")

#######################################
# 3. Effects coding for tasks matrix  #
#######################################

tasks.matrix <- as.matrix(tasks[,c("screen","RAM","processor","price","brand")])
X.mat <- efcode.attmat.f(tasks.matrix)
pricevec <- tasks$price - mean(tasks$price) # create new price variable that is centered on its mean
X.brands <- X.mat[ ,9:11] # get the columns from X.mat that represent brand
X.BrandByPrice <- X.brands*pricevec # multiply each column in X.brands by pricevec
X.matrix <- cbind(X.mat, X.BrandByPrice) # combine the two matrices

extras2.matrix <- as.matrix(extras2[,c("screen","RAM","processor","price","brand")])
X.mat2 <- efcode.attmat.f(extras2.matrix)
pricevec2 <- extras2$price - 1 # create new price variable that is centered on its mean
X.brands2 <- X.mat2[ ,9:11] # get the columns from X.mat that represent brand
X.BrandByPrice2 <- X.brands2*pricevec2 # multiply each column in X.brands by pricevec
X.matrix2 <- cbind(X.mat2, X.BrandByPrice2) # combine the two matrices

##################################
# 4. Create responses data frame #
##################################

ydata <- resp[,4:39]
names(ydata)
ydata <- na.omit(ydata)
ydata <- as.matrix(ydata)
zowner <- 1 * (!is.na(resp$vList3)) # whether the respondent owned STC

##################################
# 5. Prepare data for rhierMnlDP #
##################################

lgtdata <- NULL # a starter placeholder for the list 
for (i in 1:424) {lgtdata[[i]] <- list(y = ydata[i, ], X = X.matrix)}
length(lgtdata)
str(lgtdata)
lgtdata[[3]]

##########################################
# 5. Prepare data for rhierMnlDP and run #
##########################################

### Set parameters
resp.count <- 1:424
iterations <- 100000
post.burndown <- 5001:20000
### Set parameters

lgtdata.model <- lgtdata[resp.count]
mcmctest <- list(R=iterations, keep=5) # specify every 5th sample is kept
Data1 <- list(p=3,lgtdata=lgtdata.model) # create the "Data" list rhierMnlDP() expects; p is choice set size
testrun1 <- rhierMnlDP(Data=Data1,Mcmc=mcmctest)
names(testrun1)
betadraw1 <- testrun1$betadraw  # array that has the draws (i.e. samples from marginal posterior distributions) for the regression coefficients
dim(betadraw1)

for (j in 1:6) {
    jpeg(paste0("Beta",j,".jpeg"), width = 900, height = 600, quality = 100)
    par(mfrow = c(3,2))
    for (i in 1:6) {
        plot(1:length(betadraw1[i,j,]),betadraw1[i,j,], xlab = "Index", ylab = "Beta Est.", main = paste("Beta",j,"Respondent",i)) # "Index" is the iteration number after "thinning."
        abline(v=5000, col = "red") ## vertical line
        abline(v=10000, col = "red") ## vertical line?jpb
    }
    dev.off()
}

plot(density(betadraw1[1,1,post.burndown],width=2)) # look at distribution
summary(betadraw1[1,1,post.burndown])
mn <- mean(betadraw1[1,1,post.burndown])
sd <- sd(betadraw1[1,1,post.burndown])
mn
sd
abline(v=0, col="red") ## vertical line
abline(v=mn, col="red") ## vertical line
prob <- pnorm(0,mean=mn, sd=sd, lower.tail = FALSE)
prob # compute the probability that person 1 beta 1 > 0

for (j in 1:14) {
    jpeg(paste0("BetaDistr",j,".jpeg"), width = 900, height = 400, quality = 100)
    par(mfrow = c(2,3))
    for (i in 1:6) {
        plot(density(betadraw1[i,j,post.burndown],width=2), xlim = c(-4,4), main = paste("Beta",j,"Respondent",i)) # look at distribution
        mn <- mean(betadraw1[i,j,post.burndown])
        abline(v=mn, col="red") ## vertical line
    }
    dev.off()
}

betameansoverall <- apply(betadraw1[,,post.burndown],c(2),mean) # overall means of the coefficients across respondents
betameansoverall
betameansindividual <- apply(betadraw1[,,post.burndown],c(1,2),mean) # matrix of coefficient means by respondent
betameansindividual

screen1 <- apply((1-(betadraw1[1,1,post.burndown]+betadraw1[1,2,post.burndown])),2,quantile,probs=c(0.05,0.10,0.25,0.5 ,0.75,0.90,0.95)) # calculate a distribution for the difference between respondent 1's 1st and 2nd coefficients (these are for the screen attribute)
ram1 <- apply((1-(betadraw1[1,3,post.burndown]+betadraw1[1,4,post.burndown])),2,quantile,probs=c(0.05,0.10,0.25,0.5 ,0.75,0.90,0.95)) # calculate a distribution for the difference between respondent 1's 1st and 2nd coefficients (these are for the screen attribute)
processor1 <- apply((1-(betadraw1[1,5,post.burndown]+betadraw1[1,6,post.burndown])),2,quantile,probs=c(0.05,0.10,0.25,0.5 ,0.75,0.90,0.95)) # calculate a distribution for the difference between respondent 1's 1st and 2nd coefficients (these are for the screen attribute)
price1 <- apply((betadraw1[1,1,post.burndown]-betadraw1[1,2,post.burndown])) # calculate a distribution for the difference between respondent 1's 1st and 2nd coefficients (these are for the screen attribute)
brand1 <- apply((betadraw1[1,1,post.burndown]-betadraw1[1,2,post.burndown])) # calculate a distribution for the difference between respondent 1's 1st and 2nd coefficients (these are for the screen attribute)
brandprice1 <- apply((betadraw1[1,1,post.burndown]-betadraw1[1,2,post.burndown])) # calculate a distribution for the difference between respondent 1's 1st and 2nd coefficients (these are for the screen attribute)

percent <- apply(betadraw1[,,post.burndown],2,quantile,probs=c(0.05,0.10,0.25,0.5 ,0.75,0.90,0.95))
write.csv(percent,"table1.csv")
percent

plot(density(betadraw1[1,1,post.burndown]-betadraw1[1,2,post.burndown],width=2))


### compute the empirical probability based on sample values
p1b1 <- betadraw1[1,1,post.burndown]
p1b1df <- as.data.frame(p1b1)
str(p1b1df)
tab <- table(sign(p1b1df$p1b1))
str(tab)
tabdf <- as.data.frame(tab)
tabdf$Pct <- tabdf$Freq/sum(tabdf$Freq)
tabdf$Pct

table(testrun1$Istardraw) # 3 is the mode

#################################################
# 6. Create a model with the ownership variable #
#################################################

zownertest <- matrix(scale(zowner[resp.count],scale=FALSE),ncol=1)
Data2 <- list(p=3, lgtdata = lgtdata.model, Z = zownertest)
testrun2 <- rhierMnlDP(Data = Data2, Mcmc = mcmctest)
names(testrun2)
dim(testrun2$Deltadraw)
apply(testrun2$Deltadraw[post.burndown,],2,mean) 
deltadrawdistr <- apply(testrun2$Deltadraw[post.burndown,],2,quantile,probs=c(0.95,0.90,0.75,0.5 ,0.25,0.10,0.05))
deltadrawdistr
write.csv(deltadrawdistr, "deltadrawdistr.csv")
betadraw2 <- testrun2$betadraw
dim(betadraw2)

betameansoverall2 <- apply(betadraw2[,,post.burndown],c(2),mean) # overall means of the coefficients across respondents
betameansoverall2
betameansindividual2 <- apply(betadraw2[,,post.burndown],c(1,2),mean) # matrix of coefficient means by respondent
betameansindividual2
summary((betadraw2[1,1,post.burndown]-betadraw2[1,2,post.burndown])) # calculate a distribution for the difference between respondent 1's 1st and 2nd coefficients (these are for the screen attribute)

percent <- apply(betadraw2[,,post.burndown],2,quantile,probs=c(0.05,0.10,0.25,0.5 ,0.75,0.90,0.95))
write.csv(percent,"table2.csv")
percent


###################################################################################################
# 7a. Calculate predicted choice probabilities for Model 1 - No Covariate (Individual Beta Means) #
###################################################################################################

dim(betameansindividual) # matrix with subjects in the rows and means in the columns
dim(t(betameansindividual))
dim(X.matrix)
xbeta <- X.matrix%*%t(betameansindividual)
dim(xbeta)
xbeta2 <- matrix(xbeta,ncol=3,byrow=TRUE) # reading values column by column and putting them into rows of new matrix with 3 columns
dim(xbeta2) 
expxbeta2 <- exp(xbeta2) # exponentiate
dim(expxbeta2)

### get predicted choice probabilities
rsumvec <- rowSums(expxbeta2)
pchoicemat <- expxbeta2/rsumvec # predicted choice probabilities
head(pchoicemat)
dim(pchoicemat)
custchoice <- max.col(pchoicemat) # prediction of customer choice
head(custchoice)
str(custchoice)

### assess model fit
ydatavec <- as.vector(t(ydata))
str(ydatavec)
table(custchoice,ydatavec[1:15264]) # provides confusion matrix
roctest <- roc(ydatavec[1:15264], custchoice, plot=TRUE) # ROC curve
auc(roctest) # Area Under the Curve
logliketest <- testrun1$loglike # -2log(likelihood) test applies to nested models only
str(logliketest)
mean(logliketest)
par(mfrow = c(1,1))
hist(logliketest, main = "Model 1 - Log Likelihood Histogram")

# predict the choices for the 36 choice sets
m <- matrix(custchoice, nrow =36, byrow=F)
m2 <- t(m)
predicted1 <- apply(m2, 2, function(x){tabulate(na.omit(x))}) 
write.csv(predicted1, "predicted1.csv")

################################################################################################
# 7b. Calculate predicted choice probabilities for Model 1 - No Covariate (Overall Beta Means) #
################################################################################################

betavec=matrix(betameansoverall,ncol=1,byrow=TRUE)
xbeta=X.matrix%*%(betavec)
dim(xbeta)
xbeta2=matrix(xbeta,ncol=3,byrow=TRUE)
dim(xbeta2)
expxbeta2=exp(xbeta2)
rsumvec=rowSums(expxbeta2)
pchoicemat=expxbeta2/rsumvec
pchoicemat
### Predicted frequencies based on overall model
pchoicematoverall1 <- round(pchoicemat*424,digits=0)
pchoicematoverall1
write.csv(pchoicematoverall1, "predictoverall1.csv")

custchoice <- max.col(pchoicematoverall1) # prediction of customer choice
head(custchoice)
str(custchoice)

### assess model fit
ydatavec <- as.vector(t(ydata))
str(ydatavec)
table(c(rep(custchoice,424)),ydatavec[1:15264]) # provides confusion matrix
roctest <- roc(ydatavec[1:15264], c(rep(custchoice,424)), plot=TRUE) # ROC curve
auc(roctest) # Area Under the Curve


#################################################################################################
# 8a. Calculate predicted choice probabilities for Model 2  - Covariate (Individual Beta Means) #
#################################################################################################

dim(betameansindividual2) # matrix with subjects in the rows and means in the columns
dim(t(betameansindividual2))
dim(X.matrix)
xbeta <- X.matrix%*%t(betameansindividual2)
dim(xbeta)
xbeta2 <- matrix(xbeta,ncol=3,byrow=TRUE) # reading values column by column and putting them into rows of new matrix with 3 columns
dim(xbeta2) 
expxbeta2 <- exp(xbeta2) # exponentiate
dim(expxbeta2)

### get predicted choice probabilities
rsumvec <- rowSums(expxbeta2)
pchoicemat <- expxbeta2/rsumvec # predicted choice probabilities
head(pchoicemat)
dim(pchoicemat)
custchoice <- max.col(pchoicemat) # prediction of customer choice
head(custchoice)
str(custchoice)

### assess model fit
ydatavec <- as.vector(t(ydata))
str(ydatavec)
table(custchoice,ydatavec[1:15264]) # provides confusion matrix
roctest <- roc(ydatavec[1:15264], custchoice, plot=TRUE) # ROC curve
auc(roctest) # Area Under the Curve
logliketest <- testrun2$loglike # -2log(likelihood) test applies to nested models only
str(logliketest)
mean(logliketest)
par(mfrow = c(1,1))
hist(logliketest, main = "Model 2 - Log Likelihood Histogram")

# predict the choices for the 36 choice sets
m <- matrix(custchoice, nrow =36, byrow=F)
m2 <- t(m)
predicted2 <- apply(m2, 2, function(x){tabulate(na.omit(x))}) 
write.csv(predicted2, "predicted2.csv") 


#############################################################################################
# 8b. Calculate predicted choice probabilities for Model 2 - Covariate (Overall Beta Means) #
#############################################################################################

betavec=matrix(betameansoverall2,ncol=1,byrow=TRUE)
xbeta=X.matrix%*%(betavec)
dim(xbeta)
xbeta2=matrix(xbeta,ncol=3,byrow=TRUE)
dim(xbeta2)
expxbeta2=exp(xbeta2)
rsumvec=rowSums(expxbeta2)
pchoicemat=expxbeta2/rsumvec
pchoicemat
### Predicted frequencies based on overall model
pchoicematoverall2 <- round(pchoicemat*424,digits=0)
pchoicematoverall2
write.csv(pchoicematoverall2, "predictoverall1.csv")

custchoice <- max.col(pchoicematoverall2) # prediction of customer choice
head(custchoice)
str(custchoice)

### assess model fit
ydatavec <- as.vector(t(ydata))
str(ydatavec)
table(c(rep(custchoice,424)),ydatavec[1:15264]) # provides confusion matrix
roctest <- roc(ydatavec[1:15264], c(rep(custchoice,424)), plot=TRUE) # ROC curve
auc(roctest) # Area Under the Curve


##################################################################
# 9a. Predict extra scenarios using betas from Model 1 - Overall #
##################################################################

Xextra.matrix <- as.matrix(extras1[,c("V1","V2","V3","V4","V5","V6","V7","V8","V9",
                                      "V10","V11","V12","V13","V14")])
betavec <- matrix(betameansoverall,ncol=1,byrow=TRUE)
xextrabeta <- Xextra.matrix%*%(betavec)
xbetaextra2 <- matrix(xextrabeta,ncol=3,byrow=TRUE)
betavec
dim(betavec)
dim(Xextra.matrix)
dim(xextrabeta)
dim(xbetaextra2)

expxbetaextra2 <- exp(xbetaextra2)
rsumvec <- rowSums(expxbetaextra2)
pchoicemat1 <- expxbetaextra2/rsumvec
pchoicemat1


#####################################################################
# 9b. Predict extra scenarios using betas from Model 1 - Individual #
#####################################################################

xextrabetaind=Xextra.matrix%*%(t(betameansindividual))
dim(xextrabetaind)
xbetaextra2ind=rbind(matrix(xextrabetaind[1:3,],ncol=3,byrow=TRUE),matrix(xextrabetaind[4:6,],ncol=3,byrow=TRUE))
dim(xbetaextra2ind)

expxbetaextra2ind=exp(xbetaextra2ind)
rsumvecind=rowSums(expxbetaextra2ind)
pchoicematind=expxbetaextra2ind/rsumvecind
pchoicematind

custchoiceind <- max.col(pchoicematind)
head(custchoiceind)
str(custchoiceind)
extra1 <- custchoiceind[1:424]
extra2 <- custchoiceind[425:848]
table(extra1)
table(extra2)

#####################################################################
# 9b. Predict extra scenarios using betas from Model 1 - Individual #
#####################################################################

xextrabetaind2=X.matrix2%*%(t(betameansindividual))
dim(xextrabetaind2)
xbetaextra2ind2=rbind(matrix(xextrabetaind2[1:3,],ncol=3,byrow=TRUE),matrix(xextrabetaind2[4:6,],ncol=3,byrow=TRUE),
                         matrix(xextrabetaind2[7:9,],ncol=3,byrow=TRUE),matrix(xextrabetaind2[10:12,],ncol=3,byrow=TRUE),
                         matrix(xextrabetaind2[13:15,],ncol=3,byrow=TRUE),matrix(xextrabetaind2[16:18,],ncol=3,byrow=TRUE))

dim(xbetaextra2ind2)

expxbetaextra2ind2=exp(xbetaextra2ind2)
rsumvecind2=rowSums(expxbetaextra2ind2)
pchoicematind1=expxbetaextra2ind2/rsumvecind2
pchoicematind1

custchoiceind2 <- max.col(pchoicematind1)
head(custchoiceind2)
str(custchoiceind2)

extra1 <- custchoiceind2[1:424]
extra2 <- custchoiceind2[425:848]
extra3 <- custchoiceind2[849:1272]
extra4 <- custchoiceind2[1273:1696]
extra5 <- custchoiceind2[1697:2120]
extra6 <- custchoiceind2[2121:2544]

table(extra1)
table(extra2)
table(extra3)
table(extra4)
table(extra5)
table(extra6)

##################################################################
# 9a. Predict extra scenarios using betas from Model 1 - Overall #
##################################################################

Xextra.matrix <- as.matrix(extras1[,c("V1","V2","V3","V4","V5","V6","V7","V8","V9",
                                      "V10","V11","V12","V13","V14")])
betavec <- matrix(betameansoverall,ncol=1,byrow=TRUE)
xextrabeta <- Xextra.matrix%*%(betavec)
xbetaextra2 <- matrix(xextrabeta,ncol=3,byrow=TRUE)
betavec
dim(betavec)
dim(Xextra.matrix)
dim(xextrabeta)
dim(xbetaextra2)

expxbetaextra2 <- exp(xbetaextra2)
rsumvec <- rowSums(expxbetaextra2)
pchoicemat1 <- expxbetaextra2/rsumvec
pchoicemat1


######################################################################
# 10b. Predict extra scenarios using betas from Model 2 - Individual #
######################################################################

xextrabetaind=Xextra.matrix%*%(t(betameansindividual2))
dim(xextrabetaind)
xbetaextra2ind=matrix(xextrabetaind,ncol=3,byrow=TRUE)
dim(xbetaextra2ind)

expxbetaextra2ind=exp(xbetaextra2ind)
rsumvecind=rowSums(expxbetaextra2ind)
pchoicematind=expxbetaextra2ind/rsumvecind
pchoicematind

custchoiceind <- max.col(pchoicematind)
head(custchoiceind)
str(custchoiceind)
extra1 <- custchoiceind[1:424]
extra2 <- custchoiceind[425:848]
table(extra1)
table(extra2)


######################################################
## chi square Calculation for each respondent
#######################################################
custchoicedf <- as.data.frame(custchoice)
ydatavecdf <- as.data.frame(ydatavec)
str(ydatavecdf)
str(custchoicedf)
custydatadf <- cbind(custchoicedf,ydatavecdf)
write.csv(custydatadf, file = "custchoiceydata.csv")
idchoice <- read.csv("idchoice.csv")
head(idchoice)
head(pchoicemat)
maxprob <- apply(pchoicemat,1,max)
str(maxprob)
head(maxprob)
maxprobdf <- as.data.frame(maxprob)
pchoicematdf <- as.data.frame(pchoicemat)
choiceydatadf <- cbind(idchoice,pchoicematdf,maxprobdf,ydatavecdf)
head(choiceydatadf)
str(choiceydatadf)
choiceydatadf$ob1 <- NULL
choiceydatadf$ob2 <- NULL
choiceydatadf$ob3 <- NULL
choiceydatadf$ob1 <- 0
choiceydatadf$ob2 <- 0
choiceydatadf$ob3 <- 0
head(choiceydatadf)

choiceydatadf$ob1=ifelse(choiceydatadf$ydatavec == 1,1,0)
choiceydatadf$ob2=ifelse(choiceydatadf$ydatavec == 2,1,0)
choiceydatadf$ob3=ifelse(choiceydatadf$ydatavec == 3,1,0)
head(choiceydatadf)

choiceydatadf$dif1 <- abs(choiceydatadf$V1-choiceydatadf$maxprob)
choiceydatadf$dif2 <- abs(choiceydatadf$V2-choiceydatadf$maxprob)
choiceydatadf$dif3 <- abs(choiceydatadf$V3-choiceydatadf$maxprob)

choiceydatadf$ex1=ifelse(choiceydatadf$dif1 < 0.000001,1,0)
choiceydatadf$ex2=ifelse(choiceydatadf$dif2 < 0.000001,1,0)
choiceydatadf$ex3=ifelse(choiceydatadf$dif3 < 0.000001,1,0)
head(choiceydatadf)
choiceydatadf$chisq <- ((choiceydatadf$ob1-choiceydatadf$ex1)**2) + ((choiceydatadf$ob2-choiceydatadf$ex2)**2) + ((choiceydatadf$ob3-choiceydatadf$ex3)**2)

head(choiceydatadf)
str(choiceydatadf)

write.csv(choiceydatadf, file = "choiceydatadf.csv")

chisquare <- aggregate(chisq~id, choiceydatadf, sum)
dim(chisquare)
head(chisquare)
str(chisquare)
chisquare
hist(chisquare$chisq, main="Histogram of Chisq values", xlab="chisq")
require(ggplot2)
plot(chisquare$chisq, chisquare$id)
chi.sub <- subset(chisquare, chisq < 30.0001)
head(chi.sub)
hist(chi.sub$chisq, main="Histogram of Chisq values < 30", xlab="chisq")
chi.sub2 <- subset(chisquare, chisq > 30)
chi.sub3 <- subset(chisquare, chisq > 60)
chi.sub3
head(chi.sub2)
hist(chi.sub2$chisq, main="Histogram of Chisq values > 30", xlab="chisq")
str(chi.sub2)
chi.sub2[order(chi.sub2$chisq),]


##################################################################
#accuracy based on confusion matrix for each of the 424 respondents
##################################################################
resp_accuracy <- NULL
for (i in 1:424) {
    start <- i*36-35
    end <- i*36
    d <- table(factor(custchoice[start:end],levels = 1:3),
               factor(ydatavec[start:end], levels = 1:3))
    resp_accuracy[i] <- sum(diag(d))/sum(d)
} 
plot(resp_accuracy, main = "Model Accuracy by Respondent")
respdf <- data.frame(resp_accuracy)
head(respdf)
str(respdf)
head(ydatadf)
rn <- rownames(ydatadf)
rndf <- as.data.frame(rn)
resp_all <- cbind(rndf,respdf)
head(resp_all)

str(resp_all)
hist(resp_all$resp_accuracy)
outlier <- subset(resp_all, resp_accuracy < 0.6)
outlier[order(outlier$resp_accuracy),]
#############################################################
#########################################################
#############end#########################################

#idcdf <- read.csv("idchoice.csv")
#rndf <- read.csv("rownamedf.csv")
#total <- merge(rndf,idcdf,by="id")
#str(total)
#head(total)
#write.csv(total,"totalidchoice.csv")