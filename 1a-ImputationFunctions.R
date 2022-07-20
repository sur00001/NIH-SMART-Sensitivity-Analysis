################################################################################
# This script defines the functions for the imputation and sensitivity analysis 
# A DEMO IS AT THE BOTTOM OF THE SCRIPT
################################################################################
#library(matrixcalc)
library(mvtnorm)
library(mice)
library(broom.mixed)



#-------------------------------------------------------------------------------
# Function to impute data 

# Inputs: long.dat and wide.dat (data sets with missingness), 
#          M (# of imputations; default is 10),
#          SA (if SA=0, sensitivity parameter kij=0, otw. kij is specified;
#              default is no sensitivity analysis ie. MAR imputation),
#          kij (the sensitivity delta vector)

# Output: data frame with all M imputed datasets
#-------------------------------------------------------------------------------


impute.data = function(long.dat, wide.dat, M=10, SA=0, kij=NULL){

#-------------------------------------------------
# Add missingness Yij indicator to long data (Rij)
#-------------------------------------------------
#Rij = 1 if Yij is missing 

long.dat$rij = is.na(long.dat$y)
long.dat$rij = ifelse(long.dat$rij==T,1,0) #if missing, Rij = 1

#-------------------------------------------------
# Initialize wide and long datasets that stores the M imputed datasets
#-------------------------------------------------
lmm.imp = data.frame()
lm.imp=data.frame()

#-------------------------------------------------
# Fit model 
#-------------------------------------------------
fit=lmer(y~0+ba+a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15+a16+a17+a18+a19+a20+a21+a22+a23+a24+
           b113+b114+b115+b116+b117+b118+b119+b120+b121+b122+b123+b124+
           b013+b014+b015+b016+b017+b018+b019+b020+b021+b022+b023+b024+
           nr113+nr114+nr115+nr116+nr117+nr118+nr119+nr120+nr121+nr122+nr123+nr124+
           nr013+nr014+nr015+nr016+nr017+nr018+nr019+nr020+nr021+nr022+nr023+nr024+
           runin+first_half+time_1st_half+time2_1st_half+time3_1st_half+second_half+time_2nd_half+time2_2nd_half+time3_2nd_half+
           (runin+ time_1st_half+ time_2nd_half | id), control=lmerControl(check.nobs.vs.nRE="ignore"), data=long.dat)
# This automatically fits the model with only the observed observations

# For each imputation m in {1,....,M}...
for (m in 1:M){
set.seed(m)
#-------------------------------------------------------------------------------
# STEP A: Draw fixed effects and random effects
#-------------------------------------------------------------------------------

#-------------------
# Draw fixed effects 
#-------------------
Beta = fixef(fit)
VarBeta = matrix(vcov(fit),nrow=82)
beta.star = rmvnorm(1,mean=Beta,sigma=VarBeta)
#is.positive.semi.definite(VarBeta) #Do I need the Chomolesky decomposition? Maybe not, because covariance matrix is always positive semi-definite?
beta.star = t(beta.star) #p by 1

#-------------------
#Draw random effects 
#-------------------
#NOTE: I hard coded it for 4 random effects, will change it later 

b = ranef(fit) #n x 4; For each individual, there are 4 random effects 
b = b[[1]] #get data frame, each row is a person

#Variance covariance matrix of random effects
RE.df = as.data.frame(VarCorr(fit))
REVC = RE.df$vcov
VarB = matrix(rep(0,16),nrow=4,ncol=4) #Initialize variance-covariance matrix for random effects

diag(VarB) = c(REVC[1],REVC[2],REVC[3],REVC[4]) #first 4 elements are variances
VarB[1,2] = VarB[2,1] = REVC[5] #covariance of intercept and run-in
VarB[1,3] = VarB[3,1] = REVC[6] #covariance of intercept and time 1st half
VarB[1,4] = VarB[4,1] = REVC[7] #covariance of intercept and time 2nd half
VarB[2,3] = VarB[3,2] = REVC[8] #covariance of run-in and time 1st half
VarB[2,4] = VarB[4,2] = REVC[9] #covariance of run-in and time 2nd half
VarB[3,4] = VarB[4,3] = REVC[10]# covariance of time 1st half and time 2nd half
colnames(VarB) = rownames(VarB) = c("Intercept","run-in","time 1st half","time 2nd half")

#Draw random effects for all individuals
b.star = data.frame(rep(0,nperson),rep(0,nperson),rep(0,nperson),rep(0,nperson)) #column for each of the 4 random effects

for (i in 1:nperson){
  bi = as.numeric(b[i,])
  b.star[i,] = rmvnorm(1,mean=bi,sigma=VarB)
}

colnames(b.star) = c("Intercept b*","run-in b*","time 1st half b*","time 2nd half b*")
b.star = t(b.star) 
dim(b.star) #q by nperson matrix

#-------------------------------------------------------------------------------
# STEP B: Calculate nij for each missing observation
#-------------------------------------------------------------------------------

#---------------------------
# Set sensitivity parameter kij
# (N by 1 where N=25*nperson)
#---------------------------

#kij = 0 if user does not want a sensitivity analysis (MAR imputation)
if (SA==0){ 
  kij = matrix(rep(0,nperson*25),nrow=nperson*25) #(nperson*25 by 1)
}

# Need to finish
if (SA==1){
  kij = as.matrix(kij) #a vector kij is specified as a parameter in this function
}

#---------------------------
# Get Design matrix X 
# (N by p where N=25*nperson)
#---------------------------

X = long.dat[,names(long.dat) %in% names(fixef(fit))] #design matrix - matrix of covariates
X = as.matrix(X) #design matrix of all observations, (25*n) by p
dim(X) #3250 (25*nperson) by 82 (p)

#---------------------------
# Get Z matrix
# (N by q where q = # of random effects)
#---------------------------

#First column of Z is 1s (for the random intercept)
Z = data.frame(Intercept=rep(1,nperson*25),runin=rep(0,nperson*25),
               time_1st_half=rep(0,nperson*25),time_2nd_half=rep(0,nperson*25))

Z[,2:4] = long.dat[,names(long.dat) %in% c("runin","time_1st_half","time_2nd_half")]

Z = as.matrix(Z) #random effect design matrix of all observations, (25*n) by q
dim(Z) #3250 (25*nperson) by 4 (q)

#---------------------------
# Calculate eta_ij for missing obs
#---------------------------

misIndv = unique(long.dat$id[long.dat$rij==1]) #individuals with at least one missing observation
eta_ijs = c() #initialize the vector of eta_ij's for missing observations

#For each person with at least one missing obs...
for (i in c(misIndv)){ 
  
#Denote ri = # of missing observations for person i, p = # of covariates, q = # of random effects
  
eta_ij = X[long.dat$rij==1 & long.dat$id==i,] %*% beta.star # (ri by  p)  *   (p by 1)
        + Z[long.dat$rij==1 & long.dat$id==i,] %*% b.star[,i] # (ri by q) * (q by 1)
        + kij[long.dat$rij==1 & long.dat$id==i,] # (ri by 1)

# (ri by 1) = (ri by p)*(p by 1) + (ri by q)*(q by 1) + (ri by 1)

eta_ijs = c(eta_ijs, eta_ij) #add to eta vector
#print(i)
}

#Check that the length of eta is the number of missing observations
length(eta_ijs) == sum(long.dat$rij==1)

#-------------------------------------------------------------------------------
# STEP C: Draw error from residual variance 
#-------------------------------------------------------------------------------

# Ask James questions about calculating q (a bit confused what the A matrix is in the Bates paper)
# Also is it neccessary to draw sigma^2 versus just using the residual variance?

#For now, draw error using residual variance
RE.df = as.data.frame(VarCorr(fit)) #Variance covariance matrix of random effects
sigma = sqrt(RE.df$vcov[11]) #last element is the residual variance

e_ijs = rnorm(sum(long.dat$rij==1),0,sigma) # (total missing obs by 1) vector

#-------------------------------------------------------------------------------
# STEP D: Impute each missing Yij = eta_ij + e_ij and create imputed dataset
#-------------------------------------------------------------------------------

#NOTE: I did a simple case, where I did not re-calculate the stage 1 avg and check
#non-responder status 

yijs = eta_ijs + e_ijs

#imputed dataset in long format (all yijs are imputed)
imp.long = long.dat 
imp.long$y[imp.long$rij==1] = yijs 
imp.long$.imp = rep(m,length(imp.long$y))

#impute dataset in wide format (y_avg is imputed)
imp.wide = wide.dat 
imp.wide$y_avg = (imp.long$y[imp.long$ftime==21]+imp.long$y[imp.long$ftime==22] 
                 + imp.long$y[imp.long$ftime==23] + imp.long$y[imp.long$ftime==24])/4
imp.wide$.imp = rep(m, length(imp.wide$y_avg))

lmm.imp = rbind(lmm.imp, imp.long)
lm.imp = rbind(lm.imp, imp.wide)
}

return(list(lmm.imp, lm.imp)) #long and wide dataset with all M imputed datasets
}


#------------------------------------------------------------------------------------------------------------------
# Fit a LM or LMM model with imputed data 

# Inputs: m (imputation number), mod.type (LMM or LM), allimp (dataset of all M imputed datasets)

# Output: fitted LM or LMM model for ONE imputed dataset
#------------------------------------------------------------------------------------------------------------------

fit.imputed = function(m,mod.type,allimp) {
  
  # subset to one imputed dataset
  mydata = allimp[allimp$.imp==m,]
  
  # fit LM or LMM model 
  if (mod.type=="LM"){
    mod = lm(y_avg~ba+y0+a+b1+b0+nr1+nr0, data=mydata)
  }
  else {
    mod = lmer(y~0+ba+runin+first_half+time_1st_half+time2_1st_half+time3_1st_half+second_half+time_2nd_half+time2_2nd_half+time3_2nd_half+
                 a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15+a16+a17+a18+a19+a20+a21+a22+a23+a24+
                 b113+b114+b115+b116+b117+b118+b119+b120+b121+b122+b123+b124+
                 b013+b014+b015+b016+b017+b018+b019+b020+b021+b022+b023+b024+
                 nr113+nr114+nr115+nr116+nr117+nr118+nr119+nr120+nr121+nr122+nr123+nr124+
                 nr013+nr014+nr015+nr016+nr017+nr018+nr019+nr020+nr021+nr022+nr023+nr024+
                 (runin+ time_1st_half+ time_2nd_half | id), control=lmerControl(check.nobs.vs.nRE="ignore"), data=mydata)
  }
  
  return(mod)
}

################################################################################
# DEMO
################################################################################

M=4

#Open Toy long.dat and wide.dat
load("C:/Users/surtc/OneDrive/Documents/NIH Internship/toy_datasets_w_missingness.rda")

#Impute data
imps = impute.data(long.dat,wide.dat,M=M,SA=0) #no sensitivity analysis
lmm.imps = imps[[1]]
lm.imps = imps[[2]]

#Run a LMM model and get pooled estimate using Rubin's rules

#Obtain list of model fits for each imputed dataset
lmm.fits = lapply(1:M,fit.imputed,mod.type="LMM",allimp=lmm.imp) 
pooled.fits = pool(lmm.fits)
summary(pooled.fits)

A_trtcoefLMM = (fit$estimate[fit$term=="a21"] + fit$estimate[fit$term=="a21"] +
                  fit$estimate[fit$term=="a22"] + fit$estimate[fit$term=="a24"])/4

 

