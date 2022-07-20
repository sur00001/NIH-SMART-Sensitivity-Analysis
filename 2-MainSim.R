##### This script is the main simulation code (uses functions from 1-SimFunctions.R) #####
setwd("~/NIH Internship")
source("1-SmartSimFunctions.R")
source("1a-ImputationFunctions.R")

################################################################################
#-------------------------------------------------------------------------------
# Set parameters for simulation 
#-------------------------------------------------------------------------------
################################################################################

#Initialize seed for random number generation
seed=99803332 
seedstart=seed
set.seed(seed)

library(survival)
library(stats)
library(gtools)
library(plyr)
library(mvtnorm)
library(lme4)
library(statmod)
library(lsmeans)
library(multcomp)
library(nlme)
library(mvtnorm)
library(mice)
library(broom.mixed)

#library(tidyverse)

#------------------------------------------------------------
# NUMBER OF REPS AND BOOTSTRAP SAMPLES, AND MODEL INDICATORS  
#------------------------------------------------------------
#Number of reps and bootstraps, model indicators and zstar
nreps = 100 #number of rcts
B = 0 #number of bootstrap samples
z975=qnorm(.975,mean=0,sd=1)
LMM = 1 #indicator variables for LMM and LM
LM = 1

#------------------------------------------------------------
# SIMULATION MODEL PARAMETERS 
#------------------------------------------------------------
alpha=10000.0
beta_age=-50.
beta_ri=-4000. #run-in 
beta_a=00.
beta_1b=00.
beta_0b=00.
# correlation=exp(-|time diff in weeks|/rho)
rho=52
gam=2000 #error sd
#beta_1nr
#beta_0nr
# beta_nr121 = -2865.421
# beta_nr122 = -2905.38
# beta_nr123 = -2816.242
# beta_nr124 = -2853.776
# beta_nr021 = -2766.204
# beta_nr022 = -2794.646
# beta_nr023 = -2717.026
# beta_nr024 = -2747.798
#p_1 = .2222
#p_0 = .2241

#------------------------------------------------------------
# SIMULATION PARAMETERS RELATED TO MODEL
#------------------------------------------------------------
minage=20
maxagedel=20

#------------------------------------------------------------
# BASIC SIMULATION PARAMETERS 
#------------------------------------------------------------
# NPERSON = # of total subjects in trial 
# FRAC_TREAT is the fraction randomized to treatment arm (A=1)
#------------------------------------------------------------
nperson=130
frac_treat=.5
#------------------------------------------------------------
#------------------------------------------------------------
# MISSINGNESS PARAMETERS 
#  miss_type = NONE    Full data analysis
#  miss_type = MCAR    Missing Completely at Random
#  miss_type = MAR     Missing at Random
#  miss_type = MNAR    Missing Not at Random
#  miss_prob           Probability of Missing at each week 
#
# MISSINGNESS LOGISTIC MODELS FOR MAR and MNAR
#  alpha_m     = intercept
#  beta_m0     = effect of y0-6000
#  beta_ma     = effect of a-frac_treat
#  beta_mt     = effect of time-12
#  beta_mnr    = effect of non-response
#  beta_mr     = effect of (Y21+Y22+Y23+Y24)/4-10000
#  beta_yt     = effect of Yt (current time point)

# IMPUTATION PARAMETERS
# M            = number of imputed datasets
# impute       = indicates whether you want results after MI 
# SA           = indicators whether you want to conduct a sensitivity analysis
#------------------------------------------------------------
miss_type="MCAR"
alpha_m=-2.2
beta_m0=-.5/1000
beta_ma=1
beta_mt=1/12
beta_mnr=.75
beta_mr=-.25/1000
beta_yt =-.25/1000
miss_prob=.1 #only relevant for MCAR
M = 10
impute = 0
SA = 0

# Results dataframe will contain the rep #, model ("LMM" or "LM"), SE_type (bootstrap or model),
#treatment coef, SE, RMSE, Coverage probability and rejection (0 or 1) for the main effect of A
res_df = data.frame()

options(warn=-1)

################################################################################
#-------------------------------------------------------------------------------
# Main code 
#-------------------------------------------------------------------------------
################################################################################

start_time = Sys.time()

for (nrep in 1:nreps){
  cat('--------------------------------------','\n')
  cat('Rep#=',nrep,'\n')
  currentseed=.Random.seed[2]
  currentseed=.Random.seed[currentseed+2]
  cat('current seed=',currentseed,'\n')
  cat('--------------------------------------','\n')
  

#-------------------------------------------
# Generate data
#-------------------------------------------
cat("Generate data","\n")

#Get long and wide data (with missingness if miss_type is not "none")
sim.dats = gen.data(seed = currentseed, alpha = alpha, beta_age = beta_age, beta_ri= beta_ri, beta_a= beta_a, 
                    beta_1b= beta_1b, beta_0b= beta_0b, rho= rho, gam= gam,
         nperson= nperson, frac_treat=frac_treat, miss_type=miss_type, miss_prob=miss_prob,
         alpha_m=alpha_m, beta_m0=beta_m0,beta_ma=beta_ma,beta_mt=beta_mt,beta_mnr = beta_mnr,
         beta_mr=beta_mr, beta_yt=beta_yt, z975,minage=minage,maxagedel=maxagedel)

wide.dat = sim.dats[[1]] #for linear model
long.dat = sim.dats[[2]] #for linear mixed model

#-------------------------------------------
# Impute data
#-------------------------------------------

if (impute==1){
  imps = impute.data(long.dat,wide.dat,M=M,SA=SA)
  lmm.imps = imps[[1]]
  lm.imps = imps[[2]]
}


#-------------------------------------------
# Get and test contrasts based on model SE
#-------------------------------------------
cat("Computing and testing contrasts","\n")

p0_hat = (sum(wide.dat$a==0)-sum(wide.dat$nr0==1))/(sum(wide.dat$a==0))
p1_hat = (sum(wide.dat$a==1)-sum(wide.dat$nr1==1))/(sum(wide.dat$a==1))


# LM model
if(LM == 1){ 
  
  #Fit model 
  if (impute == 1){ #using imputed data
    
    #Obtain list of model fits for each imputed dataset
    lm.fits = lapply(1:M,fit.imputed,mod.type="LM",allimp=lm.imp) #fit.imputed() is a function in "1a-ImputationFunctions.R"
    pooled.fits = pool(lm.fits)
    fit = summary(pooled.fits)
    
    A_trtcoefLM = fit$estimate[fit$term=="a"] 
    
    #Standard error calculated with Rubin's rules
    A_seLM = fit$std.error[fit$term=="a"]
  } 
  else #using observed data 
  {
    fit=lm(y_avg~ba+y0+a+b1+b0+nr1+nr0,data=wide.dat) 

    # TEST FOR MAIN EFFECT OF A - LM
    K=rep(0,times=8) #8 elements because we took out ym1
    dim(K)=c(1,8)
    K[4] = 1   #Beta_a hat 
    c_LM=glht(fit,linfct=K)
    #z=summary(c_LM,test=Chisqtest())
    
    #Trt coef and SE
    A_trtcoefLM = coef(c_LM)[1]
    A_seLM = sqrt(vcov(c_LM)[1,1])
  }
  
  #Compute RMSE, Coverage probability for 95% CI and hypothesis test
  
  #Old Contrast: true_LM = beta_a + (1-p_1)*beta_1nr + (1-p_0)*beta_0nr
  true_LM = beta_a
  testA_LM = RMSE_CovP_Rej(truth=true_LM, trtcoef=A_trtcoefLM, SE=A_seLM, zstar= z975) #function defined in 1-SimFunctions.R
  A_RMSELM = testA_LM[1]; A_CovPLM = testA_LM[2]; A_rejLM = testA_LM[3]
  #Note:testA_LM is a 3 element vector of the RMSE value, Coverage prob and whether we rejected (0 or 1)
  
  #Add results to data frame for main effect A
  res_df = rbind(res_df,c(nrep,"LM","Ordinary SE",A_trtcoefLM,A_seLM,A_RMSELM,A_CovPLM,A_rejLM))
}

# LMM model
if(LMM == 1){ 
  
  # Fit model
  if (impute == 1){ #using imputed data
    
    #Obtain list of model fits for each imputed dataset
    lmm.fits = lapply(1:M,fit.imputed,mod.type="LMM",allimp=lmm.imp) #fit.imputed() is a function in "1a-ImputationFunctions.R"
    pooled.fits = pool(lmm.fits)
    fit = summary(pooled.fits)
    
    A_trtcoefLMM = (fit$estimate[fit$term=="a21"] + fit$estimate[fit$term=="a21"] +
                   fit$estimate[fit$term=="a22"] + fit$estimate[fit$term=="a24"])/4
    #Still working on this
    A_seLMM 
  } 
  else #using raw data
  { 
  fit=lmer(y~0+ba+runin+first_half+time_1st_half+time2_1st_half+time3_1st_half+second_half+time_2nd_half+time2_2nd_half+time3_2nd_half+
               a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15+a16+a17+a18+a19+a20+a21+a22+a23+a24+
               b113+b114+b115+b116+b117+b118+b119+b120+b121+b122+b123+b124+
               b013+b014+b015+b016+b017+b018+b019+b020+b021+b022+b023+b024+
               nr113+nr114+nr115+nr116+nr117+nr118+nr119+nr120+nr121+nr122+nr123+nr124+
               nr013+nr014+nr015+nr016+nr017+nr018+nr019+nr020+nr021+nr022+nr023+nr024+
               (runin+ time_1st_half+ time_2nd_half | id), control=lmerControl(check.nobs.vs.nRE="ignore"), data=long.dat)
 
  # TEST FOR MAIN EFFECT OF A - LMM
  K=rep(0,times=82)
  dim(K)=c(1,82)
  K[31]=1/4
  K[32]=1/4
  K[33]=1/4
  K[34]=1/4
  c_LMM=glht(fit,linfct=K)
  #z=summary(c_LMM,test=Chisqtest())
  
  #Trt coef and SE
  A_trtcoefLMM = coef(c_LMM)[1]
  A_seLMM = sqrt(vcov(c_LMM)[1,1])
  }
  
  #Compute RMSE, Coverage probability for 95% CI and hypothesis test
  
 #Old contrast: true_LMM = beta_a + ((1-p_1)*(beta_nr121+beta_nr122+beta_nr123+beta_nr124))/4 -
   # ((1-p_0)*(beta_nr021+beta_nr022+beta_nr023+beta_nr024))/4
  true_LMM = beta_a 
  testA_LMM = RMSE_CovP_Rej(truth=true_LMM, trtcoef=A_trtcoefLMM, SE=A_seLMM, zstar= z975) #function defined in 1-SimFunctions.R
  A_RMSELMM = testA_LMM[1]; A_CovPLMM = testA_LMM[2]; A_rejLMM = testA_LMM[3]
  #Note: testA_LMM is a 3 element vector of the RMSE value, Coverage prob and whether we rejected (0 or 1)
  
  #Add results to data frame for main effect A
  res_df = rbind(res_df,c(nrep,"LMM","Ordinary SE",A_trtcoefLMM,A_seLMM,A_RMSELMM,A_CovPLMM,A_rejLMM))
}

#-------------------------------------------
# Bootstrap 
#-------------------------------------------
if (B!=0){
cat("Starting bootstrap","\n")

boot.trtcoefs_LM = c() #initialize
boot.trtcoefs_LMM = c() #initialize

for (b in 1:B){ 
cat("Bootstrap #:",b,"\n")

#----------------
#Bootstrap sample
#----------------
  nperson = dim(wide.dat)[1] #For people who have ALL y21-y24 missing, I deleted them since they have no y_avg
  
  #Bootstrap wide data (for 'LM' model)
   wide.boot <- sample_n(wide.dat,nperson,replace=TRUE)

   #Bootstrap long data (for 'LMM')
   long.boot <- data.frame()
   id.n=1
   for (i in wide.boot$id){
     idat <- long.dat %>% filter(id==i)
     idat$id = id.n
     long.boot <- rbind(long.boot,idat)
     id.n= id.n+1
   }
  wide.boot$id = 1:nperson #edit wide bootstrap data id column

#----------------
#Get contrast 
#----------------

  p0_boot = (sum(wide.boot$a==0)-sum(wide.boot$nr0==1))/(sum(wide.boot$a==0))
  p1_boot = (sum(wide.boot$a==1)-sum(wide.boot$nr1==1))/(sum(wide.boot$a==1))
  
  if(LM == 1){  #Get contrast based on wide.boot data
    MainA_bootLM = LM_contrasts(fit_dat = wide.boot,p1_hat = p1_boot, p0_hat = p0_boot)
    boot.trtcoefs_LM = c(boot.trtcoefs_LM, MainA_bootLM[1])
  }
  
  if(LMM == 1){ #Get contrast based on long.boot data
    MainA_bootLMM = LMM_contrasts(fit_dat = long.boot,p1_hat = p1_boot, p0_hat = p0_boot)
    boot.trtcoefs_LMM = c(boot.trtcoefs_LMM, MainA_bootLMM[1]) #add contrast for bootstrap sample
  }
  
} #End bootstrap loop

#----------------
# Bootstrap SE, RMSE etc.
#----------------

if (LM==1){
  boot.A_seLM = sd(boot.trtcoefs_LM) #bootstrap SE
  
  #Get Coverage and do hypothesis testing based on bootstrap SE
  boot.testA_LM = RMSE_CovP_Rej(truth=true_LM, trtcoef=A_trtcoefLM, SE= boot.A_seLM, zstar= z975) #function defined in 1-SimFunctions.R
  boot.A_RMSE = boot.testA_LM[1]; boot.A_CovP = boot.testA_LM[2]; boot.A_rej = boot.testA_LM[3]
  
  #Add bootstrap results to data frame for main effect A
  res_df = rbind(res_df,c(nrep,"LM","Bootstrap SE",A_trtcoefLM,boot.A_seLM,"--",boot.A_CovP,boot.A_rej))
}

if (LMM==1){
  boot.A_seLMM = sd(boot.trtcoefs_LMM)

  #Get Coverage and do hypothesis testing based on bootstrap SE
  boot.testA_LMM = RMSE_CovP_Rej(truth=true_LMM, trtcoef=A_trtcoefLMM, SE= boot.A_seLMM, zstar= z975) #function defined in 1-SimFunctions.R
  boot.A_RMSE = boot.testA_LMM[1]; boot.A_CovP = boot.testA_LMM[2]; boot.A_rej = boot.testA_LMM[3]
  
  #Add bootstrap results to data frame for main effect A
  res_df = rbind(res_df,c(nrep,"LMM","Bootstrap SE", A_trtcoefLMM, boot.A_seLMM, "--", boot.A_CovP, boot.A_rej))
}
}

} #end rep loop
names(res_df) = c("RCT #","Model","SE_type","A_trtcoef","A_se","A_RMSE","A_CovP","A_rej")
res_df$A_trtcoef = as.numeric(res_df$A_trtcoef); res_df$A_se = as.numeric(res_df$A_se)
res_df$A_RMSE = as.numeric(res_df$A_RMSE); res_df$A_CovP=as.numeric(res_df$A_CovP)
res_df$A_rej = as.numeric(res_df$A_rej)


#-------------------------------------------------------------------------------
# Compute simulation standard error
#-------------------------------------------------------------------------------

# Empirical simulation SE - the standard error of the final trt A estimate (across all the reps)
lm_A_sq.err = (res_df$A_trtcoef[res_df$Model=="LM"] - mean(res_df$A_trtcoef[res_df$Model=="LM"],na.rm=TRUE))^2 #vector of sq.errors
lm_A_sse = sqrt(sum(lm_A_sq.err,na.rm=TRUE)/((nreps-1)*nreps))

lmm_A_sq.err = (res_df$A_trtcoef[res_df$Model=="LMM"] - mean(res_df$A_trtcoef[res_df$Model=="LMM"],na.rm=TRUE))^2 #vector of sq.errors
lmm_A_sse = sqrt(sum(lmm_A_sq.err,na.rm=TRUE)/((nreps-1)*nreps)) 


#Print Results
write_results=function() {
  cat('Simulation to mimic SMART trial','\n'
      ,'James F. Troendle and Aparajita Sur, July 2022','\n'
      ,'-------------------------------------------------------------------------','\n'
      ,'Random seed for R=                        ',seedstart,'\n'
      ,'Number of reps=                           ',nreps,'\n'
      ,'Number of subjects=                       ',nperson,'\n'
      ,'Fraction of subjects in first stage arm1= ',frac_treat,'\n'
      ,'Number of bootstraps for variance est=    ',B,'\n'
      ,'-------------------------------------------------------------------------','\n'
      ,'-------------------------------------------------------------------------','\n'
      ,'Model Alpha=                              ',alpha,'\n'
      ,'Beta for Baseline Age=                    ',beta_age,'\n'
      ,'Model Beta for 1st stage trt=             ',beta_a,'\n'
      ,'Model Beta for 2nd stage trt when a=1=    ',beta_1b,'\n'
      ,'Model Beta for 2nd stage trt when a=0=    ',beta_0b,'\n'
      #,'Model Beta for NR when a=1=               ',beta_1nr,'\n'
      #,'Model Beta for NR when a=-1=              ',beta_0nr,'\n'
      #,'prob of response when a=1=                ',p_1,'\n'
      #,'prob of response when a=0=                ',p_0,'\n'
      ,'corr parameter=                           ',rho,'\n'
      ,'error sd parameter=                       ',gam,'\n'
      ,'Missingness Mechanism=                    ',miss_type,'\n'
      ,'-------------------------------------------------------------------------','\n'
      ,'---------------------------------------------------------------------------------','\n'
      ,' MAIN A                       LM      LMM      ','\n'
      ,'---------------------------------------------------------------------------------','\n')
  cat(' ---------------------------------------------------------------------------------','\n')
  cat(sprintf("%s%8.1f%s%8.1f\n",
              ' Avg. treatment coef=    ',round(mean(res_df$A_trtcoef[res_df$Model=="LM"],na.rm=TRUE),1),
              ' ',round(mean(res_df$A_trtcoef[res_df$Model=="LMM"],na.rm=TRUE),1)))
  cat(sprintf("%s%8.1f%s%8.1f\n",
              ' Avg. SE(trt coef)=      ',round(mean(res_df$A_se[res_df$Model=="LM"&res_df$SE_type=="Ordinary SE"]),1),
              ' ',round(mean(res_df$A_se[res_df$Model=="LMM"]),1)))
  cat(sprintf("%s%8.1f%s%8.1f\n",
              ' Sim 2*SE(trt coef)=     ',round(lm_A_sse,1),
              ' ',round(lmm_A_sse,1)))
  cat(sprintf("%s%8.1f%s%8.1f\n",
              ' RMSE treatment coef=    ',round(mean(res_df$A_RMSE[res_df$Model=="LM"],na.rm=TRUE),1),
              ' ',round(mean(res_df$A_RMSE[res_df$Model=="LMM"],na.rm=TRUE),1)))
  cat(sprintf("%s%6.1f%s%6.1f\n",
              ' Power/Type I error(%)=    ',round((sum(res_df$A_rej[res_df$Model=="LM"&res_df$SE_type=="Ordinary SE"])/nreps)*100,2),
              '   ',round((sum(res_df$A_rej[res_df$Model=="LMM"&res_df$SE_type=="Ordinary SE"])/nreps)*100,2)))
  cat(sprintf("%s%6.1f%s%6.1f\n",
              ' Coverage of 95% CI(%)=    ',round((sum(res_df$A_CovP[res_df$Model=="LM"&res_df$SE_type=="Ordinary SE"])/nreps)*100,2),
              '   ',round((sum(res_df$A_CovP[res_df$Model=="LMM"&res_df$SE_type=="Ordinary SE"])/nreps)*100,2)))
  cat(sprintf("%s%8.1f%s%8.1f\n",
              ' Avg. BOOT SE(trt coef)= ',round(mean(res_df$A_se[res_df$Model=="LM"&res_df$SE_type=="Bootstrap SE"]),1),
              ' ',round(mean(res_df$A_se[res_df$Model=="LMM"&res_df$SE_type=="Bootstrap SE"]),1)))
  cat(sprintf("%s%6.1f%s%6.1f\n",
              ' BOOTPower/Type I error(%)=',round((sum(res_df$A_rej[res_df$Model=="LM"&res_df$SE_type=="Bootstrap SE"])/nreps)*100,2),
              '   ',round((sum(res_df$A_rej[res_df$Model=="LMM"&res_df$SE_type=="Bootstrap SE"])/nreps)*100,2)))
  cat(sprintf("%s%6.1f%s%6.1f\n",
              ' BOOTCoverage of 95% CI(%)=',round((sum(res_df$A_CovP[res_df$Model=="LM"&res_df$SE_type=="Bootstrap SE"])/nreps)*100,2),
              '   ',round((sum(res_df$A_CovP[res_df$Model=="LMM"&res_df$SE_type=="Bootstrap SE"])/nreps)*100,2)))
  cat(' ---------------------------------------------------------------------------------','\n')
  cat(' Tests target 5% significance level','\n')
  cat(' ---------------------------------------------------------------------------------','\n')
  cat('----------------------------------------------------------------------------------','\n')
  
  names(res_df) = c("RCT #","Model","SE Type","Treatment Coef. of A","SE","RMSE","Coverage Prob (%)","Reject")
  
  cat('INDIVIDUAL REP RESULTS','\n','\n')
  print(res_df)
}

write_results()

end_time = Sys.time()
end_time-start_time #how long the simulation takes
