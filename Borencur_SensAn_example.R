################################## DATA SIMULATION FUNCTION ########################
# This function generates longitudinal data from an LMM with drop-outs such
# that drop-out probability depends on the last observed outcome, following the
# the set-up of Part I of the simulation study presented in the next section.
simul_data<-function(b=0.2,sd_bi=1,percent=0.4)
{ require(boot) #needed for inv.logit function
  ID<-as.factor(1:1000)
  GROUP<-c(rep(0,500),rep(1,500))
  RAN_INT<-rnorm(1000,0,sd_bi)
  RAN_SLOPE<-rnorm(1000,0,sd_bi)
  data<-data.frame(ID,GROUP)
  for(i in 0:5)
  {
    VISIT<-rep(i,1000)
    OUT<-b*VISIT*GROUP+rnorm(1000,0,1)+RAN_INT+RAN_SLOPE*VISIT
    VISIT_GROUP<-VISIT*GROUP
    data1<-data.frame(data,VISIT,VISIT_GROUP,OUT)
    if(i==0) datafinal<-data1
    else datafinal<-rbind(datafinal,data1)
  }
  lambda<-0.08462 #Change to obtain different global percentage of missing data
  ID<-1:1000
  for(i in 1:5)
  {
    j<-i-1
    OUT1<-datafinal$OUT[datafinal$VISIT==j]
    OUT2<-OUT1[!is.na(ID)]
    proba<-inv.logit(logit(lambda)-0.5*log(2)*OUT2) #Change to obtain different
    #missingness mechanism
    MISS<-rbinom(sum(!is.na(ID)),1,proba)
    ID[!is.na(ID)][MISS==1]<-NA
    is.na(datafinal$OUT[datafinal$VISIT==i])<-(is.na(ID)) }
  return(datafinal)
}
######################################## EXAMPLE ###################################
# Simulate Data
datamis<-simul_data(0.2,1,0.4)
############## MAR Analysis of the group effect theta at visit 5 #############
##### (following the sleep-maintenance insomnia study, Section 4.2 of main text)
# Create imputations using functions given in Section 1
ini <- mice2(datamis, maxit=0)
pred <- ini$pred
pred["OUT",] <- c(-2, 1, 2, 1, 0)
# NOTE! -2 indicates variable with patient id
# 2 indicates covariate with fixed and random effect
imp1<- mice2(datamis, m=20, meth=c("","","","","lmm"), pred=pred,maxit=1)
# Modify mice object to keep only data and imputations corresponding to Visit 5
mis<-imp1$data[is.na(imp1$data$OUT),]
mis<-mis$VISIT==5
imp1$imp$OUT<-imp1$imp$OUT[mis,]
imp1$data<-imp1$data[imp1$data$VISIT==5,]
imp1$nmis[5]<-sum(mis)
# Fit linear model to obtain MAR estimate of group effect at visit 5
fitMAR<-with(data=imp1, exp=lm(OUT~1+GROUP))
MAR<-pool(fitMAR,"q")
summary(MAR)
# Test significance of group effect at visit 5
fitMAR2<-with(data=imp1, exp=lm(OUT~1))
pool.compare(fitMAR,fitMAR2)$pvalue
########### Sensitivity analysis on the group effect theta at visit 5 ##########
##### (following the sleep-maintenance insomnia study, Section 4.3 of main text)
# Set informativity parameters k_0 and k_1 and calculate SD of observed outcomes
# at last visit.
k_0<-1
k_1<-1
sigma5<-sd(datamis[datamis$VISIT==5,"OUT"],na.rm=TRUE)
# Create new imputations using functions given in Section 1.
# NOTE! We do not recycle imputations from the MAR analysis. We recommend that
# imputations are redrawn for each set of informativity parameters so that the
# final results reflect the variability arising from the imputation procedure.
imp2<- mice2(datamis, m=20, meth=c("","","","","lmm"), pred=pred,maxit=1)
# As before, modify mice object to keep only data and imputations corresponding
# to Visit 5
mis<-imp2$data[is.na(imp2$data$OUT),]
mis<-mis$VISIT==5
imp2$imp$OUT<-imp2$imp$OUT[mis,]
imp2$data<-imp2$data[imp2$data$VISIT==5,]
imp2$nmis[5]<-sum(mis)
# Modify imputations according to chosen MNAR pattern mixture model
G<-imp2$data[is.na(imp2$data$OUT),]
G0<-G$GROUP==0
imp2$imp$OUT[G0,]<-imp2$imp$OUT[G0,]+k_0*sigma5
G1<-G$GROUP==1
imp2$imp$OUT[G1,]<-imp2$imp$OUT[G1,]+k_1*sigma5
# Refit linear model with new imputations to obtain new estimate
fitMNAR<-with(data=imp2, exp=lm(OUT~1+GROUP))
MNAR<-pool(fitMNAR,"q")
summary(MNAR)
# Retest significance of group effect
fitMNAR2<-with(data=imp2, exp=lm(OUT~1))
pool.compare(fitMNAR,fitMNAR2)$pvalue
