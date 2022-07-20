################################################################################
# This script defines functions to generate simulation data, estimate contrasts
# RMSE, Coverage probability and conduct a hypothesis test 


#-------------------------------------------------------------------------------
# Function to generate Y-values for all people 

# Inputs: t (timepoint),z (random error terms), gam (scales the random variance),
# true parameters (alpha,beta_age,beta_ri,beta_a,beta_1b,beta_0b)
# and data (ba = baseline age, a = first arm trt, b = 2nd arm trt)

# Output: yt value
#-------------------------------------------------------------------------------
#Notes: beta_1b, beta_0b and b (2nd arm trt) are only needed for stage 2 so the default is NULL 
#       (because generating stage 1 data does not use those parameters)
Gen.yt <- function(t,z,gam,alpha,beta_age,ba,beta_ri,
                   beta_a,a,beta_1b=NULL,beta_0b=NULL,b=NULL){
  
  if (t < 1){ 
    yt = alpha +beta_age*(ba+(t/52)) +beta_ri + z[,t+2]*gam
  } 
  #Stage 1 simulated values (weeks 1-12)
  else if (t < 13){
    yt = alpha + beta_age*(ba+(t/52)) + a*beta_a + z[,t+2]*gam
  } 
  #Stage 2 simulated values (weeks 13-24)
  else {
    yt = alpha + beta_age*(ba+(t/52)) + a*beta_a +a*b*beta_1b + (1-a)*b*beta_0b + z[,t+2]*gam
  }
  # z[,t+2] is the column/vector of random errors for yt, where each element is the yt error for a person
  return(yt)
}


#-------------------------------------------------------------------------------
# Function to induce missingness

# Inputs: nperson, t (time point), yt (vector of y observations at time t),
#      miss_type, miss_prob (if mcar), wide.df (current data frame),
#      and logistic parameters to determine missingness


# Output: data frame with missing values (aka the 'observed data')
#-------------------------------------------------------------------------------

#Function to calculate missingness for a time-point based on missingness mechanism (MAR or MNAR)
induce.missYt = function(nperson,t,yt,wide.df,alpha_m=NULL, beta_m0=NULL, 
                         beta_ma=NULL,beta_mt=NULL,beta_mnr=NULL,beta_yt=NULL,beta_mr=NULL,frac_treat=NULL,miss_type, miss_prob =NULL){
  
  rand.prob = runif(nperson,min=0,max=1) #vector of uniform distributed values for each person 
  
  if (miss_type=="MCAR"){yt.missprob = rep(miss_prob,nperson)} 
    
  if(miss_type=="MAR"){ 
    if (t<13){ #stage 1
      yt.missprob = exp(alpha_m+beta_m0*(wide.df$y0-6000)+beta_ma*(wide.df$a-frac_treat)+beta_mt*(t-12))/
        (1+exp(alpha_m+beta_m0*(wide.df$y0-6000)+beta_ma*(wide.df$a-frac_treat)+beta_mt*(t-12)))
    }
    if (t>12) #stage 2
      yt.missprob = exp(alpha_m+beta_m0*(wide.df$y0-6000)+beta_ma*(wide.df$a-frac_treat)+beta_mt*(t-12)+beta_mnr*(wide.df$nr)*(t>13))/
        (1+exp(alpha_m+beta_m0*(wide.df$y0-6000)+ beta_ma*(wide.df$a-frac_treat)+beta_mt*(t-12)+beta_mnr*(wide.df$nr)*(t>13)))
  }
  
  if(miss_type=="MNAR"){ 
    if (t<13){ #stage 1
      yt.missprob = exp(alpha_m+beta_m0*(wide.df$y0-6000)+beta_ma*(wide.df$a-frac_treat)+beta_mt*(t-12)+beta_yt*(yt-mean(yt)))/
        (1+exp(alpha_m+beta_m0*(wide.df$y0-6000)+beta_ma*(wide.df$a-frac_treat)+beta_mt*(t-12)+beta_yt*(yt-mean(yt))))
    }
    if (t>12) #stage 2
      yt.missprob = exp(alpha_m+beta_m0*(wide.df$y0-6000)+beta_ma*(wide.df$a-frac_treat)+beta_mt*(t-12)+beta_mnr*(wide.df$nr)*(t>13)+
                          beta_yt*(yt-mean(yt)) + beta_mr*((wide.df$y21+wide.df$y22+wide.df$y23+wide.df$y24)/4))/
        (1+exp(alpha_m+beta_m0*(wide.df$y0-6000)+ beta_ma*(wide.df$a-frac_treat)+beta_mt*(t-12)+beta_mnr*(wide.df$nr)*(t>13)+
                 beta_yt*(yt-mean(yt)) + beta_mr*((wide.df$y21+wide.df$y22+wide.df$y23+wide.df$y24)/4)))
  }
  yt[rand.prob < yt.missprob] = NA
  return(yt)
}


#-------------------------------------------------------------------------------
# Function to generate data for simulation

#Inputs: true parameters and other sim parameters 

#Output: wide and long data frames of generated data for that rep
#-------------------------------------------------------------------------------

gen.data = function(seed,alpha, beta_age, beta_ri, beta_a, beta_1b, beta_0b, rho, gam,
                    nperson, frac_treat, miss_type, miss_prob, alpha_m, 
                    beta_m0,beta_ma,beta_mt,beta_mnr,
                    beta_mr, beta_yt, z975,minage,maxagedel){
  
  set.seed(seed)
  #Generate everyone's baseline ages
  ba=round(runif(nperson,min=minage,max=minage+maxagedel),digits=0)
  
  #Generate everyone's treatment assignment 
  a = round(rbinom(nperson,1,prob=.5), digits = 0)
  
  #Define correlation structure and determine random variance to add to expected values
  mu=rep(0,times=26)
  sigma=rep(1,times=26**2)
  dim(sigma)=c(26,26)
  ii=1
  while (ii <= 26) {
    jj=1
    while (jj <= 26) {
      if (ii != jj) {sigma[ii,jj]=exp(-1*abs(ii-jj)/rho)}
      jj=jj+1
    }
    ii=ii+1
  }
  # Calculate error terms for Ys
  z=rmvnorm(nperson,mu,sigma)
  
  #Create dataframe of simulated yt values from weeks -1 to 12 for all participants 
  wide.dat = data.frame(id=1:nperson,a=a,ba=ba)
  
  for (t in -1:12){ #Will code it to parsapply later to parallelize loop
    yt = Gen.yt(t,z,gam,alpha,beta_age,ba,beta_ri,beta_a,a) 
    wide.dat = cbind(wide.dat,yt)
  }
  names(wide.dat) = c("id","a","ba","ym1","y0","y1","y2","y3","y4","y5","y6","y7","y8","y9","y10","y11","y12")
  
  #Induce missingness in first stage (skip below if doing full data analysis)
  if (miss_type!= "NONE"){
      for (time in 1:12){ #will parallelize this later
        
        t.index = which(names(wide.dat)==paste("y",time,sep="")) #find the column where yt is in the wide data frame
        yt = wide.dat[,t.index] #extract the yt vector (everyone's y values at time t)
        
        if (miss_type=="MCAR"){yt = induce.missYt(nperson=nperson,t=time,yt = yt,
                                                  wide.df=wide.dat,miss_type="MCAR",miss_prob = miss_prob)}
        
        if (miss_type=="MAR"){yt = induce.missYt(nperson=nperson,t=time,yt = yt,
                              wide.df=wide.dat,alpha_m=alpha_m, beta_m0=beta_m0,beta_ma=beta_ma,beta_mt=beta_mt,
                              frac_treat=frac_treat,miss_type="MAR")}
        
        if (miss_type=="MNAR"){yt = induce.missYt(nperson=nperson,t=time,yt = yt,wide.df=wide.dat,
                                                alpha_m=alpha_m, beta_m0=beta_m0,beta_ma=beta_ma,beta_mt=beta_mt,
                                                beta_yt=beta_yt,frac_treat=frac_treat,miss_type="MNAR")}
        wide.dat[,t.index] = yt
    }
  }
  #NOTE: I set up this function so we could also have different missingness mechanisms at different times. Ex: miss_type could also be a vector 
  #where each element in the vector is the missingness mechanism for each time-point (miss_type[t]; it would be a 24 element vector for the 24 weeks).
  
  # Determine response status and compute nr1,nr0,y_S1avg (average of wks 9-12)
  #Initialize
  wide.dat$nr = rep(0,nperson) #non-responders 
  wide.dat$nr0 = rep(0,nperson) #non-responders in a=0
  wide.dat$nr1 = rep(0,nperson) #non-responders in a=1
  wide.dat$b = rep(0,nperson) #b = 0, -1 or 1
  wide.dat$b1 = rep(0,nperson)#b variable*I(a=1)
  wide.dat$b0 = rep(0,nperson) #b variable*I(a=0) 
  
  wide.dat$y_S1avg = rowMeans(wide.dat[,which(names(wide.dat)=="y9"):which(names(wide.dat)=="y12")],na.rm=TRUE) #Take mean of observed y9..y12 for each person (remove NAs)
  # Checked that the line of code below works
  wide.dat$nr[is.na(wide.dat$y_S1avg)| wide.dat$y_S1avg<10000] = 1 #avg is less than 10000 steps OR all 4 weeks are missing
  wide.dat$nr0[wide.dat$nr==1 & wide.dat$a==0] = 1 #0 otherwise bc vector is initialized as 0's
  wide.dat$nr1[wide.dat$nr==1 & wide.dat$a==1] = 1 #0 otherwise bc vector is initialized as 0's
  
  #Number of non-responders overall and in each first stage arm 
  n_nr = sum(wide.dat$nr) # sums of 1s (# of non-responders)
  n_nr1 = sum(wide.dat$nr1)
  n_nr0 = sum(wide.dat$nr0)
  
  #Determine 2nd stage treatment arm for non-responders and calculate b, b1, and b0
  wide.dat$b[wide.dat$nr==1] = (2*rbinom(n_nr,1,.5))-1 #b=-1or1 for the non-responders & 0 for responders (since vector is initialized to 0)
  
  #b1 and b0
  wide.dat$b1[wide.dat$a==1] = wide.dat$b[wide.dat$a==1] #otherwise 0 (bc b was a vector initialized with all 0s)
  wide.dat$b0[wide.dat$a==0] = wide.dat$b[wide.dat$a==0] #otherwise 0 
  
  #Generate Y13-24 and add to dataset
  for (t in 13:24){
    yt = Gen.yt(t,z,gam,alpha,beta_age,ba,beta_ri,beta_a,a,beta_1b,beta_0b,b=wide.dat$b) 
    wide.dat = cbind(wide.dat,yt)
  }
  names(wide.dat) = c(names(wide.dat)[1:24],"y13","y14","y15","y16","y17","y18","y19","y20","y21","y22","y23","y24") #24 variables in dataset before adding yts
  #Note: names(wide.dat)[1:24] are the names of the columns before y13-y24 columns
  
  #Induce missingness in 2nd stage (skip below if doing full data analysis)
  
  full.dat = wide.dat #For MNAR: need y21 to 24 to be not NA, otherwise we cannot compute the prob of missingness based on y21-y24
  if (miss_type!= "NONE"){
    for (time in 13:24){ #will parallelize this later
      
      t.index = which(names(wide.dat)==paste("y",time,sep="")) #find the column where yt is in the wide data
      yt = wide.dat[,t.index] #extract the yt vector (everyone's y values at time t)
      
      if (miss_type=="MCAR"){yt = induce.missYt(nperson=nperson,t=time,yt = yt,
                                                wide.df=wide.dat,miss_type="MCAR",miss_prob = miss_prob)}
      
      if (miss_type=="MAR"){yt = induce.missYt(nperson=nperson,t=time,yt = yt,
                                               wide.df=wide.dat,alpha_m=alpha_m, beta_m0=beta_m0,beta_ma=beta_ma,beta_mt=beta_mt,
                                               beta_mnr = beta_mnr,frac_treat=frac_treat,miss_type="MAR")}
      
      if (miss_type=="MNAR"){yt = induce.missYt(nperson=nperson,t=time,yt = yt,wide.df=full.dat,
                                                alpha_m=alpha_m, beta_m0=beta_m0,beta_ma=beta_ma,beta_mt=beta_mt,
                                                beta_mnr = beta_mnr, beta_mr=beta_mr, beta_yt=beta_yt,
                                                frac_treat=frac_treat,miss_type="MNAR")}
      wide.dat[,t.index] = yt
    }
  }
  
  # Calculate final outcome 
  wide.dat$y_avg = rowMeans(wide.dat[,which(names(wide.dat)=="y21"):which(names(wide.dat)=="y24")],na.rm=TRUE) #Average y21,y22,y23,y24

  # Average the run-in periods
  wide.dat$y0 = (wide.dat$ym1 + wide.dat$y0)/2 #Averaged the 2 run-ins 
  wide.dat$ym1 = NULL 
  
  # long.d = wide.dat %>% pivot_longer(wide.dat,cols=c("y0","y1","y2","y3","y4","y5","y6",
  #                                                      "y7","y8","y9","y10","y11","y12",
  #                                                      "y13","y14","y15","y16","y17","y18",
  #                                                      "y19","y20","y21","y22","y23","y24"),
  #                                      names_to= c("time"), values_to = c("y"))
  
  #Create long data - will make this faster later
  num_response=nperson*25
  almm=rep(0,times=num_response)
  b1lmm=rep(0,times=num_response)
  b0lmm=rep(0,times=num_response)
  nr1lmm=rep(0,times=num_response)
  nr0lmm=rep(0,times=num_response)
  ylmm=rep(0,times=num_response)
  idlmm=rep(0,times=num_response)
  balmm=rep(0,times=num_response)
  timelmm=rep(0,times=num_response)
  
  i=1
  newi=0
  while (i <= nperson) {
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y0[i]
    timelmm[newi]=0
    
    newi=newi+1
    almm[newi]= wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y1[i]
    timelmm[newi]=1
    
    newi=newi+1
    almm[newi]= wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y2[i]
    timelmm[newi]=2
    
    newi=newi+1
    almm[newi]= wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y3[i]
    timelmm[newi]=3
    
    newi=newi+1
    almm[newi]= wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y4[i]
    timelmm[newi]=4
    
    newi=newi+1
    almm[newi]= wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y5[i]
    timelmm[newi]=5
    
    newi=newi+1
    almm[newi]= wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y6[i]
    timelmm[newi]=6
    
    newi=newi+1
    almm[newi]= wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y7[i]
    timelmm[newi]=7
    
    newi=newi+1
    almm[newi]= wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y8[i]
    timelmm[newi]=8
    
    newi=newi+1
    almm[newi]= wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y9[i]
    timelmm[newi]=9
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y10[i]
    timelmm[newi]=10
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y11[i]
    timelmm[newi]=11
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y12[i]
    timelmm[newi]=12
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y13[i]
    timelmm[newi]=13
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y14[i]
    timelmm[newi]=14
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y15[i]
    timelmm[newi]=15
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y16[i]
    timelmm[newi]=16
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y17[i]
    timelmm[newi]=17
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y18[i]
    timelmm[newi]=18
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y19[i]
    timelmm[newi]=19
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y20[i]
    timelmm[newi]=20
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y21[i]
    timelmm[newi]=21
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y22[i]
    timelmm[newi]=22
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y23[i]
    timelmm[newi]=23
    
    newi=newi+1
    almm[newi]=wide.dat$a[i]
    b1lmm[newi]=wide.dat$b1[i]
    b0lmm[newi]=wide.dat$b0[i]
    nr1lmm[newi]=wide.dat$nr1[i]
    nr0lmm[newi]=wide.dat$nr0[i]
    idlmm[newi]=wide.dat$id[i]
    balmm[newi]=wide.dat$ba[i]
    ylmm[newi]= wide.dat$y24[i]
    timelmm[newi]=24
    
    i=i+1
  }
  
  lmm2data=data.frame(id=idlmm,ba=balmm,a=almm,b1=b1lmm,b0=b0lmm,nr1=nr1lmm,nr0=nr0lmm,y=ylmm,time=timelmm)
  lmm2data$a1=lmm2data$a*(lmm2data$time==1)
  lmm2data$a2=lmm2data$a*(lmm2data$time==2)
  lmm2data$a3=lmm2data$a*(lmm2data$time==3)
  lmm2data$a4=lmm2data$a*(lmm2data$time==4)
  lmm2data$a5=lmm2data$a*(lmm2data$time==5)
  lmm2data$a6=lmm2data$a*(lmm2data$time==6)
  lmm2data$a7=lmm2data$a*(lmm2data$time==7)
  lmm2data$a8=lmm2data$a*(lmm2data$time==8)
  lmm2data$a9=lmm2data$a*(lmm2data$time==9)
  lmm2data$a10=lmm2data$a*(lmm2data$time==10)
  lmm2data$a11=lmm2data$a*(lmm2data$time==11)
  lmm2data$a12=lmm2data$a*(lmm2data$time==12)
  lmm2data$a13=lmm2data$a*(lmm2data$time==13)
  lmm2data$a14=lmm2data$a*(lmm2data$time==14)
  lmm2data$a15=lmm2data$a*(lmm2data$time==15)
  lmm2data$a16=lmm2data$a*(lmm2data$time==16)
  lmm2data$a17=lmm2data$a*(lmm2data$time==17)
  lmm2data$a18=lmm2data$a*(lmm2data$time==18)
  lmm2data$a19=lmm2data$a*(lmm2data$time==19)
  lmm2data$a20=lmm2data$a*(lmm2data$time==20)
  lmm2data$a21=lmm2data$a*(lmm2data$time==21)
  lmm2data$a22=lmm2data$a*(lmm2data$time==22)
  lmm2data$a23=lmm2data$a*(lmm2data$time==23)
  lmm2data$a24=lmm2data$a*(lmm2data$time==24)
  lmm2data$b113=lmm2data$b1*(lmm2data$time==13)
  lmm2data$b114=lmm2data$b1*(lmm2data$time==14)
  lmm2data$b115=lmm2data$b1*(lmm2data$time==15)
  lmm2data$b116=lmm2data$b1*(lmm2data$time==16)
  lmm2data$b117=lmm2data$b1*(lmm2data$time==17)
  lmm2data$b118=lmm2data$b1*(lmm2data$time==18)
  lmm2data$b119=lmm2data$b1*(lmm2data$time==19)
  lmm2data$b120=lmm2data$b1*(lmm2data$time==20)
  lmm2data$b121=lmm2data$b1*(lmm2data$time==21)
  lmm2data$b122=lmm2data$b1*(lmm2data$time==22)
  lmm2data$b123=lmm2data$b1*(lmm2data$time==23)
  lmm2data$b124=lmm2data$b1*(lmm2data$time==24)
  lmm2data$b013=lmm2data$b0*(lmm2data$time==13)
  lmm2data$b014=lmm2data$b0*(lmm2data$time==14)
  lmm2data$b015=lmm2data$b0*(lmm2data$time==15)
  lmm2data$b016=lmm2data$b0*(lmm2data$time==16)
  lmm2data$b017=lmm2data$b0*(lmm2data$time==17)
  lmm2data$b018=lmm2data$b0*(lmm2data$time==18)
  lmm2data$b019=lmm2data$b0*(lmm2data$time==19)
  lmm2data$b020=lmm2data$b0*(lmm2data$time==20)
  lmm2data$b021=lmm2data$b0*(lmm2data$time==21)
  lmm2data$b022=lmm2data$b0*(lmm2data$time==22)
  lmm2data$b023=lmm2data$b0*(lmm2data$time==23)
  lmm2data$b024=lmm2data$b0*(lmm2data$time==24)
  lmm2data$nr113=lmm2data$nr1*(lmm2data$time==13)
  lmm2data$nr114=lmm2data$nr1*(lmm2data$time==14)
  lmm2data$nr115=lmm2data$nr1*(lmm2data$time==15)
  lmm2data$nr116=lmm2data$nr1*(lmm2data$time==16)
  lmm2data$nr117=lmm2data$nr1*(lmm2data$time==17)
  lmm2data$nr118=lmm2data$nr1*(lmm2data$time==18)
  lmm2data$nr119=lmm2data$nr1*(lmm2data$time==19)
  lmm2data$nr120=lmm2data$nr1*(lmm2data$time==20)
  lmm2data$nr121=lmm2data$nr1*(lmm2data$time==21)
  lmm2data$nr122=lmm2data$nr1*(lmm2data$time==22)
  lmm2data$nr123=lmm2data$nr1*(lmm2data$time==23)
  lmm2data$nr124=lmm2data$nr1*(lmm2data$time==24)
  lmm2data$nr013=lmm2data$nr0*(lmm2data$time==13)
  lmm2data$nr014=lmm2data$nr0*(lmm2data$time==14)
  lmm2data$nr015=lmm2data$nr0*(lmm2data$time==15)
  lmm2data$nr016=lmm2data$nr0*(lmm2data$time==16)
  lmm2data$nr017=lmm2data$nr0*(lmm2data$time==17)
  lmm2data$nr018=lmm2data$nr0*(lmm2data$time==18)
  lmm2data$nr019=lmm2data$nr0*(lmm2data$time==19)
  lmm2data$nr020=lmm2data$nr0*(lmm2data$time==20)
  lmm2data$nr021=lmm2data$nr0*(lmm2data$time==21)
  lmm2data$nr022=lmm2data$nr0*(lmm2data$time==22)
  lmm2data$nr023=lmm2data$nr0*(lmm2data$time==23)
  lmm2data$nr024=lmm2data$nr0*(lmm2data$time==24)
  lmm2data$runin=as.numeric(lmm2data$time==0)
  lmm2data$first_half=as.numeric(lmm2data$time <= 12 & lmm2data$time > 0)
  lmm2data$time_1st_half=lmm2data$time*(lmm2data$time <= 12)
  lmm2data$time2_1st_half=lmm2data$time_1st_half**2
  lmm2data$time3_1st_half=lmm2data$time_1st_half**3
  lmm2data$second_half=as.numeric(lmm2data$time > 12)
  lmm2data$time_2nd_half=lmm2data$time*(lmm2data$time > 12)
  lmm2data$time2_2nd_half=lmm2data$time_2nd_half**2
  lmm2data$time3_2nd_half=lmm2data$time_2nd_half**3

  # Change the variable types 
  wide.dat$a = as.factor(wide.dat$a)
  wide.dat$nr0 = as.factor(wide.dat$nr0)
  wide.dat$nr1 = as.factor(wide.dat$nr1)

  lmm2data$a=as.factor(lmm2data$a)
  lmm2data$nr1=as.factor(lmm2data$nr1)
  lmm2data$nr0=as.factor(lmm2data$nr0)
  lmm2data$ftime=as.factor(lmm2data$time)

  
  #return long and wide data 
  dats = list(wide.dat, lmm2data)
  return(dats)
}


#-------------------------------------------------------------------------------
# Function to get MSE, do hypothesis testing, coverage probability etc.

#Inputs: truth (value), trtcoef(value), SE (bootstrap or model SE) and zstar (critical value)

#Output: 3 element vector of MSE, Coverage Prob (binary) and Rejection (binary) 
#-------------------------------------------------------------------------------
RMSE_CovP_Rej = function(truth, trtcoef, SE,zstar){
  
  #MSE
  MSE = (trtcoef - truth)**2
  
  #Calculate coverage for 95% CI
  CovP = ifelse(truth >= trtcoef-1.96*SE & truth <= trtcoef+1.96*SE, 1, 0) 

  #Test H_0
  z= abs(trtcoef)/SE 
  Rej = ifelse(z >= zstar,1,0) 
  
  #Vector of MSE, CovP and Reg
  vec = c(sqrt(MSE),CovP,Rej) #sqrt(MSE) because we are interested in RMSE
  
  return(vec)
}



#NOTE: The main code does not use the below functions anymore
#-------------------------------------------------------------------------------
# Functions to get contrasts and SE - for LMM and LM models

#Inputs: fit_dat (data used to fit model), p1_hat and p0_hat

#Output: 2 element vector of contrast value and SE
#-------------------------------------------------------------------------------


LMM_contrasts = function(fit_dat,p1_hat, p0_hat){
  
  #Fit model
  fit=lmer(y~0+ba+runin+first_half+time_1st_half+time2_1st_half+time3_1st_half+second_half+time_2nd_half+time2_2nd_half+time3_2nd_half+
             a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15+a16+a17+a18+a19+a20+a21+a22+a23+a24+
             b113+b114+b115+b116+b117+b118+b119+b120+b121+b122+b123+b124+
             b013+b014+b015+b016+b017+b018+b019+b020+b021+b022+b023+b024+
             nr113+nr114+nr115+nr116+nr117+nr118+nr119+nr120+nr121+nr122+nr123+nr124+
             nr013+nr014+nr015+nr016+nr017+nr018+nr019+nr020+nr021+nr022+nr023+nr024+
             (runin+ time_1st_half+ time_2nd_half | id), control=lmerControl(check.nobs.vs.nRE="ignore"), data=fit_dat)
  
  #--------------------------------
  # TEST FOR MAIN EFFECT OF A - LMM
  #--------------------------------
  K=rep(0,times=82)
  dim(K)=c(1,82)
  K[31]=1/4
  K[32]=1/4
  K[33]=1/4
  K[34]=1/4
  # K[67]=(1-p1_hat)/4
  # K[68]=(1-p1_hat)/4
  # K[69]=(1-p1_hat)/4
  # K[70]=(1-p1_hat)/4
  # K[79]=-1*(1-p0_hat)/4
  # K[80]=-1*(1-p0_hat)/4
  # K[81]=-1*(1-p0_hat)/4
  # K[82]=-1*(1-p0_hat)/4
  c=glht(fit,linfct=K)
  z=summary(c,test=Chisqtest())
  
  #Trt coef and SE
  A_trtcoef = coef(c)[1]
  A_se = sqrt(vcov(c)[1,1])
  
  vec = c(A_trtcoef,A_se)
  return(vec)
}

LM_contrasts = function(fit_dat,p1_hat, p0_hat){
  
  #--------------------------------
  # TEST FOR MAIN EFFECT OF A - LM
  #--------------------------------
  
  fit=lm(y_avg~ba+y0+a+b1+b0+nr1+nr0,data=fit_dat)
  
  K=rep(0,times=8) #8 elements because we took out ym1
  dim(K)=c(1,8)
  
  #Beta_a hat 
  K[4] = 1
  
  #CONTRAST - what we used before
  # K[4]=1
  # K[7]=1-p1_hat
  # K[8]=-1*(1-p0_hat)
  
  c=glht(fit,linfct=K)
  z=summary(c,test=Chisqtest())
  
  #Trt coef and SE
  A_trtcoef = coef(c)[1]
  A_se = sqrt(vcov(c)[1,1])
  vec = c(A_trtcoef,A_se)
  return(vec)
}




