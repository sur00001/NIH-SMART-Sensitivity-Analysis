# NIH-SMART-Sensitivity-Analysis

Notes from 7/10/22:
FYI When inducing MAR and MNAR missingness, a person can be missing all y21 to y24 observations so there is no final outcome (y_avg) when doing the linear model. We cannot run this unless we impute those people. 

Also, when doing LMM with MAR - I am getting the below warning and error, for some of the bootstrap samples. I think it is because the bootstrap happened to choose people all from the same treatment group, so the model is rank deficient. When that coefficient is dropped, the contrast matrix (K) is a different dimension than the # of coefficients. I'm still investigating this and I'll need to code a "check" to make sure people are not from the same treatment group in the bootstrap sample. 

"fixed-effect model matrix is rank deficient so dropping 1 column / coefficient" is a warning and later this error appears: Error in glht.matrix(fit, linfct = K) : 
‘ncol(linfct)’ is not equal to ‘length(coef(model))’

If we are defining the "truth" to be beta_a then we don't need the bootstrap and don't need to worry about this error.
