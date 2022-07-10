require(mice)
################################ mice.impute.lmm #####################################
### THE FOLLOWING IS A FUNCTION TO PERFORM ONE SET OF IMPUTATIONS. BY SETTING
### method="lmm" IN A CALL TO THE FUNCTION mice FROM THE MICE PACKAGE THE DESIRED
### NUMBER OF IMPUTATIONS CAN BE PERFORMED. A FEW MINOR ADJUSTMENTS TO FUNCTIONS sampler
### AND mice FROM THE MICE PACKAGE ARE REQUIRED FOR THIS NEW METHOD TO BE RECOGNISED
### AND FOR THE CORRECT ARGUMENTS TO BE PASSED ON. THE MODIFIED FUNCTIONS sampler2
### AND mice2 ARE ALSO GIVEN BELOW.
mice.impute.lmm <-
  function (ry, namesX, namesZ, nameG, nam, data,weights=NULL,logistic=FALSE)
  {c
    ### INPUT PARAMETERS ###
    #ry : vector indicating whether response was observed (1) or not (0).
    #namesX : vector with names of the covariates with fixed effect
    #namesZ: vector with names of the covariates with random effect (must be a
    # subset of namesX)
    #nameG : vector with name of grouping variable for the random effects
    #(i.e. id of patients)
    #nam : vector with name of outcome variable
    #data : data set to interpret the above names
    #weights : optional precision weights vector to be used in model fitting
    #logistic: if true, residual errors drawn from a logistic distribution
    #instead of a normal distribution, as in the sleep-maintenance insomnia study
    require(lme4) # Package to fit LMM
    ### CREATE MODEL FORMULA ###
    if(length(nameG)==0)
      stop(paste("No grouping factor to include random effects in model."))
    if (length(namesZ)==0)
      formulaGLMER=as.formula(paste(nam,"~1+",paste(namesX,collapse="+"),
                                    paste("+(1|",nameG,")")))
    if (length(namesZ)!=0)
      formulaGLMER=as.formula(paste(nam,"~1+",paste(namesX,collapse="+"),
                                    paste("+(1+",paste(namesZ,collapse="+"),"|",nameG,")")))
    ### FIT LMM FROM AVAILABLE CASES ###
    suppressWarnings(fit <- glmer(formulaGLMER, data[ry==1,], family =
                                    gaussian(link = "identity"),weights=weights[ry==1]))
    ### RETRIEVE FIXED EFFECTS ESTIMATES ###
    Beta <- fit@fixef
    VarBeta<-t(chol(as(vcov(fit),"matrix")))
    ### RETRIEVE RANDOM EFFECTS PREDICTIONS AND CONDITIONAL VAR-COV MATRIX ###
    B<-fit@ranef
    temp<- attr(ranef(fit, postVar = TRUE)[[1]],"postVar",exact=TRUE)
    nRanef<-length(namesZ)+1
    
    VarB<-matrix(0,nRanef*nLevels,nRanef*nLevels)
    for(j in 1:nLevels)
    { for(i in 1:nRanef)
    {for(k in 1:nRanef)
    {VarB[j+(i-1)*nLevels,j+(k-1)*nLevels]=temp[i,k,j]}
    }
    }
    VarB<-t(chol(VarB))
    ### DRAW FIXED AND RANDOM EFFECTS FROM THEIR ASYMPTOTIC DISTRIBUTION ###
    b.star<- B + VarB %*% rnorm(ncol(VarB))
    beta.star <- Beta + VarBeta %*% rnorm(ncol(VarBeta))
    ### GET TRANSPOSE AND DIAGONAL OPERATORS FOR SPARSE MATRICES ###
    ### FROM MATRIX PACKAGE ###
    tr<-get("t",env=environment(rep2abI),mode="function")
    diago<-get("diag",env=environment(rep2abI),mode="function")
    ### CALCULATE DEGREES OF FREEDOM OF ESTIMATOR OF RESIDUAL VARIANCE ###
    A<-fit@A
    X1<-fit@X
    M<-cBind(tr(A),X1)
    I <- matrix(0,nrow=nrow(as(A,"matrix")),ncol=nrow(as(A,"matrix")))
    I[row(I)==col(I)] <- 1
    Zero<-matrix(0,nrow=nrow(I),ncol=ncol(X1))
    N<-cBind(rBind(tr(A), I),rBind(X1,Zero))
    ddl<-nrow(X1)-sum(diago(M%*%solve(tr(N)%*%N)%*%tr(M)))
    # = number of observations used to fit the model - trace of the hat matrix
    ### DRAW THE RESIDUAL VARIANCE AND A GAUSSIAN OR LOGISTIC ERROR ##
    sigma<-sqrt((attr(VarCorr(fit),"sc",exact=TRUE)^2)*ddl/rchisq(1,ddl))
    if(logistic==FALSE)
      e<-rnorm(sum(!ry),0,sigma) else
        e<-rlogis(sum(!ry),location=0,scale=sqrt((3/(pi^2))*(sigma^2)))
    ### OBTAIN MODEL MATRICES IN SPARSE REPRESENTATION FOR ENTIRE DATA SET ###
    # 1. Create a dummy COMPLETE outcome variable and call it y2
    data2<-cbind(rnorm(nrow(data)),data)
    names(data2)[1]<-"y2"
    # 2. Call the same model using the entire data to this dummy outcome
    # variable. doFit=FALSE implies model won't actually be fitted but
    # model matrices will be created.
    if (length(namesZ)==0)
      formulaGLMER2<-as.formula(paste("y2~1+",paste(namesX,collapse="+"),
                                      paste("+(1|",nameG,")"))) else
                                        formulaGLMER2<-as.formula(paste("y2~1+",paste(namesX,collapse="+"),
                                                                        paste("+(1+",paste(namesZ,collapse="+"),"|",nameG,")")))
    suppressWarnings(fit2 <- glmer(formulaGLMER2, data2,
                                   family = gaussian(link = "identity"),doFit=FALSE))
    #3. Recover model matrices
    X<-fit2$fr$X
    Z<-t(as(fit2$FL$trms[[1]]$Zt,"matrix"))
    ### VALUES TO IMPUTE MISSING OUTCOMES = LINEAR PREDICTOR + ERROR ###
    nLevels<-length(levels(as.factor(data[,nameG])))
    
    vec <- X[!ry, ] %*% beta.star + Z[!ry, ] %*% b.star + e
    return(vec)
  }
### COMMAND FOR FUNCTION TO BE RECOGNISED BY MICE
environment(mice.impute.lmm)<-environment(mice)
###################################SAMPLER2#########################################
sampler2<-
  function (p, data, m, imp, r, visitSequence, maxit, printFlag, weights,logistic)
  {
    if (maxit > 0)
      chainVar <- chainMean <- array(0, dim = c(length(visitSequence),
                                                maxit, m), dimnames = list(dimnames(data)[[2]][visitSequence],
                                                                           1:maxit, paste("Chain", 1:m)))
    else chainVar <- chainMean <- NULL
    if (maxit < 1)
      iteration <- 0
    else {
      if (printFlag)
        cat("\n iter imp variable")
      for (k in 1:maxit) {
        iteration <- k
        for (i in 1:m) {
          if (printFlag)
            cat("\n ", iteration, " ", i)
          for (j in visitSequence) p$data[!r[, j], j] <- imp[[j]][,i]
          for (j in setdiff(p$visitSequence, visitSequence)) {
            cat.columns <- p$data[, p$categories[j, 4]]
            p$data[, (j:(j + p$categories[p$categories[j,4], 2] - 1))] <-
              matrix((model.matrix(~cat.columns -1)[, -1]),
                     ncol = p$categories[p$categories[j,4], 2], nrow = nrow(p$data))
          }
          for (j in p$visitSequence) {
            theMethod <- p$method[j]
            if (printFlag & theMethod != "dummy")
              cat(" ", dimnames(p$data)[[2]][j])
            if (theMethod != "" & (!is.passive(theMethod)) &
                theMethod != "dummy") {
              if (substring(theMethod, 1, 2) != "2l" &
                  substr(theMethod,1,3) != "lmm") {
                x <- p$data[, p$predictorMatrix[j, ] ==1, drop = FALSE]
                y <- p$data[, j]
                ry <- r[, j]
                nam <- dimnames(p$data)[[2]][j]
                f <- paste("mice.impute", theMethod, sep = ".")
                keep <- remove.lindep(x, y, ry)
                x <- x[, keep, drop = FALSE]
                imp[[j]][, i] <- do.call(f, args = list(y, ry, x))
              }
              else if (substring(theMethod, 1, 2) == "2l"){
                predictors <- p$predictorMatrix[j, ] !=0
                x <- p$data[, predictors, drop = FALSE]
                y <- p$data[, j]
                ry <- r[, j]
                type <- p$predictorMatrix[j, predictors]
                nam <- dimnames(p$data)[[2]][j]
                f <- paste("mice.impute", theMethod, sep = ".")
                keep <- remove.lindep(x, y, ry)
                x <- x[, keep, drop = FALSE]
                type <- type[keep]
                imp[[j]][, i] <- do.call(f, args = list(y,ry, x, type))
              }
              else if(substring(theMethod, 1, 3) == "lmm"){
                predictors_fixed <- p$predictorMatrix[j, ]>0
                predictors_random<- p$predictorMatrix[j, ] ==2
                group <- p$predictorMatrix[j, ] ==-2
                x <- p$data[, predictors_fixed, drop = FALSE]
                z <- p$data[, predictors_random, drop = FALSE]
                y <- p$data[, j]
                ry <- r[, j]
                nam <- dimnames(p$data)[[2]][j]
                namesX1<-dimnames(p$data)[[2]][predictors_fixed]
                namesZ1<-dimnames(p$data)[[2]][predictors_random]
                nameG<-dimnames(p$data)[[2]][group]
                nameG<-strsplit(as.character(nameG[1]),split="\\.")[[1]][1]
                f <- paste("mice.impute", theMethod, sep = ".")
                keep <- remove.lindep(x, y, ry)
                x<- x[, keep, drop = FALSE]
                namesX<-namesX1[keep]
                namesZ<-intersect(namesX,namesZ1)
                z <- z[, namesZ, drop = FALSE]
                imp[[j]][, i] <- do.call(f,
                                         args = list(ry, namesX, namesZ, nameG, nam,
                                                     p$data,weights,logistic))
              }
              p$data[!r[, j], j] <- imp[[j]][, i]
            }
            else if (is.passive(theMethod)) {
              imp[[j]][, i] <- model.frame(as.formula(theMethod),
                                           p$data[!r[, j], ])
              p$data[!r[, j], j] <- imp[[j]][, i]
            }
            else if (theMethod == "dummy") {
              cat.columns <- p$data[, p$categories[j, 4]]
              p$data[, (j:(j + p$categories[p$categories[j,4], 2] - 1))] <-
                matrix((model.matrix(~cat.columns -1)[, -1]),
                       ncol = p$categories[p$categories[j,4], 2], nrow = nrow(p$data))
              remove("cat.columns")
            }
            cmd <- p$post[j]
            if (cmd != "") {
              eval(parse(text = cmd))
              p$data[!r[, j], j] <- imp[[j]][, i]
            }
          }
        }
        for (j in 1:length(visitSequence)) {
          jj <- visitSequence[j]
          if (!is.factor(data[, jj])) {
            chainVar[j, k, ] <- apply(imp[[jj]], 2, var)
            chainMean[j, k, ] <- colMeans(as.matrix(imp[[jj]]))
          }
          if (is.factor(data[, jj])) {
            for (mm in 1:m) {
              nc <- as.integer(factor(imp[[jj]][, mm],
                                      levels = levels(data[, jj])))
              chainVar[j, k, mm] <- var(nc)
              chainMean[j, k, mm] <- mean(nc)}}}}
      if (printFlag)
        cat("\n")
    }
    return(list(iteration = iteration, imp = imp, chainMean = chainMean,
                chainVar = chainVar))
  }
environment(sampler2)<-environment(mice)
######################################MICE2#########################################
mice2<-
  function (data, m = 5, method = vector("character", length = ncol(data)),
            predictorMatrix = (1 - diag(1, ncol(data))),
            visitSequence = (1:ncol(data))[apply(is.na(data),2, any)],
            post = vector("character", length = ncol(data)),
            defaultMethod = c("pmm", "logreg", "polyreg"), maxit = 5,
            diagnostics = TRUE, printFlag = TRUE, seed = NA, imputationMethod = NULL,
            defaultImputationMethod = NULL, weights=as.numeric(!vector(,nrow(data))),
            logistic=FALSE)
  {
    check.visitSequence <- function(visitSequence, nmis, nvar) {
      if (!is.numeric(visitSequence)) {
        code <- pmatch(visitSequence, c("roman", "arabic","monotone",
                                        "revmonotone"))
        if (!is.na(code) && code == 1)
          visitSequence <- (1:nvar)[nmis > 0]
        if (!is.na(code) && code == 2)
          visitSequence <- rev((1:nvar)[nmis > 0])
        if (!is.na(code) && code == 3)
          visitSequence <- order(nmis)[nmis > 0]
        if (!is.na(code) && code == 4)
          visitSequence <- rev(order(nmis)[nmis > 0])
        if (is.na(code))
          stop("Argument visitSequence not recognized.\n")
      }
      if (all(nmis[visitSequence] == 0))
        stop(paste("No missing values found."))
      flags <- nmis == 0 & is.element(1:nvar, visitSequence)
      if (any(flags))
        stop(paste("Columns ", paste((1:nvar)[flags], collapse = " ",sep = ","),
                   " requested to be imputed, but contain no missing values."))
      flags <- visitSequence > nvar
      if (any(flags))
        stop(paste("Column numbers ", paste(visitSequence[flags],collapse = " ",
                                            sep = ","), " in visitSequence too large."))
      flags <- visitSequence < 1
      if (any(flags))
        stop(paste("Column numbers ", paste(visitSequence[flags], collapse = " ",
                                            sep = ","), " in visitSequence too small."))
      return(visitSequence)
    }
    check.predictorMatrix <- function(predictorMatrix, method,varnames, nmis, nvar)
    {
      if (!is.matrix(predictorMatrix))
        stop("Argument predictorMatrix should be a square matrix.")
      if (nvar != nrow(predictorMatrix) | nvar != ncol(predictorMatrix))
        stop(paste("The predictorMatrix has", nrow(predictorMatrix),
                   "rows and", ncol(predictorMatrix), "columns. Both should be",nvar, "."))
      dimnames(predictorMatrix) <- list(varnames, varnames)
      if (sum(diag(predictorMatrix) != 0))
        stop("The diagonal of predictorMatrix may contain only zeroes.")
      for (j in 1:nvar) {
        if (method[j] == "" & any(predictorMatrix[, j] ==
                                  1) & nmis[j] > 0)
          stop(paste("Variable", dimnames(predictorMatrix)[[1]][j],
                     "is used, has missing values, but is not imputed"))
        if (nmis[j] == 0 & any(predictorMatrix[j, ] == 1)) {
          predictorMatrix[j, ] <- rep(0, nvar)
        }
      }
      return(predictorMatrix)
    }
    check.method <- function(method, defaultMethod, visitSequence,
                             data, nmis, nvar) {
      if (all(method == "")) {
        for (j in visitSequence) {
          if (is.numeric(data[, j]))
            method[j] <- defaultMethod[1]
          else if (nlevels(data[, j]) == 2)
            method[j] <- defaultMethod[2]
          else if (nlevels(data[, j]) > 2)
            method[j] <- defaultMethod[3]
          else if (is.logical(data[, j]))
            method[j] <- defaultMethod[2]
          else method[j] <- defaultMethod[1]
        }
      }
      if (length(method) == 1) {
        if (is.passive(method))
          stop("Cannot have a passive imputation method for every column.")
        method <- rep(method, nvar)
      }
      if (length(method) != nvar)
        stop(paste("The length of method (", length(method),
                   ") does not match the number of columns in the data (",
                   nvar, ").", sep = ""))
      active <- !is.passive(method) & nmis > 0 & !(method =="")
      fullNames <- paste("mice.impute", method[active], sep = ".")
      notFound <- !sapply(fullNames, exists, mode = "function",
                          inherit = TRUE)
      if (any(notFound))
        stop(paste("The following functions were not found:",
                   paste(fullNames[notFound], collapse = ", ")))
      for (j in visitSequence) {
        if (is.numeric(data[, j]) & (method[j] %in% c("logreg","polyreg", "lda")))
          warning("mice: type mismatch for variable ",
                  dimnames(data)[[2]][j], ", numeric imputation function needed.",
                  call. = FALSE)
        else if (is.factor(data[, j]) & nlevels(data[, j]) ==
                 2 & (method[j] %in% c("pmm", "norm", "norm.nob","mean", "2l.norm")))
          warning("mice: type mismatch for variable ",
                  dimnames(data)[[2]][j], ", binary imputation function needed.",
                  call. = FALSE)
        else if (is.factor(data[, j]) & nlevels(data[, j]) >
                 2 & (method[j] %in% c("pmm", "norm", "norm.nob","mean",
                                       "2l.norm", "logreg")))
          warning("mice: type mismatch for variable ",
                  dimnames(data)[[2]][j], ", categorical imputation function needed.",
                  call. = FALSE)
      }
      return(method)
    }
    check.data <- function(data, predictorMatrix, method, nmis,
                           nvar) {
      is.predictor <- colSums(predictorMatrix) > 0
      is.constant <- (unlist(lapply(data[, is.predictor, drop = FALSE],
                                    var, na.rm = TRUE)) < 1e-05) | (nmis[is.predictor] ==
                                                                      rep(nrow(data), sum(is.predictor)))
      if (any(is.constant & !is.passive(method[is.predictor])))
        warning(paste("Constant predictor(s) detected:",
                      dimnames(data)[[2]][is.predictor][is.constant]))
    }
    call <- match.call()
    if (!is.na(seed))
      set.seed(seed)
    if (!(is.matrix(data) | is.data.frame(data)))
      stop("Data should be a matrix or data frame")
    if ((nvar <- ncol(data)) < 2)
      stop("Data should contain at least two columns")
    data <- as.data.frame(data)
    nmis <- apply(is.na(data), 2, sum)
    if (sum(nmis) == 0)
      stop("No missing values found")
    varnames <- dimnames(data)[[2]]
    if (!is.null(imputationMethod))
      method <- imputationMethod
    if (!is.null(defaultImputationMethod))
      defaultMethod <- defaultImputationMethod
    visitSequence <- check.visitSequence(visitSequence, nmis,
                                         nvar)
    method <- check.method(method, defaultMethod, visitSequence,
                           data, nmis, nvar)
    predictorMatrix <- check.predictorMatrix(predictorMatrix,
                                             method, varnames, nmis, nvar)
    if (maxit > 0)
      check.data(data, predictorMatrix, method, nmis, nvar)
    p <- padModel(data, method, predictorMatrix, visitSequence,
                  post, nmis, nvar)
    if (sum(duplicated(names(p$data))) > 0)
      stop("Column names of padded data should be unique")
    r <- (!is.na(p$data))
    imp <- vector("list", ncol(p$data))
    if (m > 0) {
      for (j in visitSequence) {
        imp[[j]] <- as.data.frame(matrix(NA, nrow = sum(!r[,j]), ncol = m))
        dimnames(imp[[j]]) <- list(row.names(data)[r[, j] ==FALSE], 1:m)
        y <- data[, j]
        ry <- r[, j]
        if (method[j] != "") {
          for (i in 1:m) {
            if (nmis[j] < nrow(data) - 1) {
              imp[[j]][, i] <- mice.impute.sample(y, ry)}
            else imp[[j]][, i] <- rnorm(nrow(data))}}}}
    q <- sampler2(p, data, m, imp, r, visitSequence, maxit, printFlag,
                  weights,logistic)
    for (j in p$visitSequence) p$data[(!r[, j]), j] <- NA
    imp <- q$imp[1:nvar]
    names(imp) <- varnames
    names(method) <- varnames
    names(post) <- varnames
    names(visitSequence) <- varnames[visitSequence]
    midsobj <- list(call = call, data = as.data.frame(p$data[,1:nvar]), m = m,
                    nmis = nmis, imp = imp, method = method,
                    predictorMatrix = predictorMatrix, visitSequence = visitSequence,
                    post = post, seed = seed, iteration = q$iteration, lastSeedValue =
                      .Random.seed,
                    chainMean = q$chainMean, chainVar = q$chainVar)
    if (diagnostics)
      midsobj <- c(midsobj, list(pad = p))
    oldClass(midsobj) <- "mids"
    return(midsobj)
  }
environment(mice2)<-environment(mice)
