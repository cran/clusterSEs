#' Wild Cluster Bootstrapped p-Values for Linear Family GLM
#'
#' This software estimates p-values using wild cluster bootstrapped t-statistics for linear family GLM models (Cameron, Gelbach, and Miller 2008). Residuals are repeatedly re-sampled by cluster to form a pseudo-dependent variable, a model is estimated for each re-sampled data set, and inference is based on the sampling distribution of the pivotal (t) statistic. Users may choose whether to impose the null hypothesis for independent variables; the null is never imposed for the intercept or any model that includes factor variables. Confidence intervals are only reported when the null hypothesis is \emph{not} imposed.
#'
#' @param mod A linear (identity link) model estimated using \code{glm}.
#' @param dat The data set used to estimate \code{mod}.
#' @param cluster A formula of the clustering variable.
#' @param ci.level What confidence level should CIs reflect? (Note: only reported when \code{impose.null == FALSE}).
#' @param impose.null Should we impose the null Ho?
#' @param boot.reps The number of bootstrap samples to draw.
#' @param report Should a table of results be printed to the console?
#' @param prog.bar Show a progress bar of the bootstrap (= TRUE) or not (= FALSE).
#'
#' @return A list with the elements
#' \item{p.values}{A matrix of the estimated p-values.}
#' \item{ci}{A matrix of confidence intervals (if null not imposed).}
#' @author Justin Esarey
#' @note Code to estimate GLM clustered standard errors by Mahmood Ara: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/.
#' @examples
#' \dontrun{
#' # predict chick weight using diet, do not impose the null hypothesis
#' # because of factor variable "Diet"
#' data(ChickWeight)
#' weight.mod <- glm(formula = weight~Diet,data=ChickWeight)
#' out <-cluster.wild(weight.mod ,dat =ChickWeight,cluster = ~Chick,boot.reps = 1000)
#' 
#' # impose null
#' dum <- model.matrix(~ ChickWeight$Diet)
#' ChickWeight$Diet2 <- as.numeric(dum[,2])
#' ChickWeight$Diet3 <- as.numeric(dum[,3])
#' ChickWeight$Diet4 <- as.numeric(dum[,4])
#' 
#' weight.mod2 <- glm(formula = weight~Diet2+Diet3+Diet4,data=ChickWeight)
#' out.2 <-cluster.wild(weight.mod2, dat =ChickWeight,cluster = ~Chick,boot.reps = 1000)
#' }
#' @rdname cluster.wild
#' @importFrom lmtest coeftest
#' @importFrom sandwich estfun
#' @importFrom sandwich sandwich
#' @references Cameron, A. Colin, Jonah B. Gelbach, and Douglas L. Miller. 2008. "Bootstrap-Based Improvements for Inference with Clustered Errors." \emph{The Review of Economics and Statistics} 90(3): 414-427.
#' @export
#' 
#

cluster.wild<-function(mod, dat, cluster, ci.level = 0.95, impose.null = TRUE, boot.reps = 1000, report = TRUE, prog.bar = TRUE){
  
  form <- mod$formula                                     # what is the formula of this model?
  variables <- all.vars(form)                             # what variables are in this model?
  dat.t <- subset(dat, select = variables)                # keep only relevant variables
  dat.t$clust <- subset(dat, select = all.vars(cluster))  # add the cluster variable into dat.t (for NA omission)
  dat <- na.omit(dat.t)                                   # drop the NAs
  clust <- as.vector(unlist(dat$clust))                   # reintegrate cluster variable w/o NA obs
  G<-length(unique(clust))                                # how many clusters are in this model?
  ind.variables <- attr(mod$terms, "term.labels")         # what independent variables are in this model?
  "%w/o%" <- function(x, y) x[!x %in% y]                  # a little function to create a without function (see ?match)
  dv <- variables %w/o% ind.variables                     # what is the dependent variable?
  ind.variables.names <- names(coefficients(mod))         # printed names of all fitted coefficients (incld. intercept)

  # check to see whether any IVs are factors
  fac <- c()
  for(i in 1:length(ind.variables)){
    fac[i] <- is.factor(dat[,ind.variables[i]])
  }
  fac <- max(fac)
  
  # do not impose the null for factor variables
  if(fac == 1 & impose.null == TRUE){
    cat("\n","\n", "Note: null not imposed (factor variables are present).", "\n")
    impose.null<-FALSE
  }
    
  
  # load in a function to create clustered standard errors
  # by Mahmood Ara: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/
  cl   <- function(dat, fm, cluster){
    #require(sandwich, quietly = TRUE)
    #require(lmtest, quietly = TRUE)
    M <- length(unique(cluster))
    N <- length(cluster)
    K <- fm$rank
    dfc <- (M/(M-1))*((N-1)/(N-K))
    uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
    vcovCL <- dfc*sandwich(fm, meat.=crossprod(uj)/N)
    coeftest(fm, vcovCL) }
  
  
  se.clust <- cl(dat, mod, clust)[ind.variables.names,2]               # retrieve the clustered SEs
  beta.mod <- coefficients(mod)[ind.variables.names]                   # retrieve the estimated coefficients
  w <- beta.mod / se.clust                                             # calculate the wald test statistics
 
  
  # if the null is to be imposed, execute the following code
  if(impose.null==TRUE){
    
    p.store <- c()                                                            # store wild boostrapped p-values
    w.store <- matrix(data=NA, nrow=boot.reps, ncol=length(ind.variables))    # store bootstrapped test statistics
    
    if(attr(mod$terms, "intercept") == 1 ){offset <- 1}else{offset <- 0}
    
    if(prog.bar==TRUE){cat("\n")}
    for(j in 1:length(ind.variables)){
      
      if(prog.bar==TRUE){cat("Independent variable being bootstrapped: ", ind.variables[j], "\n")}
      
      # run model imposing the null hypothesis
      form.null <- as.formula( paste( variables[1], "~", ind.variables[1:length(ind.variables) %w/o% j] ) )
      mod.null <- glm(form.null, data = dat, family = mod$family)
      null.resid <- residuals(mod.null)
      
      boot.dat <- dat           # copy the data set into a bootstrap resampling dataset
      wald.store <- c()         # create a container for storing the test statistics    
      
      if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
      for(i in 1:boot.reps){
        
        if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
        
        weight <- c(1, -1)[rbinom(G, size=1, prob=0.5) + 1][match(clust, unique(clust))]       # assign wild bootstrap weights
        pseudo.resid <- null.resid*weight                                # create pseudo-residuals using weights
        pseudo.dv <- predict(mod.null)+ pseudo.resid                     # create pseudo-observations using pseudo-residuals
        boot.dat[,dv] <- pseudo.dv                                       # create a bootstrap replicate data
        
        boot.mod <- glm(form, data = boot.dat, family = mod$family)      # run a model on the bootstrap replicate data
        
        se.boot <- cl(boot.dat, boot.mod, clust)[offset+j,2]                  # retrieve the bootstrap clustered SE
        beta.boot <- coefficients(boot.mod)[offset+j]                         # store the bootstrap beta coefficient
        wald.store[i] <- beta.boot / se.boot                                  # store the bootstrap test statistic
        
      }
      if(prog.bar==TRUE){close(pb)}
      
      p.store[j] <- 1 - ( sum( abs(w[offset + j]) > abs(wald.store) ) / boot.reps )    # calculate the wild bootstrap p-value
      w.store[,j] <- wald.store
       
    }
    
    # calculate t-stat for intercept, if present, w/o imposing the null
    if(attr(mod$terms, "intercept") == 1 ){
      
      if(prog.bar==TRUE){cat("Independent variable being bootstrapped:  Intercept (null not imposed)", "\n")}
  
      # don't impose the null for the constant (but still call it null.resid)
      null.resid <- residuals(mod)
      
      boot.dat <- dat           # copy the data set into a bootstrap resampling dataset
      wald.store <- c()         # create a container for storing the test statistics    
      
      if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
      for(i in 1:boot.reps){
        
        if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
        
        weight <- c(1, -1)[rbinom(G, size=1, prob=0.5) + 1][match(clust, unique(clust))]    # assign wild bootstrap weights
        pseudo.resid <- null.resid*weight                                                   # create pseudo-residuals using weights
        pseudo.dv <- predict(mod)+ pseudo.resid                                             # create pseudo-observations using pseudo-residuals
        boot.dat[,dv] <- pseudo.dv                                                          # create a bootstrap replicate data
        
        boot.mod <- glm(form, data = boot.dat, family = mod$family)      # run a model on the bootstrap replicate data
        
        se.boot <- cl(boot.dat, boot.mod, clust)[1,2]                  # retrieve the bootstrap clustered SE
        beta.boot <- coefficients(boot.mod)[1]                         # store the bootstrap beta coefficient
        wald.store[i] <- (beta.boot - beta.mod[1]) / se.boot           # store the bootstrap test statistic
        
      }
      if(prog.bar==TRUE){close(pb)}
      
      p.store <- c( 1 - ( sum( abs(w[1]) > abs(wald.store) ) / boot.reps ), p.store)    # calculate the wild bootstrap p-value
      w.store <- cbind(wald.store, w.store)
      
      ci.lo = NULL
      ci.hi = NULL
      print.ci = NULL
      out.ci = NULL
      
    }

  # if the null is NOT to be imposed...
  }else{
    
    if(prog.bar==TRUE){cat("Wild Cluster bootstrapping w/o imposing null...", "\n")}

    resid <- residuals(mod)
    
    boot.dat <- dat                                                                 # copy the data set into a bootstrap resampling dataset
    w.store <- matrix(data=NA, nrow=boot.reps, ncol=length(ind.variables.names))    # store bootstrapped test statistics
    
    resid <- residuals(mod)                                                         # get the residuals for the model
    
    if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
    for(i in 1:boot.reps){
      
      if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
      
      weight <- c(1, -1)[rbinom(G, size=1, prob=0.5) + 1][match(clust, unique(clust))]       # assign wild bootstrap weights
      pseudo.resid <- resid*weight                                                           # create pseudo-residuals using weights
      pseudo.dv <- predict(mod)+ pseudo.resid                                                # create pseudo-observations using pseudo-residuals
      boot.dat[,dv] <- pseudo.dv                                                             # create a bootstrap replicate data
      
      boot.mod <- glm(form, data = boot.dat, family = mod$family)           # run a model on the bootstrap replicate data
      
      se.boot <- cl(boot.dat, boot.mod, clust)[,2]                          # retrieve the bootstrap clustered SE
      beta.boot <- coefficients(boot.mod)                                   # store the bootstrap beta coefficient
      w.store[i,] <- (beta.boot-beta.mod) / se.boot                         # store the bootstrap test statistic
      
    }
    if(prog.bar==TRUE){close(pb)}
    
    comp.fun<-function(vec2, vec1){as.numeric(vec1>vec2)}                              # a simple function comparing v1 to v2
    p.store.s <- t(apply(X = abs(w.store), FUN=comp.fun, MARGIN = 1, vec1 = abs(w)))   # compare the BS test stats to orig. result
    p.store <- 1 - ( colSums(p.store.s) / dim(w.store)[1] )                            # calculate the cluster bootstrap p-value
    

    # compute critical t-statistics for CIs
    crit.t <- apply(X=abs(w.store), MARGIN=2, FUN=quantile, probs=ci.level )
    ci.lo <- beta.mod - crit.t*se.clust
    ci.hi <- beta.mod + crit.t*se.clust

    

    print.ci <- cbind(ind.variables.names, ci.lo, ci.hi)
    print.ci <- rbind(c("variable name", "CI lower", "CI higher"), print.ci)
    
    out.ci <- cbind(ci.lo, ci.hi)
    rownames(out.ci) <- ind.variables.names
    colnames(out.ci) <- c("CI lower", "CI higher")
        
  }
  
  out <- matrix(p.store, ncol=1)
  colnames(out) <- c("cluster-adjusted p-value")
  rownames(out) <- ind.variables.names
  out.p <- cbind(ind.variables.names, round(out, 3))
  out.p <- rbind(c("variable name", "cluster-adjusted p-value"), out.p)
  
  printmat <- function(m){
    write.table(format(m, justify="right"), row.names=F, col.names=F, quote=F, sep = "   ")
  }
  
  if(report==T){
    
    cat("\n", "\n", "Wild Cluster Bootstrapped p-values: ", "\n", "\n")
    printmat(out.p)
    if(is.null(print.ci) == FALSE){
      cat("\n", "Confidence Intervals (derived from bootstrapped t-statistics): ", "\n", "\n")
      printmat(print.ci)
    }

   }
  
  out.list<-list()
  out.list[["p.values"]]<-out
  out.list[["ci"]] <- out.ci
  return(invisible(out.list))
  
  
   
}

