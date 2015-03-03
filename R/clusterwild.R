#' Wild Cluster Bootstrapped p-Values for Linear Family GLM
#'
#' This software estimates p-values using wild cluster bootstrapped t-statistics for linear family GLM models (Cameron, Gelbach, and Miller 2008). Residuals are repeatedly re-sampled by cluster to form a pseudo-dependent variable, a model is estimated for each re-sampled data set, and inference is based on the sampling distribution of the pivotal (t) statistic. 
#'
#' @param mod A linear (identity link) model estimated using \code{glm}.
#' @param dat The data set used to estimate \code{mod}.
#' @param cluster A formula of the clustering variable.
#' @param boot.reps The number of bootstrap samples to draw.
#' @param report Should a table of results be printed to the console?
#' @param prog.bar Show a progress bar of the bootstrap (= TRUE) or not (= FALSE).
#'
#' @return A list with the elements
#' \item{p.values}{A matrix of the estimated p-values.}
#' @author Justin Esarey
#' @note Code to estimate GLM clustered standard errors by Mahmood Ara: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/.
#' @examples
#' \dontrun{
#' # predict prestige score of occupation
#' require(effects)
#' data(BEPS)
#' linear.model <- glm(Europe ~ age + gender + economic.cond.national, data=BEPS)
#' summary(linear.model)
#' 
#' # compute wild cluster bootstrapped p-values
#' clust.wd.p <- cluster.wild(linear.model, BEPS, ~ vote, report = T)
#' }
#' @rdname cluster.wild
#' @importFrom lmtest coeftest
#' @importFrom sandwich estfun
#' @importFrom sandwich sandwich
#' @references Cameron, A. Colin, Jonah B. Gelbach, and Douglas L. Miller. 2008. "Bootstrap-Based Improvements for Inference with Clustered Errors." \emph{The Review of Economics and Statistics} 90(3): 414-427.
#' @export
#' 
#

cluster.wild<-function(mod, dat, cluster, boot.reps = 1000, report = TRUE, prog.bar = TRUE){
  
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
  ind.variables.names <- names(coefficients(mod))
    
  
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
  w <- beta.mod / se.clust                                       # calculate the wald test statistics
  
  
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
  
  # calculate SE for intercept, if present
  if(attr(mod$terms, "intercept") == 1 ){
    
    if(prog.bar==TRUE){cat("Independent variable being bootstrapped:  Intercept", "\n")}

    # don't impose the null for the constant (but still call it null.resid)
    null.resid <- residuals(mod)
    
    boot.dat <- dat           # copy the data set into a bootstrap resampling dataset
    wald.store <- c()         # create a container for storing the test statistics    
    
    if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
    for(i in 1:boot.reps){
      
      if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
      
      weight <- c(1, -1)[rbinom(G, size=1, prob=0.5) + 1][match(clust, unique(clust))]       # assign wild bootstrap weights
      pseudo.resid <- null.resid*weight                                # create pseudo-residuals using weights
      pseudo.dv <- predict(mod)+ pseudo.resid                     # create pseudo-observations using pseudo-residuals
      boot.dat[,dv] <- pseudo.dv                                       # create a bootstrap replicate data
      
      boot.mod <- glm(form, data = boot.dat, family = mod$family)      # run a model on the bootstrap replicate data
      
      se.boot <- cl(boot.dat, boot.mod, clust)[1,2]                  # retrieve the bootstrap clustered SE
      beta.boot <- coefficients(boot.mod)[1]                         # store the bootstrap beta coefficient
      wald.store[i] <- (beta.boot - beta.mod[1]) / se.boot                             # store the bootstrap test statistic
      
    }
    if(prog.bar==TRUE){close(pb)}
    
    p.store <- c( 1 - ( sum( abs(w[1]) > abs(wald.store) ) / boot.reps ), p.store)    # calculate the wild bootstrap p-value
    w.store <- cbind(wald.store, w.store)
    
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
   }
  
  return(invisible(out))
   
}

