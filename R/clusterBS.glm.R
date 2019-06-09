#' Pairs Cluster Bootstrapped p-Values For GLM
#'
#' This software estimates p-values using pairs cluster bootstrapped t-statistics for GLM models (Cameron, Gelbach, and Miller 2008). The data set is repeatedly re-sampled by cluster, a model is estimated, and inference is based on the sampling distribution of the pivotal (t) statistic. 
#'
#' @param mod A model estimated using \code{glm}.
#' @param dat The data set used to estimate \code{mod}.
#' @param cluster A formula of the clustering variable.
#' @param ci.level What confidence level should CIs reflect?
#' @param boot.reps The number of bootstrap samples to draw.
#' @param stratify Sample clusters only (= FALSE) or clusters and observations by cluster (= TRUE).
#' @param cluster.se Use clustered standard errors (= TRUE) or ordinary SEs (= FALSE) for bootstrap replicates.
#' @param report Should a table of results be printed to the console?
#' @param prog.bar Show a progress bar of the bootstrap (= TRUE) or not (= FALSE).
#' @param output.replicates Should the cluster bootstrap coefficient replicates be output (= TRUE) or not (= FALSE)?
#' @param seed Random number seed for replicability (default is NULL).
#'
#' @return A list with the elements
#' \item{p.values}{A matrix of the estimated p-values.}
#' \item{ci}{A matrix of confidence intervals.}
#' \item{replicates}{Optional: A matrix of the coefficient estimates from each cluster bootstrap replicate.}
#' @author Justin Esarey
#' @note Code to estimate GLM clustered standard errors by Mahmood Arai: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/. Cluster SE degrees of freedom correction = (M/(M-1)) with M = the number of clusters.
#' @examples
#' \dontrun{
#' 
#' ##################################################################
#' # example one: predict whether respondent has a university degree
#' ##################################################################
#' require(effects)
#' data(WVS)
#' logit.model <- glm(degree ~ religion + gender + age, data=WVS, family=binomial(link="logit"))
#' summary(logit.model)
#' 
#' # compute pairs cluster bootstrapped p-values
#' clust.bs.p <- cluster.bs.glm(logit.model, WVS, ~ country, report = T)
#' 
#'   
#' ######################################
#' # example two: predict chicken weight
#' ######################################
#' rm(list=ls())
#' data(ChickWeight)
#' 
#' dum <- model.matrix(~ ChickWeight$Diet)
#' ChickWeight$Diet2 <- as.numeric(dum[,2])
#' ChickWeight$Diet3 <- as.numeric(dum[,3])
#' ChickWeight$Diet4 <- as.numeric(dum[,4])
#' 
#' weight.mod2 <- glm(formula = weight~Diet2+Diet3+Diet4+log(Time+1),data=ChickWeight)
#' 
#' # compute pairs cluster bootstrapped p-values
#' clust.bs.w <- cluster.bs.glm(weight.mod2, ChickWeight, ~ Chick, report = T)
#' 
#' 
#' ###################################################################
#' # example three: murder rate by U.S. state, with interaction term
#' ###################################################################
#' rm(list=ls())
#' require(datasets)
#' 
#' state.x77.dat <- data.frame(state.x77)
#' state.x77.dat$Region <- state.region
#' state.x77.dat$IncomeXHS <- state.x77.dat$Income * state.x77.dat$HS.Grad
#' income.mod <- glm( Murder ~ Income + HS.Grad + IncomeXHS , data=state.x77.dat)
#' 
#' # compute pairs cluster bootstrapped p-values
#' clust.bs.inc <- cluster.bs.glm(income.mod, state.x77.dat, ~ Region, 
#'                                report = T, output.replicates=T, boot.reps=10000)
#' 
#' # compute effect of income on murder rate, by percentage of HS graduates
#' # using conventional standard errors
#' HS.grad.vec <- seq(from=38, to=67, by=1)
#' me.income <- coefficients(income.mod)[2] + coefficients(income.mod)[4]*HS.grad.vec
#' plot(me.income ~ HS.grad.vec, type="l", ylim=c(-0.0125, 0.0125), 
#'      xlab="% HS graduates", ylab="ME of income on murder rate")
#' se.income <- sqrt( vcov(income.mod)[2,2] + vcov(income.mod)[4,4]*(HS.grad.vec)^2 +
#'                    2*vcov(income.mod)[2,4]*HS.grad.vec )
#' ci.h <- me.income + qt(0.975, lower.tail=T, df=46) * se.income
#' ci.l <- me.income - qt(0.975, lower.tail=T, df=46) * se.income
#' lines(ci.h ~ HS.grad.vec, lty=2)
#' lines(ci.l ~ HS.grad.vec, lty=2)
#' 
#' # use pairs cluster bootstrap to compute CIs, including bootstrap bias-correction factor
#' # including bootstrap bias correction factor
#' # cluster on Region
#' ################################################
#' # marginal effect replicates =
#' me.boot <- matrix(data = clust.bs.inc$replicates[,2], nrow=10000, ncol=30, byrow=F) +
#'            as.matrix(clust.bs.inc$replicates[,4]) %*% t(HS.grad.vec)
#' # compute bias-corrected MEs
#' me.income.bias.cor <- 2*me.income - apply(X=me.boot, FUN=mean, MARGIN=2)
#' # adjust bootstrap replicates for bias
#' me.boot.bias.cor <- me.boot + matrix(data = 2*(me.income - 
#'                                      apply(X=me.boot, FUN=mean, MARGIN=2)),
#'                                      ncol=30, nrow=10000, byrow=T)
#' # compute pairs cluster bootstrap 95% CIs, including bias correction
#' me.boot.plot <- apply(X = me.boot.bias.cor, FUN=quantile, MARGIN=2, probs=c(0.025, 0.975))
#' # plot bootstrap bias-corrected marginal effects
#' lines(me.income.bias.cor ~ HS.grad.vec, lwd=2)
#' # plot 95% Cis
#' # a little lowess smoothing applied to compensate for discontinuities 
#' # arising from shifting between replicates
#' lines(lowess(me.boot.plot[1,] ~ HS.grad.vec), lwd=2, lty=2)
#' lines(lowess(me.boot.plot[2,] ~ HS.grad.vec), lwd=2, lty=2)
#' 
#' # finishing touches to plot
#' legend(lty=c(1,2,1,2), lwd=c(1,1,2,2), "topleft", 
#'        legend=c("Model Marginal Effect", "Conventional 95% CI", 
#'                 "BS Bias-Corrected Marginal Effect", "Cluster Bootstrap 95% CI"))
#' 
#' }
#' @rdname cluster.bs.glm
#' @import stats
#' @importFrom utils write.table
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom lmtest coeftest
#' @importFrom sandwich estfun
#' @importFrom sandwich sandwich
#' @references Esarey, Justin, and Andrew Menger. 2017. "Practical and Effective Approaches to Dealing with Clustered Data." \emph{Political Science Research and Methods} forthcoming: 1-35. <URL:http://jee3.web.rice.edu/cluster-paper.pdf>.
#' @references Cameron, A. Colin, Jonah B. Gelbach, and Douglas L. Miller. 2008. "Bootstrap-Based Improvements for Inference with Clustered Errors." \emph{The Review of Economics and Statistics} 90(3): 414-427. <DOI:10.1162/rest.90.3.414>.
#' @export

cluster.bs.glm<-function(mod, dat, cluster, ci.level = 0.95, boot.reps = 1000, stratify = FALSE, 
                         cluster.se = TRUE, report = TRUE, prog.bar = TRUE, output.replicates = FALSE,
                         seed = NULL){
  
  if(is.null(seed)==F){                                               # if user supplies a seed, set it
    
    tryCatch(set.seed(seed),
             error = function(e){return("seed must be a valid integer")}, 
             warning = function(w){return(NA)}) 
    
  }
  
  form <- mod$formula                                                 # what is the formula of this model?  
  variables <- all.vars(form)                                         # what variables are in this model?
  clust.name <- all.vars(cluster)                                     # what is the name of the clustering variable?
  used.idx <- which(rownames(dat) %in% rownames(mod$model))           # what were the actively used observations in the model?
  dat <- dat[used.idx,]                                               # keep only active observations (drop the missing)
  clust <- as.vector(unlist(dat[[clust.name]]))                       # store cluster index in convenient vector
  G<-length(unique(clust))                                            # how many clusters are in this model?
  ind.variables.full <- names(coefficients(mod))                      # what independent variables are in this model?
  ind.variables <- rownames(summary(mod)$coefficients)                # what non-dropped independent variables in this model?
  
  
  # load in a function to create clustered standard errors
  # by Mahmood Arai: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/
  cl   <- function(dat, fm, cluster){
    #require(sandwich, quietly = TRUE)
    #require(lmtest, quietly = TRUE)
    M <- length(unique(cluster))
    N <- length(cluster)
    K <- fm$rank
    dfc <- (M/(M-1))
    uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
    vcovCL <- dfc*sandwich(fm, meat.=crossprod(uj)/N)
    coeftest(fm, vcovCL) }
  
   if(cluster.se == T){
     
     se.clust <- cl(dat, mod, clust)[ind.variables,2]               # retrieve the clustered SEs
     beta.mod <- coefficients(mod)[ind.variables]                   # retrieve the estimated coefficients
     w <- beta.mod / se.clust                                       # calculate the t-test statistic
     
   }else{

     se.beta <- summary(mod)$coefficients[ind.variables,2]          # retrieve the vanilla SEs
     beta.mod <- coefficients(mod)[ind.variables]                   # retrieve the estimated coefficients
     w <- beta.mod / se.beta                                        # calculate the t-test statistic

   }
  
  w.store <- matrix(data=NA, nrow=boot.reps, ncol=length(ind.variables))  # store bootstrapped test statistics
  
  # keep track of the beta bootstrap replicates for possible output
  rep.store <- matrix(data=NA, nrow=boot.reps, ncol=length(beta.mod))
  colnames(rep.store) <- ind.variables
  
  if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
  for(i in 1:boot.reps){
    
    if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
    
    boot.sel <- sample(1:G, size=G, replace=T)                            # randomly select clusters
    
    # pick the observations corresponding to the randomly selected clusters
    boot.ind <- c()                                                       # where the selected obs will be stored
    boot.clust <- c()                                                     # create new cluster index for the bootstrap data
    
    for(k in 1:G){

      obs.sel <- which(clust == unique(clust)[boot.sel[k]])               # which observations are in the sampled cluster?
      if(stratify==T){
        
        obs.samp <- sample(obs.sel, size = length(obs.sel), replace=T)    # sample randomly from the selected cluster
        boot.ind <- c(boot.ind, obs.samp)                                 # append the selected obs index to existing index
        
      }else{
        
        boot.ind <- c(boot.ind, obs.sel)                                  # append the selected obs index to existing index
                  
      }
      boot.clust <- c(boot.clust, rep(k, length(obs.sel)))                # store the new bootstrap cluster index
      
    }
    
    boot.dat <- dat[boot.ind,]                                            # create the bootstrapped data

    # run a model on the bootstrap replicate data
    boot.mod <- suppressWarnings(tryCatch(glm(form, data = boot.dat, family = mod$family), 
                error = function(e){return(NULL)}))                                    

    if(is.null(boot.mod) == FALSE ){
      if(boot.mod$converged == 0){boot.mod <- NULL}                    # judge GLM as failure if convergence not achieved
    }
    fail <- is.null(boot.mod)                                          # determine whether the GLM process created an error
    
    if(fail==0){                                                     # proceed if the GLM model was not in error

      if(cluster.se == T){
        
        se.boot <- tryCatch(cl(boot.dat, boot.mod, boot.clust)[ind.variables,2],
                   error = function(e){return(NA)}, 
                   warning = function(w){return(NA)})                              # retrieve the bootstrap clustered SE
        beta.boot <- tryCatch(coefficients(boot.mod)[ind.variables],
                     error = function(e){return(NA)}, 
                     warning = function(w){return(NA)})                            # store the bootstrap beta coefficient
        w.store[i,] <- (beta.boot-beta.mod) / se.boot                              # store the bootstrap test statistic
        
        rep.store[i,] <- beta.boot                                                 # store the bootstrap beta for output
        
      }else{
        
        se.boot <- tryCatch(summary(boot.mod)$coefficients[ind.variables,2],
                   error = function(e){return(NA)}, 
                   warning = function(w){return(NA)})                               # retrieve the bootstrap vanilla SE
        beta.boot <- tryCatch(coefficients(boot.mod)[ind.variables],
                     error = function(e){return(NA)}, 
                     warning = function(w){return(NA)})                             # retrieve the bootstrap beta coefficient
        w.store[i,] <- (beta.boot-beta.mod) / se.boot                               # calculate the t-test statistic
        
        rep.store[i,] <- beta.boot                                                  # store the bootstrap beta for output
                
      }
    
    }else{
      w.store[i,] <- NA                                                  # if model didn't converge, store NA as a result 
      rep.store[i,] <- NA
    }
  
  }
  if(prog.bar==TRUE){close(pb)}
  
  num.fail <- length(attr(na.omit(w.store), "na.action"))         # count the number of times something went wrong
  w.store <- na.omit(w.store)                                     # drop the erroneous bootstrap replicates
  
  
  comp.fun<-function(vec2, vec1){as.numeric(vec1>vec2)}                              # a simple function comparing v1 to v2
  p.store.s <- t(apply(X = abs(w.store), FUN=comp.fun, MARGIN = 1, vec1 = abs(w)))   # compare the BS test stats to orig. result
  p.store <- 1 - ( colSums(p.store.s) / dim(w.store)[1] )                            # calculate the cluster bootstrap p-value

  # compute critical t-statistics for CIs
  crit.t <- apply(X=abs(w.store), MARGIN=2, FUN=quantile, probs=ci.level )
  if(cluster.se == TRUE){
    ci.lo <- beta.mod - crit.t*se.clust
    ci.hi <- beta.mod + crit.t*se.clust
  }else{
    ci.lo <- beta.mod - crit.t*se.beta
    ci.hi <- beta.mod + crit.t*se.beta
  }
  
  print.ci <- cbind(ind.variables, ci.lo, ci.hi)
  print.ci <- rbind(c("variable name", "CI lower", "CI higher"), print.ci)
  
  out.ci <- cbind(ci.lo, ci.hi)
  rownames(out.ci) <- ind.variables
  colnames(out.ci) <- c("CI lower", "CI higher")
  
  out <- matrix(p.store, ncol=1)
  colnames(out) <- c("clustered bootstrap p-value")
  rownames(out) <- ind.variables
  out.p <- cbind(ind.variables, out)
  out.p <- rbind(c("variable name", "cluster bootstrap p-value"), out.p)
  
  
  printmat <- function(m){
    write.table(format(m, justify="right"), row.names=F, col.names=F, quote=F, sep = "   ")
  }
  
  if(report==T){
    
    if(num.fail!=0){
    cat("\n", "\n", "\n", "****", "Warning: ", num.fail, " out of ", boot.reps, " bootstrap replicate models failed to estimate.", "****", "\n", sep="")
    }
    
    cat("\n", "Cluster Bootstrap p-values: ", "\n", "\n")
    printmat(out.p)

    cat("\n", "Confidence Intervals (derived from bootstrapped t-statistics): ", "\n", "\n")
    printmat(print.ci)

    if(length(ind.variables) < length(ind.variables.full)){
    cat("\n", "\n", "****", "Note: ", length(ind.variables.full) - length(ind.variables), " variables were unidentified in the model and are not reported.", "****", "\n", sep="")
    cat("Variables not reported:", "\n", sep="")
    cat(ind.variables.full[!ind.variables.full %in% ind.variables], sep=", ")
    cat("\n", "\n")
    }
    
  }
  
  out.list<-list()
  out.list[["p.values"]]<-out
  out.list[["ci"]] <- out.ci
  if(output.replicates == TRUE){out.list[["replicates"]] <- rep.store}
  return(invisible(out.list))
  
}

