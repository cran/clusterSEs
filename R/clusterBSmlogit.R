#' Pairs Cluster Bootstrapped p-Values For mlogit
#'
#' This software estimates p-values using pairs cluster bootstrapped t-statistics for multinomial logit models (Cameron, Gelbach, and Miller 2008). The data set is repeatedly re-sampled by cluster, a model is estimated, and inference is based on the sampling distribution of the pivotal (t) statistic. 
#'
#' @param mod A model estimated using \code{mlogit}.
#' @param dat The data set used to estimate \code{mod}, but in standard (not \code{mlogit.data}) form..
#' @param cluster A formula of the clustering variable.
#' @param ci.level What confidence level should CIs reflect?
#' @param boot.reps The number of bootstrap samples to draw.
#' @param stratify Sample clusters only (= FALSE) or clusters and observations by cluster (= TRUE).
#' @param cluster.se Use clustered standard errors (= TRUE) or ordinary SEs (= FALSE) for bootstrap replicates.
#' @param report Should a table of results be printed to the console?
#' @param prog.bar Show a progress bar of the bootstrap (= TRUE) or not (= FALSE).
#'
#' @return A list with the elements
#' \item{p.values}{A matrix of the estimated p-values.}
#' \item{ci}{A matrix of confidence intervals.}
#' @author Justin Esarey
#' @note Code to estimate GLM clustered standard errors by Mahmood Ara: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/, although modified slightly to work for \code{mlogit} models.
#' @examples
#' \dontrun{
#' # predict method of hospital admission
#' require(VGAMdata)
#' data(vtinpat)
#' 
#' # to save time, take a sample of this data
#' set.seed(32149)
#' idx <- sample(1:dim(vtinpat)[1], 5000, replace=F)
#' vtinpat <- vtinpat[idx,]
#' 
#' vtinpat$hos.num <- as.numeric(vtinpat$hospital)
#' vtinpat$age <- as.numeric(vtinpat$age.group)
#' vtinpat.mlogit <- mlogit.data(vtinpat, choice = "admit", shape="wide")
#' vt.mod <- mlogit(admit ~ 0 | age + sex, data = vtinpat.mlogit)
#' summary(vt.mod)
#' 
#' # compute cluster bootstrapped p-values (takes a while)
#' clust.p <- cluster.bs.mlogit(vt.mod, dat=vtinpat, cluster = ~ hos.num, report=TRUE)
#' }
#' @rdname cluster.bs.mlogit
#' @importFrom lmtest coeftest
#' @importFrom sandwich estfun
#' @importFrom sandwich sandwich
#' @import mlogit
#' @references Cameron, A. Colin, Jonah B. Gelbach, and Douglas L. Miller. 2008. "Bootstrap-Based Improvements for Inference with Clustered Errors." \emph{The Review of Economics and Statistics} 90(3): 414-427.
#' @export

cluster.bs.mlogit<-function(mod, dat, cluster, ci.level = 0.95, boot.reps = 1000, stratify = FALSE, cluster.se = TRUE, report = TRUE, prog.bar = TRUE){
  
  form <- mod$formula                                                    # what is the formula of this model?  
  variables <- all.vars(form)                                            # what variables are in this model?
  dat.t <- subset(dat, select = variables)                               # keep only relevant variables
  dat.t$clust <- subset(dat, select = all.vars(cluster))                 # add the cluster variable into dat.t (for NA omission)
  dat <- na.omit(dat.t)                                                  # drop the NAs
  clust <- as.vector(unlist(dat$clust))                                  # reintegrate cluster variable w/o NA obs
  G<-length(unique(clust))                                               # how many clusters are in this model?
  ind.variables <- names(coefficients(mod))                              # what independent variables are in this model?
  
  
  # load in a function to create clustered standard errors for mlogit models
  # initial code by Mahmood Ara: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/
  # slightly modified for mlogit models by Justin Esarey on 3/3/2015
  
  cl.mlogit   <- function(fm, cluster){
    
    # fm: a fitted mlogit model
    # cluster: a data vector with the cluster
    #          identity of each observation in fm
    
    #require(sandwich, quietly = TRUE)
    #require(lmtest, quietly = TRUE)
    M <- length(unique(cluster))
    N <- length(cluster)
    K <- length(coefficients(fm))
    dfc <- (M/(M-1))*((N-1)/(N-K))
    uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
    vcovCL <- dfc*sandwich(fm, meat.=crossprod(uj)/N)
    coeftest(fm, vcovCL) 
  }
  
   if(cluster.se == T){
     
     se.clust <- cl.mlogit(mod, clust)[ind.variables,2]               # retrieve the clustered SEs
     beta.mod <- coefficients(mod)[ind.variables]                   # retrieve the estimated coefficients
     w <- beta.mod / se.clust                                       # calculate the t-test statistic
     
   }else{

     se.beta <- summary(mod)$coefficients[ind.variables,2]          # retrieve the vanilla SEs
     beta.mod <- coefficients(mod)[ind.variables]                   # retrieve the estimated coefficients
     w <- beta.mod / se.beta                                        # calculate the t-test statistic

   }
  
  w.store <- matrix(data=NA, nrow=boot.reps, ncol=length(ind.variables))    # store bootstrapped test statistics
  
  if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
  for(i in 1:boot.reps){
    
    if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
    
    boot.sel <- sample(1:G, size=G, replace=T)                         # randomly select clusters
    
    # pick the observations corresponding to the randomly selected clusters
    boot.ind <- c()                                                    # where the selected obs will be stored
    boot.clust <- c()                                                  # create + store a new cluster index for the bootstrap data
    
    for(k in 1:G){

      obs.sel <- which(clust == unique(clust)[boot.sel[k]])                           # which observations are in the sampled cluster?
      if(stratify==T){
        
        obs.samp <- sample(obs.sel, size = length(obs.sel), replace=T)    # sample randomly from the selected cluster
        boot.ind <- c(boot.ind, obs.samp)                                 # append the selected obs index to existing index
        
      }else{
        
        boot.ind <- c(boot.ind, obs.sel)                                  # append the selected obs index to existing index
                  
      }
      boot.clust <- c(boot.clust, rep(k, length(obs.sel)))             # store the new bootstrap cluster index
      
    }
    
    boot.dat <- dat[boot.ind,]                                                         # create the bootstrapped data
    boot.dat <- mlogit.data(boot.dat, choice = "admit", shape="wide")

    # run a model on the bootstrap replicate data
    boot.mod <- suppressWarnings(tryCatch(mlogit(form, data = boot.dat), 
                error = function(e){return(NA)}))                                    

    fail <- max(is.na(boot.mod))                                     # determine whether the mlogit process created an error
    
    if(fail==0){                                                     # proceed if the mlogit model was not in error

      if(cluster.se == T){
        
        se.boot <- tryCatch(cl.mlogit(boot.mod, boot.clust)[ind.variables,2],
                   error = function(e){return(NA)}, 
                   warning = function(w){return(NA)})                              # retrieve the bootstrap clustered SE
        beta.boot <- tryCatch(coefficients(boot.mod)[ind.variables],
                     error = function(e){return(NA)}, 
                     warning = function(w){return(NA)})                            # store the bootstrap beta coefficient
        w.store[i,] <- (beta.boot-beta.mod) / se.boot                              # store the bootstrap test statistic
        
      }else{
        
        se.boot <- tryCatch(summary(boot.mod)$coefficients[ind.variables,2],
                   error = function(e){return(NA)}, 
                   warning = function(w){return(NA)})                               # retrieve the bootstrap vanilla SE
        beta.boot <- tryCatch(coefficients(boot.mod)[ind.variables],
                     error = function(e){return(NA)}, 
                     warning = function(w){return(NA)})                             # retrieve the bootstrap beta coefficient
        w.store[i,] <- (beta.boot-beta.mod) / se.boot                               # calculate the t-test statistic
                
      }
    
    }else{
      w.store[i,] <- NA                                                  # if model didn't converge, store NA as a result 
    }
  
  }
  if(prog.bar==TRUE){close(pb)}
  
  num.fail <- length(attr(na.omit(w.store), "na.action"))         # count the number of times something went wrong
  w.store <- na.omit(w.store)                                     # drop the erroneous bootstrap replicates
  
  
  comp.fun<-function(vec2, vec1){as.numeric(vec1>vec2)}                              # a simple function comparing v1 to v2
  p.store.s <- t(apply(X = abs(w.store), FUN=comp.fun, MARGIN = 1, vec1 = abs(w)))   # compare the BS test stats to orig. result
  p.store <- 1 - ( colSums(p.store.s) / dim(w.store)[1] )                                       # calculate the cluster bootstrap p-value

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
  out.p <- rbind(c("variable name", "clustered bootstrap p-value"), out.p)
  

  printmat <- function(m){
    write.table(format(m, justify="right"), row.names=F, col.names=F, quote=F, sep = "   ")
  }
  
  if(report==T){
    
    if(num.fail!=0){
    cat("\n", "\n", "\n", "****", "Warning: ", num.fail, " out of ", boot.reps, "bootstrap replicate models failed to estimate.", "****", "\n")
    }
    
    cat("\n", "Cluster Bootstrap p-values: ", "\n", "\n")
    printmat(out.p)

    cat("\n", "Confidence Intervals (derived from bootstrapped t-statistics): ", "\n", "\n")
    printmat(print.ci)
    
    
  }
  
  out.list<-list()
  out.list[["p.values"]]<-out
  out.list[["ci"]] <- out.ci
  return(invisible(out.list))
  
}

