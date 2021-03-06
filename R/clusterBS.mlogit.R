#' Pairs Cluster Bootstrapped p-Values For mlogit
#'
#' This software estimates p-values using pairs cluster bootstrapped t-statistics for multinomial logit models (Cameron, Gelbach, and Miller 2008). The data set is repeatedly re-sampled by cluster, a model is estimated, and inference is based on the sampling distribution of the pivotal (t) statistic. 
#'
#' @param mod A model estimated using \code{mlogit}.
#' @param dat The data set used to estimate \code{mod}.
#' @param cluster A formula of the clustering variable.
#' @param ci.level What confidence level should CIs reflect?
#' @param boot.reps The number of bootstrap samples to draw.
#' @param cluster.se Use clustered standard errors (= TRUE) or ordinary SEs (= FALSE) for bootstrap replicates.
#' @param report Should a table of results be printed to the console?
#' @param prog.bar Show a progress bar of the bootstrap (= TRUE) or not (= FALSE).
#' @param output.replicates Should the cluster bootstrap coefficient replicates be output (= TRUE) or not (= FALSE)?
#' @param seed Random number seed for replicability (default is NULL).
#'
#' @return A list with the elements
#' \item{p.values}{A matrix of the estimated p-values.}
#' \item{ci}{A matrix of confidence intervals.}
#' @author Justin Esarey
#' @note Code to estimate GLM clustered standard errors by Mahmood Arai: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/, although modified slightly to work for \code{mlogit} models. Cluster SE degrees of freedom correction = (M/(M-1)) with M = the number of clusters.
#' @examples
#' \dontrun{
#' 
#' #######################################
#' # example one: train ticket selection
#' #######################################
#' require(mlogit)
#' data("Train", package="mlogit")
#' Train$choiceid <- 1:nrow(Train)
#' 
#' Tr <- dfidx(Train, shape = "wide", varying = 4:11, sep = "_", 
#'           choice = "choice", idx = list(c("choiceid", "id")), 
#'           idnames = c(NA, "alt"))
#' Tr$price <- Tr$price/100 * 2.20371
#' Tr$time <- Tr$time/60
#' 
#' ml.Train <- mlogit(choice ~ price + time + change + comfort | -1, Tr)
#' 
#' # compute pairs cluster bootstrapped p-values
#' # note: few reps to speed up example
#' cluster.bs.tr <- cluster.bs.mlogit(ml.Train, Tr, ~ id, boot.reps=100)
#' 
#' 
#' ##################################################################
#' # example two: predict type of heating system installed in house
#' ##################################################################
#' require(mlogit)
#' data("Heating", package = "mlogit")
#' H <- Heating
#' H$region <- as.numeric(H$region)
#' H.ml <- dfidx(H, shape="wide", choice="depvar", varying=c(3:12),
#'          idx = list(c("idcase", "region")))
#' m <- mlogit(depvar~ic+oc, H.ml)
#' 
#' # compute pairs cluster bootstrapped p-values
#' cluster.bs.h <- cluster.bs.mlogit(m, H.ml, ~ region, boot.reps=1000)
#' 
#' }
#' @rdname cluster.bs.mlogit
#' @import stats
#' @importFrom dfidx dfidx idx idx_name
#' @importFrom utils write.table
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom lmtest coeftest
#' @importFrom sandwich estfun
#' @importFrom sandwich sandwich
#' @importFrom mlogit mlogit hmftest mFormula is.mFormula mlogit.optim cov.mlogit cor.mlogit rpar scoretest med rg stdev qrpar prpar drpar
#' @references Esarey, Justin, and Andrew Menger. 2017. "Practical and Effective Approaches to Dealing with Clustered Data." \emph{Political Science Research and Methods} forthcoming: 1-35. <URL:http://jee3.web.rice.edu/cluster-paper.pdf>.
#' @references Cameron, A. Colin, Jonah B. Gelbach, and Douglas L. Miller. 2008. "Bootstrap-Based Improvements for Inference with Clustered Errors." \emph{The Review of Economics and Statistics} 90(3): 414-427. <DOI:10.1162/rest.90.3.414>.
#' @export

cluster.bs.mlogit<-function(mod, dat, cluster, ci.level = 0.95, boot.reps = 1000, cluster.se = TRUE, 
                            report = TRUE, prog.bar = TRUE, output.replicates = FALSE,
                            seed = NULL){
  
  # compensate for bizarre R formula updating bug
  # thanks to Jason Thorpe for reporting!
  form.old <- update(mod$formula, 1 ~ 1 )
  while(form.old != mod$formula){
    form.old <- mod$formula
    invisible(mod <- update(mod, new = .~.))
  }
  
  if(is.null(seed)==F){                                               # if user supplies a seed, set it
    
    tryCatch(set.seed(seed),
             error = function(e){return("seed must be a valid integer")}, 
             warning = function(w){return(NA)}) 
    
  }
    
  form <- mod$formula                                                    # what is the formula of this model?  
  variables <- all.vars(form)                                            # what variables are in this model?
  used.idx <- which(rownames(dat) %in% rownames(mod$mod))                # what observations are used?
  dat <- dat[used.idx,]                                                  # keep only active observations
  ind.variables <- names(coefficients(mod))                              # what independent variables are in this model?
  "%w/o%" <- function(x, y) x[!x %in% y]                                 # create a without function (see ?match)
  dv <- variables %w/o% all.vars(update(form, 1 ~ .))                    # what is the dependent variable?
  
  # obtain the clustering variable
  clust.name <- all.vars(cluster)                                        # name of the cluster variable
  dat.rs <- as.data.frame(subset(idx(dat), select = clust.name ))        # select cluster variable from data set
  dat.rs$id.zz <- idx(dat, n=1)                                          # choice index
  dat.rs$ti.zz <- idx(dat, n=2)                                          # alternative index
  clust <- reshape(dat.rs, timevar="ti.zz",                              # reshape long to wide, store as clust
        idvar=c("id.zz", clust.name), direction="wide")[[clust.name]]  
  if(sum(is.na(clust)>0)){stop("missing cluster indices")}               # check for missing cluster indices
  G <- length(unique(clust))                                             # how many clusters are in this model?
  
  
  # load in a function to create clustered standard errors for mlogit models
  # initial code by Mahmood Arai: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/
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
    dfc <- (M/(M-1))
    uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
    vcovCL <- dfc*sandwich(fm, meat.=crossprod(uj)/N)
    coeftest(fm, vcovCL) 
  }
  
   if(cluster.se == T){
     
     se.clust <- cl.mlogit(mod, clust)[ind.variables,2]             # retrieve the clustered SEs
     beta.mod <- coefficients(mod)[ind.variables]                   # retrieve the estimated coefficients
     w <- beta.mod / se.clust                                       # calculate the t-test statistic
     
   }else{

     se.beta <- summary(mod)$CoefTable[ind.variables,2]             # retrieve the vanilla SEs
     beta.mod <- coefficients(mod)[ind.variables]                   # retrieve the estimated coefficients
     w <- beta.mod / se.beta                                        # calculate the t-test statistic

   }
  
  # keep track of the beta bootstrap replicates for possible output
  rep.store <- matrix(data=NA, nrow=boot.reps, ncol=length(beta.mod))
  colnames(rep.store) <- ind.variables
  
  w.store <- matrix(data=NA, nrow=boot.reps, ncol=length(ind.variables))      # store bootstrapped test statistics
  
  if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
  for(i in 1:boot.reps){
    
    if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
    
    boot.sel <- sample(1:G, size=G, replace=T)                                 # randomly select clusters
    
    # pick the observations corresponding to the randomly selected clusters
    boot.ind <- c()                                                            # where the selected obs will be stored
    boot.clust <- c()                                                          # create + store a new cluster index for the bootstrap data
    for(k in 1:G){

      obs.sel <- which(dat[[clust.name]] ==                                    # which observations are in the sampled cluster?
                  unlist(unique(clust))[boot.sel[k]])
      boot.ind <- c(boot.ind, obs.sel)                                         # append the selected obs index to existing index
      boot.clust <- c(boot.clust, rep(k, length(obs.sel)))                     # store the new bootstrap cluster index
      
    }
    
    boot.dat <- dat[boot.ind,]                                                 # create the bootstrapped data
    eval(parse(text=paste("boot.dat$idx[[\"", clust.name,                      # add boot-specific cluster variable
                          "\"]]=boot.clust", sep="")))      
    alt.num <- length(unique(idx(dat, n=2)))                                   # how many alternatives are there?
    ch.num <- round( dim(boot.dat)[1] / alt.num )                              # how many choices are there? (rounding for imprecision)
    eval(parse(text=paste("boot.dat$idx[[\"",idx_name(dat, n=1),               # create new choice index
                          "\"]]=rep(1:ch.num, each=alt.num)", sep="")))    
#    eval(parse(text=paste("boot.dat$idx[[\"", idx_name(dat, n=2),             # create new alternative index
#                          "\"]]=value=rep(1:alt.num, ch.num)", sep="")))
    rownames(boot.dat) <- NULL                                                 # purge old (duplicated) row names
    
    
    boot.mod.call <- mod$call                                                  # get original model call
    boot.mod.call[[3]] <- quote(boot.dat)                                      # modify call for bootstrap data set
    
    boot.mod <- suppressWarnings(tryCatch(eval(boot.mod.call),                 # estimate model on bootstrap dataset
                            error = function(e){return(NULL)}))    
    
    fail <- is.null(boot.mod)                                                  # determine whether the mlogit process created an error
    
    # obtain the bootstrap clustering variable
    boot.dat.rs <- as.data.frame(
            subset(idx(boot.dat), select = clust.name ))                       # select cluster variable from BS data set
    boot.dat.rs$id.zz <- idx(boot.dat, n=1)                                    # choice index
    boot.dat.rs$ti.zz <- idx(boot.dat, n=2)                                    # alternative index
    boot.clust.n <- reshape(boot.dat.rs, timevar="ti.zz",                      # reshape long to wide, store as clust
              idvar=c("id.zz", clust.name),
              direction="wide")[[clust.name]]  
    
    if(fail==0){                                                     # proceed if the mlogit model was not in error

      if(cluster.se == T){
        
        se.boot <- tryCatch(cl.mlogit(boot.mod, boot.clust.n)[ind.variables,2],
                   error = function(e){return(NA)}, 
                   warning = function(w){return(NA)})                              # retrieve the bootstrap clustered SE
        beta.boot <- tryCatch(coefficients(boot.mod)[ind.variables],
                     error = function(e){return(NA)}, 
                     warning = function(w){return(NA)})                            # store the bootstrap beta coefficient
        w.store[i,] <- (beta.boot-beta.mod) / se.boot                              # store the bootstrap test statistic
        
        rep.store[i,] <- beta.boot                                                 # store the bootstrap beta for output
        
        
      }else{
        
        se.boot <- tryCatch(summary(boot.mod)$CoefTable[ind.variables,2],
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
  if(output.replicates == TRUE){out.list[["replicates"]] <- rep.store}
  return(invisible(out.list))
  
}

