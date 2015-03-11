#' Cluster-Adjusted Confidence Intervals and p-Values for GLM
#'
#' Computes p-values and confidence intervals for GLM models based on cluster-specific model estimation (Ibragimov and Muller 2010). A separate model is estimated in each cluster, and then p-values and confidence intervals are computed based on a t/normal distribution of the cluster-specific estimates.
#'
#' @param mod A model estimated using \code{glm}.
#' @param dat The data set used to estimate \code{mod}.
#' @param cluster A formula of the clustering variable.
#' @param ci.level What confidence level should CIs reflect?
#' @param report Should a table of results be printed to the console?
#' @param drop Should clusters within which a model cannot be estimated be dropped?
#'
#' @return A list with the elements
#' \item{p.values}{A matrix of the estimated p-values.}
#' \item{ci}{A matrix of confidence intervals.}
#' @author Justin Esarey
#' @note Confidence intervals are centered on the cluster averaged estimate, which can diverge from original model estimates if clusters have different numbers of observations. Consequently, confidence intervals may not be centered on original model estimates.
#' @examples
#' \dontrun{
#' # predict whether respondent has a university degree
#' require(effects)
#' data(WVS)
#' logit.model <- glm(degree ~ religion + gender + age, data=WVS, family=binomial(link="logit"))
#' summary(logit.model)
#' 
#' # compute cluster-adjusted p-values
#' clust.p <- cluster.im(logit.model, WVS, ~ country, ci.level = 0.95, report = T, drop = FALSE) 
#' }
#' @rdname cluster.im
#' @references Ibragimov, Rustam, and Ulrich K. Muller. 2010. "t-Statistic Based Correlation and Heterogeneity Robust Inference." \emph{Journal of Business & Economic Statistics} 28(4): 453-468. 
#' @export

cluster.im<-function(mod, dat, cluster, ci.level = 0.95, report = TRUE, drop = FALSE){
  
  form <- mod$formula                                     # what is the formula of this model?
  
  variables <- all.vars(form)                             # what variables are in this model?
  dat.t <- subset(dat, select = variables)                # keep only relevant variables
  dat.t$clust <- subset(dat, select = all.vars(cluster))  # add the cluster variable into dat.t (for NA omission)
  dat <- na.omit(dat.t)                                   # drop the NAs
  clust <- as.vector(unlist(dat$clust))                   # reintegrate cluster variable w/o NA obs
  G<-length(unique(clust))                                # how many clusters are in this model?
  ind.variables <- names(coefficients(mod))               # what independent variables are in this model?
 
  
  b.clust <- matrix(data = NA, nrow = G, ncol = length(ind.variables))     # a matrix to store the betas
  n.clust <- c() 
  
  G.o <- G
  for(i in 1:G){
     
    clust.ind <- which(clust == unique(clust)[i])                         # select obs in cluster i
    
    clust.dat <- dat[clust.ind,]                                                         # create the cluster i data set
    clust.mod <- suppressWarnings(glm(form, data = clust.dat, family = mod$family))      # run a model on the cluster i data
    
    # should we stop if one cluster-specific model does not converge?
    if(drop==FALSE){
      if(clust.mod$converged == F){stop("cluster-specific model did not converge", call.=FALSE)}
      b.clust[i,] <- coefficients(clust.mod)                                                 # store the cluster i beta coefficient
      
    }else{
      if(clust.mod$converged == T){
        b.clust[i,] <- coefficients(clust.mod)                                               # store the cluster i beta coefficient
      }else{
        b.clust[i,] <- NA
      }
    }
  
  }

  if(drop==TRUE){
    b.clust <- na.omit(b.clust)
    G <- dim(b.clust)[1]
  }
  
  b.hat <- colMeans(b.clust)                                # calculate the avg beta across clusters
  b.dev <- sweep(b.clust, MARGIN = 2, STATS = b.hat)        # sweep out the avg betas
  #s.hat <- sqrt( (1 / (G-1)) * colSums(b.dev^2) )          # manually calculate the SE of beta estimates across clusters (deprecated)
  vcv.hat <- cov(b.dev)                                     # calculate VCV matrix
  rownames(vcv.hat) <- ind.variables
  colnames(vcv.hat) <- ind.variables
  s.hat <- sqrt(diag(vcv.hat))                              # calculate standard error

  t.hat <- sqrt(G) * (b.hat / s.hat)                        # calculate t-statistic
  
  # compute p-val based on # of clusters
  p.out <- 2*pmin( pt(t.hat, df = G-1, lower.tail = TRUE), pt(t.hat, df = G-1, lower.tail = FALSE) )
  
  
  # compute CIs
  ci.lo <- b.hat - qt((1-ci.level)/2, df=(G-1), lower.tail=FALSE)*(s.hat/sqrt(G))
  ci.hi <- b.hat + qt((1-ci.level)/2, df=(G-1), lower.tail=FALSE)*(s.hat/sqrt(G))
  
  out <- matrix(p.out, ncol=1)
  rownames(out) <- ind.variables

  out.p <- cbind( ind.variables, round(out, 3))
  out.p <- rbind(c("variable name", "cluster-adjusted p-value"), out.p)
  
  out.ci <- cbind(ci.lo, ci.hi)
  rownames(out.ci) <- ind.variables
  colnames(out.ci) <- c("CI lower", "CI higher")
  
  print.ci <- cbind(ind.variables, ci.lo, ci.hi)
  print.ci <- rbind(c("variable name", "CI lower", "CI higher"), print.ci)
  

  printmat <- function(m){
    write.table(format(m, justify="right"), row.names=F, col.names=F, quote=F, sep = "   ")
  }
  
  if(report==T){
    cat("\n", "Cluster-Adjusted p-values: ", "\n", "\n")
    printmat(out.p)
    
    cat("\n", "Confidence Intervals (centered on cluster-averaged results):", "\n", "\n")
    printmat(print.ci)
        
    if(G.o > G){
      cat("\n", "Note:", G.o - G, "clusters were dropped due to non-convergence.", "\n", "\n")
    }
    
  }
 
  
  out.list<-list()
  out.list[["p.values"]]<-out
  out.list[["ci"]]<-out.ci
  return(invisible(out.list))
  
}

