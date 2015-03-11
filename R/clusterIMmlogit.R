#' Cluster-Adjusted Confidence Intervals and p-Values for mlogit
#'
#' Computes p-values and confidence intervals for multinomial logit models based on cluster-specific model estimation (Ibragimov and Muller 2010). A separate model is estimated in each cluster, and then p-values and confidence intervals are computed based on a t/normal distribution of the cluster-specific estimates.
#'
#' @param mod A model estimated using \code{mlogit}.
#' @param dat The data set used to estimate \code{mod}.
#' @param cluster A formula of the clustering variable.
#' @param ci.level What confidence level should CIs reflect?
#' @param report Should a table of results be printed to the console?
#' @param truncate Should outlying cluster-specific beta estimates be excluded?
#'
#' @return A list with the elements
#' \item{p.values}{A matrix of the estimated p-values.}
#' \item{ci}{A matrix of confidence intervals.}
#' @author Justin Esarey
#' @note Confidence intervals are centered on the cluster averaged estimate, which can diverge from original model estimates if clusters have different numbers of observations. Consequently, confidence intervals may not be centered on original model estimates.
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
#' # compute cluster-adjusted p-values (takes a while)
#' clust.p <- cluster.im.mlogit(vt.mod, dat=vtinpat.mlogit, cluster = ~ hos.num, 
#'            ci.level = 0.95, report=TRUE, truncate=TRUE)          
#' }
#' @rdname cluster.im.mlogit
#' @import mlogit
#' @references Ibragimov, Rustam, and Ulrich K. Muller. 2010. "t-Statistic Based Correlation and Heterogeneity Robust Inference." \emph{Journal of Business & Economic Statistics} 28(4): 453-468. 
#' @export

cluster.im.mlogit<-function(mod, dat, cluster, ci.level = 0.95, report = TRUE, truncate = FALSE){
  
  form <- mod$formula                                     # what is the formula of this model?
  
  variables <- all.vars(form)                             # what variables are in this model?
  clust <- subset(dat, select = all.vars(cluster))    # add the cluster variable into dat.t
  G<-dim(unique(clust))[1]                                # how many clusters are in this model?
  
  b.clust <- NULL
  
  for(i in 1:G){
    
    clust.ind <- which(clust == unlist(unique(clust))[i])                         # select obs in cluster i
    
    clust.dat <- dat[clust.ind,]                                          # create the cluster i data set
    clust.mod <- tryCatch( mlogit(form, data = clust.dat), error = function(e){NA})                  # run a model on the cluster i data

    if(class(clust.mod)=="mlogit"){
      
      if(is.null(b.clust)){b.clust<-matrix(data=NA, nrow=G, ncol=length(coefficients(clust.mod)))}
      b.clust[i,] <- coefficients(clust.mod)                              # store the cluster i beta coefficient

    }else{
      stop("model does not estimate in at least one cluster")
    }
  }
  
  # remove clusters with outlying betas
  dropped <- 0
  if(truncate==TRUE){

    IQR <- apply(FUN=quantile, MARGIN=2, X=b.clust, probs=c(0.25, 0.75))
    
    b.clust.save <- b.clust
    for(i in 1:dim(b.clust)[2]){
      b.clust.save[,i] <- ifelse( abs(b.clust[,i]) > (abs(mean(b.clust[,i])) + 6*abs(IQR[2,i] - IQR[1,i])), 0, 1)
    }
    
    save.clust <- apply(X=b.clust.save, MARGIN=1, FUN=min)
    dropped <- dim(b.clust)[1] - sum(save.clust)
    
    b.clust.adj <- cbind(b.clust, save.clust)
    
    b.clust <- subset(b.clust, subset=(save.clust==1), select=1:dim(b.clust)[2])
    
  }
  
  b.hat <- colMeans(b.clust)                                # calculate the avg beta across clusters
  b.dev <- sweep(b.clust, MARGIN = 2, STATS = b.hat)        # sweep out the avg betas
  #s.hat <- sqrt( (1 / (G-1)) * colSums(b.dev^2) )          # manually calculate the SE of beta estimates across clusters (deprecated)
  vcv.hat <- cov(b.dev)                                     # calculate VCV matrix
  s.hat <- sqrt(diag(vcv.hat))                              # calculate standard error
  
  t.hat <- sqrt(G) * (b.hat / s.hat)                        # calculate t-statistic
  
  # compute p-val based on # of clusters
  p.out <- 2*pmin( pt(t.hat, df = G-1, lower.tail = TRUE), pt(t.hat, df = G-1, lower.tail = FALSE) )
  
  # compute CIs
  ci.lo <- b.hat - qt((1-ci.level)/2, df=(G-1), lower.tail=FALSE)*(s.hat/sqrt(G))
  ci.hi <- b.hat + qt((1-ci.level)/2, df=(G-1), lower.tail=FALSE)*(s.hat/sqrt(G))
    
  out <- matrix(p.out, ncol=1)
  out.p <- cbind( names(coefficients(summary(mod))), round(out,3) )
  
  out.ci <- cbind(ci.lo, ci.hi)
  rownames(out.ci) <- names(coefficients(summary(mod)))
  colnames(out.ci) <- c("CI lower", "CI higher")
  
  print.ci <- cbind(names(coefficients(summary(mod))), ci.lo, ci.hi)
  print.ci <- rbind(c("variable name", "CI lower", "CI higher"), print.ci)
  

  printmat <- function(m){
    write.table(format(m, justify="right"), row.names=F, col.names=F, quote=F, sep = "   ")
  }
    
  if(report==T){
    
    cat("\n", "Cluster-Adjusted p-values: ", "\n", "\n")
    printmat(out.p)
    
    cat("\n", "Confidence Intervals (centered on cluster-averaged results):", "\n", "\n")
    printmat(print.ci)
        
    if(dropped > 0){cat("\n", "Note:", dropped, "clusters were dropped as outliers.", "\n", "\n")}
    
  }
  

  out.list<-list()
  out.list[["p.values"]]<-out
  out.list[["ci"]]<-out.ci
  return(invisible(out.list))
  
}
