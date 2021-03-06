% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterBS.mlogit.R
\name{cluster.bs.mlogit}
\alias{cluster.bs.mlogit}
\title{Pairs Cluster Bootstrapped p-Values For mlogit}
\usage{
cluster.bs.mlogit(
  mod,
  dat,
  cluster,
  ci.level = 0.95,
  boot.reps = 1000,
  cluster.se = TRUE,
  report = TRUE,
  prog.bar = TRUE,
  output.replicates = FALSE,
  seed = NULL
)
}
\arguments{
\item{mod}{A model estimated using \code{mlogit}.}

\item{dat}{The data set used to estimate \code{mod}.}

\item{cluster}{A formula of the clustering variable.}

\item{ci.level}{What confidence level should CIs reflect?}

\item{boot.reps}{The number of bootstrap samples to draw.}

\item{cluster.se}{Use clustered standard errors (= TRUE) or ordinary SEs (= FALSE) for bootstrap replicates.}

\item{report}{Should a table of results be printed to the console?}

\item{prog.bar}{Show a progress bar of the bootstrap (= TRUE) or not (= FALSE).}

\item{output.replicates}{Should the cluster bootstrap coefficient replicates be output (= TRUE) or not (= FALSE)?}

\item{seed}{Random number seed for replicability (default is NULL).}
}
\value{
A list with the elements
\item{p.values}{A matrix of the estimated p-values.}
\item{ci}{A matrix of confidence intervals.}
}
\description{
This software estimates p-values using pairs cluster bootstrapped t-statistics for multinomial logit models (Cameron, Gelbach, and Miller 2008). The data set is repeatedly re-sampled by cluster, a model is estimated, and inference is based on the sampling distribution of the pivotal (t) statistic.
}
\note{
Code to estimate GLM clustered standard errors by Mahmood Arai: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/, although modified slightly to work for \code{mlogit} models. Cluster SE degrees of freedom correction = (M/(M-1)) with M = the number of clusters.
}
\examples{
\dontrun{

#######################################
# example one: train ticket selection
#######################################
require(mlogit)
data("Train", package="mlogit")
Train$choiceid <- 1:nrow(Train)

Tr <- dfidx(Train, shape = "wide", varying = 4:11, sep = "_", 
          choice = "choice", idx = list(c("choiceid", "id")), 
          idnames = c(NA, "alt"))
Tr$price <- Tr$price/100 * 2.20371
Tr$time <- Tr$time/60

ml.Train <- mlogit(choice ~ price + time + change + comfort | -1, Tr)

# compute pairs cluster bootstrapped p-values
# note: few reps to speed up example
cluster.bs.tr <- cluster.bs.mlogit(ml.Train, Tr, ~ id, boot.reps=100)


##################################################################
# example two: predict type of heating system installed in house
##################################################################
require(mlogit)
data("Heating", package = "mlogit")
H <- Heating
H$region <- as.numeric(H$region)
H.ml <- dfidx(H, shape="wide", choice="depvar", varying=c(3:12),
         idx = list(c("idcase", "region")))
m <- mlogit(depvar~ic+oc, H.ml)

# compute pairs cluster bootstrapped p-values
cluster.bs.h <- cluster.bs.mlogit(m, H.ml, ~ region, boot.reps=1000)

}
}
\references{
Esarey, Justin, and Andrew Menger. 2017. "Practical and Effective Approaches to Dealing with Clustered Data." \emph{Political Science Research and Methods} forthcoming: 1-35. <URL:http://jee3.web.rice.edu/cluster-paper.pdf>.

Cameron, A. Colin, Jonah B. Gelbach, and Douglas L. Miller. 2008. "Bootstrap-Based Improvements for Inference with Clustered Errors." \emph{The Review of Economics and Statistics} 90(3): 414-427. <DOI:10.1162/rest.90.3.414>.
}
\author{
Justin Esarey
}
