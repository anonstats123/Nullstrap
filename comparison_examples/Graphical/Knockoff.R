library(base)  # symmetrize matrices
library(Matrix)  # symmetrize matrices
library(matrixcalc)  # to access triangular matrix conviniently
library(ppcor)  # calculate sample partial correlations
options(digits = 12)  # to avoid rounding issues in the thresholding

GenerateKnockoffs <- function(p, n){
  # Construct knock-off counterparts to the sample partial correlations. 
  #
  # Args:  
  #   p: number of variables of interest.
  #   n: number of samples.
  # 
  # Returns:
  #   p * p matrix of knock-off counterparts.
  
  elements.temp <- rt( p*(p-1)/2 , n - p )
  elements <- elements.temp/sqrt(n - p  + elements.temp^2)
  knockoffs.temp <- diag(p)
  knockoffs.temp[ lower.tri(knockoffs.temp) ] <- elements
  knockoffs <- as.matrix( forceSymmetric(knockoffs.temp, uplo = "L") )
  return( knockoffs ) 
}

SoftTholdStatistics <- function(sample.partial.cor.matrix, knockoffs){
  # Calculate the test statistics Wij.
  #
  # Args:
  #   sample.partial.cor.matrix: sample partial correaltion matrix.
  #   knockoffs: knock-off counterparts for sample partial correaltion matrix.
  #
  # Returns:
  #   matrix: p * p test statistics matrix.
  #   vector: vector of absolute values of the test statistics Wij with 1 <= i < j <= p.
  
  diff <- abs(sample.partial.cor.matrix) - abs(knockoffs)
  maximum <- (abs(sample.partial.cor.matrix) + abs(knockoffs) + abs(diff))/2
  W <- maximum * sign(diff)
  w <- abs(W[lower.tri(W)]) 
  return(list("matrix" = W,  "vector" = w))
}

FDPCheck <- function(W.vec, t, q, plus=FALSE){
  # Working function to select data-dependent threshold for statistics Wij.
  #
  # Args:
  #   W.vec: vector of test statistics Wij with 1 <= i < j <= p,
  #          where p is the number of variables of interest.
  #   t: threshold.
  #   q: target FDR level.
  #   plus: switch for regular knockoff/knockoff plus.
  #
  # Returns:
  #   Return 1 if the estimated FDP at threshold t is not greater than q. 
  #   Otherwise, return 0.
  
  return(ifelse((sum(W.vec <= -t) + as.numeric(plus))/max(sum(W.vec >= t), 1) <= q , 1, 0))
}

GraphEstimation <- function(dat, FDRtarget, plus){
  # Estimate a connectivity graph by using knock-offs.
  #
  # Args: 
  #   dat: n * p data matrix, where n is sample size and p is the number of variables.
  #   FDRtarget: target FDR level.
  #   plus: switch for regular knockoff/knockoff plus.
  #
  # Returns:
  #   edge.set: two-column matrix with each row corresponding to an 
  #             estimated edge. Empty if estimation returns no edge.
  #   adjacency.matrix: corresponding adjacency matrix.
  #   continuous.raw: estimated matrix of conditional partial correlations.
  #   continuous.thresh: thresholded estimated matrix of conditional partial correlations.
  #
  # Raises:
  #   Error if n <= p.
  
  n <- dim(dat)[1]
  p <- dim(dat)[2]
  
  if(n <= p){
    stop("Sample size must be greater than the number of variables.")
  }
  
  sample.partial.cor.matrix <- pcor(dat)$estimate
  #sample.partial.cor.matrix <- round(sample.partial.cor.matrix, digits = 10)
  knockoffs <- GenerateKnockoffs(p, n)
  
  W <- SoftTholdStatistics(sample.partial.cor.matrix, knockoffs)$matrix
  W.vec <- W[lower.tri(W)]  # to avoid double counts
  
  # Construct grid for threshonding ranging between min(|Wij|) and max(|Wij|).
  t.grid <- seq(min(abs(W.vec)), max(abs(W.vec)), length.out = 1e4)
  
  # Calculate threshold t. 
  w <- sort(SoftTholdStatistics(sample.partial.cor.matrix, knockoffs)$vector)
  thold.approx <- min(sapply(t.grid, function(x) ifelse(FDPCheck(W.vec, x, FDRtarget, plus)==1, x, Inf)))
  thold <- ifelse(length(w[w >= thold.approx])==0, Inf, min(w[w >= thold.approx])) 
  
  # Get estimated edge set.
  est.edge <- which(upper.triangle(round(W, digits = 10)) >= round(thold, digits = 10), arr.ind = T)
  rownames(est.edge) <- c()
  
  # Construct adjacency matrix based on the estimated edge set.
  est.adj <- abs(round(W, digits = 10) >= round(thold, digits = 10))
  
  # Construct a continuous estimates
  est.cont.raw <- sample.partial.cor.matrix
  est.cont.tresh <- est.cont.raw * est.adj  # hard-thresholding
  
  return(list("edge.set" = data.frame(est.edge), "adjacency.matrix" = est.adj, "continuous.raw" = est.cont.raw, "continuous.thresh" = est.cont.tresh))
}

example <- function()
{
  ### Example
  FDRtarget <- 0.2
  brain.dat <- read.csv("Data/1.csv", header = FALSE)
  
  set.seed(1)
  edge.set <- GraphEstimation(brain.dat, FDRtarget, FALSE)$edge.set
  adjacency.matrix <- GraphEstimation(brain.dat, FDRtarget, FALSE)$adjacency.matrix
  
  dim(edge.set)[1]  # number of estimated edges
  edge.set[1:3,]  # displays the first three edges (2 3; 1 7; 4 8)
}
