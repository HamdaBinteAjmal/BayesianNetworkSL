## This code has functions to compare the G1DBN, shrinkage and lasso methods


library(G1DBN)
library(bnlearn)
library(caret)
library(longitudinal)
library(fdrtool)
library("GeneNet") # version 1.2 or later
library(lars)

# load functions for estimating shrinkage VAR model
source("shrink-VAR-R-code\\mvr.shrink.R") # define mvr.shrink function
source("shrink-VAR-R-code\\shrink.VAR.R") 

# functions for plotting the resulting VAR network
source("shrink-VAR-R-code\\makedotfile.VAR.R")

library(lars)

## This function calculates the lasso scores for the genes. This is an inner function called by 
## ApplyLars function. X = predictors and Y = response variables.
LASSO <- function(X,Y)
{  cv <- lars::cv.lars(x = X, y = Y, plot.it = FALSE)
  ideal_l1_ratio <- cv$index[which.max(cv$cv - cv$cv.error <= min(cv$cv))]
  obj <- lars::lars(x = X, y = Y, type = "lasso")
  scaled_coefs <- scale(obj$beta, FALSE, 1 / obj$normx)
  l1 <- apply(X = scaled_coefs, MARGIN = 1, FUN = function(x) sum(abs(x)))
  c1 = coef(obj)[which.max(l1 / tail(l1, 1) > ideal_l1_ratio),]
  return(c1)
}

## Outer function to calculate lasso scores for the variables. Require data frame as the input. 
## Automatically shifts the data over a time slice. Predictors are variables of previous time slice. 
## Response are variables of next time slice. 
ApplyLars <- function(data)
{
  predictors  = as.matrix(data[1:(nrow(data)-1),])
  response = data[2:nrow(data),]
  nR = ncol(response)
  arcs = lapply(1:nR, function(y) 
  {
    print(y)
    Y = response[,y]
    arcs = LASSO(X = predictors, Y = Y)
    
  })
  scores = abs(do.call(rbind, arcs))
  return(scores)
}

## Function for the shrinkage application.
## Inner function
Shrinkage <- function(out.Var, df=10, plot.locfdr=1)
{
  
  B_est<-out.Var$B_est
  r.mat<-out.Var$r.mat
  p<-ncol(B_est)
  value<-as.vector(r.mat)
  
  node.list <- value
  ## calculate significance
  
  pc <- node.list
  
  # locfdr
  require("fdrtool")
  if (any(pc > 1) || any(pc < -1))
    stop("Data out of range: input correlations must be in [-1; 1]")
  out <-fdrtool(z.transform(pc),"correlation", plot = FALSE)#, cutoff.method = "locfdr")   # 
  
  prob.nonzero <- 1 - out$lfdr
  node.sig <- cbind(node.list, out$pval, out$qval, out$lfdr, prob.nonzero)#
  pval = matrix(out$pval, p)
  qval = matrix(out$qval, p)
  lfdr = matrix(out$lfdr, p)
  pval = matrix(out$pval, p)
   return(pval)
}
GetScores <- function(data, index )
{
  source("shrink-VAR-R-code/mvr.shrink.R") # define mvr.shrink function
  source("shrink-VAR-R-code/shrink.VAR.R") 
   source("shrink-VAR-R-code/makedotfile.VAR.R")
  
  print(paste0("We are at :", index))
  #################################################################
  
  # load Arabidopsis thaliana data set
  data = as.longitudinal(data)
  VAR.coeff <- shrink.VAR(data)
  
  names(VAR.coeff)
  # "B_est"      "r.mat"      "lambda"     "lambda.var"
  
  # assign local fdr values to each edge (prob=1-local fdr)
  results.VAR <- Shrinkage(VAR.coeff)
  return(results.VAR)
  
  
}

## Same as G1DBN PRCurve with slightly different methods to calculate the Precision and Recall
PRCurve_GeneNet <- function (score, validMat, dec = TRUE, index = 0) 
{
  # print(index)
  p <- dim(score)[1]
  q <- dim(score)[2]
  
  nbPos = sum(validMat)
  tri <- sort(score, decreasing = dec)
  precision = array(0, length(tri))
  recall = array(0, length(tri))
  for (j in 1:length(tri)) {
    if (length(tri) > 10) {
      if ((j%%(round(length(tri)/10))) == 0) {
        cat(round(10 * j/(round(length(tri)/10))), "% ")
      }
    }
    
    if(dec)
    {
      learnedNet = (score >= tri[j])
    }
    else
    {
      learnedNet = (score <= tri[j])
    }
    learnedNet[which(is.na(learnedNet))] = FALSE
    TP = sum((learnedNet==TRUE) * (validMat == TRUE), na.rm = TRUE)
    FP = sum((learnedNet == TRUE)  * (validMat == FALSE ), na.rm = TRUE)
    FN = sum((learnedNet == FALSE) * (validMat == TRUE), na.rm = TRUE)
    recall[j] =  TP / (TP + FN)
    precision[j] = TP / (TP + FP)
    
    
  }
  cat("\n")
  return(list(recall = recall, precision = precision))
}

ExecuteG1DBNS2 <- function(S1, data, alpha1 ,alpha2, euler = FALSE)
{
  S2 <- DBNScoreStep2(S1$S1ls, data, method='ls', alpha1)
  S2 <- DBNScoreStep2Euler(S1$S1ls, data, method='ls', alpha1)
  G2 <- BuildEdges(S2,threshold=alpha2,dec=FALSE)
  G2
}

ExecuteG1DBNS1 <- function(data, alpha1, euler = FALSE)
{
  S1 <- DBNScoreStep1(data, method='ls')
  G1 <- BuildEdges(S1$S1ls,threshold=alpha1,dec=FALSE)
  return(list("mat" = S1, "graph" = G1))
}
SimulateData <- function(Genes, timepoints, noiseLevel = 0, seed = 1)
{
  n = timepoints
  p = Genes
  set.seed(seed)
  ## the network - adjacency Matrix
  prop = 0.05
  MyNet <- SimulNetworkAdjMatrix(p,prop,c(-0.95,-0.05,0.05,0.95))
  ## initializing the B vector
  ## initializing the variance of the noise
  sigmaEps <- runif(p,0.03, 0.08)
  B1 <- runif(p,-0.95, -0.05)
  B2 <- runif(p,0.05, 0.95)
  B = unlist(lapply(1:p, function (x) sample(c(B1[x], B2[x]), 1))) 
  ## initializing the process Xt
  X0 <- B +  runif(p,0.03, 0.08)#rnorm(1,0,sigmaEps*noiseLevel)
  ## number of time points
  ## the AR(1) time series process
  #sigmaEps <- c(0.03,0.08)
  Xn <- SimulGeneExpressionAR1(MyNet$A,B,X0,sigmaEps,n)
  plot(1:n, Xn[,1],type="l", xlab="Time t", ylab="X(t)",
       main="Simulated AR(1) time series", ylim=c(min(Xn),max(Xn)))
  for (i in 2:p){
    lines(1:n,Xn[,i],col=i)
  }
  ml = list("data" = Xn, "RealNet"  = MyNet)
  return (ml)
}
ConvertToBN <-function(MyNet)
{
  adjMat = MyNet$AdjMatrix
  p = ncol(adjMat)
  newAdjMat = matrix(0, p*2,p*2)
  names = paste0("V", 1:p)
  names = c(names, paste0(names,"_1"))
  colnames(newAdjMat) = names
  rownames(newAdjMat) =  colnames(newAdjMat)
  newAdjMat[1:p, 1:p+p] = t(MyNet$AdjMatrix)
  new_bn = bnlearn::empty.graph(colnames(newAdjMat))
  amat(new_bn) = newAdjMat
  new_bn
}
ConvertToBNFromG <- function(G, geneCount)
{
  p = geneCount
  newAdjMat = matrix(0, p*2,p*2)
  names = paste0("V", 1:p)
  names = c(names, paste0(names,"_1"))
  colnames(newAdjMat) = names
  rownames(newAdjMat) =  colnames(newAdjMat)
  for(i in 1:nrow(G))
  {
    pred = G[i,1]
    tar = G[i,2]
    newAdjMat[pred,p+tar] = 1
  }
  new_bn = bnlearn::empty.graph(colnames(newAdjMat))
  amat(new_bn) = newAdjMat
  plot(new_bn)
  new_bn
}

ShiftData <- function(data)
{
  
  n = nrow(data)
  p = ncol(data)
  data2 = as.data.frame(matrix(NA, ncol = 2*p, nrow = n-1))
  data2[1:n-1, 1:p] = data[1:n-1,]
  new_names = vector(length = p)
  colnames(data) = paste0("V", 1:p)
  for(i in 1:ncol(data))
  {
    new_names[i] = paste0(colnames(data)[i] , "_1")
    
    
    data2[,i+ncol(data)] = data[2:n, i]
  }
  colnames(data2) =  c(colnames(data), new_names)
  data2
}
CreateBlackList <- function(data)
{
  Slice1names = paste0("V", 1:ncol(data), "_1")
  Slice0names = paste0("V", 1:ncol(data))
  bl1 = expand.grid(Slice1names, Slice0names)
  bl2 = expand.grid(Slice0names, Slice0names)
  bl3 = expand.grid(Slice1names, Slice1names)
  blacklist = rbind(bl1,bl2,bl3)
  blacklist <- as.matrix(blacklist)
  colnames(blacklist) = c("from", "to")
  blacklist
}

GetAMAT <- function(mat)
{
  genes = nrow(mat) / 2
  newMat = mat[1:genes, (genes+1):ncol(mat)] 
  
}

## Generate Results Table ##
GenerateResults <- function(dataS, DBNList, RealDBN)
{
  
  conf_mats = lapply(DBNList, function(x) confusionMatrix( reference = factor(GetAMAT(amat(RealDBN))), data = factor(GetAMAT(amat(x))), positive = "1"))
  
  precision = unlist(lapply(conf_mats, function(x) x$byClass[3]))
  recall = unlist(lapply(conf_mats, function(x) x$byClass[1]))
  hamming = unlist(lapply(DBNList, function(x) hamming(learned = x, true = RealDBN)))
  loglik = unlist(lapply(DBNList, function(x) score(x = x, data = dataS, type = "bge") ))
  aic = unlist(lapply(DBNList, function(x) score(x = x, data = dataS, type = "loglik-g") ))
  bic = unlist(lapply(DBNList, function(x) score(x = x, data = dataS, type = "aic-g") ))
  bge = unlist(lapply(DBNList, function(x) score(x = x, data = dataS, type = "bic-g") ))
  accuracy= unlist(lapply(conf_mats, function(x) x$overall[1]))
  balanced_accuracy = unlist(lapply(conf_mats, function(x) x$byClass[11]))
  f1 = unlist(lapply(conf_mats, function(x) x$byClass[7]))
  TP = unlist(lapply(conf_mats, function(x) x$table[2,2]))
  TN =  unlist(lapply(conf_mats, function(x) x$table[1,1]))
  FP =  unlist(lapply(conf_mats, function(x) x$table[2,1]))
  FN =  unlist(lapply(conf_mats, function(x) x$table[1,2]))
  results_Noise1 = data.frame(row.names = c("F1", "Balanced Accuracy" ,"Accuracy", "PPV", "SENS", "HAMMING", "TP", "TN", "FP", "FN", "LOGLIK", "AIC", "BIC", "BGE"))
  results_Noise1 = rbind(f1, balanced_accuracy, accuracy, precision, recall, hamming, TP, TN, FP, FN, loglik,aic, bic, bge )
  colnames(results_Noise1) = names(DBNList)
  results_Noise1
}
CalculateMeans <- function(curves)
{
  max = max(unlist(lapply(curves, function(x) length(x$recall))))
  
  for(i in 1:length(curves))
  {
    length(curves[[i]]$recall) = max
    length(curves[[i]]$precision) = max
  }
  newdt = lapply(curves, as.data.frame)
  ss = do.call(cbind, newdt)
  all_recalls = ss[, grep("recall", colnames(ss))]
  all_precisions = ss[, grep("precision", colnames(ss))]
  mean_recalls = rowMeans(all_recalls, na.rm = TRUE)
  mean_precisions = rowMeans(all_precisions, na.rm = TRUE)
  return(list("recall"  = mean_recalls, "precisions" = mean_precisions))
}
