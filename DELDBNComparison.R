## DELDBN


library(bnlearn)
data= dataset$data
data2= dataS
dataD <- DifferentialData(data)
DifferentialData <- function(data)
{
  n = nrow(data)
  p = ncol(data)
  data2 = as.data.frame(matrix(NA, ncol = 2*p, nrow = n-1))
  new_names = vector(length = p)
  names = vector(length = p)
  for(i in 1:p)
  {
    new_names[i] = paste0("V", i , "_1")
    names[i] = paste0("V",i)
  }
  
  dataD <- data[2:nrow(data),] - data[1:(nrow(data)-1),]
  data2[,1:p] <- data[-1,]
  data2[,(p+1):(p*2)] <- dataD
  names(data2) <- c(names, new_names)
  return(data2)
    
}

DELDBN <- function(dataS, alpha)
{
  all_parents <- names(dataS)[which(!grepl("_", names(dataS)))]
  all_targets <- names(dataS)[which(grepl("_", names(dataS)))]

  bns <- lapply(all_targets, function(target)
  {
    data_temp <- dataS[, c(all_parents, target )]
    learn.mb(data_temp, target, alpha = alpha, "gs")
  })

  bn <- empty.graph(names(dataS))
  for (i in 1:length(bns))
  {
    target <- paste0("V", i ,"_1")
    parents <- bns[[i]]
    if (length(parents) > 0)
    {
      for (p in 1:length(parents))
      {
        bn <- set.arc(bn, from = parents[p], to = target)
      }
    }
  }
  return (bn)
}
delDBN_0.001 <- DELDBN(dataS,0.001)
delDBN_0.005 <- DELDBN(dataS,0.005)
delDBN_0.01 <- DELDBN(dataS,0.01)
delDBN_0.05 <- DELDBN(dataS, 0.05)
delDBN_0.1 <- DELDBN(dataS, 0.1)


## Not working with tDelta values, because, maybe, the data isnt simulated thatt way. 
delDBN_D_0.001 <- DELDBN(dataD, 0.001)
delDBN_D_0.005 <- DELDBN(dataD, 0.005)
delDBN_D_0.01 <- DELDBN(dataD, 0.01)
delDBN_D_0.05 <- DELDBN(dataD, 0.05)
delDBN_D_0.1 <- DELDBN(dataD, 0.1)
DBNList <- list("DELDBN_0.001" = delDBN_0.001, "DELDBN_0.005" = delDBN_0.005, "DELDBN_0.01" = delDBN_0.01 ,"DELDBN_0.05" = delDBN_0.05, "DELDBN_0.1" = delDBN_0.1, "real" = realDBN)
GenerateResults(dataS, DBNList, realDBN)

DBNList <- list("DELDBN_D_0.001" = delDBN_D_0.001, "DELDBN_D_0.005" = delDBN_D_0.005, "DELDBN_D_0.01" = delDBN_D_0.01 ,"DELDBN_D_0.05" = delDBN_D_0.05, "DELDBN_D_0.1" = delDBN_D_0.1, "real" = realDBN)
GenerateResults(dataD, DBNList, realDBN)


delDBN_0.001_s_100 <- mclapply(datasets, function(x) DELDBN(ShiftData(x$data),0.001))
CalculatePrecisionAndRecallForMultiple(delDBN_0.001_s_100, datasets)
delDBN_0.005_s_100 <- mclapply(datasets, function(x) DELDBN(ShiftData(x$data),0.005))
CalculatePrecisionAndRecallForMultiple(delDBN_0.005_s_100, datasets)

delDBN_0.01_s_100 <- mclapply(datasets, function(x) DELDBN(ShiftData(x$data),0.01))
CalculatePrecisionAndRecallForMultiple(delDBN_0.01_s_100, datasets)

delDBN_0.05_s_100 <- mclapply(datasets, function(x) DELDBN(ShiftData(x$data),0.05))
CalculatePrecisionAndRecallForMultiple(delDBN_0.05_s_100, datasets)

delDBN_0.1_s_100 <- mclapply(datasets, function(x) DELDBN(ShiftData(x$data),0.1))
CalculatePrecisionAndRecallForMultiple(delDBN_0.1_s_100, datasets)



delDBN_0.001_d_100 <- mclapply(datasets, function(x) DELDBN(DifferentialData(x$data),0.001))
CalculatePrecisionAndRecallForMultiple(delDBN_0.001_d_100, datasets)

delDBN_0.005_d_100 <- mclapply(datasets, function(x) DELDBN(DifferentialData(x$data),0.005))
CalculatePrecisionAndRecallForMultiple(delDBN_0.005_d_100, datasets)

delDBN_0.01_d_100 <- mclapply(datasets, function(x) DELDBN(DifferentialData(x$data),0.01))
CalculatePrecisionAndRecallForMultiple(delDBN_0.01_d_100, datasets)

delDBN_0.05_d_100 <- mclapply(datasets, function(x) DELDBN(DifferentialData(x$data),0.05))
CalculatePrecisionAndRecallForMultiple(delDBN_0.05_d_100, datasets)

delDBN_0.1_d_100 <- mclapply(datasets, function(x) DELDBN(DifferentialData(x$data),0.1))
CalculatePrecisionAndRecallForMultiple(delDBN_0.1_d_100, datasets)


