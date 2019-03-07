## DELDBN


library(bnlearn)
data= dataset$data#ead.table("data_1.txt")
data2= dataS

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

DBNList <- list("DELDBN_0.001" = delDBN_0.001, "DELDBN_0.005" = delDBN_0.005, "DELDBN_0.01" = delDBN_0.01 ,"DELDBN_0.05" = delDBN_0.05, "DELDBN_0.1" = delDBN_0.1, "real" = realDBN)
GenerateResults(dataS, DBNList, realDBN)

