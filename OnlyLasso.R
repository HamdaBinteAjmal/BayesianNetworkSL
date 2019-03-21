## OnlyLasso 
lasso = ApplyLars(data$data)
lasso_inv <- 1-lasso

lasso_BN_0.9 <- ConvertToBNFromG(G1DBN::BuildEdges(lasso_inv, threshold = 0.9), 50)
lasso_BN_0.95 <- ConvertToBNFromG(G1DBN::BuildEdges(lasso_inv, threshold = 0.95), 50)
lasso_BN_0.85 <- ConvertToBNFromG(G1DBN::BuildEdges(lasso_inv, threshold = 0.85), 50)

DBNList <- list("lasso_0.9" = lasso_BN_0.9, "lasso_0.95" = lasso_BN_0.95, "lasso_0.85" = lasso_BN_0.85, "real" = realDBN)
GenerateResults(dataS, DBNList, realDBN)

## For 100
load(file = "lasso_100.rds")
lasso_BN_0.9_100 <- mclapply(lasso_100, function(x) ConvertToBNFromG(G1DBN::BuildEdges(1-x, threshold = 0.9), 50))
lasso_BN_0.95_100 <-  mclapply(lasso_100, function(x) ConvertToBNFromG(G1DBN::BuildEdges(1-x, threshold = 0.95), 50))
lasso_BN_0.85_100 <-  mclapply(lasso_100, function(x) ConvertToBNFromG(G1DBN::BuildEdges(1-x, threshold = 0.85), 50))

CalculatePrecisionAndRecallForMultiple(lasso_BN_0.9_100, datasets)
CalculatePrecisionAndRecallForMultiple(lasso_BN_0.95_100, datasets)
CalculatePrecisionAndRecallForMultiple(lasso_BN_0.85_100, datasets)