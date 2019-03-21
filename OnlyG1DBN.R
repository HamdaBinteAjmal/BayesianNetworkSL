## Generate results for G1DBN at alpha1 = 0.7 and alpha2 = 0.01, 0.05, 0.1



alpha1=0.7
out<-DBNScoreStep1(data$data,method='ls')
S2<-DBNScoreStep2(out$S1ls,data = data$data,method='ls',alpha1=alpha1)



S2_0.01 <- G1DBN::BuildEdges(S2, threshold = 0.01)
S2_0.01 <- ConvertToBNFromG(S2_0.01, ncol(data$data))


S2_0.05 <- G1DBN::BuildEdges(S2, threshold = 0.05)
S2_0.05 <- ConvertToBNFromG(S2_0.05, ncol(data$data))

S2_0.1 <- G1DBN::BuildEdges(S2, threshold = 0.1)
S2_0.1 <- ConvertToBNFromG(S2_0.1, ncol(data$data))

DBNList <- list("S2_0.01" = S2_0.01, "S2_0.05" = S2_0.05, "S2_0.1" = S2_0.1, "real" = realDBN)
GenerateResults(dataS, DBNList, realDBN)


load(file = "G1s.rds")
G1DBN_0.7_100 = lapply(G1s, function(x) ConvertToBNFromG(BuildEdges(x$mat$S1ls, 0.7),50))
CalculatePrecisionAndRecallForMultiple(G1DBN_0.7_100, datasets)

G1DBN_0.3_100 = lapply(G1s, function(x) ConvertToBNFromG(BuildEdges(x$mat$S1ls, 0.3),50))
CalculatePrecisionAndRecallForMultiple(G1DBN_0.3_100, datasets)

G1DBN_0.1_100 = lapply(G1s, function(x) ConvertToBNFromG(BuildEdges(x$mat$S1ls, 0.1),50))
CalculatePrecisionAndRecallForMultiple(G1DBN_0.1_100, datasets)