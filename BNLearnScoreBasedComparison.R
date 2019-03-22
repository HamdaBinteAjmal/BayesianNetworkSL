## Compare with different score based algorithms in bnlearn
bl <- CreateBlackList(df)
dataS <- ShiftData(df)
#hc_bic <- hc(x = dataS, blacklist = bl, score = "bic-g"  , restart = 100 )
hc_bic_maxP_3 <- hc(x = dataS, blacklist = bl, score = "bic-g"  , restart = 10, maxp = 3)


tabu_bic_maxP_3 <- tabu(x = dataS, start = NULL, whitelist = NULL, blacklist = bl, score = "bic-g",debug = FALSE, tabu = 100, max.tabu = 100, maxp = 3, optimized = TRUE)

DBNList <- list("hc" = hc_bic_maxP_3, "tabu" = tabu_bic_maxP_3, "real" = realDBN)
GenerateResults(dataS, DBNList, realDBN)



hc_maxP_3_100 <- mclapply(datasets, function(x) hc(x = ShiftData(x$data), blacklist = bl, score = "bic-g"  , restart = 1, maxp = 3))

tabu_maxP_3_100 <- mclapply(datasets, function(x) tabu(x = ShiftData(x$data), start = NULL, whitelist = NULL, blacklist = bl, score = "bic-g",debug = FALSE, tabu = 100, max.tabu = 100, maxp = 3, optimized = TRUE))

CalculatePrecisionAndRecallForMultiple(hc_maxP_3_100, datasets)

CalculatePrecisionAndRecallForMultiple(tabu_maxP_3_100, datasets)

