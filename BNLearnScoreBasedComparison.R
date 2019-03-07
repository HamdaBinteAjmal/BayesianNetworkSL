## Compare with different score based algorithms in bnlearn
bl <- CreateBlackList(df)
dataS <- ShiftData(df)
#hc_bic <- hc(x = dataS, blacklist = bl, score = "bic-g"  , restart = 100 )
hc_bic_maxP_5 <- hc(x = dataS, blacklist = bl, score = "bic-g"  , restart = 100, maxp = 5)


tabu_bic_maxP_5 <- tabu(x = dataS, start = NULL, whitelist = NULL, blacklist = bl, score = "bic-g",debug = FALSE, tabu = 100, max.tabu = 100, maxp = 5, optimized = TRUE)

DBNList <- list("hc" = hc_bic_maxP_5, "tabu" = tabu_bic_maxP_5, "real" = realDBN)
GenerateResults(dataS, DBNList, realDBN)