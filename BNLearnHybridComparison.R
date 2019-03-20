## Compare againt hybrid algorithms in BNLearn
library(bnlearn)
library(G1DBN)
seeds = 1
p = 50
n = 20
noiseLevel = 10

#dataset = SimulateData(p,n,noiseLevel,1)

#dataS = ShiftData(dataset$data)
bl <- CreateBlackList(dataset$data)
mmhc_bic_hc_0.01_max5 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.01), maximize.args = list(restart = 100, maxp = 5) )
mmhc_bic_hc_0.05_max5 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.05), maximize.args = list(restart = 100, maxp = 5) )
mmhc_bic_hc_0.1_max5 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.1), maximize.args = list(restart = 100, maxp = 5) )
mmhc_bic_hc_0.5_max5 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.5), maximize.args = list(restart = 100, maxp = 5) )
mmhc_bic_hc_0.7_max5 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.7), maximize.args = list(restart = 100, maxp = 5) )
mmhc_bic_hc_0.9_max5 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.9), maximize.args = list(restart = 100, maxp = 5) )
mmhc_bic_hc_0.99_max5 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.99), maximize.args = list(restart = 100, maxp = 5) )

DBNList <- list("mmhc_0.01_max5" = mmhc_bic_hc_0.01_max5, "mmhc_0.05_max5" = mmhc_bic_hc_0.05_max5, "mmhc_0.1_max5" = mmhc_bic_hc_0.1_max5,
                "mmhc_0.5_max5" = mmhc_bic_hc_0.5_max5, "mmhc_0.7_max5" = mmhc_bic_hc_0.7_max5,"mmhc_0.9_max5" = mmhc_bic_hc_0.9_max5,
                "mmhc_0.99_max5" = mmhc_bic_hc_0.99_max5,
                "real" = realDBN)
realDBN <- ConvertToBN(dataset$RealNet)
GenerateResults(dataS, DBNList, realDBN)


mmhc_bic_hc_0.01 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.01), maximize.args = list(restart = 100) )
mmhc_bic_hc_0.05 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.05), maximize.args = list(restart = 100) )
mmhc_bic_hc_0.1 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.1), maximize.args = list(restart = 100) )
mmhc_bic_hc_0.5 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.5), maximize.args = list(restart = 100) )
mmhc_bic_hc_0.7 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.7), maximize.args = list(restart = 100) )
mmhc_bic_hc_0.9 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.9), maximize.args = list(restart = 100) )
mmhc_bic_hc_0.99 <- mmhc(x = dataS, blacklist = bl, restrict.args = list(alpha = 0.99), maximize.args = list(restart = 100) )



DBNList <- list("mmhc_0.01" = mmhc_bic_hc_0.01, "mmhc_0.05" = mmhc_bic_hc_0.05, "mmhc_0.1" = mmhc_bic_hc_0.1,
                "mmhc_0.5" = mmhc_bic_hc_0.5, "mmhc_0.7" = mmhc_bic_hc_0.7,"mmhc_0.9" = mmhc_bic_hc_0.9,"mmhc_0.99" = mmhc_bic_hc_0.99,
                "real" = realDBN)
realDBN <- ConvertToBN(data$RealNet)
GenerateResults(dataS, DBNList, realDBN)


rs_hc_0.01 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
       maximize = "hc", restrict.args = list(alpha = 0.01), maximize.args = list(restart = 100), debug = FALSE)

rs_hc_0.05 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                     maximize = "hc", restrict.args = list(alpha = 0.05), maximize.args = list(restart = 100), debug = FALSE)

rs_hc_0.1 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                     maximize = "hc", restrict.args = list(alpha = 0.1), maximize.args = list(restart = 100), debug = FALSE)

rs_hc_0.5 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                     maximize = "hc", restrict.args = list(alpha = 0.5), maximize.args = list(restart = 100), debug = FALSE)

rs_hc_0.7 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                     maximize = "hc", restrict.args = list(alpha = 0.7), maximize.args = list(restart = 100), debug = FALSE)

rs_hc_0.9 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                    maximize = "hc", restrict.args = list(alpha = 0.9), maximize.args = list(restart = 100), debug = FALSE)

rs_hc_0.99 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                    maximize = "hc", restrict.args = list(alpha = 0.99), maximize.args = list(restart = 100), debug = FALSE)
DBNList <- list("rs_hc_0.01" = rs_hc_0.01, "rs_hc_0.05" = rs_hc_0.05, "rs_hc_0.1" = rs_hc_0.1,
                "rs_hc_0.5" = rs_hc_0.5, "rs_hc_0.7" = rs_hc_0.7,"rs_hc_0.9" = rs_hc_0.9,"rs_hc_0.99" = rs_hc_0.99,
                "real" = realDBN)
GenerateResults(dataS, DBNList, realDBN)


rs_tabu_0.01 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                     maximize = "tabu", restrict.args = list(alpha = 0.01), maximize.args = list(max.tabu = 100, tabu = 100), debug = FALSE)

rs_tabu_0.05 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                     maximize = "tabu", restrict.args = list(alpha = 0.05), maximize.args = list(max.tabu = 100, tabu = 100), debug = FALSE)

rs_tabu_0.1 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                    maximize = "tabu", restrict.args = list(alpha = 0.1), maximize.args = list(max.tabu = 100, tabu = 100), debug = FALSE)

rs_tabu_0.5 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                    maximize = "tabu", restrict.args = list(alpha = 0.5), maximize.args = list(max.tabu = 100, tabu = 100), debug = FALSE)

rs_tabu_0.7 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                    maximize = "tabu", restrict.args = list(alpha = 0.7), maximize.args = list(max.tabu = 100, tabu = 100), debug = FALSE)

rs_tabu_0.9 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                    maximize = "tabu", restrict.args = list(alpha = 0.9), maximize.args = list(max.tabu = 100, tabu = 100), debug = FALSE)

rs_tabu_0.99 <- rsmax2(dataS,blacklist = bl, restrict = "si.hiton.pc",
                     maximize = "tabu", restrict.args = list(alpha = 0.99), maximize.args = list(max.tabu = 100, tabu = 100), debug = FALSE)

DBNList <- list("rs_tabu_0.01" = rs_tabu_0.01, "rs_tabu_0.05" = rs_tabu_0.05, "rs_tabu_0.1" = rs_tabu_0.1,
                "rs_tabu_0.5" = rs_tabu_0.5, "rs_tabu_0.7" = rs_tabu_0.7,"rs_tabu_0.9" = rs_tabu_0.9,"rs_tabu_0.99" = rs_tabu_0.99,
                "real" = realDBN)
GenerateResults(dataS, DBNList, realDBN)


### 100 runs ##
#datasets


## mmhc 3 parents
mmhc_bic_hc_0.7_max3_100 <- mclapply(datasets, function(x) mmhc(x = ShiftData(x$data), blacklist = bl, restrict.args = list(alpha = 0.7), maximize.args = list(restart = 100, maxp = 3) ))
mmhc_bic_hc_0.9_max3_100 <- mclapply(datasets, function(x) mmhc(x = ShiftData(x$data), blacklist = bl, restrict.args = list(alpha = 0.9), maximize.args = list(restart = 100, maxp = 3) ))
#mmhc_bic_hc_0.99_max3_100 <- mclapply(datasets, function(x) mmhc(x = ShiftData(x$data), blacklist = bl, restrict.args = list(alpha = 0.99), maximize.args = list(restart = 100, maxp = 3) ))

##mmhc 
mmhc_bic_hc_0.1_100 <- mclapply(datasets, function(x) mmhc(x = ShiftData(x$data), blacklist = bl, restrict.args = list(alpha = 0.1), maximize.args = list(restart = 1) ))
mmhc_bic_hc_0.7_100 <- mclapply(datasets, function(x) mmhc(x = ShiftData(x$data), blacklist = bl, restrict.args = list(alpha = 0.7), maximize.args = list(restart = 1) ))
mmhc_bic_hc_0.9_100 <- mclapply(datasets, function(x) mmhc(x = ShiftData(x$data), blacklist = bl, restrict.args = list(alpha = 0.9), maximize.args = list(restart = 1) ))
#mmhc_bic_hc_0.99_100 <- mclapply(datasets, function(x) mmhc(x = ShiftData(x$data), blacklist = bl, restrict.args = list(alpha = 0.99), maximize.args = list(restart = 100) ))


## rsmax2
rs_hc_0.7_100 <- mclapply(datasets,function(x) rsmax2(ShiftData(x$data), blacklist = bl, restrict = "si.hiton.pc",
                    maximize = "hc", restrict.args = list(alpha = 0.7), maximize.args = list(restart = 1), debug = FALSE))

rs_hc_0.9_100 <-  mclapply(datasets,function(x) rsmax2(ShiftData(x$data), blacklist = bl, restrict = "si.hiton.pc",
                    maximize = "hc", restrict.args = list(alpha = 0.9), maximize.args = list(restart = 1), debug = FALSE))

rs_hc_0.1_100 <- mclapply(datasets,function(x) rsmax2(ShiftData(x$data), blacklist = bl, restrict = "si.hiton.pc",
                     maximize = "hc", restrict.args = list(alpha = 0.1), maximize.args = list(restart = 1), debug = FALSE))

rs_tabu_0.9_100 <- mclapply(datasets,function(x) rsmax2(ShiftData(x$data),blacklist = bl, restrict = "si.hiton.pc",
                      maximize = "tabu", restrict.args = list(alpha = 0.9), maximize.args = list(max.tabu = 100, tabu = 100), debug = FALSE))

rs_tabu_0.01_100 <- mclapply(datasets,function(x) rsmax2(ShiftData(x$data),blacklist = bl, restrict = "si.hiton.pc",
                                             maximize = "tabu", restrict.args = list(alpha = 0.01), maximize.args = list(max.tabu = 100, tabu = 100), debug = FALSE))
