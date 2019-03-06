## Same problem try to solve by different algos  in bnlearn
##Comparison with bnlearn constraint base algos

bl <- CreateBlackList(df)

## constrained Based
pc_0.05 <- pc.stable(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl,
          alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
pc_0.01 <- pc.stable(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl,
                     alpha = 0.01, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
pc_0.1 <- pc.stable(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl,
                     alpha = 0.1, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)

gs_0.05 <- gs(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, 
   alpha = 0.05, debug = FALSE, optimized = FALSE, strict = FALSE, undirected = FALSE)


gs_0.01 <- gs(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, 
              alpha = 0.01, debug = FALSE, optimized = FALSE, strict = FALSE, undirected = FALSE)
gs_0.1 <- gs(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, 
              alpha = 0.1, debug = FALSE, optimized = FALSE, strict = FALSE, undirected = FALSE)


iamb_0.05 <- iamb(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
     alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
     strict = FALSE, undirected = FALSE)
iamb_0.01 <- iamb(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
                  alpha = 0.01, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
                  strict = FALSE, undirected = FALSE)
iamb_0.1 <- iamb(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
                  alpha = 0.1, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
                  strict = FALSE, undirected = FALSE)
iamb_0.2 <- iamb(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
                 alpha = 0.2, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
                 strict = FALSE, undirected = FALSE)

fiamb_0.05 <- fast.iamb(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
          alpha = 0.05, B = NULL, max.sx = NULL, debug = TRUE, optimized = FALSE,
          strict = FALSE, undirected = FALSE)
fiamb_0.01 <- fast.iamb(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
                        alpha = 0.01, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
                        strict = FALSE, undirected = FALSE)
fiamb_0.1 <- fast.iamb(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
                        alpha = 0.1, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
                        strict = FALSE, undirected = FALSE)

iiamb_0.05 <- inter.iamb(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
           alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
           strict = FALSE, undirected = FALSE)
iiamb_0.01 <- inter.iamb(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
                         alpha = 0.01, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
                         strict = FALSE, undirected = FALSE)
iiamb_0.1 <- inter.iamb(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
                         alpha = 0.1, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
                         strict = FALSE, undirected = FALSE)

mmpc_0.05 <- mmpc(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
     alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
     strict = FALSE, undirected = FALSE)

mmpc_0.01 <- mmpc(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
                  alpha = 0.01, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
                  strict = FALSE, undirected = FALSE)

mmpc_0.1 <- mmpc(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
                  alpha = 0.1, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
                  strict = FALSE, undirected = FALSE)

hiton_0.05 <- si.hiton.pc(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
            alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
            strict = FALSE, undirected = FALSE)
hiton_0.01 <- si.hiton.pc(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
                          alpha = 0.01, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
                          strict = FALSE, undirected = FALSE)
hiton_0.1 <- si.hiton.pc(x = dataS, cluster = NULL, whitelist = NULL, blacklist = bl, test = NULL,
                          alpha = 0.1, B = NULL, max.sx = NULL, debug = FALSE, optimized = FALSE,
                          strict = FALSE, undirected = FALSE)



DBNList <- list("pc_0.01" = pc_0.01, "pc_0.05" = pc_0.05, "pc_0.1" = pc_0.1,
                "gs_0.01" = gs_0.01, "gs_0.05" = gs_0.05, "gs_0.1" = gs_0.1,
                "iamb_0.01" = iamb_0.01, "iamb_0.05" = iamb_0.05, "iamb_0.1" = iamb_0.1,
                "iiamb_0.01" = iiamb_0.01, "iiamb_0.05" = iiamb_0.05, "iiamb_0.1" = iiamb_0.1,
                "mmpc_0.01" = mmpc_0.01, "mmpc_0.05" = mmpc_0.05, "mmpc_0.1" = mmpc_0.1,
                "real" = realDBN)
GenerateResults(dataS, DBNList, realDBN)


