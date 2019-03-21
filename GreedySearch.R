## Hill Climbinf for max BIC + weight
library(parallel)
## Searchs greedily for the max scoring parent. If it increases the current score, add it, 
## then repeat. Stop when
## no other parents increases score
GreedySearchParentsWithAdjustedWeight <- function(data, p = 50, n = 20, adjustment = 0, multiplier = -1, score_mat, maxP = Inf)
{
  # score_mat <- G1DBN::DBNScoreStep1(data)
  dataS <- ShiftData(data)
  targets <- 1:p
  all_parents <-  names(dataS)[1:p]
  local_bns_arcs <- lapply(targets, function(x)
    
    GreedySearchWithAdjustedWeight_OneTarget(score_mat = score_mat, targetIdx = x, all_parents = all_parents, multiplier = multiplier, dataS = dataS, adjustment = adjustment, maxP = maxP))
  
  fullBN <- do.call(what = rbind, args = local_bns_arcs)
  nodes <- names(dataS)
  bn<- bnlearn::empty.graph(nodes)
  arcs(bn) <- fullBN
  return (bn)
  
  
}

## Inner function
GreedySearchWithAdjustedWeight_OneTarget<- function(score_mat, targetIdx, all_parents, adjustment = 0, multiplier = -1, dataS, maxP = maxP)
{
  keep_adding = 1
  targetNode <- paste0("V", targetIdx , "_1")
  nodes <- c(all_parents, targetNode)
  
 # print(targetNode)
  data_subset <- dataS[,nodes]
  multi  = multiplier
  weights <- multi * (log(score_mat) + adjustment)
  weight_vector <-  weights[targetIdx,]
  
  
  bn <- bnlearn::empty.graph(nodes)
 # empty_bn_score <-  score(x  =bn, data = data_subset, type = "bic-g", by.node = TRUE )
  old_bic <-  score(x  =bn, data = data_subset, type = "bic-g",by.node = TRUE )[targetNode]
  TotalParents <- 0
  while (keep_adding > 0 && TotalParents < maxP)
  {
    
   
     new_bics <- unlist(lapply(all_parents, function(x)
    {
      if(!is.na(x))
      {
  
        bn_1 <- bnlearn::set.arc(bn, from = x, to = targetNode)
        score <- score(x  = bn_1, data = data_subset, type = "bic-g", by.node = TRUE )[targetNode]
      }
      else
      {
        score <- NA
      }

      score
    })) 
    bic_and_weight <- new_bics + weight_vector
    difference <- bic_and_weight - old_bic
    keep_adding <- length(which(difference > 0)) > 0
    
    if(keep_adding)
    {
      max_idx <- which.max(difference)
      best_parent <- all_parents[max_idx]
      bn <- set.arc(bn, from = best_parent, to = targetNode)
      TotalParents <- TotalParents + 1
      #print(paste0("Adding Parent: ", best_parent ))
      all_parents[max_idx] <- NA
      old_bic <- bic_and_weight[max_idx]
    }
    
  }
  #print(TotalParents)
  return (arcs(bn))
  
}


main <- function()
{
  source("ComparisonFunctions.R")
  set.seed(1)
  p = 50
  n = 20
  
  data <- SimulateData(Genes = 50, timepoints = 20,seed = 1)
  G1_weights = ExecuteG1DBNS1(data$data, 0.7)
  G1_mat = G1_weights$mat$S1ls
  
  blacklist = CreateBlackList(data = data$data)
  dataS = ShiftData(data$data)
  df = data$data
  empty_weights <- matrix(data = 1, nrow = 50, ncol = 50)
  gs_bic_only <- GreedySearchParentsWithAdjustedWeight(data = df, p = 50, n = 20, score_mat = empty_weights)
  
  gs_weight_adj_0 <- GreedySearchParentsWithAdjustedWeight(data = df, p = 50, n = 20, score_mat = G1_mat, multiplier = -1, adjustment = 0)
  gs_weight_adj_2.0 <- GreedySearchParentsWithAdjustedWeight(data = df, p = 50, n = 20, score_mat = G1_mat, multiplier = -1, adjustment = 2.0)
  
  gs_weight_adj_2.5 <- GreedySearchParentsWithAdjustedWeight(data = df, p = 50, n = 20, score_mat = G1_mat, multiplier = -1, adjustment = 2.5)
  
  gs_lasso_adj_0 <- GreedySearchParentsWithAdjustedWeight(data = df, p = 50, n = 20, score_mat = lasso, multiplier = 1, adjustment = 0)
  gs_lasso_adj_2.5 <- GreedySearchParentsWithAdjustedWeight(data = df, p = 50, n = 20, score_mat = lasso, multiplier = 1, adjustment = 2.5)
  gs_lasso_adj_3.5 <- GreedySearchParentsWithAdjustedWeight(data = df, p = 50, n = 20, score_mat = lasso, multiplier = 1, adjustment = 3.5)
  
  DBNList <-  list("gs_bic" = gs_bic_only, "gs_G1DBN_0" = gs_weight_adj_0, "gs_G1DBN_2.0" = gs_weight_adj_2.0,
                   "gs_G1DBN_2.5" = gs_weight_adj_2.5,
               "gs_lasso_0" = gs_lasso_adj_0
               , "gs_lasso_2.5" = gs_lasso_adj_2.5, "gs_lasso_3.5" = gs_lasso_adj_3.5)
  #DBNList <- append(DBNList, DBNList1)
  GenerateResults(dataS, DBNList, RealDBN = realDBN)
  
  ## Do it 4 100 datasets, gonna be very time consuming.
  seeds = 1:100
  p = 50
  n = 20
  noiseLevel = 10
 # load(file = "scors_100.rds")
  
  #
  #load(file = "G1DBNScoresFor100.rds")
  #G1s = lapply(datasets, function(x) ExecuteG1DBNS1(x$data, alpha1))
  #save(G1s, file = "G1s.rds")
  load("G1s.rds")
  datasets = lapply(seeds, function(x) SimulateData(p,n,noiseLevel, x ))
  ## not using mapply here becasue empty_weights is essentially a matrix with all 1s so not needed different for 100 bns
  gs_bic_100 <- mclapply(datasets, function(x)  GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = empty_weights) )
  CalculatePrecisionAndRecallForMultiple(gs_bic_100, datasets)
  

  gs_weight_0_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y$mat$S1ls, multiplier = -1, adjustment = 0), datasets, G1s, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple(gs_weight_0_100, datasets)
  
  
  
  gs_weight_2_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y$mat$S1ls, multiplier = -1, adjustment = 2), datasets, G1s, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple(gs_weight_2_100, datasets)
  
  
  
  gs_weight_2.5_100 <-  mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y$mat$S1ls, multiplier = -1, adjustment = 2.5), datasets, G1s, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple(gs_weight_2.5_100, datasets)
  
  load(file = "lasso_100.rds")
  
  gs_lasso_adj_0_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y, multiplier = 1, adjustment = 0), datasets, lasso_100, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple( gs_lasso_adj_0_100, datasets)
  
  gs_lasso_adj_2.5_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y, multiplier = 1, adjustment = 2.5), datasets, lasso_100, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple( gs_lasso_adj_2.5_100, datasets)
  
  gs_lasso_adj_3.5_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y, multiplier = 1, adjustment = 3.5), datasets, lasso_100, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple(gs_lasso_adj_3.5_100, datasets)
  
  
  ### Very important, limit max number of parents in GreedySearch to improve precision. Lets say 5
  gs_bic_maxP_100 <- mclapply(datasets, function(x)  GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = empty_weights,maxP = 5) )
  CalculatePrecisionAndRecallForMultiple(gs_bic_maxP_100, datasets)
  
  gs_weight_maxP_0_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y$mat$S1ls, multiplier = -1, adjustment = 0, maxP = 5), datasets, G1s, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple(gs_weight_maxP_0_100, datasets)
  
  
  
  gs_weight_maxP_2_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y$mat$S1ls, multiplier = -1, adjustment = 2, maxP = 5), datasets, G1s, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple(gs_weight_maxP_2_100, datasets)
  
  
  
  gs_weight_maxP_2.5_100 <-  mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y$mat$S1ls, multiplier = -1, adjustment = 2.5, maxP = 5), datasets, G1s, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple(gs_weight_maxP_2.5_100, datasets)
  
  
  ## MaxP = 3
  gs_bic_maxP3_100 <- mclapply(datasets, function(x)  GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = empty_weights,maxP = 3) )
  CalculatePrecisionAndRecallForMultiple(gs_bic_maxP3_100, datasets)
  
  gs_weight_maxP3_0_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y$mat$S1ls, multiplier = -1, adjustment = 0, maxP = 3), datasets, G1s, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple(gs_weight_maxP3_0_100, datasets)
  
  
  
  gs_weight_maxP3_2_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y$mat$S1ls, multiplier = -1, adjustment = 2, maxP = 3), datasets, G1s, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple(gs_weight_maxP3_2_100, datasets)
  
  
  gs_weight_maxP3_2.5_100 <-  mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y$mat$S1ls, multiplier = -1, adjustment = 2.5, maxP = 3), datasets, G1s, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple(gs_weight_maxP3_2.5_100, datasets)
  
  ## MAx p  = 3 for lasss
  
  load(file = "lasso_100.rds")
  
  gs_lasso_adj_MaxP3_0_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y, multiplier = 1, adjustment = 0, maxP = 3), datasets, lasso_100, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple( gs_lasso_adj_MaxP3_0_100, datasets)
  
  gs_lasso_adj_MaxP3_2.5_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y, multiplier = 1, adjustment = 2.5, maxP = 3), datasets, lasso_100, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple( gs_lasso_adj_MaxP3_2.5_100, datasets)
  
  gs_lasso_adj_MaxP3_3.5_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y, multiplier = 1, adjustment = 3.5, maxP = 3), datasets, lasso_100, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple(gs_lasso_adj_MaxP3_3.5_100, datasets)
  
  ## Maxp = 5 for lasso
  
  gs_lasso_adj_MaxP5_0_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y, multiplier = 1, adjustment = 0, maxP = 5), datasets, lasso_100, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple( gs_lasso_adj_MaxP5_0_100, datasets)
  
  gs_lasso_adj_MaxP5_2.5_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y, multiplier = 1, adjustment = 2.5, maxP = 5), datasets, lasso_100, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple( gs_lasso_adj_MaxP5_2.5_100, datasets)
  
  gs_lasso_adj_MaxP5_3.5_100 <- mcmapply(function(x,y) GreedySearchParentsWithAdjustedWeight(data = x$data, p = 50, n = 20, score_mat = y, multiplier = 1, adjustment = 3.5, maxP = 5), datasets, lasso_100, SIMPLIFY = FALSE)
  CalculatePrecisionAndRecallForMultiple(gs_lasso_adj_MaxP5_3.5_100, datasets)
  
  
}