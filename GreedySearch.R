## Greedy Search for max BIC + weight

## Searchs greedily for the max scoring parent. If it increases the current score, add it, 
## then repeat. Stop when
## no other parents increases score
GreedySearchParentsWithAdjustedWeight <- function(data, p = 50, n = 20, adjustment = 0, negative = TRUE, score_mat)
{
  # score_mat <- G1DBN::DBNScoreStep1(data)
  dataS <- ShiftData(data)
  targets <- 1:p
  all_parents <-  names(dataS)[1:p]
  local_bns_arcs <- lapply(targets, function(x)
    
    GreedySearchWithAdjustedWeight_OneTarget(score_mat = score_mat, targetIdx = x, all_parents = all_parents, negative = negative, dataS = dataS, adjustment = adjustment))
  
  fullBN <- do.call(what = rbind, args = local_bns_arcs)
  nodes <- names(dataS)
  bn<- bnlearn::empty.graph(nodes)
  arcs(bn) <- fullBN
  return (bn)
  
  
}

## Inner function
GreedySearchWithAdjustedWeight_OneTarget<- function(score_mat, targetIdx, all_parents, adjustment = 0, negative = TRUE,dataS)
{
  keep_adding = 1
  targetNode <- paste0("V", targetIdx , "_1")
  nodes <- c(all_parents, targetNode)
  
  print(targetNode)
  data_subset <- dataS[,nodes]
  multi  = 1
  if(negative == TRUE)
  {
    multi = -1
  }
  weights <- multi * (log(score_mat) + adjustment)
  weight_vector <-  weights[targetIdx,]
  
  
  bn <- bnlearn::empty.graph(nodes)
 # empty_bn_score <-  score(x  =bn, data = data_subset, type = "bic-g", by.node = TRUE )
  old_bic <-  score(x  =bn, data = data_subset, type = "bic-g",by.node = TRUE )[targetNode]
  while (keep_adding != 0)
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
      print(paste0("Adding Parent: ", best_parent ))
      all_parents[max_idx] <- NA
      old_bic <- bic_and_weight[max_idx]
    }
  }
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
  
  gs_weight_adj_0 <- GreedySearchParentsWithAdjustedWeight(data = df, p = 50, n = 20, score_mat = G1_mat,negative = TRUE, adjustment = 0)
  
  gs_weight_adj_2.5 <- GreedySearchParentsWithAdjustedWeight(data = df, p = 50, n = 20, score_mat = G1_mat,negative = TRUE, adjustment = 2.5)
  
  gs_lasso_adj_0 <- GreedySearchParentsWithAdjustedWeight(data = df, p = 50, n = 20, score_mat = lasso,negative = FALSE, adjustment = 0)
  
  gs_lasso_adj_3 <- GreedySearchParentsWithAdjustedWeight(data = df, p = 50, n = 20, score_mat = lasso,negative = FALSE, adjustment = 3)
  
  DBNList <-  list("gs_bic" = gs_bic_only, "gs_G1DBN_0" = gs_weight_adj_0, "gs_G1DBN_2.5" = gs_weight_adj_2.5,
               "gs_lasso_0" = gs_lasso_adj_0, "gs_lasso_3" = gs_lasso_adj_3)
  #DBNList <- append(DBNList, DBNList1)
  GenerateResults(dataS, DBNList, RealDBN = realDBN)
  
}