## Greedy Search for max BIC + weight

## Searchs greedily for the max scoring parent. If it increases the current score, add it, 
## then repeat. Stop when
## no other parents increases score
GreedySearchParentsWithAdjustedWeight <- function(data, p = 50, n = 20, adjustment = 0, negative = TRUE, score_mat)
{
  # score_mat <- G1DBN::DBNScoreStep1(data)
  dataS <- ShiftData(data)
  targets <- 1:p
  all_parents <-  names(dataS)[1:genes]
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
  
  while (keep_adding != 0)
  {
    old_bic <-  score(x  =bn, data = data_subset, type = "bic-g",by.node = TRUE )[targetNode]
    new_bics <- unlist(lapply(all_parents, function(x)
    {
      
      if(!is.na(x))
      {
        bn_1 <- bnlearn::set.arc(bn, from = x, to = targetNode)
        score <- score(x  =bn_1, data = data_subset, type = "bic-g", by.node = TRUE )[targetNode]
      }
      else
      {
        score <- NA
      }
      scor
      
      bic_and_weight <- new_bics+weight_vector
      differencee
    })) <- bic_and_weight - old_bic
    keep_adding <- length(which(difference > 0))
    
    if(keep_adding > 0)
    {
      max_idx <- which.max(difference)
      best_parent <- all_parents[max_idx]
      bn <- set.arc(bn, from = best_parent, to = targetNode)
      print(paste0("Adding: ", best_parent ))
      all_parents[max_idx] <- NA
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
  #scores = FullExhaustiveSearch(dataS) ##Expensive
  save(scores, file = "scores.RData")
  load(file = "scores.RData")
  empty_weights <- matrix(data = 0, nrow = 50, ncol = 50)
  gs_bic_only <- GreedySearchParentsWithAdjustedWeight(data = data, p = 50, n = 20, score_mat = empty_weights)
  
  gs_weight_adj_0 <- GreedySearchParentsWithAdjustedWeight(data = data, p = 50, n = 20, score_mat = score_mat,negative = TRUE, adjustment = 0))
  
  gs_weight_adj_2.5 <- GreedySearchParentsWithAdjustedWeight(data = data, p = 50, n = 20, score_mat = score_mat,negative = TRUE, adjustment = 2.5)
  
  gs_lasso_adj_0 <- GreedySearchParentsWithAdjustedWeight(data = data, p = 50, n = 20, score_mat = lasso,negative = FALSE, adjustment = 0)
  
  gs_lasso_adj_3 <- GreedySearchParentsWithAdjustedWeight(data = data, p = 50, n = 20, score_mat = lasso,negative = FALSE, adjustment = 3)
  
  
  
  
}