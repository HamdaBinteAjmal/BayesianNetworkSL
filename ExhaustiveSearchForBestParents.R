
## Full Exhaustive Search to calculate scores of all the all possible sets of parents (max parents = 3)


## This coded combines BIC score with the weights calculated using the G1DBN score and the lasso score.
## For G1DBN score, do BIC_score -log(G1DBN_score)-adjustment
## Adjustment is calculated by looking at the graphFullExhaustiveSearch <- function(dataS, beta = NULL)
{
  allNodes = names(dataS)
  All_possible_parents = allNodes[!grepl("_1", allNodes)]
  Allcombinations = MakeAllPossibleCombination(All_possible_parents)
  All_possible_targets = allNodes[grepl("_1", allNodes)]
  ScoreAllTargets = lapply(All_possible_targets, function(x) ScoreEntirepossibleSetsOfParents(dataS, x, Allcombinations, beta))
  
  
}
## Makes all possible combinates of parents (max = 3)
MakeAllPossibleCombination <- function(all_parents)
{
  combinations0 = combn(all_parents, 0)
  combinations1 = combn(all_parents, 1) 
  combinations2 = combn(all_parents,2)
  combinations3 = combn(all_parents, 3)
  
  parents_0 <- c(NA)
  parents_1 <- lapply(1:ncol(combinations1), function(x) return(c(combinations1[x])))
  parents_2 <- lapply(1:ncol(combinations2), function(x) return(c(combinations2[1,x], combinations2[2,x])))
  parents_3 <- lapply(1:ncol(combinations3), function(x) return(c(combinations3[1,x], combinations3[2,x], combinations3[3,x])))
  
  all_parents <- c(parents_0, parents_1, parents_2, parents_3)
  return (all_parents)
}
## Score all possible sets of parents for one target Node and returns a vector of the scores. 
# This is an expensive method
ScoreEntirepossibleSetsOfParents <- function(dataS, tarNode, all_parents, beta)
{
  print(tarNode)
  
  tarIndex = sub("V", "", tarNode)
  tarIndex = as.numeric(sub("_1", "", tarIndex))
  allScores <- array(data = NA, dim = length(all_parents))
  #  beta = beta[which(beta$to == tarNode),]
  allNodes = names(dataS)
  All_possible_parents = allNodes[!grepl("_1", allNodes)]
  
  bn = empty.graph(nodes = c(All_possible_parents, tarNode))
  ScoreOfGraph = score(bn, dataS[,c(nodes(bn))], type = "bic-g", by.node = TRUE)[tarNode]
  
  allScores[1] = ScoreOfGraph
  
  for (i in 2:length(all_parents))
  {
    bn = empty.graph(nodes = c(All_possible_parents, tarNode))
    
    parents = all_parents[[i]]
    for (j in 1:length(parents))
    {
      if(!is.na(parents[j]))
      {
        bn <- set.arc(bn, from = parents[j], to = tarNode)
      }
      
    }
    
    newScore = score(bn, dataS[,c(nodes(bn))], type = "bic-g", by.node = TRUE)[tarNode]
    if(newScore == -Inf)
      print("!blah")
    
    allScores[i] = newScore 
    
  }
  return (allScores)
}


## Inner function , negative = TRUE if G1DBN, else for LASSO, negative = false
## This takes all parent combinations as input and returns a list of the calculated weights for each parent combination

CalculateWeightsOfOneTarget <- function(G1_mat,all_parent_combinations, target_idx, adjustment = 0, negative = TRUE)
{
  multiplier = 1
  if (negative)
  {
    multiplier = -1
    
  }
  all_weights <- lapply(all_parent_combinations, function(parents)
  {
    parents_idx <- as.numeric(substring(parents,2))
    weights <- lapply(parents_idx, function(parent) 
      return (G1_mat[target_idx, parent]))
    
    sum(unlist(lapply(weights, function(x) multiplier * (log(x) + adjustment))))
  })
  return(all_weights)
}

CalculateWeightsOfAllTargets <- function(all_parent_combinations, genes = 50, G1_mat, adjustment = 0, negative = TRUE)
{
  all_targets_weights <- lapply(1:genes, function(x) unlist(CalculateWeightsOfOneTarget(G1_mat, all_parent_combinations, x , adjustment, negative )))
}
## Takes input a list of all parent combinatios, and the id of the max scoring parents sets for
## each target node. Combine all into one big BN
ConstructBNUsingMaxScoringParents <- function(BICandWeight,  all_parent_combinations, allNodes)
{
  max_ids <- apply(BICandWeight, 1, which.max)
  #allNodes = names(dataS)
  All_possible_targets = allNodes[grepl("_1", allNodes)]
  max_parents <- lapply(max_ids, function(x) unlist(all_parent_combinations[x]))
  bn <- empty.graph(allNodes)
  for (i in 1:length(All_possible_targets))
    
  {
    target <- All_possible_targets[[i]]
    parents <- max_parents[[i]]
    
    for (par in parents)
    {
      bn <- set.arc(bn, from = par, to = target)
    }
  }
  return (bn)
}
## Combines BIC with weights.
## If only interested in BIC, just set all_combintations_weights = 0
## bic_scores is the return value of FullExhaustiveSearch. Its the BIC scores of 
## each local bayesian network with all the combinations of parent sets. Its the output of calculate weights
## of all target function.
CombineBICandWeight <- function(all_targets_weights, bic_scores)
{
  BICandWeight <- t(mapply(function(x,y) x+y, bic_scores, all_targets_weights, SIMPLIFY = TRUE))
  
}

Main_ExhaustiveSearchForBestParentSets <- function()
{
  source("ScriptForComparison.R")
  set.seed(1)
  p = 50
  n = 20
  
  data <- SimulateData(Genes = 50, timepoints = 20,seed = 1)
  G1_weights = ExecuteG1DBNS1(data$data, 0.7)
  G1_mat = G1_weights$mat$S1ls
  
  blacklist = CreateBlackList(data = data$data)
  dataS = ShiftData(data$data)
  # scores = FullExhaustiveSearch(dataS) ## Expensive
  save(scores, file = "scores.RData")
  load(file = "scores.RData")
  
  
  all_parents <-  names(dataS)[1:p]
  all_parent_combinations <- MakeAllPossibleCombination(all_parents)
  G1_weights = ExecuteG1DBNS1(data$data, 0.7)
  G1_mat = G1_weights$mat$S1ls
  lasso = ApplyLars(data$data)
  
  ## Only max BIC
  Max_bic <- ConstructBNUsingMaxScoringParents(scores, all_parent_combinations, names(dataS))
  
  ## BIC + G1DBN. Adjustment = 0.
  all_targets_weights <- CalculateWeightsOfAllTargets( all_parent_combinations, genes, G1_mat,  adjustment = 0, negative = TRUE)
  bic_and_weights <- CombineBICandWeight(all_targets_weights, scores)
  Max_BN <- ConstructBNUsingMaxScoringParents(bic_and_weights, all_parent_combinations, names(dataS))
  

  ## BIC + G1DBN. adjustment = 2.5
  all_targets_weights_2.5 <- CalculateWeightsOfAllTargets( all_parent_combinations, genes, G1_mat,  adjustment = 2.5, negative = TRUE)
  bic_and_weights_2.5 <- CombineBICandWeight(all_targets_weights_2.5, scores)
  Max_BN_2.5 <- ConstructBNUsingMaxScoringParents(bic_and_weights_2.5, all_parent_combinations, names(dataS))
  
  ## BIC + G1DBN. Adjustment = 2.0
  all_targets_weights_2.0 <- CalculateWeightsOfAllTargets( all_parent_combinations, genes, G1_mat,  adjustment = 2.0, negative = TRUE)
  bic_and_weights_2.0 <- CombineBICandWeight(all_targets_weights_2.0, scores)
  Max_BN_2.0 <- ConstructBNUsingMaxScoringParents(bic_and_weights_2.0, all_parent_combinations, names(dataS))
  
  ## BIC + lasso,, Adjustment = 0
  all_target_lasso <- CalculateWeightsOfAllTargets( all_parent_combinations, p, lasso,  adjustment = 0, negative = FALSE)
  bic_and_lasso <- CombineBICandWeight(all_target_lasso, scores)
  max_BN_lasso <- ConstructBNUsingMaxScoringParents(bic_and_lasso, all_parent_combinations, names(dataS))

  ## BIC + lasso, Adjustment = 2.5
  all_target_lasso_2.5 <- CalculateWeightsOfAllTargets( all_parent_combinations, p, lasso,  adjustment = 2.5, negative = FALSE)
  bic_and_lasso <- CombineBICandWeight(all_target_lasso_2.5, scores)
  max_BN_lasso_2.5 <- ConstructBNUsingMaxScoringParents(bic_and_lasso_2.5, all_parent_combinations, names(dataS))
  
  
  
  ## BIC + lasso, Adjustment = 3.5
  all_target_lasso_3.5 <- CalculateWeightsOfAllTargets( all_parent_combinations, p, lasso,  adjustment = 3.5, negative = FALSE)
  bic_and_lasso <- CombineBICandWeight(all_target_lasso_3.5, scores)
  max_BN_lasso_3.5 <- ConstructBNUsingMaxScoringParents(bic_and_lasso_2.5, all_parent_combinations, names(dataS))
  
  
  
  }
