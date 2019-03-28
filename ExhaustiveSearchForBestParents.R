
## Full Exhaustive Search to calculate scores of all the all possible sets of parents (max parents = 3)


## This coded combines BIC score with the weights calculated using the G1DBN score and the lasso score.
## For G1DBN score, do BIC_score -log(G1DBN_score)-adjustment
## Adjustment is calculated by looking at the graph
ApplyLars_selectedParents <- function(data, pred_pos)
{
  predictors  = as.matrix(data[1:(nrow(data)-1), pred_pos])
  response = data[2:nrow(data),]
  nR = ncol(response)
  arcs = lapply(1:nR, function(y) 
  {
    print(y)
    Y = response[,y]
    arcs = LASSO(X = predictors, Y = Y)
    
  })
  scores = abs(do.call(rbind, arcs))
  return(scores)
}
FullExhaustiveSearch_selectedParents <- function(dataS, parents, beta = NULL)
{
  allNodes = names(dataS)
  All_possible_parents = parents
  #All_possible_parents = allParents[!grepl("_1", allNodes)]
  Allcombinations = MakeAllPossibleCombination(All_possible_parents)
  All_possible_targets = allNodes[grepl("_1", allNodes)]
  ScoreAllTargets = lapply(All_possible_targets, function(x) ScoreEntirepossibleSetsOfParents(dataS, x, Allcombinations, beta))
  
  
}
FullExhaustiveSearch <- function(dataS, beta = NULL)
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

CalculateWeightsOfOneTarget <- function(G1_mat,all_parent_combinations, target_idx)#, adjustment = 0, negative = TRUE)
{
  # multiplier = 1
  # if (negative)
  # {
  #   multiplier = -1
  #   
  # }
  all_weights <- lapply(all_parent_combinations, function(parents)
  {
    parents_idx <- as.numeric(substring(parents,2))
    weights <- lapply(parents_idx, function(parent) 
      return (G1_mat[target_idx, parent]))
    
    #sum(unlist(lapply(weights, function(x) multiplier * (log(x) + adjustment))))
    sum(unlist(lapply(weights, function(x) log(x))))
  })
  return(all_weights)
}

CalculateWeightsOfAllTargets <- function(all_parent_combinations, genes = 50, G1_mat)#, adjustment = 0, negative = TRUE)
{
  all_targets_weights <- lapply(1:genes, function(x) unlist(CalculateWeightsOfOneTarget(G1_mat, all_parent_combinations, x )))#, adjustment, negative )))
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
    if(!is.na(parents))
    {
      for (par in parents)
      {
        bn <- set.arc(bn, from = par, to = target)
      }
    }
  }
  return (bn)
}
## Combines BIC with weights.
## If only interested in BIC, just set all_combintations_weights = 0
## bic_scores is the return value of FullExhaustiveSearch. Its the BIC scores of 
## each local bayesian network with all the combinations of parent sets. Its the output of calculate weights
## of all target function.
CombineBICandWeight <- function(all_targets_weights, bic_scores, adjustment = 0, multiplier = -1, all_parent_combinations)
{
  lengths <- unlist(lapply(all_parent_combinations,length))
  adjusted <- lapply(all_targets_weights, function(x)
    {
      mapply(function(wt,lts) 
      {
        
      multiplier * (wt + (adjustment * lts))
      }, x,lengths, SIMPLIFY = FALSE)
  })
  adjusted <- lapply(adjusted, unlist)
  
  
  BICandWeight <- t(mapply(function(x,y) {
    x +y
    }
    , bic_scores, adjusted, SIMPLIFY = TRUE))
  BICandWeight 
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
   
  dataS = ShiftData(data$data)
 # scores = FullExhaustiveSearch(dataS) ## Expensive
  save(scores, file = "scores.Rda")
    load(file = "scores.RData")
  
  
  all_parents <-  names(dataS)[1:p]
  all_parent_combinations <- MakeAllPossibleCombination(all_parents)
 
   lasso = ApplyLars(data$data)
  
  ## Only max BIC

  
  bic_and_zeros <- do.call(rbind,scores)
  Max_bic <- ConstructBNUsingMaxScoringParents(bic_and_zeros, all_parent_combinations, names(dataS))
  
  ## BIC + G1DBN. Adjustment = 0.
  all_targets_weights <- CalculateWeightsOfAllTargets( all_parent_combinations, p, G1_mat)
  bic_and_weights <- CombineBICandWeight(all_targets_weights, scores,adjustment = 0, multiplier = -1, all_parent_combinations = all_parent_combinations)
  Max_BN_0 <- ConstructBNUsingMaxScoringParents(bic_and_weights, all_parent_combinations, names(dataS))
  

  ## BIC + G1DBN. adjustment = 2.5
 # all_targets_weights_2.5 <- CalculateWeightsOfAllTargets( all_parent_combinations, p, G1_mat)
  bic_and_weights_2.5 <- CombineBICandWeight(all_targets_weights, scores, adjustment = 2.5, multiplier = -1, all_parent_combinations = all_parent_combinations)
  Max_BN_2.5 <- ConstructBNUsingMaxScoringParents(bic_and_weights_2.5, all_parent_combinations, names(dataS))
  
  ## BIC + G1DBN. Adjustment = 2.0
 # all_targets_weights_2.0 <- CalculateWeightsOfAllTargets( all_parent_combinations, p, G1_mat)
  bic_and_weights_2.0 <- CombineBICandWeight(all_targets_weights, scores,  adjustment = 2.0, multiplier = -1, all_parent_combinations = all_parent_combinations)
  Max_BN_2.0 <- ConstructBNUsingMaxScoringParents(bic_and_weights_2.0, all_parent_combinations, names(dataS))
  
  #
  # BIC + lasso,, Adjustment = 0
  all_target_lasso <- CalculateWeightsOfAllTargets( all_parent_combinations, p, lasso)#,  adjustment = 0, negative = FALSE)
  bic_and_lasso <- CombineBICandWeight(all_target_lasso, scores, adjustment = 0, multiplier = 1, all_parent_combinations)
  max_BN_lasso <- ConstructBNUsingMaxScoringParents(bic_and_lasso, all_parent_combinations, names(dataS))

  ## BIC + lasso, Adjustment = 2.5
 # all_target_lasso_2.5 <- CalculateWeightsOfAllTargets( all_parent_combinations, p, lasso)#,  adjustment = 2.5, negative = FALSE)
  bic_and_lasso_2.5 <- CombineBICandWeight(all_target_lasso, scores, adjustment = 2.5, multiplier = 1, all_parent_combinations)
  max_BN_lasso_2.5 <- ConstructBNUsingMaxScoringParents(bic_and_lasso_2.5, all_parent_combinations, names(dataS))
  
  
  
  ## BIC + lasso, Adjustment = 3.5
 # all_target_lasso_3.5 <- CalculateWeightsOfAllTargets( all_parent_combinations, p, lasso)#,  adjustment = 3.5, negative = FALSE)
  bic_and_lasso_3.5 <- CombineBICandWeight(all_target_lasso, scores, adjustment = 3.5, multiplier = -1, all_parent_combinations)
  max_BN_lasso_3.5 <- ConstructBNUsingMaxScoringParents(bic_and_lasso_2.5, all_parent_combinations, names(dataS))
  
  realDBN <- ConvertToBN(data$RealNet)   
  DBNList <- list("BIC" = Max_bic, "BIC_G1DBN" = Max_BN_0, "BIC_G1DBN_2." = Max_BN_2.0, "BIC_G1DBN_2.5" = Max_BN_2.5,
                  "BIC_lasso" = max_BN_lasso, "BIC_lasso_2.5" = max_BN_lasso_2.5, "BIC_lasso_3.5" = max_BN_lasso_3.5, "Real" = realDBN)
   
   
  GenerateResults(dataS, DBNList, realDBN)
}

## Do it 4 100 datasets, gonna be very time consuming.
seeds = 1:100
p = 50
n = 20
noiseLevel = 10

datasets = lapply(seeds, function(x) SimulateData(p,n,noiseLevel, x ))


dataS_100 <- lapply(datasets, function(x) ShiftData(x$data))
#scores_100 <-  lapply(dataS_100, function(x) FullExhaustiveSearch(x)) ## Expensive
#save(scores_100,file = "scors_100.rds")
load(file = "scors_100.rds")

## Only BIC score 
## This has a problem that it mostly favours graphs with most nodes
bic_and_zeros_100 <- lapply(scores_100, function(x) do.call(rbind,x))
nodes <- names(dataS_100[[1]])
all_parents <-  nodes[1:p]
all_parent_combinations <- MakeAllPossibleCombination(all_parents)

Max_bics_100 <- lapply(bic_and_zeros_100, function(x)
  {
    ConstructBNUsingMaxScoringParents(x, all_parent_combinations, nodes)
  })

results_max_bics_100 <- CalculatePrecisionAndRecallForMultiple(Max_bics_100, datasets)

## BIC + G1DBN
# expensive
#all_targets_weights_100 <- lapply(G1s, function(x) CalculateWeightsOfAllTargets( all_parent_combinations, p, x$mat$S1ls))
#save(file = "G1DBNScoresFor100.rds",all_targets_weights_100)
# BIC+G1DBN+adjustment zero adjustment
load(file = "G1DBNScoresFor100.rds")
bic_and_weights_100 <- mapply(function(x,y) CombineBICandWeight(x, y, adjustment = 0, multiplier = -1, all_parent_combinations), all_targets_weights_100, scores_100, SIMPLIFY = FALSE)
Max_BN_0_100 <- lapply(bic_and_weights_100 , function(x) ConstructBNUsingMaxScoringParents(x, all_parent_combinations, nodes))
CalculatePrecisionAndRecallForMultiple(Max_BN_0_100, datasets)


## BIC + G1DBN, adjustment = 2.5
#all_targets_weights_2.5_100 <- CalculateWeightsOfAllTargets( all_parent_combinations, p, G1_mat,  adjustment = 2.5, negative = TRUE)
bic_and_weights_2.5_100 <- mapply(function(x,y) CombineBICandWeight(x, y, adjustment = 2.5, multiplier = -1, all_parent_combinations), all_targets_weights_100, scores_100, SIMPLIFY = FALSE)
Max_BN_2.5_100 <- lapply(bic_and_weights_2.5_100 , function(x) ConstructBNUsingMaxScoringParents(x, all_parent_combinations, nodes))
CalculatePrecisionAndRecallForMultiple(Max_BN_2.5_100, datasets)

## BIC + G1DBN, adjustment = 2.0
bic_and_weights_2.0_100 <- mapply(function(x,y) CombineBICandWeight(x, y, adjustment = 2.0, multiplier = -1, all_parent_combinations), all_targets_weights_100, scores_100, SIMPLIFY = FALSE)
Max_BN_2.0_100 <- lapply(bic_and_weights_2.0_100 , function(x) ConstructBNUsingMaxScoringParents(x, all_parent_combinations, nodes))
DBNList <- Max_BN_2.0_100
CalculatePrecisionAndRecallForMultiple(DBNList, datasets)


## lasso 
lasso_100 = lapply(datasets, function(x) ApplyLars(x$data))
save(lasso_100, file = "lasso_100.rds")
## BIC + lasso,, Adjustment = 0

all_targets_lasso_100 <- lapply(lasso_100, function(x) CalculateWeightsOfAllTargets( all_parent_combinations, p, x))
save(file = "all_targets_lasso.rds", all_targets_lasso_100)
load(file = "all_targets_lasso.rds")
## BIC + lasso,, Adjustment = 0
bic_and_lasso_100 <- mapply(function(x,y) CombineBICandWeight(x, y, adjustment = 0, multiplier = 1, all_parent_combinations), all_targets_lasso_100, scores_100, SIMPLIFY = FALSE)
max_BN_lasso_100 <- lapply(bic_and_lasso_100 , function(x) ConstructBNUsingMaxScoringParents(x, all_parent_combinations, nodes))
CalculatePrecisionAndRecallForMultiple(max_BN_lasso_100, datasets)
## BIC + lasso,, Adjustment = 2.5

bic_and_lasso_2.5_100 <- mapply(function(x,y) CombineBICandWeight(x, y, adjustment = 2.5, multiplier = 1, all_parent_combinations), all_targets_lasso_100, scores_100, SIMPLIFY = FALSE)
max_BN_lasso_2.5_100 <- lapply(bic_and_lasso_2.5_100 , function(x) ConstructBNUsingMaxScoringParents(x, all_parent_combinations, nodes))
CalculatePrecisionAndRecallForMultiple(max_BN_lasso_2.5_100, datasets)
## BIC + lasso,, Adjustment = 3.5

bic_and_lasso_3.5_100 <- mapply(function(x,y) CombineBICandWeight(x, y, adjustment = 3.5, multiplier = 1, all_parent_combinations), all_targets_lasso_100, scores_100, SIMPLIFY = FALSE)
max_BN_lasso_3.5_100 <- lapply(bic_and_lasso_3.5_100 , function(x) ConstructBNUsingMaxScoringParents(x, all_parent_combinations, nodes))
CalculatePrecisionAndRecallForMultiple(max_BN_lasso_3.5_100, datasets)



#