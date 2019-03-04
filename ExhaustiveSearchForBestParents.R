
## Full Exhaustive Search to calculate scores of all the all possible sets of parents (max parents = 3)
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
