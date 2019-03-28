# Yeast cycle data
yeast_data <- readxl::read_xls("Yeast_data.xls")

library(dplyr)
library(readr)
yeast_data <- yeast_data %>% 
  t %>%
  as_tibble %>%
  type_convert

colnames(yeast_data) <- yeast_data[1,]
yeast_data <- yeast_data[-1,]
yeast_data <- as.data.frame(yeast_data)
yeast_data[] <- lapply(yeast_data, as.numeric)

##Apply G1DBN as per the authors instructions...
alpha1 = 0.1 
alpha2 = 0.0059
predictors = c("ACE2","FKH1","FKH2", "GAT3" ,"MBP1", "MCM1", "MIG2", "NDD1", "PHD1", "RAP1", "RME1", "STB1", "SUT1", "SWI4", 
               "SWI5", "SWI6", "TEC1", "YOX1")
pred_pos <- which(colnames(yeast_data) %in% predictors) 
G1_yeast  = G1DBN::DBNScoreStep1(yeast_data, predPosition = pred_pos)
G1_str <- BuildEdges(G1_yeast$S1ls,threshold = alpha1, predNames = as.factor(predictors), targetNames = as.factor(colnames(yeast_data)))


G2_yeast <- G1DBN::DBNScoreStep2(S1 = G1_yeast$S1ls, alpha1 = alpha1, predPosition = pred_pos,data = yeast_data)
G2_str <- BuildEdges(G2_yeast,threshold = alpha2, predNames = predictors, targetNames = colnames(yeast_data))
G2_str <- as.data.frame(G2_str, stringAsFactors = TRUE)
G2_arcs <- G2_str[,c(1,2)]
names(G2_arcs) <- c("from", "to")
common <- inner_join(G2_arcs, yeast_arcs_validated)


## add more arcs atleast 500
G2_str_500 <- BuildEdges(G2_yeast,threshold = 1, predNames = predictors, targetNames = colnames(yeast_data))
G2_str_500 <- as.data.frame(G2_str_500, stringAsFactors = TRUE)

G2_arcs_500 <- G2_str_500[,c(1,2)]
names(G2_arcs_500) <- c("from", "to")
common <- inner_join(G2_arcs_500, yeast_arcs_validated)

##Some tweaking because I have to fit it in my own code :(

data <- yeast_data
colnames(data) <- paste0("V", 1:length(data))

# Score_mat is the G1DBN score put into full 786*786 matrix
score_mat <- matrix(NA, nrow = ncol(yeast_data), ncol = ncol(yeast_data))
for(col in 1:length(pred_pos))
{
  score_mat[, pred_pos[col]] <- G1_yeast$S1ls[,col]
  
}
bl <- CreateBlackList(data)
pred_names <- paste0("V", pred_pos)
targetNames <- paste0("V", 1:ncol(data))
onlyTargets <- targetNames[-pred_pos]
bl_ <- expand.grid(onlyTargets, targetNames)
colnames(bl_) <- colnames(bl)
bl <- rbind(bl,bl_)

# impute the data, by using the value in previous time slice.
yeast_data_imputed <- yeast_data
  for(i in 1:nrow(yeast_data_imputed))
  {
    for(j in 1:ncol(yeast_data_imputed))
    {
      if(is.na(yeast_data_imputed[i,j]) && i == 1)
      {
        yeast_data_imputed[i,j] <- yeast_data_imputed[i+1,j]
      }
      if(is.na(yeast_data_imputed[i,j]) && i!=1)
      {
        yeast_data_imputed[i,j] <- yeast_data_imputed[i-1,j]
      }
    }
  }
dataS = ShiftData(yeast_data_imputed)
yeast_arcs_validated <- xlsx::read.xlsx2("ValidattedEdges.xls", sheetName = "ValidatedEdges")

## Exhaustive search using BIC only
empty_weights <- matrix(data = 1, nrow = nrow(score_mat), ncol = ncol(score_mat))
all_parents <-  pred_names
full_yeast_score_bics <- FullExhaustiveSearch_selectedParents(dataS, all_parents)
yeast_bic_and_zeros <- do.call(rbind,full_yeast_score_bics)
yeast_all_parent_combinations <- MakeAllPossibleCombination(all_parents)
yeast_Max_bic <- ConstructBNUsingMaxScoringParents(yeast_bic_and_zeros, yeast_all_parent_combinations, names(dataS))

yeast_max_bic_arcs <- MakeArcSet(yeast_Max_bic,yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_max_bic_arcs, yeast_arcs_validated)

## BIC + G1DBN. Adjustment = 0.
yeast_all_targets_weights <- CalculateWeightsOfAllTargets( yeast_all_parent_combinations, genes = nrow(score_mat), score_mat)
yeast_bic_and_weights <- CombineBICandWeight(yeast_all_targets_weights, full_yeast_score_bics ,adjustment = 0, multiplier = -1, all_parent_combinations = yeast_all_parent_combinations)
yeast_Max_BN_0 <- ConstructBNUsingMaxScoringParents(yeast_bic_and_weights, yeast_all_parent_combinations, names(dataS))

yeast_Max_BN_0_arcs <- MakeArcSet(yeast_Max_BN_0, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_Max_BN_0_arcs, yeast_arcs_validated)


## BIC + G1DBN. Adjustment = 2.0.
#yeast_all_targets_weights <- CalculateWeightsOfAllTargets( yeast_all_parent_combinations, genes = nrow(score_mat), score_mat)
yeast_bic_and_weights_2 <- CombineBICandWeight(yeast_all_targets_weights, full_yeast_score_bics ,adjustment = 2, multiplier = -1, all_parent_combinations = yeast_all_parent_combinations)
yeast_Max_BN_2 <- ConstructBNUsingMaxScoringParents(yeast_bic_and_weights_2, yeast_all_parent_combinations, names(dataS))

yeast_Max_BN_2_arcs <- MakeArcSet(yeast_Max_BN_2, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_Max_BN_2_arcs, yeast_arcs_validated)

## BIC + Lasso  
lasso <-   ApplyLars_selectedParents(yeast_data_imputed, pred_pos)
lasso_big <- matrix(NA, nrow = ncol(yeast_data), ncol = ncol(yeast_data))
for(col in 1:length(pred_pos))
{
  lasso_big[, pred_pos[col]] <- lasso[,col]
  
}

yeast_all_target_lasso <- CalculateWeightsOfAllTargets( yeast_all_parent_combinations, genes = nrow(lasso_big), lasso_big)#,  adjustment = 0, negative = FALSE)
yeast_bic_and_lasso <- CombineBICandWeight(yeast_all_target_lasso, full_yeast_score_bics, adjustment = 0, multiplier = 1, yeast_all_parent_combinations)
yeast_max_BN_lasso <- ConstructBNUsingMaxScoringParents(yeast_bic_and_lasso, yeast_all_parent_combinations, names(dataS))


yeast_max_BN_lasso_arcs <-  MakeArcSet(yeast_max_BN_lasso, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_max_BN_lasso_arcs, yeast_arcs_validated)

### Greedy Search BIC only
yeast_empty_weights <- matrix(data = 1, nrow = ncol(yeast_data), ncol = ncol(yeast_data))
yeast_gs_bic_only <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = yeast_empty_weights, pred_pos = pred_pos)
yeast_gs_bic_only_arcs <-  MakeArcSet(yeast_gs_bic_only, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_bic_only_arcs, yeast_arcs_validated)

### Greedy Search BIC + G1DBN weight
yeast_gs_weight_adj_0 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = score_mat, multiplier = -1, adjustment = 0, pred_pos = pred_pos)
yeast_gs_weight_adj_arcs <- MakeArcSet(yeast_gs_weight_adj_0, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_weight_adj_arcs, yeast_arcs_validated)

### Greedy Search BIC + G1DBN weight adjustment = 2
yeast_gs_weight_adj_2 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = score_mat, multiplier = -1, adjustment = 2, pred_pos = pred_pos)

yeast_gs_weight_adj_2_arcs <- MakeArcSet(yeast_gs_weight_adj_2, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_weight_adj_2_arcs, yeast_arcs_validated)

###Greedy Search + lasso adj = 0
yeast_gs_lasso_adj_0 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = lasso_big, multiplier = 1, adjustment = 0, pred_pos = pred_pos)
yeast_gs_lasso_adj_0_arcs <- MakeArcSet(yeast_gs_lasso_adj_0, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_lasso_adj_0_arcs, yeast_arcs_validated)

###Greedy Search + lasso adj = 2.5
yeast_gs_lasso_adj_2.5 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = lasso_big, multiplier = 1, adjustment = 2.5, pred_pos = pred_pos)
yeast_gs_lasso_adj_2.5_arcs <- MakeArcSet(yeast_gs_lasso_adj_2.5, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_lasso_adj_2.5_arcs, yeast_arcs_validated)

MakeArcSet <- function(bn, data)
{
  arc_set <- arcs(bn)
  froms <- unlist(lapply(arc_set[,1], function(from) as.numeric(gsub("V", from, replacement = ""))))
  tos <- unlist(lapply(arc_set[,2], function(to) {
    
    to <- gsub("V", to, replacement = "")
    to <- as.numeric(gsub("_1",to, replacement = "" ))}))
  
  froms <- names(data)[froms]
  tos <- names(data)[tos]
  arc_set <- as.data.frame(cbind("from" = froms, "to" = tos ))
  return(arc_set)
  
  
}

PercentageValidatedFromRecovered <- function(arc_set, validated_set)
{
  common <- inner_join(arc_set, validated_set)
  return(nrow(common)/nrow(arc_set) * 100)
  
}