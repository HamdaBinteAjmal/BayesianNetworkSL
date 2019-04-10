## Yeastt Spellman data for all the transcripts now only 18 TFs.
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


validated_edges_full <- xlsx::read.xlsx(file = "Yeast_Spellman_ValidatedEdges_Full.xlsx", sheetIndex = 1)
names(validated_edges_full) <- c("from", "to")
yeast_arcs_validated_35 <- validated_edges_full

yeast_arcs_potential_35 <- xlsx::read.xlsx(file = "Yeast_spellman_full_potential.xlsx", sheetIndex = 1)

## Just checked the validated edges and there are only 37 out of 786 that can be regulators ## so why not 
## use only these 37
potential_parents <- levels(validated_edges_full$from)
potential_parents <- potential_parents[c(-14,-15)]
pred_pos_35 <- which(names(yeast_data) %in% potential_parents)

## Only getting 35 because MSA1 and MSA2 are not present in yeasttract database rather MSB1 and MSB2 are present. 
alpha1 = 0.1 
#alpha2 = 0.0059
G1_yeast_35  = G1DBN::DBNScoreStep1(yeast_data, predPosition = pred_pos_35)
G1_str_35 <- BuildEdges(G1_yeast_35$S1ls,threshold = alpha1, predNames = potential_parents, targetNames = colnames(yeast_data))


G2_yeast_35 <- G1DBN::DBNScoreStep2(S1 = G1_yeast_35$S1ls, alpha1 = alpha1, predPosition = pred_pos_35,data = yeast_data)
G2_str_35 <- BuildEdges(G2_yeast_35,threshold = alpha2, predNames = potential_parents, targetNames = colnames(yeast_data))
G2_str_35 <- as.data.frame(G2_str_35, stringAsFactors = TRUE)
G2_arcs_35 <- G2_str_35[,c(1,2)]

names(G2_arcs_35) <- c("from", "to")
result <- PercentageValidatedFromRecovered(G2_arcs_35, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(G2_arcs_35, yeast_arcs_potential_35)
MakeRandomAndEvaluate(nrow(G2_arcs_35),targetNames, potential_parents, yeast_arcs_validated_35, yeast_arcs_potential_35, result, result_potential)

## just a test ##
G1_temp<- BuildEdges(G1_yeast_35$S1ls,threshold = 1, predNames = potential_parents, targetNames = colnames(yeast_data))
G1_temp <- as.data.frame(G1_temp, stringAsFactors = TRUE)
G1_temp <- G1_temp[1:2500,c(1,2)]
names(G1_temp) <- c("from", "to")
result <- PercentageValidatedFromRecovered(G1_temp, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(G1_temp, yeast_arcs_potential_35)

###

score_mat_35 <- matrix(NA, nrow = ncol(yeast_data), ncol = ncol(yeast_data))
for(col in 1:length(pred_pos_35))
{
  score_mat_35[, pred_pos_35[col]] <- G1_yeast_35$S1ls[,col]
  
}
## Not do the drill :o 
## Am I doing exhaustive search???? lets try with HC first
yeast_empty_weights <- matrix(data = 1, nrow = ncol(yeast_data), ncol = ncol(yeast_data))
yeast_gs_bic_only_35 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = yeast_empty_weights, pred_pos = pred_pos_35)
yeast_gs_bic_only_arcs_35 <-  MakeArcSet(yeast_gs_bic_only_35, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_bic_only_arcs_35, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(yeast_gs_bic_only_arcs_35, yeast_arcs_potential_35)
MakeRandomAndEvaluate(nrow(yeast_gs_bic_only_arcs_35),targetNames, potential_parents, yeast_arcs_validated_35, yeast_arcs_potential_35, result, result_potential)


## GS + maxP
yeast_gs_bic_only_35_maxP3 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = yeast_empty_weights, pred_pos = pred_pos_35, maxP = 4)
yeast_gs_bic_only_arcs_35_maxP3 <-  MakeArcSet(yeast_gs_bic_only_35_maxP3, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_bic_only_arcs_35_maxP3, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(yeast_gs_bic_only_arcs_35_maxP3, yeast_arcs_potential_35)
MakeRandomAndEvaluate(nrow(yeast_gs_bic_only_arcs_35_maxP3),targetNames, potential_parents, yeast_arcs_validated_35, yeast_arcs_potential_35, result, result_potential)

## GS + G1DBN 
yeast_gs_weight_adj_0_35 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = score_mat_35, multiplier = -1, adjustment = 0, pred_pos = pred_pos_35)
yeast_gs_weight_adj_arcs_35 <- MakeArcSet(yeast_gs_weight_adj_0_35, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_weight_adj_arcs_35, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(yeast_gs_weight_adj_arcs_35, yeast_arcs_potential_35)
MakeRandomAndEvaluate(nrow(yeast_gs_weight_adj_arcs_35),targetNames, potential_parents, yeast_arcs_validated_35, yeast_arcs_potential_35, result, result_potential)

## GS + G1DBN + adj = 2.5
yeast_gs_weight_adj_2.5_35 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = score_mat_35, multiplier = -1, adjustment = 2.5, pred_pos = pred_pos_35)
yeast_gs_weight_adj_2.5_arcs_35 <- MakeArcSet(yeast_gs_weight_adj_2.5_35, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_weight_adj_2.5_arcs_35, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(yeast_gs_weight_adj_2.5_arcs_35, yeast_arcs_potential_35)
MakeRandomAndEvaluate(nrow(yeast_gs_weight_adj_2.5_arcs_35),targetNames, potential_parents, yeast_arcs_validated_35, yeast_arcs_potential_35, result, result_potential)


## GS + G1DBN + adj = 3.5
yeast_gs_weight_adj_3.5_35 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = score_mat_35, multiplier = -1, adjustment = 3.5, pred_pos = pred_pos_35)
yeast_gs_weight_adj_3.5_arcs_35 <- MakeArcSet(yeast_gs_weight_adj_3.5_35, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_weight_adj_3.5_arcs_35, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(yeast_gs_weight_adj_3.5_arcs_35, yeast_arcs_potential_35)


## GS + G1DBN + maxP3
yeast_gs_weight_adj_0_35_maxP3 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = score_mat_35, multiplier = -1, adjustment = 0, pred_pos = pred_pos_35, maxP = 4)
yeast_gs_weight_adj_arcs_35_maxP3 <- MakeArcSet(yeast_gs_weight_adj_0_35_maxP3, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_weight_adj_arcs_35_maxP3, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(yeast_gs_weight_adj_arcs_35_maxP3, yeast_arcs_potential_35)
MakeRandomAndEvaluate(nrow(yeast_gs_weight_adj_2.5_arcs_35_maxP3),targetNames, potential_parents, yeast_arcs_validated_35, yeast_arcs_potential_35, result, result_potential)


## GS + G1DBN + adjustmnet 2.5 + maxP3
yeast_gs_weight_adj_2.5_35_maxP3 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = score_mat_35, multiplier = -1, adjustment = 2.5, pred_pos = pred_pos_35, maxP = 4)
yeast_gs_weight_adj_arcs_2.5_35_maxP3 <- MakeArcSet(yeast_gs_weight_adj_2.5_35_maxP3, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_weight_adj_arcs_2.5_35_maxP3, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(yeast_gs_weight_adj_arcs_2.5_35_maxP3, yeast_arcs_potential_35)





lasso <-   ApplyLars_selectedParents(yeast_data_imputed, pred_pos_35)
lasso_big_35 <- matrix(NA, nrow = ncol(yeast_data), ncol = ncol(yeast_data))
for(col in 1:length(pred_pos_35))
{
  lasso_big_35[, pred_pos_35[col]] <- lasso[,col]
  
}

## GS + lasso
yeast_gs_lasso_adj_0_35 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = lasso_big_35, multiplier = 1, adjustment = 0, pred_pos = pred_pos_35)
yeast_gs_lasso_adj_arcs_35 <- MakeArcSet(yeast_gs_lasso_adj_0_35, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_lasso_adj_arcs_35, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(yeast_gs_lasso_adj_arcs_35, yeast_arcs_potential_35)
MakeRandomAndEvaluate(nrow(yeast_gs_lasso_adj_arcs_35),targetNames, potential_parents, yeast_arcs_validated_35, yeast_arcs_potential_35, result, result_potential)

## GS + lasso + maxP = 3
yeast_gs_lasso_adj_0_35_maxP3 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = lasso_big_35, multiplier = 1, adjustment = 0, pred_pos = pred_pos_35, maxP = 4)
yeast_gs_lasso_adj_arcs_35_maxP3 <- MakeArcSet(yeast_gs_lasso_adj_0_35_maxP3, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_lasso_adj_arcs_35_maxP3, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(yeast_gs_lasso_adj_arcs_35_maxP3, yeast_arcs_potential_35)
MakeRandomAndEvaluate(nrow(yeast_gs_lasso_adj_arcs_35_maxP3),targetNames, potential_parents, yeast_arcs_validated_35, yeast_arcs_potential_35, result, result_potential)


## GS + lasso + adjustmet
yeast_gs_lasso_adj_2_35 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = lasso_big_35, multiplier = 1, adjustment = 2, pred_pos = pred_pos_35)
yeast_gs_lasso_adj_2_arcs_35 <- MakeArcSet(yeast_gs_lasso_adj_2_35, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_lasso_adj_2_arcs_35, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(yeast_gs_lasso_adj_2_arcs_35, yeast_arcs_potential_35)

## GS + lasso + adjustmet = maxP3
yeast_gs_lasso_adj_2_35_maxP3 <- GreedySearchParentsWithAdjustedWeight_selectedParents(data = yeast_data_imputed, p = ncol(yeast_data_imputed), n = nrow(yeast_data_imputed), score_mat = lasso_big_35, multiplier = 1, adjustment = 2, pred_pos = pred_pos_35, maxP = 4)
yeast_gs_lasso_adj_2_arcs_35_maxP3 <- MakeArcSet(yeast_gs_lasso_adj_2_35_maxP3, yeast_data_imputed)
result <- PercentageValidatedFromRecovered(yeast_gs_lasso_adj_2_arcs_35_maxP3, yeast_arcs_validated_35)
result_potential <-  PercentageValidatedFromRecovered(yeast_gs_lasso_adj_2_arcs_35_maxP3, yeast_arcs_potential_35)
