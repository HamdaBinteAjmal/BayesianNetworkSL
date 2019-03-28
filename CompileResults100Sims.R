## Compile Results
library(e1071)
CalculatePrecisionAndRecallForMultiple <- function(DBNList, datasets)
{
  dbns <- lapply(DBNList, function(x) GetAMAT(amat(x)))
  reals <- lapply(datasets, function(x) t(x$RealNet$AdjMatrix))
  realdbns <- lapply(datasets, function (x) ConvertToBN(x$RealNet))
  confs_mats <- mapply(function(dbn, real)
  {
    confusionMatrix( reference = factor(real), data = factor(dbn), positive = "1")
    
  }, dbns, reals, SIMPLIFY = FALSE)
  precisions = unlist(lapply(confs_mats, function(x) x$byClass[3]))
  recalls = unlist(lapply(confs_mats, function(x) x$byClass[1]))
  balanced_Acc <- unlist(lapply(confs_mats, function(x) x$byClass[11]))
  hammings = unlist(mapply(function(x,y) hamming(learned = x, true = y), DBNList, realdbns, SIMPLIFY = FALSE))
  # hammings = unlist(mapply(function(dbn, real) 
  #   {
  #     length(which(dbn&real) == FALSE)
  #   
  #   }, dbns, reals, SIMPLIFY = FALSE))
  # 
  TN <- unlist(lapply(confs_mats, function(x) x$table[1,1]))
  FN <- unlist(lapply(confs_mats, function(x) x$table[1,2]))
  FP <- unlist(lapply(confs_mats, function(x) x$table[2,1]))
  TP <- unlist(lapply(confs_mats, function(x) x$table[2,2]))
  
  precision_mean <- mean(precisions, na.rm = TRUE)
  precision_sd <- sd(precisions,  na.rm = TRUE)
  recall_mean <- mean(recalls, na.rm = TRUE)
  recall_sd <- sd(recalls, na.rm = TRUE)
  hamming_mean <- mean(hammings, na.rm = TRUE)
  hamming_sd <- sd(hammings, na.rm = TRUE)
  TN_mean <- mean(TN, na.rm = TRUE)
  TN_sd <- sd(TN, na.rm = TRUE)
  FN_mean <- mean(FN, na.rm = TRUE)
  FN_sd <- sd(FN, na.rm = TRUE)
  FP_mean <- mean(FP, na.rm = TRUE)
  FP_sd <- sd(FP, na.rm = TRUE)
  TP_mean <- mean(TP, na.rm = TRUE)
  TP_sd <- sd(TP, na.rm = TRUE)
  acc_mean <-  mean(balanced_Acc, na.rm = TRUE)
  acc_sd <- sd(balanced_Acc, na.rm = TRUE)
  return (list("precision_mean" = precision_mean, "precision_sd" = precision_sd,
               "recall_mean" = recall_mean, "recall_sd" = recall_sd,
               "hamming_mean" = hamming_mean, "hamming_sd" = hamming_sd,
               "TN_mean" = TN_mean, "TN_sd" = TN_sd,
               "FN_mean" = FN_mean, "FN_sd" = FN_sd,
               "FP_mean" = FP_mean, "FP_sd" = FP_sd,
               "TP_mean" = TP_mean, "TP_sd" = TP_sd,
               "Accuracy" = acc_mean, "Accuracy_sd" = acc_sd))
  
  
}
names <- c("G1DBN_0.1", "G1DBN_0.3", "G1DBN_0.7",
           "Lasso_0.85", "Lasso_0.90", "Lasso_0.95",
           
           "Ext_BIC", "Ext_BIC_G1DBN_0", "Ext_BIC_G1DBN_2.0", "Ext_BIC_G1DBN_2.5",
           "Ext_BIC_Lasso_0", "Ext_BIC_Lasso_2.5", "Ext_BIC_Lasso_3.5",
           
           "HC_BIC", "HC_BIC_G1DBN_0", "HC_BIC_G1DBN_2.0", "HC_BIC_G1DBN_2.5",
           "HC_BIC_Lasso_0", "HC_BIC_Lasso_2.5", "HC_BIC_Lasso_3.5",
           
           "HC_BIC_MaxP3", "HC_BIC_G1DBN_0_MaxP3", "HC_BIC_G1DBN_2.0_MaxP3", "HC_BIC_G1DBN_2.5_MaxP3",
           "HC_BIC_Lasso_0_MaxP3", "HC_BIC_Lasso_2.5_MaxP3", "HC_BIC_Lasso_3.5_MaxP3",
           
           "HC_BIC_MaxP5", "HC_BIC_G1DBN_0_MaxP5", "HC_BIC_G1DBN_2.0_MaxP5", "HC_BIC_G1DBN_2.5_MaxP5",
           "HC_BIC_Lasso_0_MaxP5", "HC_BIC_Lasso_2.5_MaxP5", "HC_BIC_Lasso_3.5_MaxP5",
           
           "PC_0.01", "PC_0.05", "PC_0.1",
           
           "GrowShrink_0.01", "GrowShrink_0.05", "GrowShrink_0.1",
           
           "IAMB_0.01", "IAMB_0.05", "IAMB_0.1",
           
           "Fast_IAMB_0.01", "Fast_IAMB_0.05", "Fast_IAMB_0.1",
           
           "MMPC_0.01", "MMPC_0.05", "MMPC_0.1",
           
           "SI_HITON_0.01", "SI_HITON_0.05", "SI_HITON_0.1",
           
           "HC_BIC_MaxP3", "TABU_BIC_MaxP3", 
           
           "MMHC_0.9_MaxP3", "MMHC_0.7_MaxP3", "MMHC_0.1_MaxP3", "MMHC_0.01_MaxP3",
           "RsMax2_0.9", "RsMax2_0.7", "RsMax2_0.1", "RsMax2_0.01",
           
           "Local_Markov_GrowShrink_0.001", "Local_Markov_GrowShrink_0.005","Local_Markov_GrowShrink_0.01", "Local_Markov_GrowShrink_0.05", "Local_Markov_GrowShrink_0.1",
           
           "DelDBN_GrowShrink_0.001", "DelDBN_GrowShrink_0.005", "DelDBN_GrowShrink_0.01", "DelDBN_GrowShrink_0.05", "DelDBN_GrowShrink_0.1" 
         )
OneBigList <- list("G1DBN_0.1" = G1DBN_0.1_100,
                   "G1DBN_0.3" = G1DBN_0.3_100,
                   "G1DBN_0.7" = G1DBN_0.7_100,
                  
                   "Lasso_0.85" = lasso_BN_0.85_100,
                   "Lasso_0.90" = lasso_BN_0.9_100,
                   "Lasso_0.95" = lasso_BN_0.95_100,
     
                   
                   "Ext_BIC" = Max_bics_100,
                   "Ext_BIC_G1DBN_0" = Max_BN_0_100,
                   "Ext_BIC_G1DBN_2.0" = Max_BN_2.0_100, 
                   "Ext_BIC_G1DBN_2.5"= Max_BN_2.5_100,
                  "Ext_BIC_Lasso_0" = max_BN_lasso_100,
                  "Ext_BIC_Lasso_2.5" = max_BN_lasso_2.5_100,
                  "Ext_BIC_Lasso_3.5" = max_BN_lasso_3.5_100,
                  
                  
                  "HC_BIC" = gs_bic_100,
                  "HC_BIC_G1DBN_0" = gs_weight_0_100,
                  "HC_BIC_G1DBN_2.0" = gs_weight_2_100,
                  "HC_BIC_G1DBN_2.5" = gs_weight_2.5_100,
                  
                  "HC_BIC_Lasso_0" = gs_lasso_adj_0_100,
                  "HC_BIC_Lasso_2.5" = gs_lasso_adj_2.5_100,
                  "HC_BIC_Lasso_3.5"= gs_lasso_adj_3.5_100,
                  
                  "HC_BIC_MaxP3"= gs_bic_maxP3_100,
                  "HC_BIC_G1DBN_0_MaxP3" = gs_weight_maxP3_0_100,
                  "HC_BIC_G1DBN_2.0_MaxP3" = gs_weight_maxP3_2_100,
                  "HC_BIC_G1DBN_2.5_MaxP3" = gs_weight_maxP3_2.5_100,
                  
                  "HC_BIC_Lasso_0_MaxP3" = gs_lasso_adj_MaxP3_0_100,
                  "HC_BIC_Lasso_2.5_MaxP3" = gs_lasso_adj_MaxP3_2.5_100,
                  "HC_BIC_Lasso_3.5_MaxP3" = gs_lasso_adj_MaxP3_3.5_100,
                   
                  "HC_BIC_MaxP5" = gs_bic_maxP_100,
                  "HC_BIC_G1DBN_0_MaxP5" = gs_weight_maxP_0_100,
                  "HC_BIC_G1DBN_2.0_MaxP5" = gs_weight_maxP_2_100,
                  "HC_BIC_G1DBN_2.5_MaxP5" = gs_weight_maxP_2.5_100,
                  
                  "HC_BIC_Lasso_0_MaxP5" =gs_lasso_adj_MaxP5_0_100,
                  "HC_BIC_Lasso_2.5_MaxP5" = gs_lasso_adj_MaxP5_2.5_100,
                  "HC_BIC_Lasso_3.5_MaxP5" = gs_lasso_adj_MaxP5_3.5_100,
                   
                  "PC_0.01" = pc_0.01_100,
                  "PC_0.05" = pc_0.05_100,
                  "PC_0.1" = pc_0.1_100,
                   
                  "GrowShrink_0.01" = gs_0.01_100,
                  "GrowShrink_0.05" = gs_0.05_100,
                  "GrowShrink_0.1" = gs_0.1_100,
                   
                  "IAMB_0.01" = iamb_0.01_100,
                  "IAMB_0.05" = iamb_0.05_100,
                  "IAMB_0.1" = iamb_0.1_100,
                   
                  "Fast_IAMB_0.01" = fiamb_0.01_100,
                  "Fast_IAMB_0.05" = fiamb_0.05_100,
                  "Fast_IAMB_0.1" = fiamb_0.1_100,
                  
                  "MMPC_0.01" = mmpc_0.01_100,
                  "MMPC_0.05" = mmpc_0.05_100,
                  "MMPC_0.1" = mmpc_0.1_100,
                  
                  "SI_HITON_0.01" = hiton_0.01_100,
                  "SI_HITON_0.05" = hiton_0.05_100,
                  "SI_HITON_0.1" = hiton_0.1_100,
                  
                  "HC_BIC_MaxP3_bnlearn" = hc_maxP_3_100,
                  "TABU_BIC_MaxP3" = tabu_maxP_3_100,
                  
                  "MMHC_0.9_MaxP3" = mmhc_bic_hc_0.9_max3_100,
                  "MMHC_0.7_MaxP3" = mmhc_bic_hc_0.7_max3_100,
                  "MMHC_0.1_MaxP3" = mmhc_bic_hc_0.1_max3_100,
                  "MMHC_0.01_MaxP3" = mmhc_bic_hc_0.01_max3_100,
                  
                  "RsMax2_0.9" = rs_hc_0.9_100,
                  "RsMax2_0.7" = rs_hc_0.7_100,
                  "RsMax2_0.1" = rs_hc_0.1_100,
                  "RsMax2_0.01" = rs_hc_0.01_100,
                  
                  "Local_Markov_GrowShrink_0.001" = delDBN_0.001_s_100,
                  "Local_Markov_GrowShrink_0.005" = delDBN_0.005_s_100,
                  "Local_Markov_GrowShrink_0.01" = delDBN_0.01_s_100,
                  "Local_Markov_GrowShrink_0.05" = delDBN_0.05_s_100,
                  "Local_Markov_GrowShrink_0.1"= delDBN_0.1_s_100,
                  
                  "DelDBN_GrowShrink_0.001" = delDBN_0.001_d_100,
                  "DelDBN_GrowShrink_0.005" = delDBN_0.005_d_100,
                  "DelDBN_GrowShrink_0.01" = delDBN_0.01_d_100,
                  "DelDBN_GrowShrink_0.05" = delDBN_0.05_d_100,
                  "DelDBN_GrowShrink_0.1" = delDBN_0.1_d_100
                  
                   )

lapply(OneBigList, function(x) save(x, file = "data.RData"))

results <- mclapply(1:length(OneBigList) ,function(x) unlist(CalculatePrecisionAndRecallForMultiple(OneBigList[[x]], datasets)))
df <- do.call(rbind, results)
# df <- data.frame(matrix(unlist(results), nrow=length(OneBigList), byrow=T))
# colnames(df) <- c("precision_mean", "precision_sd" ,
#                       "recall_mean" , "recall_sd",
#                       "hamming_mean", "hamming_sd",
#                       "TN_mean", "TN_sd" ,
#                       "FN_mean", "FN_sd",
#                       "FP_mean", "FP_sd" ,
#                       "TP_mean", "TP_sd" )

#row.names(df) <- names(OneBigList)
df <- as.data.frame(df, row.names = names(OneBigList))
library(ggplot2)
ggplot(data = df) +
  geom_point(aes(x = recall_mean, y  = precision_mean))


mapply(function(x,y){
       x = amat(x)
       y =amat(y)
       sum(x==y & x ==1)
       }, Max_bics_100, Max_BN_0_100)