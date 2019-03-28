# Put G1DBN, Shrinkage on G2DBN on same graph
## For p = 50 and n = 20,50
seeds = 1:100
p = 50
n = 20
noiseLevel = 10

datasets = lapply(seeds, function(x) SimulateData(p,n,noiseLevel, x ))

alpha1 = 0.7
G1s = lapply(datasets, function(x) ExecuteG1DBNS1(x$data, alpha1))
curves_G1<- mapply(function(G1, ds) PRCurve_GeneNet(score=G1$mat$S1ls,validMat=abs(ds$RealNet$AdjMatrix )>0,dec=FALSE), G1s, datasets, SIMPLIFY = FALSE)

G1_means = CalculateMeans(curves_G1)

save(G1s, file = "G1s.rds")


## DO G2 ## 

G2s = mapply(function(x,y,z) {
  print(z)
  DBNScoreStep2(S1= x$mat$S1ls,data=y$data,alpha1=alpha1)
}, G1s, datasets, 1:length(G1s), SIMPLIFY = FALSE)
badbad = which(unlist(lapply(G2s, function(x) sum(x, na.rm = TRUE))) == 0)
## last point here hamda  ###
curves_G2 = mapply(function(G2, ds) PRCurve_GeneNet(score=G2,validMat=abs(ds$RealNet$AdjMatrix)>0,dec=FALSE)
                   , G2s[-badbad], datasets[-badbad], SIMPLIFY = FALSE)

G2_means = CalculateMeans(curves_G2)



## Do Shrinkage 

shrinks = mapply(function(ds, ind) tryCatch(GetScores(ds$data,ind  ), error = function(e) NULL), datasets, 1:length(datasets), SIMPLIFY = FALSE)
badbad = which(unlist(lapply(shrinks, is.null)))

curves_shrink<-  mapply(function(s, ds) PRCurve_GeneNet(score=s,validMat=t(ds$RealNet$AdjMatrix ), dec = FALSE), shrinks, datasets, SIMPLIFY = FALSE)

shrink_means = CalculateMeans(curves_shrink)




lassos = lapply(datasets, function(x) tryCatch(ApplyLars(x$data), error = function(e) NULL))
badbad = which(unlist(lapply(lassos, is.null)))


lasso1 <- lapply(lassos, function(x) 1 - x)

curves_lasso <- mapply(function(lasso, ds) PRCurve_GeneNet(score=lasso,validMat=abs(ds$RealNet$AdjMatrix )>0,dec=TRUE), lassos, datasets, SIMPLIFY = FALSE)
lasso_means <- CalculateMeans(curves_lasso)


# Put G1DBN, Shrinkage on G2DBN on same graph
seeds = 1:100
p = 50
n = 50
noiseLevel = 10

datasets_50 = lapply(seeds, function(x) SimulateData(p,n,noiseLevel, x ))

alpha1 = 0.7
G1s_50 = lapply(datasets_50, function(x) ExecuteG1DBNS1(x$data, alpha1))
curves_G1_50<- mapply(function(G1, ds) PRCurve_GeneNet(score=G1$mat$S1ls,validMat=abs(ds$RealNet$AdjMatrix )>0,dec=FALSE), G1s_50, datasets_50, SIMPLIFY = FALSE)

G1_50_means = CalculateMeans(curves_G1_50)




## DO G2 ## 

G2s_50 = mapply(function(x,y,z) {
  print(z)
  DBNScoreStep2(S1= x$mat$S1ls,data=y$data,alpha1=alpha1)
}, G1s_50, datasets_50, 1:length(G1s_50), SIMPLIFY = FALSE)
badbad = 46 # some pb regression 
curves_G2_50 = mapply(function(G2, ds) PRCurve_GeneNet(score=G2,validMat=abs(ds$RealNet$AdjMatrix)>0,dec=FALSE)
                      , G2s_50, datasets_50, SIMPLIFY = FALSE)
G2_50_means = CalculateMeans(curves_G2_50)


## Do Shrinkage 

shrinks_50 = mapply(function(ds, ind) tryCatch(GetScores(ds$data,ind  ), error = function(e) NULL), datasets_50, 1:length(datasets_50), SIMPLIFY = FALSE)
badbad = which(unlist(lapply(shrinks_50, is.null)))
#shrinks_50 = shrinks_50[-badbad]
#datasets_ = datasets[-badbad]
curves_shrink_50<-  mapply(function(s, ds) PRCurve_GeneNet(score=s,validMat=t(ds$RealNet$AdjMatrix>0 ), dec = FALSE), shrinks_50[-badbad], datasets_50[-badbad], SIMPLIFY = FALSE)
shrink_50_means = CalculateMeans(curves_shrink_50)



lassos_50 = lapply(datasets_50, function(x) tryCatch(ApplyLars(x$data), error = function(e) NULL))
badbad = which(unlist(lapply(lassos_50, is.null)))

curves_lasso_50<- mapply(function(lasso, ds) PRCurve_GeneNet(score=lasso,validMat=abs(ds$RealNet$AdjMatrix )>0,dec=TRUE), lassos_50[-badbad], datasets_50[-badbad], SIMPLIFY = FALSE)
lasso_means_50 = CalculateMeans(curves_lasso_50)


plot(x = G1_50_means$recall, y = G1_50_means$precisions, type = "l", col = "red",lty=2,lwd=4, ylim=c( 0,1),xlim = c(0,1), xlab = "recall", ylab = "precision", main = "n = 50, p = 50") 
lines(x = G2_50_means$recall, y = G2_50_means$precisions, type = "l", col = "red",lty=1,lwd=3, ylim=c( 0,1),xlim = c(0,1))
lines(x = shrink_50_means$recall, y = shrink_50_means$precisions, type = "l", col = "green",lty=2,lwd=4, ylim=c( 0,1),xlim = c(0,1)) 
lines(x = lasso_means_50$recall, y = lasso_means_50$precisions, type = "l", col = "blue",lty=2,lwd=4, ylim=c( 0,1),xlim = c(0,1)) 
temp <- legend("topright", legend = c(" ", " ", " ", " "),
               text.width = strwidth("1,000,000,000"),
               lty = c(2,1,2,2), lwd = c(2,1,2,2), xjust = 1, yjust = 1,
               title = "", col = c("red", "red", "green", "blue"))
text(temp$rect$left + temp$rect$w, temp$text$y,
     c("Step 1", "Step 2", "shrinkage", "lasso"), pos = 2)



plot(x = G1_means$recall, y = G1_means$precisions, type = "l", col = "red",lty=2,lwd=4, ylim=c( 0,1),xlim = c(0,1), xlab = "recall", ylab = "precision", main = "n = 20, p = 50") 
lines(x = G2_means$recall, y = G2_means$precisions, type = "l", col = "red",lty=1,lwd=3, ylim=c( 0,1),xlim = c(0,1))
lines(x = shrink_means$recall, y = shrink_means$precisions, type = "l", col = "green",lty=2,lwd=4, ylim=c( 0,1),xlim = c(0,1)) 
lines(x = lasso_means$recall, y = lasso_means$precisions, type = "l", col = "blue",lty=2,lwd=4, ylim=c( 0,1),xlim = c(0,1)) 
temp <- legend("topright", legend = c(" ", " ", " ", " "),
               text.width = strwidth("1,000,000,000"),
               lty = c(2,1,2,2), lwd = c(2,1,2,2), xjust = 1, yjust = 1,
               title = "", col = c("red", "red", "green", "blue"))
text(temp$rect$left + temp$rect$w, temp$text$y,
     c("Step 1", "Step 2", "shrinkage", "lasso"), pos = 2)
