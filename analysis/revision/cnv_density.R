library("dbscan");
setwd('/Users/sujunchen/Document/labs/2020/scRNA');
to.plot <- readRDS('./2020-05-26_cnv_score_pri_0518.rds');
x <- to.plot[, 1:2]
kNNdistplot(x, k = 3)
abline(h=.01, col = "red", lty=2);
res <- dbscan(x, eps = 0.01, minPts = 3);
###
res <- optics(x, eps = 10, minPts = 10)