#Scripts 1-3 are adapted from LTLA at https://github.com/LTLA/MatrixEval2017/blob/master/real/10X/
#The scripts should be executed in the following order:
#preprocess.R filter cells, check some qc metrics
#normalize.R normalize data, fast cluster for each sample(batch) first then calculate normalization factor
#dimred.R run MNN to removal batch effect, for the current data with >30k cells, need request 125G memory to run
#findclust.R run SNN to find clusters, do some tuning for parameter k. We find k=20 works well for the current data
#miscellaneous/tunetsne.R using different max_iter and perplexity, check the output visually, we find the default works just well for the current analysis
#analyzeclust_tsne.R run TSNE on the data, generates tsne plot to inspect batch removal effect, clusters, etc. 
#annotype.R run QuSage to annotate cell type for each cluster
#analyzeclust_plot.R 
#reeval_vcaf.R 
#
####
