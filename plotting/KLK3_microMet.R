library(BoutrosLab.plotting.general);
library(readxl);
library(pROC);
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/tcell/recurrence/')
#dat <- data.frame(read_excel('/Users/sujunchen/Dropbox/me/singleCell/data_revision/Tcell_KLK3/FFPE_RNA_combined.xlsx'));
#dat <- data.frame(read_excel('FFPE_RNA_combined.xlsx'));
dat <- read.table('FFPE_combined_28samples.txt', header = TRUE);
to.plot <- dat[dat$gene=='PSA_3', ];
pval <- scientific.notation(wilcox.test(to.plot$mean~to.plot$Status)$p.value);
t_col <- function(color, percent = 20, name = NULL) {
  rgb.val <- col2rgb(color)
## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  return(t.col)
  };

col_b <- c('blue', 'red');
col_f <- unlist(lapply(col_b, function(x) t_col(x, 40)));

create.boxplot(
	data = to.plot,
	formula = mean~Status,
    xlab.label = 'Recurrence after surgery',
    ylab.label = 'Relative KLK3 expression',
    filename = NULL,
    #filename = paste0(Sys.Date(), '_tcell_KLK3_recurrence.pdf'),
    xaxis.lab = c('No', 'Yes'),
    add.stripplot = TRUE,
    xlab.cex = 1.5,
    ylab.cex = 1.5,
    style = 'Nature',
    #ylimits = c(0, 0.75),
    #yat = seq(0, 0.7, 0.2),
    lwd = 2,
    border.col = col_b,
    #col = col_f,
    height = 5,
    width = 4.5,
    add.text = TRUE,
    text.labels = pval,
    text.x = 1.5,
    text.y = max(to.plot$mean) + 0.02
    );
###
roc.mytest <- roc(factor(to.plot$Status), to.plot$mean);
plot.data <- data.frame(sensitivity = roc.mytest$sensitivities, specificity = roc.mytest$specificities);

create.scatterplot(
  formula = sensitivity~1-specificity,
  data = plot.data,
  type = c('l', 'p'),
  lwd = 2,
  xlimits = c(0, 1),
  ylimits = c(0, 1),
  add.xyline = TRUE, xyline.col = 'black', xyline.lty = 2,
  #filename = generate.filename('AUC_klk3', 'microMet', 'pdf'),
  style = 'Nature',
  cex = .1, 
  key = list(
    text = list(
      lab = paste0('AUC: ', round(roc.mytest$auc, 2)), cex = 1, col = 'black'
      ),
    x = 0.6, y = 0.1
    )
  );
