library(BoutrosLab.plotting.general);
to.plot <- read.table('/Users/sujunchen/Dropbox/me/singleCell/figures/F3/elements/bob/Figure_Dif_Again.txt', header = TRUE, row.names =1);
rownames(to.plot) <- gsub('tumor', 'T', rownames(to.plot));
to.plot <- to.plot[grep('_T', rownames(to.plot)), ]
to.plot <- to.plot[!grepl('_N_T', rownames(to.plot)), ]
colours <- colorRampPalette(c('grey', 'red'))(8);
to.plot.f <- data.frame(t(to.plot));
to.plot.f <- to.plot.f[rev(rownames(to.plot.f)), ];
to.plot <- data.frame(type = rep(seq(nrow(to.plot.f)), times = 5), samp = rep(seq(5), each = nrow(to.plot.f)));

p1 <- create.scatterplot(
	formula = type ~ samp,
	data = to.plot,
	cex = unlist(to.plot.f)*5 + 0.1,
	col = colours[round(10*unlist(to.plot.f))+1],
	xlab.label = '',
	ylab.label = '',
	xaxis.lab = gsub('_T', '', colnames(to.plot.f)),
	yaxis.lab = rownames(to.plot.f),
	yat = seq(19),
	xat = seq(5),
	xaxis.rot = 45,
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	#xlimits = c(0.5, 3.5),
	width = 3,
	style = 'BoutrosLab'
	);
pdf(generate.filename('tcell', 'TumorMarker', 'pdf'), width = 5);
print(p1);
dev.off();

p2 <- create.scatterplot(
	formula = type ~ samp,
	data = to.plot,
	cex = rev(c(0, 1, 2, 3, 4) + 0.1),
	col = colours[rev(round(2*c(0, 1, 2, 3, 4))+1)]
	);
pdf(generate.filename('tcell', 'TumorMarker_legend', 'pdf'), width = 5	);
print(p2);
dev.off();
