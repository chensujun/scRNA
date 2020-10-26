library(mclust);
data(Hidalgo1872, package = "MMST");
Thickness <- Hidalgo1872$thickness
Year <- rep(c("1872", "1873-74"), c(289, 196)) 
dens <- densityMclust(Thickness);
br <- seq(min(Thickness), max(Thickness), length = 21)
plot(dens, what = "density", data = Thickness, breaks = br);
h1 <- hist(Thickness[Year == "1872"], breaks = br, plot = FALSE)
h1$density <- h1$density*prop.table(table(Year))[1]
h2 <- hist(Thickness[Year == "1873-74"], breaks = br, plot = FALSE) 
h2$density <- h2$density*prop.table(table(Year))[2]
x <- seq(min(Thickness)-diff(range(Thickness))/10,
	max(Thickness)+diff(range(Thickness))/10, length = 200) 
cdens <- predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
col <- adjustcolor(mclust.options("classPlotColors")[1:2], alpha = 0.3)
plot(h1, xlab = "Thickness", freq = FALSE, main = "", border = FALSE, col = col[1],
                   xlim = range(x), ylim = range(h1$density, h2$density, cdens))
plot(h2, add = TRUE, freq = FALSE, border = FALSE, col = col[2])
matplot(x, cdens, type = "l", lwd = 1, add = TRUE, lty = 1:3, col = 1)
###
plot_dens_default <- function(ddat, status, name = ''){
	dens <- densityMclust(ddat,  modelNames = 'V');
	br <- seq(min(ddat), max(ddat), length = 50);
	plot(dens, what = 'density', data = ddat, breaks = br);
	h1 <- hist(ddat[status=='N'], breaks = br, plot = FALSE);
	h1$density <- h1$density*prop.table(table(status))['N'];
	h2 <- hist(ddat[status=='M'], breaks = br, plot = FALSE);
	h2$density <- h2$density*prop.table(table(status))['M'];
	x <- seq(min(ddat)-diff(range(ddat))/10, max(ddat)+diff(range(ddat))/10, length = 200);
	cdens <- predict(dens, x, what = 'cdens');
	cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro));
	col <- adjustcolor(mclust.options("classPlotColors")[1:2], alpha = 0.3)
	plot(h1, xlab = 'Score', freq = FALSE, main = name, border = FALSE, col = col[1], 
		xlim = range(x), ylim = c(0, max(h1$density, h2$density, cdens)));
	plot(h2, add = TRUE, freq = FALSE, border = FALSE, col = col[2])
	matplot(x, cdens, type = "l", lwd = 1, add = TRUE, lty = 1:3, col = 1);
	box(bty = 'l')
};

ddat <- mydata$score;
status <- ifelse(ddat<0.03, 'N', 'M');
plot_dens_default(ddat, status);

ddat <- mydata$corr;
status <- ifelse(ddat<0.3, 'N', 'M');
plot_dens_default(ddat, status);
annot <- readRDS('2020-05-15_metadata_13tumor.rds');
rownames(annot) <- gsub('-', '.', rownames(annot))
to.plot$samp <- annot[rownames(to.plot), ]$orig.ident;
to.plot$type <- annot[rownames(to.plot), ]$type;
for(samp in unique(to.plot$samp)){
	ddat <- to.plot[to.plot$samp==samp, ]$score;
	status <- ifelse(ddat<0.03, 'N', 'M');
	plot_dens_default(ddat, status, samp);
};

for(type in unique(to.plot$type)){
	ddat <- to.plot[to.plot$type==type, ]$score;
	status <- ifelse(ddat<0.03, 'N', 'M');
	plot_dens_default(ddat, status, type);
};
#####

#####
cnv.mel <- readRDS('2020-05-26_cnv_score_mel_0521.rds');
ddat <- cnv.mel$score;

ddat <- cnv.mel$cor;
status <- ifelse(ddat<0.3, 'N', 'M');
plot_dens_default(ddat, status);

