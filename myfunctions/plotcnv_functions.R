calc_cnv_sep <- function(mycnv.all, annot, otfile_name, cmethod = 'pearson'){
	cnv.all <- colMeans((mycnv.all-1)^2);
	saveRDS(cnv.all, file = generate.filename('infercnv_meansquare', otfile_name, 'rds'));
	cnv.epi <- data.frame(cnv = cnv.all, type = annot[match(names(cnv.all), rownames(annot)), ]$type);
	cnv.epi <- cnv.epi[order(-cnv.epi$cnv), ];
	cnv.top <- cnv.epi[1:(nrow(cnv.epi)/20), ];
	avg.top <- rowMeans(mycnv.all[, rownames(cnv.top)]);
	mycor.all <- apply(mycnv.all, 2, function(x) cor(x, avg.top, method = cmethod));
	saveRDS(mycor.all, file = generate.filename('infercnv_correlation', otfile_name, 'rds'));
};

plot_distr_sep <- function(to.plot, otfile_name, text.labels, 	epikey = FALSE, l2 = FALSE){
	mypch <- c(3,17)
	names(mypch) <- c("epithelia", "non");
	to.plot$pch <- mypch[to.plot$type];
	to.plot$col <- 'black';
	if(nrow(to.plot[to.plot$score>0.04&to.plot$cor>0.4, ])>0){
	to.plot[to.plot$score>0.04&to.plot$cor>0.4, ]$col <- 'red';		
	}
	if(nrow(to.plot[to.plot$score<0.04&to.plot$cor<0.4, ])>0){
	to.plot[to.plot$score<0.04&to.plot$cor<0.4, ]$col <- 'blue';
	};
	q1 <- round(100*nrow(to.plot[to.plot$col=='red'&to.plot$type=='epithelia', ])/nrow(to.plot[to.plot$col=='red', ]));
	q2 <- round(100*nrow(to.plot[to.plot$col=='black'&to.plot$type=='epithelia', ])/nrow(to.plot[to.plot$col=='black', ]));
	q3 <- round(100*nrow(to.plot[to.plot$col=='blue'&to.plot$type=='epithelia', ])/nrow(to.plot[to.plot$col=='blue', ]));
	if(epikey){
		mykey <- list(
	            text = list(
	            lab = paste0(c(q1, q2, q3), '%'),
	            cex = 1, 
	            col = c('red', 'black', 'blue')
	            ),
	    	points = list(
	            pch = 20, 
	            col = c('red', 'black', 'blue'), 
	            fill = c('red', 'black', 'blue'), 
	            cex = 1
	            ),
	    x = 0.75, 
	    y = 0.95
	    )} else{
			mykey <- NULL
	    };
	if(l2){
		vline = c(0.04, 0.02);
		hline = c(0.4, 0.2);
		tline = c(1, 2);
	}else{
		vline = 0.04;
		hline = 0.4;
		tline = 1;
	};

	xlimit <- max(0.1, max(to.plot$score));
	create.scatterplot(
		formula = cor~score,
		data = to.plot,
		col = to.plot$col,
		cex = .5,
		#pch = to.plot$pch,
		xlab.label = '',
		ylab.label = '',
		xaxis.fontface = 'plain', 
		yaxis.fontface = 'plain',
		abline.h = hline,
		abline.v = vline,
		abline.lty = tline,
		xlimits = c(0, xlimit),
		ylimits = c(-0.5, 1),
		style = 'Nature',
		yat = seq(-0.5, 1, 0.5),
		xat = seq(0, round(max(to.plot$score)*10)/10, round(max(to.plot$score)*10)/20),
		xlimit = c(0, max(to.plot$score)),
		filename = generate.filename('cnvplot', otfile_name, 'pdf'),
		width = 4,
		height = 4,
		add.text = TRUE,
		text.x = xlimit*3/4,
		text.y = -0.1,
		text.labels = text.labels,
		text.cex = 1.5,
		covariate.legend = cov.legend,
		print.colour.key = color.key,
		key = mykey
		);
};

plot_hm <- function(to.plot, mycnv, mygene, name, otfile_name){
	mypch <- c(3,17)
	names(mypch) <- c("epithelia", "non");
	to.plot$pch <- mypch[to.plot$type];
	#to.plot$cor <- to.plot$cor.5;
	#to.plot[to.plot$type=='epithelia', ]$cor <- apply(to.plot[to.plot$type=='epithelia', 2:4], 1, max);
	to.plot$col <- 'black';
	if(nrow(to.plot[to.plot$score>0.04&to.plot$cor>0.4, ])>0){
	to.plot[to.plot$score>0.04&to.plot$cor>0.4, ]$col <- 'red';		
	}
	if(nrow(to.plot[to.plot$score<0.04&to.plot$cor<0.4, ])>0){
	to.plot[to.plot$score<0.04&to.plot$cor<0.4, ]$col <- 'blue';
	};
	dat.plot <- to.plot;
	####
	mygene <- mygene[mygene$V1%in%rownames(mycnv), ];
	#mygene <- mygene[mixedorder(mygene$V2), ];
	mychr <- table(mygene$V2);
	mychr <- mychr[mixedorder(names(mychr))];
	gaps_col <- sapply(seq(length(mychr)), function(x) sum(mychr[1:x]));
	#mybreaks <- unique(c(seq(0, 0.5, 0.1), seq(0.5, 1.5, 0.01), seq(1.5, 2, 0.1)));
	mybreaks <- seq(0, 2, length.out = 16);
	#mybreaks <- seq(0, 2, length.out = 16);
	#mycol <- colorRampPalette(c('blue', 'white', 'red'))(length(mybreaks));
	custom_pal <- color.palette(c("darkblue", "white", "darkred"),c(2, 2));
	mycol <- custom_pal(length(mybreaks) + 1);
	to.plot <- mycnv[match(mygene$V1, rownames(mycnv)), ];
	mygrid <- rep(NA, nrow(to.plot));
	mygrid[gaps_col] <- 'black';
	### create top covariate for chromosome annotation
	col.chr <- colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(length(unique(names(mychr))));
	names(col.chr) <- unique(names(mychr));
	col.type <- c('#66C2A5', '#FFED6F', 'cornflowerblue');
	names(col.type) <- c('red', 'blue', 'black');
	if(name %in%c('JD1800153SL','JD1800172SL')){
		chr.covariate <- list(
			rect = list(
				col = rep(col.chr, times = mychr),
				fill = rep(col.chr, times = mychr)
				)
			)
		color.key <- FALSE
	}else if(name == 'JD1800177SL'){
		chr.covariate <- NULL
		color.key <- FALSE

	}else {
		chr.covariate <- NULL
		color.key <- FALSE
		cov.legend <- NULL
	};
	chr.covariate <- list(
	rect = list(
		col = rep(col.chr, times = mychr),
		fill = rep(col.chr, times = mychr)
		)
	);
	mytype <- dat.plot[colnames(to.plot), ]$col;
	to.plot <- cbind(to.plot[, mytype=='red'], to.plot[, mytype=='black'], to.plot[, mytype=='blue']);
	mytype <- dat.plot[colnames(to.plot), ]$col;
	gaps_row <- c(table(mytype)['red'], sum(table(mytype)[c('black', 'red')]));
	mygrid.row <- rep(NA, ncol(to.plot));
	mygrid.row[gaps_row] <- 'black';

	cell.covariate <- list(
		rect = list(
			col = col.type[mytype],
			fill = col.type[mytype]
			)
		);
	print('start ploting')
	p <- create.heatmap(
		x = to.plot,
		filename = NULL,
		colour.scheme = mycol,
		total.colours = length(mybreaks) + 1,
		at = mybreaks,
		cluster.dimensions = 'none',
		row.colour = 'black',
		col.lines = gaps_col,
		col.lwd = 0.5,
		row.lwd = 0.5,
		row.lines = gaps_row,
		covariates.top = chr.covariate,
		covariates = cell.covariate,
		covariate.legend = NULL,
		grid.col = TRUE,
		force.grid.col = TRUE,
		grid.row = TRUE,
		force.grid.row = TRUE,
		print.colour.key = color.key
		);
	png(generate.filename(otfile_name, name, 'png'), width = 10, height = 3.5, units = 'in', res = 300);
	print(p);
	dev.off();
};
###
plot_hm_type <- function(to.plot, mycnv, mygene, name, otfile_name){
	mypch <- c(3,17)
	names(mypch) <- c("epithelia", "non");
	to.plot$pch <- mypch[to.plot$type];
	#to.plot$cor <- to.plot$cor.5;
	#to.plot[to.plot$type=='epithelia', ]$cor <- apply(to.plot[to.plot$type=='epithelia', 2:4], 1, max);
	to.plot$col <- to.plot$type;
	dat.plot <- to.plot;
	####
	mygene <- mygene[mygene$V1%in%rownames(mycnv), ];
	#mygene <- mygene[mixedorder(mygene$V2), ];
	mychr <- table(mygene$V2);
	mychr <- mychr[mixedorder(names(mychr))];
	gaps_col <- sapply(seq(length(mychr)), function(x) sum(mychr[1:x]));
	#mybreaks <- unique(c(seq(0, 0.5, 0.1), seq(0.5, 1.5, 0.01), seq(1.5, 2, 0.1)));
	mybreaks <- seq(0, 2, length.out = 16);
	#mybreaks <- seq(0, 2, length.out = 16);
	#mycol <- colorRampPalette(c('blue', 'white', 'red'))(length(mybreaks));
	custom_pal <- color.palette(c("darkblue", "white", "darkred"),c(2, 2));
	mycol <- custom_pal(length(mybreaks) + 1);
	to.plot <- mycnv[match(mygene$V1, rownames(mycnv)), ];
	mygrid <- rep(NA, nrow(to.plot));
	mygrid[gaps_col] <- 'black';
	### create top covariate for chromosome annotation
	col.chr <- colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(length(unique(names(mychr))));
	names(col.chr) <- unique(names(mychr));
	col.type <- default.colours(2);
	names(col.type) <- c('epithelia', 'non');
	if(name %in%c('JD1800153SL','JD1800172SL')){
		chr.covariate <- list(
			rect = list(
				col = rep(col.chr, times = mychr),
				fill = rep(col.chr, times = mychr)
				)
			)
		color.key <- FALSE
	}else if(name == 'JD1800177SL'){
		chr.covariate <- NULL
		color.key <- FALSE

	}else {
		chr.covariate <- NULL
		color.key <- FALSE
		cov.legend <- NULL
	};
	chr.covariate <- list(
	rect = list(
		col = rep(col.chr, times = mychr),
		fill = rep(col.chr, times = mychr)
		)
	);
	mytype <- dat.plot[colnames(to.plot), ]$col;
	to.plot <- cbind(to.plot[, mytype=='epithelia'], to.plot[, mytype=='non']);
	mytype <- dat.plot[colnames(to.plot), ]$col;
	gaps_row <- c(table(mytype)['epithelia']);
	mygrid.row <- rep(NA, ncol(to.plot));
	mygrid.row[gaps_row] <- 'black';

	cell.covariate <- list(
		rect = list(
			col = col.type[mytype],
			fill = col.type[mytype]
			)
		);
	print('start ploting')
	p <- create.heatmap(
		x = to.plot,
		filename = NULL,
		colour.scheme = mycol,
		total.colours = length(mybreaks) + 1,
		at = mybreaks,
		cluster.dimensions = 'none',
		row.colour = 'black',
		col.lines = gaps_col,
		col.lwd = 0.5,
		row.lwd = 0.5,
		row.lines = gaps_row,
		covariates.top = chr.covariate,
		covariates = cell.covariate,
		covariate.legend = NULL,
		grid.col = TRUE,
		force.grid.col = TRUE,
		grid.row = TRUE,
		force.grid.row = TRUE,
		print.colour.key = color.key
		);
	png(generate.filename(otfile_name, name, 'png'), width = 10, height = 3.5, units = 'in', res = 300);
	print(p);
	dev.off();
};

plot_distr_class <- function(to.plot, otfile_name, text.labels, epikey = FALSE, l2 = FALSE){
	mypch <- c(3,17)
	names(mypch) <- c("epithelia", "non");
	to.plot$pch <- mypch[to.plot$type];

	mycol <- c('red', 'black', 'blue');
	names(mycol) <- c('2', '3', '1');
	to.plot$col <- mycol[as.character(as.vector(to.plot$class))];

	q1 <- round(100*nrow(to.plot[to.plot$col=='red'&to.plot$type=='epithelia', ])/nrow(to.plot[to.plot$col=='red', ]));
	q2 <- round(100*nrow(to.plot[to.plot$col=='black'&to.plot$type=='epithelia', ])/nrow(to.plot[to.plot$col=='black', ]));
	q3 <- round(100*nrow(to.plot[to.plot$col=='blue'&to.plot$type=='epithelia', ])/nrow(to.plot[to.plot$col=='blue', ]));
	if(epikey){
		mykey <- list(
	            text = list(
	            lab = paste0(c(q1, q2, q3), '%'),
	            cex = 1, 
	            col = c('red', 'black', 'blue')
	            ),
	    	points = list(
	            pch = 20, 
	            col = c('red', 'black', 'blue'), 
	            fill = c('red', 'black', 'blue'), 
	            cex = 1
	            ),
	    x = 0.75, 
	    y = 0.95
	    )} else{
			mykey <- NULL
	    };
	if(l2){
		vline = c(0.04, 0.02);
		hline = c(0.4, 0.2);
		tline = c(1, 2);
	}else{
		vline = 0.04;
		hline = 0.4;
		tline = 1;
	};

	xlimit <- max(0.1, max(to.plot$score));
	create.scatterplot(
		formula = corr~score,
		data = to.plot,
		col = to.plot$col,
		cex = .5,
		#pch = to.plot$pch,
		xlab.label = '',
		ylab.label = '',
		xaxis.fontface = 'plain', 
		yaxis.fontface = 'plain',
		abline.h = hline,
		abline.v = vline,
		abline.lty = tline,
		xlimits = c(0, xlimit),
		ylimits = c(-0.5, 1),
		style = 'Nature',
		yat = seq(-0.5, 1, 0.5),
		xat = seq(0, round(max(to.plot$score)*10)/10, round(max(to.plot$score)*10)/20),
		xlimit = c(0, max(to.plot$score)),
		filename = generate.filename('cnvclust', otfile_name, 'pdf'),
		width = 4,
		height = 4,
		add.text = TRUE,
		text.x = xlimit*3/4,
		text.y = -0.1,
		text.labels = text.labels,
		text.cex = 1.5,
		covariate.legend = cov.legend,
		print.colour.key = color.key,
		key = mykey
		);
};

