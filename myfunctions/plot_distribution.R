
plot_distribution <- function(percent.tab, name, annot, ncell.max = 610, mytype = 'tiff', width = 13, height = 3, mycol = mycol, legend = NULL, right.padding = 2){
	#percent.tab <- as.data.frame.matrix(table(test.data@meta.data$orig.ident, test.data@meta.data$res.1));
	percent.tab$idc <- annot[match(rownames(percent.tab), annot$sample), ]$IDCP;
	percent.tab <- percent.tab[order(percent.tab$idc), ];
	percent.tab$gs <- annot[match(rownames(percent.tab), annot$sample), ]$pGS;
	percent.tab <- percent.tab[order(percent.tab$gs), ];
	percent.tab <- percent.tab[, colnames(percent.tab)!='gs'];
	annot <- annot[match(rownames(percent.tab), annot$sample), ];
	cov.annot <- create.covariate(annot[seq(nrow(annot), 1), ]);
	to.plot <- apply(percent.tab[, grep('idc', colnames(percent.tab), invert = TRUE)], 1, function(x) 100*x/sum(x));
	colnames(to.plot) <- gsub('JD1800|SL', '', colnames(to.plot));
	to.plot <- melt(to.plot);
	#to.plot$Var1 <- gsub('Normal epithelial', 'Epithelial', to.plot$Var1);
	to.plot$Var2 <- factor(to.plot$Var2, levels = gsub('JD1800|SL', '', rownames(percent.tab)));
	to.plot$Var1 <- factor(to.plot$Var1);
	#mycol <- rainbow(2*ncol(percent.tab[, grep('idc', colnames(percent.tab), invert = TRUE)]))[seq(1, 2*ncol(percent.tab[, grep('idc', colnames(percent.tab), invert = TRUE)]), 2)]
	#names(mycol) <- levels(to.plot$Var1);

	mycol.t <- c("#B3EEB4", "#2F6D60", "#5EC2AA");
   	names(mycol.t) <- c("T2c", "T3a", "T3b");
    mycol.gs <- c('yellow', 'orange', '#B33187', '#6E1573')
    names(mycol.gs) <- c('7(3+4)', '7(4+3)', 9, 10);

	covariates.legend <- legend.grob(
		    list(
		        legend = list(
		                colours = mycol.gs,
		                    labels = names(mycol.gs),
		                    title = expression(underline("Gleason score")),
		                    border = 'white'
	                    ),
	            legend = list(
	                    colours = c("#FEE6CE","#FDAE6B","#E6550D"),
	                    labels = c("0 - 9.9", "10 - 19.9", expression(''>=20)),
	                    title = expression(underline('PSA (ng/mL)')),
	                    border = 'white'
	                    ),
	            legend = list(
	                    colours = mycol.t,
	                    labels = names(mycol.t),
	                    title = expression(underline('T category')),
	                    border = 'white'
	                    ),
				legend = list(
	                    colours = c("firebrick3","lightskyblue"),
	                    labels = c('TRUE', 'FALSE'),
	                    title = expression(underline('IDC')),
	                    border = 'white'
	                    ),
				legend = list(
	                    colours = c("darkblue","lightskyblue"),
	                    labels = c('TRUE', 'FALSE'),
	                    title = expression(underline('Node Mets')),
	                    border = 'white'
	                    ),
				legend = list(
	                    colours = c("darkgreen","lightskyblue"),
	                    labels = c('TRUE', 'FALSE'),
	                    title = expression(underline('Bone Mets')),
	                    border = 'white'
	                    ),
	            legend = list(
					colours = mycol,
					labels = names(mycol),
					title = expression(underline('Cell type')),
					border = 'white',
					cex = 0.5
					),
	            legend = list(
					colours = default.colours(2),
					labels = c('TRUE', 'FALSE'),
					title = expression(underline('Epithelial')),
					border = 'white',
					cex = 0.5
					),
	            legend = list(
					colours = default.colours(4)[3:4],
					labels = c('immune', 'stroma'),
					title = expression(underline('TME')),
					border = 'white',
					cex = 0.5
					)

		        ),
		    size = 2.6,
		    label.cex = 1.5,
		    title.cex = 1.5,
		    title.just = 'left',
		    between.row = 0.15,
		    layout = c(2, 5)
		    );

	bar.plot <- create.barplot(
		#file = filename,
		formula = Var2~value,
		data = to.plot,
		groups = Var1,
		col = mycol,
		stack = TRUE,
		plot.horizontal = TRUE,
		ylab.label = 'Patient',
		xlab.label = 'Percentage',
		xlimits = c(0, 100),
		xat = seq(0, 100, 20),
		yaxis.rot = 90,
		style = 'Nature',
		width = 8,
		height = 6,
		y.spacing = -1,
		border.col = 'white',
		box.ratio = 8
		);
	bar.ncell.data <- data.frame(ncell = rowSums(percent.tab[, grep('idc', colnames(percent.tab), invert = TRUE)]), patient = gsub('JD1800|SL', '', names(rowSums(percent.tab))));
	bar.ncell.data$patient <- factor(bar.ncell.data$patient, levels = bar.ncell.data$patient);
	bar.ncell.data$idc <- percent.tab$idc;
	bar.ncell <- create.barplot(
		formula = patient~ncell,
		data = bar.ncell.data,
		plot.horizontal = TRUE,
		ylab.label = 'Patient',
		xlab.label = '# cells',
		yaxis.lab = rep('', nrow(bar.ncell.data)),
		yaxis.rot = 90,
		xlimits = c(0, ncell.max),
		xat = c(0, floor(ncell.max/100)*100),
		style = 'Nature',
		width = 8,
		height = 6,
		y.spacing = -1,
		border.col = 'white',
		box.ratio = 8
		);	

	bar.idc <- create.barplot(
		formula = patient~idc,
		data = bar.ncell.data,
		plot.horizontal = TRUE,
		ylab.label = 'Patient',
		xlab.label = '% IDC',
		xlimits = c(0, 85),
		xat = c(0, 80),
		yaxis.lab = rep('', nrow(bar.ncell.data)),
		yaxis.rot = 90,
		style = 'Nature',
		width = 8,
		height = 6,
		y.spacing = -1,
		border.col = 'white',
		box.ratio = 8
		);	
	percent.tab.e <- data.frame(epi = rowSums(percent.tab[, c('Basal/intermediate', 'Luminal')]), 
		non = rowSums(percent.tab[, c('Endothelial', 'Fibroblast', 'Mast', 'Monocytic', 'T')]))
	to.plot <- apply(percent.tab.e[, grep('idc', colnames(percent.tab.e), invert = TRUE)], 1, function(x) 100*x/sum(x));
	colnames(to.plot) <- gsub('JD1800|SL', '', colnames(to.plot));
	to.plot <- melt(to.plot);
	#to.plot$Var1 <- gsub('Normal epithelial', 'Epithelial', to.plot$Var1);
	to.plot$Var2 <- factor(to.plot$Var2, levels = gsub('JD1800|SL', '', rownames(percent.tab)));
	to.plot$Var1 <- factor(to.plot$Var1);
	bar.epi <- create.barplot(
		#file = 'test.pdf',
		formula = Var2~value,
		data = to.plot,
		groups = Var1,
		col = default.colours(2),
		stack = TRUE,
		plot.horizontal = TRUE,
		yaxis.lab = rep('', nrow(bar.ncell.data)),
		ylab.label = 'Patient',
		xlab.label = 'Percentage',
		xlimits = c(0, 100),
		xat = seq(0, 100, 50),
		yaxis.rot = 90,
		style = 'Nature',
		width = 8,
		height = 6,
		y.spacing = -1,
		border.col = 'white',
		box.ratio = 8
		);

	percent.tab.t <- data.frame(immune = rowSums(percent.tab[, c('Mast', 'Monocytic', 'T')]), 
		stroma = rowSums(percent.tab[, c('Endothelial', 'Fibroblast')]))
	to.plot <- apply(percent.tab.t[, grep('idc', colnames(percent.tab.t), invert = TRUE)], 1, function(x) 100*x/sum(x));
	colnames(to.plot) <- gsub('JD1800|SL', '', colnames(to.plot));
	to.plot <- melt(to.plot);
	to.plot$Var2 <- factor(to.plot$Var2, levels = gsub('JD1800|SL', '', rownames(percent.tab)));
	to.plot$Var1 <- factor(to.plot$Var1);
	bar.tme <- create.barplot(
		#file = 'test.pdf',
		formula = Var2~value,
		data = to.plot,
		groups = Var1,
		col = default.colours(4)[3:4],
		stack = TRUE,
		plot.horizontal = TRUE,
		yaxis.lab = rep('', nrow(bar.ncell.data)),
		ylab.label = 'Patient',
		xlab.label = 'Percentage',
		xlimits = c(0, 100),
		xat = seq(0, 100, 50),
		yaxis.rot = 90,
		style = 'Nature',
		width = 8,
		height = 6,
		y.spacing = -1,
		border.col = 'white',
		box.ratio = 8
		);
	if(is.null(legend)){
		mylegend <- NULL
		pdf(generate.filename('percent_patient', paste0(name, '_legend'), 'pdf'));
		grid.draw(covariates.legend)
		dev.off()
	}else{
		mylegend <- list(
		        inside = list(
		            x = 1.05,
		            y = 1,
		            fun = covariates.legend
		            )
		        );

	};

	if(mytype=='tiff'){
		filename = generate.filename('percent_patient', name, 'tiff');
		create.multiplot(
	    plot.objects = list(
	    	cov.annot,
	    	bar.idc,
	        bar.plot,
	        bar.epi,
	        bar.ncell,
	        bar.tme
	        ),
	    file = filename,
	    plot.layout = c(6,1),
	    panel.heights = 1, 
	    panel.widths = c(0.2, 0.15, 0.5, 0.2, 0.2, 0.2),
	    y.spacing = c(0,0, 0, 0, 1),
	    x.spacing = c(1, 1, 1, 1, 2),    
	    x.relation = 'free',
	    y.relation = 'free',
	    yaxis.tck = 0, 
	    xaxis.tck = 0.75,
	    print.new.legend = TRUE,
	    right.padding = right.padding,
	    left.padding = 0.01,
	    yaxis.cex = 1,
	    legend = mylegend,
	    style = 'Nature',
	    height = height,
	    width = width
	    );
	}else{
		filename = generate.filename('percent_patient', name, 'pdf');
		create.multiplot(
	    plot.objects = list(
	    	cov.annot,
	    	bar.idc,
	        bar.plot,
	        bar.epi,
	        bar.ncell,
	        bar.tme
	        ),
	    file = filename,
	    plot.layout = c(6,1),
	    panel.heights = 1, 
	    panel.widths = c(0.2, 0.15, 0.5, 0.2, 0.2, 0.2),
	    y.spacing = c(0,0, 0, 0, 1),
	    x.spacing = c(1, 1, 1, 1, 2),    
	    x.relation = 'free',
	    y.relation = 'free',
	    xaxis.lab = list(
	    	NULL, NULL, NULL,
	    	seq(0, 100, 50),
	    	NULL,
	    	seq(0, 100, 50)
	    	),
	    yaxis.tck = 0, 
	    xaxis.tck = 0.75,
	    print.new.legend = TRUE,
	    right.padding = right.padding,
	    left.padding = 0.01,
	    yaxis.cex = 1,
	    legend = mylegend,
	    style = 'Nature',
	    height = height,
	    width = width
	    );
	}
};
####
create.covariate = function(annot) {
    annot.numeric = annot[, c('pGS', 'pre.treatment.psa', 'pathological_t', 'idc', 'pN', 'Clinical_M')];
    #mycol.gs <- force.colour.scheme(levels(as.factor(annot$pGS)), 'gleason.score');
    mycol.gs <- c('yellow', 'orange', '#B33187', '#6E1573')
    names(mycol.gs) <- c('7(3+4)', '7(4+3)', 9, 10);
    mycol <- force.colour.scheme(levels(as.factor(annot$pathological_t)), 'clinicalt9');
   	names(mycol) <- levels(as.factor(annot$pathological_t));
   	mycol <- mycol[levels(as.factor(annot$pathological_t))];

    annot.numeric$pGS <- as.numeric(as.factor(annot.numeric$pGS));
    annot.numeric[which(annot$pre.treatment.psa < 10), 'pre.treatment.psa'] = length(mycol.gs) + 1
    annot.numeric[which(annot$pre.treatment.psa >= 10 & annot$pre.treatment.psa < 20), 'pre.treatment.psa'] = length(mycol.gs) + 2
    annot.numeric[which(annot$pre.treatment.psa >= 20), 'pre.treatment.psa'] = length(mycol.gs) + 3
   	annot.numeric[,'pathological_t'] = length(mycol.gs) + 3 + as.numeric(as.factor(annot.numeric$pathological_t));
    annot.numeric[, 'idc'] <- length(mycol.gs) + 3 + length(mycol) + as.numeric(as.factor(annot.numeric$idc));
    annot.numeric[, 'pN'] <- max(annot.numeric[, 1:4]) + ifelse(annot.numeric$pN=='N0', 2, 1);
    annot.numeric[, 'Clinical_M'] <- max(annot.numeric[, 1:5]) + ifelse(annot.numeric$Clinical_M=='M1b', 1, 2);
    ####
    
	colour.sch <- c(
			mycol.gs,
            "#FEE6CE","#FDAE6B","#E6550D",
            mycol,
            'lightskyblue','firebrick3',
            'darkblue', 'lightskyblue',
            'darkgreen', 'lightskyblue'
            );
	cov.heatmap = create.heatmap(
	        x = annot.numeric,
	        clustering.method = 'none',
	        same.as.matrix = TRUE,
	        colour.scheme = colour.sch,
	        grid.col = TRUE,
	        grid.row = TRUE,
	        row.colour = 'white',
	        row.lwd = 3,
	        col.colour = 'white',
	        col.lwd = 3,
	        total.colours = length(colour.sch) + 1,
	        fill.colour = 'slategray',
	        print.colour.key = FALSE,
	        xat = rep('', 6),
	       	xat.top = TRUE,
	        );
};

