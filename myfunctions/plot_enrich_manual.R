plot_enrich <- function(myego, name, width = 9, height = 10, xrot = 0,
        spot.size.function = function(x) {abs(x)/2}, mybg = rownames(pros13@scale.data), nterms = 15, 	key.sizes = c(1, 5, 10)){
	grep_data <- function(myego){
		my.bg <- data.frame(matrix(NA, nrow = 0, ncol = length(myego)));
	    colnames(my.bg) <- names(myego);
	    my.or <- my.bg;
		for(i in seq(length(myego))){
			ego <- myego[[i]];
			ipatient <- names(myego);
			print(ipatient);
			go <- ego[order(ego$FDR), ];
			bg.dat <- -log10(ego[1:nterms, 'FDR', drop = FALSE]);
			OR.dat <- ego[1:min(nterms, nrow(go)), 'OR', drop = FALSE];
			rownames(bg.dat) <- rownames(OR.dat) <- go$TermName[1:nrow(bg.dat)];
			my.bg[rownames(bg.dat), ipatient] <- bg.dat[, 1];
			my.or[rownames(OR.dat), ipatient] <- log2(OR.dat[, 1]);
		}
    return(list(bg = my.bg, or = my.or))
	};
	to.plot <- grep_data(myego);
	my.bg <- to.plot[['bg']];
    my.or <- to.plot[['or']];

    spot.colour.function <- function(x) {
        colours <- rep("white", length(x));
        colours[sign(x) == -1] <- default.colours(2, palette.type = "dotmap")[1]; 
        colours[sign(x) == 1] <- default.colours(2, palette.type = "dotmap")[2]; 
        return(colours);
    };

	dot.key <- list(
	    # indicate which side of the plot the legend will appear
	        space = "right",
	    points = list(
	        cex = spot.size.function(key.sizes),
	        col = spot.colour.function(key.sizes),
	        pch = 19
	        ),
	    # dot labels
	    text = list(
	        lab = as.character(key.sizes),
	        cex = 1,
	        adj = 1.5,
	        fontface = "bold"
	        ),
	    title = expression(underline('log'[2]*'OR')),
	    x = 3, 
	    y =7
	    );


    create.dotmap(
            file = generate.filename('dotmap_go', name, 'pdf'),
            x = my.or,
            xaxis.cex = 1,
            yaxis.cex = 1.2,
            left.padding = 0,
            bottom.padding = 4,
            # use specified spot size and colour functions
            spot.size.function = spot.size.function,
            spot.colour.function = spot.colour.function,
            # create a legend matching the dot sizes
            key = dot.key,
            key.top = 1,
            xaxis.lab = gsub('JD1800|SL', '', colnames(my.bg)),
            yaxis.lab = rownames(my.or),
            xaxis.rot = xrot,
            pch = 21,
            pch.border.col = 'transparent',
            # add the background
            bg.data = my.bg,
            # add a colourkey
            colourkey = TRUE,
            colour.scheme = c("white", "black"),
            total.colour = 5,
            bg.alpha = 1,
            at = c(0, -log10(0.05), 5, 10),
            colourkey.labels.at = c(0, -log10(0.05), 10, 50),
            colourkey.labels = c(1, expression(0.05), expression(10^-10), expression(''<=10^-50)),
            width = width,
            height = height,
            na.spot.size = 3,
            add.grid = TRUE,
            col.lwd = 1,
            style = 'Nature',
            col.colour = 'black', 
            row.colour = 'black', 
            );
};
