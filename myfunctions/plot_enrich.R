####
# this script takes gprofiler output and generate dotplot for the enrichment results
# this script is built upon the function create.dotmap from R package BoutrosLab.plotting.general (BPG)
# note the function filter for terms with size >=10 and <=500, modify the script directly to disable this
# myego: a list of gprofiler result(s) generated using the 'gprofiler' function
# gene.go: reference for the ontology/pathway that you used for the enrichemnt, you can download the reference gmt file from gprofiler website and load it with read.gmt function from clusterProfiler
# mydomain:the type of ontology that you want to plot, e.g. BP/CC/MF, for a full list, see'src_filter' parameter from 'gprofiler' function
# name: name stem for output plot
# width & height: adjust for optimal visualization effect
# mybg: background that you used for gene enrichment analysis, default to all genes in gene.go
# spot.size.function: modify this to control for spot size
# noOR: if TRUE show the ratio of overlap.size/query.size instead of odds ratio, default to FALSE
# the following show an example for how to use the function
# load('2020-07-07_test_data_plot_enrich.rda');
# myego <- list(up = ego.up, dn = ego.dn);
# plot_enrich(myego, gene.go, 'keg', 
#   'test', 10, 8, mybg = gene.go$gene, 
#   xrot = 90, noOR = TRUE, spot.size.function = function(x) {abs(x)*2}
# );
# sometimes local R can have problem output plot files
# you can cheat solve the problem by doing the following:
# annotate out line: 89, 90, 124
# de-annotate lines: 91,92, 128-130
# then you'll get a output, the font can be wierd but you can adjust it manually afterwards
# 
plot_enrich <- function(myego, gene.go, mydomain, name, width = 9, height = 10, xrot = 0,
        spot.size.function = function(x) {abs(x)/2}, mybg = gene.go$gene, noOR = FALSE, 
        plot.legend = FALSE, key.sizes = c(1, 5, 10)){
    #gene.go <- readRDS('~/circRNA/star-circ/circRNA_landscape/rnaseq_landscape/scRNA/2019-01-05_gp_hsapiens.GO.Name.rds');
    my.ego <- myego;
    my.bg <- data.frame(matrix(NA, nrow = 0, ncol = length(my.ego)));
    colnames(my.bg) <- names(my.ego);
    my.or <- my.bg;
    calc_data <- function(my.ego){
        for(i in seq(length(my.ego))){
        	ego <- my.ego[[i]];
            ipatient <- names(my.ego)[i];
            print(ipatient)
        #	go <- ego[ego$domain%in%c('rea')&ego$significant==TRUE&ego$relative.depth == 4, 1:13];
        	go <- ego[ego$domain%in%mydomain&ego$significant==TRUE&ego$term.size>=10&ego$term.size<=500, 1:13];
            if(nrow(go)>0){
            	go <- go[order(go$p.value), ];
            	nbg <- length(intersect(mybg, gene.go$gene))
            	go$log2OR <- apply(go[, 4:6],1, function(x) log2((x[3]/(x[2]-x[3]))/(x[1]/(nbg-x[1]))));
            	bg.dat <- -log10(go[1:min(10, nrow(go)), 'p.value', drop = FALSE]);
            	OR.dat <- go[1:min(10, nrow(go)), 'log2OR', drop = FALSE];
            	rownames(bg.dat) <- go$term.name[1:min(10, nrow(bg.dat))];
            	rownames(OR.dat) <- go$term.name[1:min(10, nrow(bg.dat))];
                my.bg[rownames(bg.dat), ipatient] <- bg.dat[, 1];
                my.or[rownames(OR.dat), ipatient] <- OR.dat[, 1];
            }
        };
    return(list(bg = my.bg, or = my.or))
    };
    calc_data_noOR <- function(my.ego){
        for(i in seq(length(my.ego))){
            ego <- my.ego[[i]];
            ipatient <- names(my.ego)[i];
            print(ipatient)
            go <- ego[ego$domain%in%mydomain&ego$significant==TRUE&ego$term.size>=10&ego$term.size<=500, 1:13];
            if(nrow(go)>0){
                go <- go[order(go$p.value), ];
                go$log2OR <- go$overlap.size/go$query.size;
                bg.dat <- -log10(go[1:min(10, nrow(go)), 'p.value', drop = FALSE]);
                OR.dat <- go[1:min(10, nrow(go)), 'log2OR', drop = FALSE];
                rownames(bg.dat) <- go$term.name[1:min(10, nrow(bg.dat))];
                rownames(OR.dat) <- go$term.name[1:min(10, nrow(bg.dat))];
                my.bg[rownames(bg.dat), ipatient] <- bg.dat[, 1];
                my.or[rownames(OR.dat), ipatient] <- OR.dat[, 1];
            }
        };
    return(list(bg = my.bg, or = my.or))
    }; 
    if(noOR == TRUE){
        to.plot <- calc_data_noOR(my.ego);
        print('NO OR calculation')

    }else{
        to.plot <- calc_data(my.ego)
    };
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
    y =7,
    height = 10
    );

    if(plot.legend){
        legend.col <- TRUE;
            create.dotmap(
            file = generate.filename(paste0('dotmap_go_', mydomain), name, 'pdf'),
#    p <-    create.dotmap(
            #file = NULL,
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
            colourkey = FALSE,
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

    }else{
        legend.col = FALSE;
        create.dotmap(
            file = generate.filename(paste0('dotmap_go_', mydomain), name, 'pdf'),
#    p <-    create.dotmap(
            #file = NULL,
            x = my.or,
            xaxis.cex = 1,
            yaxis.cex = 1.2,
            left.padding = 0,
            bottom.padding = 4,
            # use specified spot size and colour functions
            spot.size.function = spot.size.function,
            spot.colour.function = spot.colour.function,
            # create a legend matching the dot sizes
            #key = dot.key,
            key.top = 1,
            xaxis.lab = gsub('JD1800|SL', '', colnames(my.bg)),
            yaxis.lab = rownames(my.or),
            xaxis.rot = xrot,
            pch = 21,
            pch.border.col = 'transparent',
            # add the background
            bg.data = my.bg,
            # add a colourkey
            colourkey = FALSE,
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
    }
#    pdf(generate.filename(paste0('dotmap_go_', mydomain), name, 'pdf'));
#    print(p);
#    dev.off()
};