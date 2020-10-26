grep_data <- function(myego, npath = 5, cut.fdr = 0.25, Nmin = 0, Nmax = 10000){
        my.bg <- data.frame(matrix(NA, nrow = 0, ncol = length(myego)));
        colnames(my.bg) <- names(myego);
        my.or <- my.bg;
        for(i in seq(length(myego))){
                ego <- myego[[i]];
        ego <- ego[ego$FDR<cut.fdr&ego$nPath>Nmin&ego$nPath<Nmax, ]
                #ego <- ego[ego$FDR<0.05, ];
                if(nrow(ego)>0){
                        ego <- ego[order(-ego$Enrichment), ];
                        bg.dat <- -log10(ego[1:min(npath, nrow(ego)), 'FDR', drop = FALSE]);
                        or.dat <- ego[1:min(npath, nrow(ego)), 'Enrichment', drop = FALSE];
                        rownames(bg.dat) <- ego[1:min(npath, nrow(ego)), ]$TermName;
                        rownames(or.dat) <- rownames(bg.dat);
                        my.bg[rownames(bg.dat), names(myego)[i]] <- bg.dat[, 1];
                        my.or[rownames(or.dat), names(myego)[i]] <- or.dat[, 1];

                };
        };
      for(i in seq(length(myego))){
            ego <- myego[[i]];
            my.bg[, names(myego)[i]] <- -log10(ego[match(rownames(my.bg), ego$TermName), 'FDR'])
            my.or[, names(myego)[i]] <- ego[match(rownames(my.or), ego$TermName), 'Enrichment'];
      }
        return(list(bg = my.bg, or = my.or))
};

grep_data_go <- function(myego, npath = 5, cut.fdr = 0.25, Nmin = 0, Nmax = 10000){
        my.bg <- data.frame(matrix(NA, nrow = 0, ncol = length(myego)));
        colnames(my.bg) <- names(myego);
        my.or <- my.bg;
        for(i in seq(length(myego))){
                ego <- myego[[i]];
		ego <- ego[ego$FDR<cut.fdr&ego$GeneInGO>Nmin&ego$GeneInGO<Nmax, ]
                #ego <- ego[ego$FDR<0.05, ];
                if(nrow(ego)>0){
                        ego <- ego[order(-ego$Enrichment), ];
                        bg.dat <- -log10(ego[1:min(npath, nrow(ego)), 'FDR', drop = FALSE]);
                        or.dat <- ego[1:min(npath, nrow(ego)), 'Enrichment', drop = FALSE];
                        rownames(bg.dat) <- ego[1:min(npath, nrow(ego)), ]$GOTerm;
                        rownames(or.dat) <- rownames(bg.dat);
                        my.bg[rownames(bg.dat), names(myego)[i]] <- bg.dat[, 1];
                        my.or[rownames(or.dat), names(myego)[i]] <- or.dat[, 1];

                };
        };
      for(i in seq(length(myego))){
            ego <- myego[[i]];
            my.bg[, names(myego)[i]] <- -log10(ego[match(rownames(my.bg), ego$GOTerm), 'FDR'])
            my.or[, names(myego)[i]] <- ego[match(rownames(my.or), ego$GOTerm), 'Enrichment'];
      }
        return(list(bg = my.bg, or = my.or))
};

grep_data_keg <- function(myego, npath = 5, cut.fdr = 0.25, Nmin = 0, Nmax = 10000){
        my.bg <- data.frame(matrix(NA, nrow = 0, ncol = length(myego)));
        colnames(my.bg) <- names(myego);
        my.or <- my.bg;
        for(i in seq(length(myego))){
                ego <- myego[[i]];
                ego <- ego[ego$FDR<cut.fdr&ego$GeneInPathway>Nmin&ego$GeneInPathway<Nmax, ]
                #ego <- ego[ego$FDR<0.05, ];
                if(nrow(ego)>0){
                        ego <- ego[order(-ego$Enrichment), ];
                        bg.dat <- -log10(ego[1:min(npath, nrow(ego)), 'FDR', drop = FALSE]);
                        or.dat <- ego[1:min(npath, nrow(ego)), 'Enrichment', drop = FALSE];
                        rownames(bg.dat) <- ego[1:min(npath, nrow(ego)), ]$PathwayTerm;
                        rownames(or.dat) <- rownames(bg.dat);
                        my.bg[rownames(bg.dat), names(myego)[i]] <- bg.dat[, 1];
                        my.or[rownames(or.dat), names(myego)[i]] <- or.dat[, 1];

                };
        };
      for(i in seq(length(myego))){
            ego <- myego[[i]];
            my.bg[, names(myego)[i]] <- -log10(ego[match(rownames(my.bg), ego$PathwayTerm), 'FDR'])
            my.or[, names(myego)[i]] <- ego[match(rownames(my.or), ego$PathwayTerm), 'Enrichment'];
      }
        return(list(bg = my.bg, or = my.or))
};


