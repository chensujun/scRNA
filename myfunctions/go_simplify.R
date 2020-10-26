library(GOSemSim);
library(tidyr);
go_simplify <- function(res, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", ontology, semData) {
    res$ID <- as.vector(res$ID);
    if (missing(semData) || is.null(semData)) {
        if (measure == "Wang") {
            semData <- godata(ont = ontology)
        } else {
            stop("godata should be provided for IC-based methods...")
        }
    } else {
        if (ontology != semData@ont) {
            msg <- paste("semData is for", semData@ont, "ontology, while enrichment result is for", ontology)
            stop(msg)
        }
    }

    sim <- mgoSim(res$ID, res$ID,
                  semData = semData,
                  measure=measure,
                  combine=NULL)

    ## to satisfy codetools for calling gather
    go1 <- go2 <- similarity <- NULL

    sim.df <- as.data.frame(sim)
    sim.df$go1 <- row.names(sim.df)
    sim.df <- gather(sim.df, go2, similarity, -go1)

    sim.df <- sim.df[!is.na(sim.df$similarity),]

    ## feature 'by' is attached to 'go1'
    sim.df <- merge(sim.df, res[, c("ID", by)], by.x="go1", by.y="ID")
    sim.df$go2 <- as.character(sim.df$go2)

    ID <- res$ID

    GO_to_remove <- character()
    for (i in seq_along(ID)) {
        ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
        ## if length(ii) == 1, then go1 == go2
        if (length(ii) < 2)
            next

        sim_subset <- sim.df[ii,]

        jj <- which(sim_subset[, by] == select_fun(sim_subset[, by]))

        ## sim.df <- sim.df[-ii[-jj]]
        GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
    }

    res[!res$ID %in% GO_to_remove, ]
};
