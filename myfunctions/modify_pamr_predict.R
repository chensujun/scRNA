softmax <-function(x, gap = FALSE) {
  d <- dim(x)
  maxdist <- x[, 1]
  pclass <- rep(1, d[1])
  for(i in seq(2, d[2])) {
    l <- x[, i] > maxdist
    pclass[l] <- i
    maxdist[l] <- x[l, i]
  }
  dd <- dimnames(x)[[2]]
  if(gap) {
    x <- abs(maxdist - x)
    x[cbind(seq(d[1]), pclass)] <- drop(x %*% rep(1, d[2]))
    gaps <- do.call("pmin", data.frame(x))
  }
  pclass <- if(is.null(dd) || !length(dd))
    pclass
  else
    factor(pclass, levels = seq(d[2]), labels = dd)
  if(gap)
    list(class = pclass, gaps = gaps)
  else
    pclass
};

diag.disc <-function(x, centroids, prior, weight) {
### Computes the class discriminant functions assuming scaled x and centroids
  if(! missing(weight)) {
    posid <- (weight > 0)
    if(any(posid)) {
      weight <- sqrt(weight[posid])
      centroids <- centroids[posid,  , drop = FALSE] * weight
      x <- x[posid,  , drop = FALSE] * weight
    }
    else {
      mat <- outer(rep(1, ncol(x)), log(prior), "*")
      dimnames(mat) <- list(NULL, dimnames(centroids)[[2]])
      return(mat)
    }
  }
  dd <- t(x) %*% centroids
  dd0 <- drop(rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
  names(dd0) <- NULL
  scale(dd, dd0, FALSE)
};

soft.shrink <-function(delta, threshold) {
  dif <- abs(delta) - threshold
  delta <- sign(delta) * dif * (dif > 0)
  nonzero <- sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
  attr(delta, "nonzero") <- nonzero
  delta
}

pamr.predict <-  function(fit, newx, threshold, type = c("class", "posterior", "centroid", "nonzero"), 
                          prior = fit$prior,  threshold.scale = fit$
                          threshold.scale) {
  norm.cen <- fit$norm.cen
  if(!is.null(norm.cen)) {
    newx <- abs(t(scale(t(newx), center = norm.cen, scale = FALSE)))
  }
  type <- match.arg(type)
  sd <- fit$sd
  centroid.overall <- fit$centroid.overall
  centroids <- fit$centroids
  centroid.overall <- centroid.overall[names(centroid.overall)%in%rownames(newx)]
  centroids <- centroids[names(centroid.overall), ];
  sd <- sd[names(centroid.overall)]
  newx <- newx[rownames(centroids), ]
  se.scale <- fit$se.scale
  delta <- scale((centroids - centroid.overall)/sd, FALSE, threshold.scale * 
                 se.scale)

  if(fit$sign.contrast=="positive"){delta <- delta*(delta>0)}
  if(fit$sign.contrast=="negative"){delta <- delta*(delta<0)}


  delta.shrunk <- scale(soft.shrink(delta, threshold), FALSE,
                        1/(  threshold.scale * se.scale))
  posid <- drop(abs(delta.shrunk) %*% rep(1, length(prior)))
                
  if(!match(type, c("centroid", "nonzero"), FALSE))
    dd <- diag.disc((newx - centroid.overall)/sd, delta.shrunk, 
                    prior, posid)
  switch(type,
         class = softmax(dd),
         posterior = {
           dd <- exp(dd)
           dd/drop(dd %*% rep(1, length(prior)))
         }
         ,
         centroid = centroid.overall + delta.shrunk * sd,
         nonzero = {
           nz <- drop(abs(delta.shrunk) %*% rep(1, ncol(centroids)
                                                )) > 0
           seq(nz)[nz]
         }
         )
}


