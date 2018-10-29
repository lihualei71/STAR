##############################################################################
## STAR for Hierarchical Testing
##############################################################################

source("generic_STAR.R")

#### Find_Candidate_Fun
find.leaves <- function(covar, mask, prop = NULL){
    tree <- induced_subgraph(covar, mask)
    leaves <- V(tree)$name[degree(tree, mode = "out") < 1]
    leaves <- which(V(covar)$name %in% leaves)
    return(as.list(leaves))
}

#### Update_Mask_Fun
## Naive Update
tree.mask.update <- function(candid, score, mask, prop){
    num.update <- max(ceiling(sum(mask) * prop), 1)
    candid.vals <- sapply(candid, function(set){
        mean(score[set])
    })
    thresh <- sort(candid.vals, decreasing = TRUE)[num.update]
    reveal.inds <- unique(unlist(candid[candid.vals >= thresh]))
    mask[reveal.inds] <- FALSE
    return(mask)
}

#### Score_Fun
## Isotonic tree regression
isoreg.censor <- function(covar, pvals, mask, fun.list,
                          score0 = NULL, num.steps = 10){
    edges <- get.edgelist(covar)[, c(2, 1)]
    n <- length(pvals)
    weights <- rep(1, n)
    h <- fun.list$h
    g <- fun.list$g
    sinv <- fun.list$sinv
    s.deriv <- fun.list$s.deriv
    Const <- fun.list$Const
    init.tdpvals <- pmin(pvals, g(pvals))
    init.ref.tdpvals <- pmax(pvals, g(pvals))
    init.tdpvals <- pmin(pvals, g(pvals))
    init.ref.tdpvals <- pmax(pvals, g(pvals))

    if (is.null(score0)){
        mu0 <- mean(-log(pvals))
        mu <- rep(mu0, n)
    } else {
        mu <- score0
    }
    tdpvals <- ifelse(mask, init.tdpvals, pvals)
    ref.tdpvals <- ifelse(mask, init.ref.tdpvals, pvals)
    
    for (i in 1:num.steps){
        temp <- -1 / s.deriv(ref.tdpvals)
        imputed.logpvals <- ifelse(
            mask,
            (tdpvals^(1/mu-1)*(-log(tdpvals))+ref.tdpvals^(1/mu-1)*temp*(-log(ref.tdpvals))) / (tdpvals^(1/mu-1)+ref.tdpvals^(1/mu-1)*temp),
            -log(tdpvals))
        fit <- activeSet(edges, "poisson", y = imputed.logpvals,
                         weights = weights, x0 = mu, ups = 1e-2)
        mu.new <- fit$x

        if (max(abs(mu.new - mu)) < 0.01) {
            converge <- TRUE
            break
        }
        mu <- pmax(mu.new, 1)
    }
    return(mu)
}


#### Plot_Fun
plot.rej.tree <- function(tree, mask, layout, main,
                          col.bg = "#000000", col.fg = "#FFB6C1",
                          vertex.size = 8,
                          edge.arrow.size = 0,
                          ...){
    color <- ifelse(mask, col.fg, col.bg)
    plot(tree, main = main, layout = layout,
         vertex.size = vertex.size,
         edge.arrow.size = edge.arrow.size,
         vertex.color = color, ...)
}

## plot.strength.tree <- function(tree, strength, layout, main,
##                                col.fg = "#000080",
##                                col.bg = "#ADD8E6",
##                                vertex.size = 8,
##                                edge.arrow.size = 0,
##                                ...){
##     rbPal <- colorRampPalette(c(col.bg, col.fg))
##     color <- rbPal(100)[as.numeric(cut(strength, breaks=100))]
##     plot(tree, main = main, layout = layout,
##          vertex.size = vertex.size,
##          edge.arrow.size = edge.arrow.size,
##          vertex.color = color, ...)
## }

tree.plot <- function(covar, mask, score, main,
                      add.plot = NULL,
                      cex.main = 2, ...){
    par(mar = c(2, 3, 3, 2), cex.main = cex.main)
    if (is.null(add.plot)){
        par(mfrow = c(1, 2))
    } else {
        par(mfrow = c(2, 2))
        for (plot in add.plot){
            eval(plot)
        }
    }
    layout1 <- layout_as_tree(covar, circular = TRUE)
    plot.rej.tree(covar, mask, layout1, main = main, ...)
    layout2 <- layout_as_tree(covar, circular = FALSE)
    plot.rej.tree(covar, mask, layout2, main = main, ...)    
}

####
STAR.tree <- function(
    pvals, tree,
    fun.list = create.fun("SeqStep", pstar = 0.5),
    alpha.list = seq(0.05, 0.3, 0.01),
    type = c("naive", "model-assist"),
    ...){

    type <- type[1]
    find.candid.fun <- find.leaves
    plot.fun <- tree.plot
    update.mask.fun <- tree.mask.update    
    if (type == "naive"){
        score.fun <- NULL
    } else if (type == "model-assist"){
        score.fun <- isoreg.censor
    }

    generic.STAR(pvals = pvals, covar = tree, type = type,
                fun.list = fun.list, alpha.list = alpha.list,
                score.fun = score.fun,
                find.candid.fun = find.candid.fun,
                update.mask.fun = update.mask.fun,
                plot.fun = plot.fun,
                ...)
}
