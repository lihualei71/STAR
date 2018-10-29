##############################################################################
## STAR for Testing on DAGs
##################################################################
############

source("generic_STAR.R")

#### Find_Candidate_Fun
find.DAG.strong <- function(covar, mask, prop = NULL,
                            order = 1){
    DAG <- induced_subgraph(covar, mask)
    n <- sum(mask)
    leaves <- V(DAG)$name[degree(DAG, mode = "out") < 1]    
    if (order == 1){
        return(as.list(as.numeric(leaves)))
    }
    ## candids <- neighborhood(DAG, order = n,
    ##                         nodes = V(DAG), mode = "out")
    ## candids <- lapply(candids, function(candid){
    ##     as.numeric(names(candid))
    ## })
    return(candids)
}

find.DAG.weak <- function(covar, mask, prop = NULL){
    DAG <- induced_subgraph(covar, mask)
    n <- sum(mask)
    child.degree <- neighborhood(DAG, order = 1, nodes = V(DAG),
                                 mode = "in", mindist = 1)
    child.degree <- sapply(child.degree, length)
    invalid.child <- V(DAG)$name[child.degree == 1]
    invalid.node <- neighborhood(DAG, order = 1,
                                 nodes = invalid.child,
                                 mode = "in", mindist = 1)
    invalid.node <- sapply(invalid.node, function(node){
        names(node)
    })
    candids <- setdiff(V(DAG)$name, invalid.node)
    candids <- as.list(as.numeric(candids))
    ## candids <- lapply(V(DAG), function(vertex){
    ##     subtree <- neighborhood(DAG, order = n,
    ##                             nodes = vertex, mode = "out")[[1]]
    ##     if (length(subtree) == 1){
    ##         return(subtree)
    ##     } else {
    ##         valid <- sapply(subtree[-1], function(node){
    ##             parents <-
    ##                 neighborhood(DAG, order = 1,
    ##                              nodes = node, mode = "in")[[1]]
    ##             all(parents %in% subtree)
    ##         })
    ##     }
    ##     subtree[c(TRUE, valid)]
    ## })
    ## candids <- lapply(candids, function(candid){
    ##     as.numeric(names(candid))
    ## })
    return(candids)
}

#### Update_Mask_Fun
## Naive Update
DAG.mask.update <- function(candid, score, mask, prop){
    num.update <- min(max(ceiling(sum(mask) * prop), 1),
                      length(candid))
    candid.vals <- sapply(candid, function(set){
        mean(score[set])
    })
    thresh <- sort(candid.vals, decreasing = TRUE)[num.update]
    reveal.inds <- unique(unlist(candid[candid.vals >= thresh]))
    mask[reveal.inds] <- FALSE
    return(mask)
}

#### Plot_Fun
plot.rej.DAG <- function(DAG, mask, main, layers = NULL, 
                         col.bg = "#000000", col.fg = "#FFB6C1",
                         vertex.size = 7,
                         edge.arrow.size = 0.5,
                         ...){
    color <- ifelse(mask, col.fg, col.bg)
    layout <- layout_with_sugiyama(DAG, layers,
                                   attributes="all")
    plot(layout$extd_graph, main = main,
         vertex.size = vertex.size,
         edge.arrow.size = edge.arrow.size,
         vertex.color = color, ...)
}

DAG.plot <- function(covar, mask, score, main,
                     add.plot = NULL,
                     cex.main = 2, ...){
    par(cex.main = cex.main)    
    if (!is.null(add.plot)){
        par(mfrow = c(1, 2))
        eval(add.plot)
    }
    plot.rej.DAG(covar, mask, main = main, ...)
}

####
STAR.DAG <- function(
    pvals, DAG,
    fun.list = create.fun("SeqStep", pstar = 0.5),
    alpha.list = seq(0.05, 0.3, 0.01),
    type = c("naive", "knockoff"),
    score = NULL,
    criterion = c("strong", "weak"),
    plot.fun = DAG.plot,
    update.mask.fun = DAG.mask.update,
    prop.high = 0.05,
    prop.low = 0.01,
    ...){

    type <- type[1]
    criterion <- criterion[1]
    if (criterion == "strong"){
        find.candid.fun <- find.DAG.strong
    } else if (criterion == "weak"){
        find.candid.fun <- find.DAG.weak
        prop.high = 1 / length(pvals)
        prop.low = 1 / length(pvals)
    }

    if (type == "naive"){
        score.fun <- NULL
    } else if (type == "knockoff"){
        score.fun <- function(covar, pvals, mask,
                              fun.list, score0){
            score
        }
    }

    generic.STAR(pvals = pvals, covar = DAG, type = type,
                fun.list = fun.list, alpha.list = alpha.list,
                score.fun = score.fun,
                find.candid.fun = find.candid.fun,
                update.mask.fun = update.mask.fun,
                plot.fun = plot.fun,
                prop.high = prop.high,
                prop.low = prop.low,
                ...)
}
