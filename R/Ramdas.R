#### DAG
DAG.hierarchy <- function(DAG){
    n <- length(V(DAG))
    nodes <- as.numeric(which(degree(DAG, mode = "in") < 1))
    new.nodes <- nodes
    hierarchy <- list()
    k <- 1
    hierarchy[[k]] <- nodes
    while (length(nodes) < n){
        new.nodes <- neighborhood(DAG, 1, new.nodes, "out",
                                  mindist = 1)
        new.nodes <- as.numeric(unique(unlist(new.nodes)))
        valid <- sapply(new.nodes, function(j){
            parents <- neighborhood(DAG, 1, j, "in",
                                    mindist = 1)
            parents <- as.numeric(unlist(parents))
            all(parents %in% nodes)
        })
        new.nodes <- new.nodes[valid]
        k <- k + 1
        hierarchy[[k]] <- new.nodes
        nodes <- c(nodes, new.nodes)
    }
    return(hierarchy)
}

eff.nleaves.nnodes <- function(DAG){
    hierarchy <- DAG.hierarchy(DAG)
    m <- length(hierarchy)
    n <- length(V(DAG))
    ell <- rep(1, n)
    nnodes <- rep(1, n)
    total.ell <- length(hierarchy[[m]])
    parents.list <- lapply(V(DAG), function(node){
        parents <- neighborhood(DAG, 1, nodes = node, mode = "in",
                                mindist = 1)[[1]]
        parents <- as.numeric(parents)
        return(parents)
    })
    nparents <- sapply(parents.list, length)
    children.list <- lapply(V(DAG), function(node){
        children <- neighborhood(DAG, 1, nodes = node,
                                 mode = "out",
                                 mindist = 1)[[1]]
        children <- as.numeric(children)
        return(children)
    })
    if (m == 1){
        return(list(ell = ell, nnodes = nnodes, 
                    total.ell = total.ell,
                    hierarchy = hierarchy,
                    parents.list = parents.list))
    }
    for (k in (m-1):1){
        nodes <- hierarchy[[k]]
        tmp.ell <- sapply(nodes, function(node){
            children <- children.list[[node]]
            val <- sum(ell[children] / nparents[children])
        })
        tmp.nnodes <- sapply(nodes, function(node){
            children <- children.list[[node]]
            val <- 1 + sum(nnodes[children] / nparents[children])
        })
        ell[nodes] <- tmp.ell
        nnodes[nodes] <- tmp.nnodes
    }
    return(list(ell = ell, nnodes = nnodes,
                total.ell = total.ell,
                hierarchy = hierarchy,
                parents.list = parents.list))
}

chen.ramdas <- function(pvals, DAG, 
                        alpha.list = seq(0.01, 0.3, 0.01)){
    n <- length(pvals)
    info <- eff.nleaves.nnodes(DAG)
    ell <- info$ell
    nnodes <- info$nnodes
    L <- info$total.ell
    hierarchy <- info$hierarchy
    parents.list <- info$parents.list    
    m <- length(hierarchy)
    rejs <- sapply(alpha.list, function(alpha){
        R <- 0
        rej <- rep(FALSE, n)
        for (k in 1:m){
            nodes <- hierarchy[[k]]
            valid <- sapply(nodes, function(node){
                parents <- parents.list[[node]]
                (length(parents) == 0) ||
                    (all(rej[parents]))
            })
            candids <- nodes[valid]
            alpha.mat <- sapply(candids, function(candid){
                nleaves <- ell[candid]
                nnodes <- nnodes[candid]
                alpha.ir <- nleaves / L * alpha *
                    (nnodes + R + 0: (length(candids) - 1)) / nnodes
                return(alpha.ir)
            })
            if (class(alpha.mat) != "matrix"){
                alpha.mat <- as.matrix(alpha.mat, ncol = 1)
            }
            tmp.pvals <- pvals[candids]
            tmp.rej <- stepup(tmp.pvals, alpha.mat)
            rej[candids] <- tmp.rej
            R <- R + sum(tmp.rej)
        }
        return(rej)
    })
    return(list(rej = rejs))
}
