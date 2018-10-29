source('STAR_tree.R')
source('Lynch.R')
library("isotone")
library("ape")
library("structSSI")
library("igraph")
library("MASS")
source('summarize_methods.R')
source("expr_func.R")

rho <- as.numeric(Sys.getenv("rho"))
type <- as.numeric(Sys.getenv("type"))
repeat.times <- as.numeric(Sys.getenv("times"))
seed <- as.numeric(Sys.getenv("seed"))

tree.expr <- function(tree, mu, H0, rho, type,
                      alpha.list,
                      repeat.times){
    n <- length(mu)
    m <- length(alpha.list)
    summary.FDP <- list()
    summary.power <- list()
    for (j in 1:14){
        summary.FDP[[j]] <- matrix(rep(NA, repeat.times * m),
                                   ncol = repeat.times)
        summary.power[[j]] <- matrix(rep(NA, repeat.times * m),
                                     ncol = repeat.times)
    }
    STAR.fun.list.set <- list(
        create.fun("SeqStep", 0.5),
        create.fun("SeqStep", 0.6),
        create.fun("SeqStep", 0.7),
        create.fun("SeqStep", 0.8),
        create.fun("SeqStep", 0.9),
        create.fun("ForwardStop"),
        create.fun("HingeExp", 0.5),
        create.fun("HingeExp", 0.7),
        create.fun("HingeExp", 0.9)
        )
    for (i in 1:repeat.times){
        pvals.flag <- FALSE
        for (ii in 1:10){
            pvals <- pvals.gen(n, mu, rho, type)
            pvals.flag <- (length(pvals) == n) &&
                (all(!is.na(pvals))) &&
                (all(pvals >= 0))
            if (pvals.flag){
                break
            }
        }
        if (!pvals.flag){
            next
        }

        for (k in 1:9){
            fun <- STAR.fun.list.set[[k]]
            STARtree <- try(STAR.tree(pvals, tree, print.quiet = TRUE,
                                     alpha.list = alpha.list,
                                     fun.list = fun))
            if (class(STARtree) != "try-error"){
                STAR.tree.result <- summary.STAR(STARtree, H0)
                summary.FDP[[k]][, i] <- STAR.tree.result[, 2]
                summary.power[[k]][, i] <- STAR.tree.result[, 3]
            }
        }
    
        names(pvals) <- 1:n
        V(tree.backup)$name <- names(pvals)
        tree.el <- get.edgelist(tree.backup)
        hyptree.result <- try(summary.hyptree(pvals, tree.el, H0,
                                              alpha.list))
        if (class(hyptree.result) != "try-error"){
            summary.FDP[[10]][, i] <- hyptree.result[, 2]
            summary.power[[10]][, i] <- hyptree.result[, 3]
        }

        BH.result <- summary.BH(pvals, H0, alpha.list)
        summary.FDP[[11]][, i] <- BH.result[, 2]
        summary.power[[11]][, i] <- BH.result[, 3]

        storey.result <- summary.storey(pvals, H0,
                                        alpha.list = alpha.list)
        summary.FDP[[12]][, i] <- storey.result[, 2]
        summary.power[[12]][, i] <- storey.result[, 3]

        lg.result1 <- try(summary.lg(pvals, tree, H0,
                                     alpha.list = alpha.list,
                                     type = 1))
        if (class(lg.result1) != "try-error"){
            summary.FDP[[13]][, i] <- lg.result1[, 2]
            summary.power[[13]][, i] <- lg.result1[, 3]
        }

        lg.result2 <- try(summary.lg(pvals, tree, H0,
                                     alpha.list = alpha.list,
                                     type = 2))
        if (class(lg.result2) != "try-error"){
            summary.FDP[[14]][, i] <- lg.result2[, 2]
            summary.power[[14]][, i] <- lg.result2[, 3]
        }
        
        print(i)
    }
    avg.FDP <- lapply(summary.FDP, function(FDP){
        apply(FDP, 1, function(x){mean(x, na.rm = TRUE)})
    })
    avg.power <- lapply(summary.power, function(power){
        apply(power, 1, function(x){mean(x, na.rm = TRUE)})
    })
    return(list(FDP = avg.FDP, power = avg.power))
}

####### Generate Tree 
set.seed(seed)
n <- 1000
ns <- 50
tree <- make_tree(n)
tree.backup <- tree
V(tree)$name <- 1:n
alpha.list <- seq(0.01, 0.3, 0.01)
output.filename <- paste0("../data/data_tree_rho_", floor(rho*10), "_type_", floor(type), "_seed_", seed, ".RData")

## Breath-First Tree
H0.ind <- 1:ns
H0 <- rep(FALSE, n)
H0[H0.ind] <- TRUE
mu <- rep(0, n)
mu[H0.ind] <- 2

result.BFS.tree <- tree.expr(tree, mu, H0, rho, type, alpha.list,
                             repeat.times)
data.tree <- list(result.BFS.tree)
save(file = output.filename, data.tree)


## Depth First Tree
H0.ind <- sort(dfs(tree, root = 1)$order[1:ns])
H0 <- rep(FALSE, n)
H0[H0.ind] <- TRUE
mu <- rep(0, n)
mu[H0.ind] <- 2

result.DFS.tree <- tree.expr(tree, mu, H0, rho, type, alpha.list,
                             repeat.times)
data.tree <- list(result.BFS.tree, result.DFS.tree)
save(file = output.filename, data.tree)

## Breath-First IsoTree
H0.ind <- 1:ns
H0 <- rep(FALSE, n)
H0[H0.ind] <- TRUE
mu <- rep(0, n)
mu[H0.ind[1:floor(ns*0.5)]] <- 2.5
mu[H0.ind[(floor(ns*0.5)+1):ns]] <- 1.5

result.BFS.isotree <- tree.expr(tree, mu, H0, rho, type, alpha.list,
                                repeat.times)
data.tree <- list(result.BFS.tree, result.DFS.tree,
                  result.BFS.isotree)
save(file = output.filename, data.tree)

## Depth First IsoTree
H0.ind <- sort(dfs(tree, root = 1)$order[1:ns])
H0 <- rep(FALSE, n)
H0[H0.ind] <- TRUE
mu <- rep(0, n)
mu[H0.ind[1:floor(ns*0.5)]] <- 2.5
mu[H0.ind[(floor(ns*0.5)+1):ns]] <- 1.5

result.DFS.isotree <- tree.expr(tree, mu, H0, rho, type, alpha.list,
                             repeat.times)
data.tree <- list(result.BFS.tree, result.DFS.tree,
                  result.BFS.isotree, result.DFS.isotree)
save(file = output.filename, data.tree)

## Breath-First AniIsoTree
H0.ind <- 1:ns
H0 <- rep(FALSE, n)
H0[H0.ind] <- TRUE
mu <- rep(0, n)
mu[H0.ind[1:floor(ns*0.5)]] <- 1.5
mu[H0.ind[(floor(ns*0.5)+1):ns]] <- 2.5

result.BFS.aniisotree <- tree.expr(tree, mu, H0, rho, type, alpha.list,
                                   repeat.times)
data.tree <- list(result.BFS.tree, result.DFS.tree,
                  result.BFS.isotree, result.DFS.isotree,
                  result.BFS.aniisotree)
save(file = output.filename, data.tree)

## Depth First AniIsoTree
H0.ind <- sort(dfs(tree, root = 1)$order[1:ns])
H0 <- rep(FALSE, n)
H0[H0.ind] <- TRUE
mu <- rep(0, n)
mu[H0.ind[1:floor(ns*0.5)]] <- 1.5
mu[H0.ind[(floor(ns*0.5)+1):ns]] <- 2.5

result.DFS.aniisotree <- tree.expr(tree, mu, H0, rho, type, alpha.list,
                                   repeat.times)

data.tree <- list(result.BFS.tree, result.DFS.tree,
                  result.BFS.isotree, result.DFS.isotree,
                  result.BFS.aniisotree, result.DFS.aniisotree)
save(file = output.filename, data.tree)

