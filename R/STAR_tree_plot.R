source("FDR_power_plot.R")
source("expr_func.R")
source("STAR_tree.R")

library("igraph")

## Plot of truth
pdf("../figs/tree_expr_truth.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4, 2, 2, 2))
n <- 1000
ns <- 50
tree <- make_tree(n)
V(tree)$name <- 1:n
layout <- layout_as_tree(tree, circular = FALSE)
mains <- c("BFS", "DFS")
for (i in 1:2){
    if (i == 1){
        H0.ind <- 1:ns
    } else {
        H0.ind <- sort(dfs(tree, root = 1)$order[1:ns])
    }
    H0 <- rep(FALSE, n)
    H0[H0.ind] <- TRUE
    plot.rej.tree(tree, H0, layout,
                  main = mains[i],
                  vertex.size = 5, vertex.label = NA)
}
dev.off()

load("../data/data_tree_0.RData")
result.BFS.tree <- data.tree[[1]]
result.DFS.tree <- data.tree[[2]]
result.BFS.isotree <- data.tree[[3]]
result.DFS.isotree <- data.tree[[4]]
result.BFS.aniisotree <- data.tree[[5]]
result.DFS.aniisotree <- data.tree[[6]]


## Compare Methods
fdp.BFS.expr1 <- list(result.BFS.tree[[1]][c(1, 10, 13, 14, 11)],
                      result.BFS.isotree[[1]][c(1, 10, 13, 14, 11)],
                      result.BFS.aniisotree[[1]][c(1, 10, 13, 14, 11)])
plot.FDR("../figs/rej_tree_BFS_methods_FDR.pdf",
         fdp.BFS.expr1,
         mains = c("BFS (Case 1)",
             "BFS (Case 2)",
             "BFS (Case 3)"),
         methods = c("STAR", "Yekutieli", "LG", "LG2", "BH"),
         cols = c("black", "red", "blue", "magenta", "orange"),
         )

power.BFS.expr1 <- list(result.BFS.tree[[2]][c(1, 10, 13, 14, 11)],
                        result.BFS.isotree[[2]][c(1, 10, 13, 14, 11)],
                        result.BFS.aniisotree[[2]][c(1, 10, 13, 14, 11)])
plot.power("../figs/rej_tree_BFS_methods_power.pdf",
           power.BFS.expr1,
           mains = c("BFS (Case 1)",
               "BFS (Case 2)",
               "BFS (Case 3)"),
           methods = c("STAR", "Yekutieli", "LG", "LG2", "BH"),
           cols = c("black", "red", "blue", "magenta", "orange"),
           )

fdp.DFS.expr1 <- list(result.DFS.tree[[1]][c(1, 10, 13, 14, 11)],
                      result.DFS.isotree[[1]][c(1, 10, 13, 14, 11)],
                      result.DFS.aniisotree[[1]][c(1, 10, 13, 14, 11)])
plot.FDR("../figs/rej_tree_DFS_methods_FDR.pdf",
         fdp.DFS.expr1,
         mains = c("DFS (Case 1)",
             "DFS (Case 2)",
             "DFS (Case 3)"),
         methods = c("STAR", "Yekutieli", "LG", "LG2", "BH"),
         cols = c("black", "red", "blue", "magenta", "orange"),
         )

power.DFS.expr1 <- list(result.DFS.tree[[2]][c(1, 10, 13, 14, 11)],
                  result.DFS.isotree[[2]][c(1, 10, 13, 14, 11)],
                  result.DFS.aniisotree[[2]][c(1, 10, 13, 14, 11)])
plot.power("../figs/rej_tree_DFS_methods_power.pdf",
           power.DFS.expr1,
           mains = c("DFS (Case 1)",
               "DFS (Case 2)",
               "DFS (Case 3)"),
           methods = c("STAR", "Yekutieli", "LG", "LG2", "BH"),
           cols = c("black", "red", "blue", "magenta", "orange"),
           )


## Compare parameters of ISS
fdp.BFS.expr2 <- list(result.BFS.tree[[1]][1:5],
                      result.BFS.isotree[[1]][1:5],
                      result.BFS.aniisotree[[1]][1:5])
plot.FDR("../figs/rej_tree_BFS_STARSS_FDR.pdf",
         fdp.BFS.expr2, 
         mains = c("BFS (Case 1)",
             "BFS (Case 2)",
             "BFS (Case 3)"),
         methods = c("p*=0.5", "p*=0.6", "p*=0.7", "p*=0.8", "p*=0.9"),
         cols = c("black", "red", "blue", "orange", "magenta")
         )

power.BFS.expr2 <- list(result.BFS.tree[[2]][1:5],
                        result.BFS.isotree[[2]][1:5],
                        result.BFS.aniisotree[[2]][1:5])
plot.power("../figs/rej_tree_BFS_STARSS_power.pdf",
           power.BFS.expr2, 
           mains = c("BFS (Case 1)",
               "BFS (Case 2)",
               "BFS (Case 3)"),
           methods = c("p*=0.5", "p*=0.6", "p*=0.7", "p*=0.8", "p*=0.9"),
           cols = c("black", "red", "blue", "orange", "magenta")
           )

fdp.DFS.expr2 <- list(result.DFS.tree[[1]][1:5],
                      result.DFS.isotree[[1]][1:5],
                      result.DFS.aniisotree[[1]][1:5])
plot.FDR("../figs/rej_tree_DFS_STARSS_FDR.pdf",
         fdp.DFS.expr2, 
         mains = c("DFS (Case 1)",
             "DFS (Case 2)",
             "DFS (Case 3)"),
         methods = c("p*=0.5", "p*=0.6", "p*=0.7", "p*=0.8", "p*=0.9"),
         cols = c("black", "red", "blue", "orange", "magenta")
         )

power.DFS.expr2 <- list(result.DFS.tree[[2]][1:5],
                        result.DFS.isotree[[2]][1:5],
                        result.DFS.aniisotree[[2]][1:5])
plot.power("../figs/rej_tree_DFS_STARSS_power.pdf",
           power.DFS.expr2, 
           mains = c("DFS (Case 1)",
               "DFS (Case 2)",
               "DFS (Case 3)"),
           methods = c("p*=0.5", "p*=0.6", "p*=0.7", "p*=0.8", "p*=0.9"),
           cols = c("black", "red", "blue", "orange", "magenta")
           )


## Compare h-function 
fdp.BFS.expr3 <- list(result.BFS.tree[[1]][c(1, 6:9)],
                      result.BFS.isotree[[1]][c(1, 6:9)],
                      result.BFS.aniisotree[[1]][c(1, 6:9)])
plot.FDR("../figs/rej_tree_BFS_h_FDR.pdf",
         fdp.BFS.expr3,
         mains = c("BFS (Case 1)",
             "BFS (Case 2)",
             "BFS (Case 3)"),
         methods = c("SS (p*=0.5)", "FS", "HE (p*=0.5)",
             "HE (p*=0.7)", "HE (p*=0.9)"),
         cols = c("black", "red", "blue", "orange",
             "magenta")
         )

power.BFS.expr3 <- list(result.BFS.tree[[2]][c(1, 6:9)],
                        result.BFS.isotree[[2]][c(1, 6:9)],
                        result.BFS.aniisotree[[2]][c(1, 6:9)])
plot.power("../figs/rej_tree_BFS_h_power.pdf",
         power.BFS.expr3,
         mains = c("BFS (Case 1)",
             "BFS (Case 2)",
             "BFS (Case 3)"),
         methods = c("SS (p*=0.5)", "FS", "HE (p*=0.5)",
             "HE (p*=0.7)", "HE (p*=0.9)"),
         cols = c("black", "red", "blue", "orange",
             "magenta")
         )

fdp.DFS.expr3 <- list(result.DFS.tree[[1]][c(1, 6:9)],
                      result.DFS.isotree[[1]][c(1, 6:9)],
                      result.DFS.aniisotree[[1]][c(1, 6:9)])
plot.FDR("../figs/rej_tree_DFS_h_FDR.pdf",
         fdp.DFS.expr3,
         mains = c("DFS (Case 1)",
             "DFS (Case 2)",
             "DFS (Case 3)"),
         methods = c("SS (p*=0.5)", "FS", "HE (p*=0.5)",
             "HE (p*=0.7)", "HE (p*=0.9)"),
         cols = c("black", "red", "blue", "orange",
             "magenta")
         )

power.DFS.expr3 <- list(result.DFS.tree[[2]][c(1, 6:9)],
                      result.DFS.isotree[[2]][c(1, 6:9)],
                      result.DFS.aniisotree[[2]][c(1, 6:9)])
plot.power("../figs/rej_tree_DFS_h_power.pdf",
         power.DFS.expr3,
         mains = c("DFS (Case 1)",
             "DFS (Case 2)",
             "DFS (Case 3)"),
         methods = c("SS (p*=0.5)", "FS", "HE (p*=0.5)",
             "HE (p*=0.7)", "HE (p*=0.9)"),
         cols = c("black", "red", "blue", "orange",
             "magenta")
         )

#### Sensitivity Analysis
for (rho in c(5, -5)){
    load(paste0("../data/data_tree_", rho, ".RData"))
    result.BFS.tree <- data.tree[[1]]
    result.DFS.tree <- data.tree[[2]]
    result.BFS.isotree <- data.tree[[3]]
    result.DFS.isotree <- data.tree[[4]]
    result.BFS.aniisotree <- data.tree[[5]]
    result.DFS.aniisotree <- data.tree[[6]]


    ## Compare Methods
    fdp.BFS.expr1 <- list(result.BFS.tree[[1]][c(1, 10, 13, 14, 11)],
                          result.BFS.isotree[[1]][c(1, 10, 13, 14, 11)],
                          result.BFS.aniisotree[[1]][c(1, 10, 13, 14, 11)])
    plot.FDR(paste0("../figs/rej_tree_BFS_methods_FDR_", rho, ".pdf"),
             fdp.BFS.expr1,
             mains = c("BFS (Case 1)",
                 "BFS (Case 2)",
                 "BFS (Case 3)"),
             methods = c("STAR", "Yekutieli", "LG", "LG2", "BH"),
             cols = c("black", "red", "blue", "magenta", "orange"),
             )

    power.BFS.expr1 <- list(result.BFS.tree[[2]][c(1, 10, 13, 14, 11)],
                            result.BFS.isotree[[2]][c(1, 10, 13, 14, 11)],
                            result.BFS.aniisotree[[2]][c(1, 10, 13, 14, 11)])
    plot.power(paste0("../figs/rej_tree_BFS_methods_power_", rho, ".pdf"),
               power.BFS.expr1,
               mains = c("BFS (Case 1)",
                   "BFS (Case 2)",
                   "BFS (Case 3)"),
               methods = c("STAR", "Yekutieli", "LG", "LG2", "BH"),
               cols = c("black", "red", "blue", "magenta", "orange"),
               )

    fdp.DFS.expr1 <- list(result.DFS.tree[[1]][c(1, 10, 13, 14, 11)],
                          result.DFS.isotree[[1]][c(1, 10, 13, 14, 11)],
                          result.DFS.aniisotree[[1]][c(1, 10, 13, 14, 11)])
    plot.FDR(paste0("../figs/rej_tree_DFS_methods_FDR_", rho, ".pdf"),
             fdp.DFS.expr1,
             mains = c("DFS (Case 1)",
                 "DFS (Case 2)",
                 "DFS (Case 3)"),
             methods = c("STAR", "Yekutieli", "LG", "LG2", "BH"),
             cols = c("black", "red", "blue", "magenta", "orange"),
             )

    power.DFS.expr1 <- list(result.DFS.tree[[2]][c(1, 10, 13, 14, 11)],
                            result.DFS.isotree[[2]][c(1, 10, 13, 14, 11)],
                            result.DFS.aniisotree[[2]][c(1, 10, 13, 14, 11)])
    plot.power(paste0("../figs/rej_tree_DFS_methods_power_", rho, ".pdf"),
               power.DFS.expr1,
               mains = c("DFS (Case 1)",
                   "DFS (Case 2)",
                   "DFS (Case 3)"),
               methods = c("STAR", "Yekutieli", "LG", "LG2", "BH"),
               cols = c("black", "red", "blue", "magenta", "orange"),
               )
}
