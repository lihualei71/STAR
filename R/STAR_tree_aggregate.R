## Read Patch Data obtained from Clusters
for (rho in c(0, 5, -5)){
    type <- ifelse(rho < 0, 1, 0)
    data.tree.sum <- list()
    for (k in 1:6){
        data.tree.sum[[k]] <- list(FDP = list(), power = list())
    }
    for (k in 1:6){
        for (j in 1:14){
            data.tree.sum[[k]]$FDP[[j]] <- rep(0, 30)
            data.tree.sum[[k]]$power[[j]] <- rep(0, 30)    
        }
    }

    for (seed in 1:20){
        filename <- paste0("../data/data_tree_rho_", rho, "_type_", type, "_seed_", seed, ".RData")
        load(filename)
        for (k in 1:6){
            for (j in 1:14){
                data.tree.sum[[k]]$FDP[[j]] <-
                    data.tree.sum[[k]]$FDP[[j]] +
                        data.tree[[k]]$FDP[[j]]
                data.tree.sum[[k]]$power[[j]] <-
                    data.tree.sum[[k]]$power[[j]] +
                        data.tree[[k]]$power[[j]]
            }
        }
    }
    for (k in 1:6){
        for (j in 1:14){
            data.tree.sum[[k]]$FDP[[j]] <-
                data.tree.sum[[k]]$FDP[[j]] / 20
            data.tree.sum[[k]]$power[[j]] <-
                data.tree.sum[[k]]$power[[j]] / 20
        }
    }
    data.tree <- data.tree.sum
    save(file = paste0("../data/data_tree_", rho, ".RData"), data.tree)
}
