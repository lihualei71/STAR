## Read Patch Data from Clusters
for (rho in c(0, 5, -5)){
    type <- ifelse(rho < 0, 1, 0)
    data.DAG.sum <- list()
    for (k in 1:5){
        data.DAG.sum[[k]] <- list(FDP = list(), power = list())
    }
    for (k in 1:5){
        for (j in 1:14){
            data.DAG.sum[[k]]$FDP[[j]] <- rep(0, 30)
            data.DAG.sum[[k]]$power[[j]] <- rep(0, 30)    
        }
    }

    
    for (seed in 1:20){
        filename <- paste0("../data/data_DAG_rho_", rho, "_type_", type, "_seed_", seed, ".RData")
        load(filename)
        for (k in 1:5){
            for (j in 1:14){
                data.DAG.sum[[k]]$FDP[[j]] <-
                    data.DAG.sum[[k]]$FDP[[j]] +
                        data.DAG[[k]]$FDP[[j]]
                data.DAG.sum[[k]]$power[[j]] <-
                    data.DAG.sum[[k]]$power[[j]] +
                        data.DAG[[k]]$power[[j]]
            }
        }
    }
    for (k in 1:5){
        for (j in 1:14){
            data.DAG.sum[[k]]$FDP[[j]] <-
                data.DAG.sum[[k]]$FDP[[j]] / 20
            data.DAG.sum[[k]]$power[[j]] <-
                data.DAG.sum[[k]]$power[[j]] / 20
        }
    }
    data.DAG <- data.DAG.sum
    save(file = paste0("../data/data_DAG_", rho, ".RData"), data.DAG)
}
