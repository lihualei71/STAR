## Read Patch Data obtained from Clusters
for (rho in c(0, 5, -5)){
    type <- ifelse(rho < 0, 1, 0)
    data.convex.sum <- list()
    for (k in 1:3){
        data.convex.sum[[k]] <- list(FDP = list(), power = list())
    }
    for (k in 1:3){
        for (j in 1:13){
            data.convex.sum[[k]]$FDP[[j]] <- rep(0, 30)
            data.convex.sum[[k]]$power[[j]] <- rep(0, 30)    
        }
    }

    for (seed in 1:20){
        filename <- paste0("../data/data_convex_rho_", rho, "_type_", type, "_seed_", seed, ".RData")
        load(filename)
        for (k in 1:3){
            for (j in 1:13){
                data.convex.sum[[k]]$FDP[[j]] <-
                    data.convex.sum[[k]]$FDP[[j]] +
                        data.convex[[k]]$FDP[[j]]
                data.convex.sum[[k]]$power[[j]] <-
                    data.convex.sum[[k]]$power[[j]] +
                        data.convex[[k]]$power[[j]]
            }
        }
    }
    for (k in 1:3){
        for (j in 1:13){
            data.convex.sum[[k]]$FDP[[j]] <-
                data.convex.sum[[k]]$FDP[[j]] / 20
            data.convex.sum[[k]]$power[[j]] <-
                data.convex.sum[[k]]$power[[j]] / 20
        }
    }
    data.convex <- data.convex.sum
    save(file = paste0("../data/data_convex_", rho, ".RData"), data.convex)
}
