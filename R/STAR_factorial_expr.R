##############################################################################
## STAR for Interaction Selection in Factorial Experiments
##################################################################
############

source("STAR_DAG.R")

effect.DAG <- function(col.names){
    split.list <- strsplit(col.names, ":")
    order <- sapply(split.list, length)
    main.ind <- which(order == 1)
    main.names <- unlist(split.list[main.ind])
    interact.ind <- setdiff(1:length(col.names), main.ind)
    edges <- lapply(interact.ind, function(i){
        effect <- split.list[[i]]
        ind1 <- main.ind[which(main.names == effect[1])]
        ind2 <- main.ind[which(main.names == effect[2])]
        if (ind1 != ind2){
            df <- data.frame(from = c(ind1, ind2), to = c(i, i))
        } else {
            df <- data.frame(from = ind1, to = i)
        }
        return(df)
    })
    edges <- Reduce(rbind, edges)
    DAG <- graph_from_data_frame(edges)
    return(DAG)
}

#### 2-levels analysis
data <- read.table("../data/two_levels.txt", stringsAsFactor = FALSE)
data <- data[, -1]
for (i in 1:6){
    data[, i] <- as.numeric(data[, i])
}
names(data) <- c("A", "B", "C", "D", "E", "F", "y")
data[["y"]] <- log(data[["y"]])
mod <- summary(lm(y ~ (A + B + C + D + E + F)^2., data = data))

nodes <- names(mod$coefficients[-1, 4])
pvals <- as.numeric(mod$coefficients[-1, 4])
DAG <- effect.DAG(nodes)

result <- STAR.DAG(pvals, DAG, criterion = "strong",
                   fun.list = create.fun(pstar = 0.5))
nodes[result$mask[, 26]]

refit <- lm(y ~ A + B + C + D + E + F + A*B + A*D + C*D, data = data)
