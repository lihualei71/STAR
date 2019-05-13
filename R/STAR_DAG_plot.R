source("FDR_power_plot.R")
source("expr_func.R")
source("STAR_DAG.R")
library("igraph")

## Plot of truth
pdf("../figs/DAG_expr_truth.pdf", width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(4, 2, 2, 2))
mains <- c("Shallow Regular DAG", "Deep Regular DAG",
           "Triangular DAG")
layers <- list()
layers[[1]] <- rep(1:4, each = 100)
layers[[2]] <- rep(1:10, each = 100)
layers[[3]] <- c(rep(1, 50/2), rep(2, 100/2), rep(3, 200/2),
                 rep(4, 300/2), rep(5, 350/2))
DAGs <- list()
DAGs[[1]] <- DAG.gen(100, type = "reg_layers")
DAGs[[2]] <- DAG.gen(100, type = "reg_layers", nlayers = 10)
DAGs[[3]] <- DAG.gen(c(50, 100, 200, 300, 350) / 2,
                     type = "nreg")
nums <- c(15, 20, 9)
for (i in 1:3){
    DAG <- DAGs[[i]]
    H0 <- rep(FALSE, length(V(DAG)))
    H0.ind <- DAG.strong.null(DAG, 1:nums[i])
    H0[H0.ind] <- TRUE
    plot.rej.DAG(DAG, H0, layers[[i]],
                 main = mains[i],
                 vertex.size = 5, vertex.label = NA)
}
dev.off()

load("../data/data_DAG_0.RData")

result.shallow.reg <- data.DAG[[1]]
result.deep.reg <- data.DAG[[2]]
result.tri.reg <- data.DAG[[3]]
#### Ignore the results for other two types of DAGs
## result.invtri.reg <- data.DAG[[4]]
## result.spindle.reg <- data.DAG[[5]]


## Compare Methods
fdp.expr1 <- list(result.shallow.reg[[1]][c(1, 10, 12, 13, 14)],
                  result.deep.reg[[1]][c(1, 10, 12, 13, 14)],
                  result.tri.reg[[1]][c(1, 10, 12, 13, 14)])
plot.FDR("../figs/rej_DAG_methods_FDR.pdf",
         fdp.expr1, 
         mains = c("Shallow Regular DAG",
                  "Deep Regular DAG",
                  "Triangular DAG"),
         methods = c("STAR", "BH", "SCR-DAG", "BH-DAG", "DAGGER"),
         cols = c("black", "red", "blue", "orange", "magenta"))

power.expr1 <- list(result.shallow.reg[[2]][c(1, 10, 12, 13, 14)],
                    result.deep.reg[[2]][c(1, 10, 12, 13, 14)],
                    result.tri.reg[[2]][c(1, 10, 12, 13, 14)])
plot.power("../figs/rej_DAG_methods_power.pdf",
         power.expr1, 
         mains = c("Shallow Regular DAG",
                  "Deep Regular DAG",
                  "Triangular DAG"),
         methods = c("STAR", "BH", "SCR-DAG", "BH-DAG", "DAGGER"),
         cols = c("black", "red", "blue", "orange", "magenta"))

## Compare parameters of ISS
fdp.expr2 <- list(result.shallow.reg[[1]][c(1:5)],
                  result.deep.reg[[1]][c(1:5)],
                  result.tri.reg[[1]][c(1:5)])
plot.FDR("../figs/rej_DAG_STARSS_FDR.pdf",
         fdp.expr2, 
         mains = c("Shallow Regular DAG",
             "Deep Regular DAG",
             "Triangular DAG"),
         methods = c("p*=0.5", "p*=0.6", "p*=0.7", "p*=0.8", "p*=0.9"),
         cols = c("black", "red", "blue", "orange", "magenta")
         )

power.expr2 <- list(result.shallow.reg[[2]][c(1:5)],
                  result.deep.reg[[2]][c(1:5)],
                  result.tri.reg[[2]][c(1:5)])
plot.power("../figs/rej_DAG_STARSS_power.pdf",
         power.expr2, 
         mains = c("Shallow Regular DAG",
             "Deep Regular DAG",
             "Triangular DAG"),
         methods = c("p*=0.5", "p*=0.6", "p*=0.7", "p*=0.8", "p*=0.9"),
         cols = c("black", "red", "blue", "orange", "magenta")
         )

## Compare h-function
fdp.expr3 <- list(result.shallow.reg[[1]][c(1, 6:9)],
                  result.deep.reg[[1]][c(1, 6:9)],
                  result.tri.reg[[1]][c(1, 6:9)])
plot.FDR("../figs/rej_DAG_h_FDR.pdf",
              fdp.expr3, 
              mains = c("Shallow Regular DAG",
                  "Deep Regular DAG",
                  "Triangular DAG"),
              methods = c("SS (p*=0.5)", "FS", "HE (p*=0.5)",
                  "HE (p*=0.7)", "HE (p*=0.9)"),
              cols = c("black", "red", "blue", "orange",
                  "magenta")
              )

power.expr3 <- list(result.shallow.reg[[2]][c(1, 6:9)],
                  result.deep.reg[[2]][c(1, 6:9)],
                  result.tri.reg[[2]][c(1, 6:9)])
plot.power("../figs/rej_DAG_h_power.pdf",
              power.expr3, 
              mains = c("Shallow Regular DAG",
                  "Deep Regular DAG",
                  "Triangular DAG"),
              methods = c("SS (p*=0.5)", "FS", "HE (p*=0.5)",
                  "HE (p*=0.7)", "HE (p*=0.9)"),
              cols = c("black", "red", "blue", "orange",
                  "magenta")
              )


#### Sensitivity Analysis
for (rho in c(5, -5)){
    load(paste0("../data/data_DAG_", rho, ".RData"))

    result.shallow.reg <- data.DAG[[1]]
    result.deep.reg <- data.DAG[[2]]
    result.tri.reg <- data.DAG[[3]]
    #### Ignore the results for other two types of DAGs
    ## result.invtri.reg <- data.DAG[[4]]
    ## result.spindle.reg <- data.DAG[[5]]


    ## Compare Methods
    fdp.expr1 <- list(result.shallow.reg[[1]][c(1, 10, 12, 13, 14)],
                      result.deep.reg[[1]][c(1, 10, 12, 13, 14)],
                      result.tri.reg[[1]][c(1, 10, 12, 13, 14)])
    plot.FDR(paste0("../figs/rej_DAG_methods_FDR_", rho, ".pdf"),
             fdp.expr1, 
             mains = c("Shallow Regular DAG",
                 "Deep Regular DAG",
                 "Triangular DAG"),
             methods = c("STAR", "BH", "SCR-DAG", "BH-DAG", "DAGGER"),
             cols = c("black", "red", "blue", "orange", "magenta"))

    power.expr1 <- list(result.shallow.reg[[2]][c(1, 10, 12, 13, 14)],
                        result.deep.reg[[2]][c(1, 10, 12, 13, 14)],
                        result.tri.reg[[2]][c(1, 10, 12, 13, 14)])
    plot.power(paste0("../figs/rej_DAG_methods_power_", rho, ".pdf"),
               power.expr1, 
               mains = c("Shallow Regular DAG",
                   "Deep Regular DAG",
                   "Triangular DAG"),
               methods = c("STAR", "BH", "SCR-DAG", "BH-DAG", "DAGGER"),
               cols = c("black", "red", "blue", "orange", "magenta"))
}
