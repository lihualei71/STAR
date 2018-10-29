source('STAR_DAG.R')
source('Lynch.R')
source('Ramdas.R')
library("igraph")
library("MASS")
source("summarize_methods.R")
source("expr_func.R")

rho <- as.numeric(Sys.getenv("rho"))
type <- as.numeric(Sys.getenv("type"))
repeat.times <- as.numeric(Sys.getenv("times"))
seed <- as.numeric(Sys.getenv("seed"))

DAG.expr <- function(DAG, mu, H0, rho,
                     alpha.list,
                     repeat.times){
    n <- length(mu)
    m <- length(alpha.list)
    summary.FDP <- list()
    summary.power <- list()
    for (j in 1:14){
        summary.FDP[[j]] <- matrix(rep(0, repeat.times * m),
                                   ncol = repeat.times)
        summary.power[[j]] <- matrix(rep(0, repeat.times * m),
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
            STARDAG <- STAR.DAG(pvals, DAG, print.quiet = TRUE,
                              alpha.list = alpha.list,
                              fun.list = fun)
            STARDAG.result <- summary.STAR(STARDAG, H0)
            summary.FDP[[k]][, i] <- STARDAG.result[, 2]
            summary.power[[k]][, i] <- STARDAG.result[, 3]
            print(paste0("************** k = ", k, "****************"))
        }

        BH.result <- summary.BH(pvals, H0, alpha.list)
        summary.FDP[[10]][, i] <- BH.result[, 2]
        summary.power[[10]][, i] <- BH.result[, 3]

        storey.result <- summary.storey(pvals, H0,
                                        alpha.list = alpha.list)
        summary.FDP[[11]][, i] <- storey.result[, 2]
        summary.power[[11]][, i] <- storey.result[, 3]

        SCR.DAG.result <- SCR.DAG(DAG, pvals,
                                  alpha.list = alpha.list)
        Lynch.result1 <- summary.Lynch(SCR.DAG.result, H0)
        summary.FDP[[12]][, i] <- Lynch.result1[, 2]
        summary.power[[12]][, i] <- Lynch.result1[, 3]
        
        BH.DAG.result <- BH.DAG(DAG, pvals,
                                alpha.list = alpha.list)
        Lynch.result2 <- summary.Lynch(BH.DAG.result, H0)
        summary.FDP[[13]][, i] <- Lynch.result2[, 2]
        summary.power[[13]][, i] <- Lynch.result2[, 3]

        CR.DAG <- chen.ramdas(pvals, DAG,
                              alpha.list = alpha.list)
        CR.result <- summary.cr(CR.DAG, H0)
        summary.FDP[[14]][, i] <- CR.result[, 2]
        summary.power[[14]][, i] <- CR.result[, 3]
        
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

### Hyper parameters
set.seed(seed)
alpha.list <- seq(0.01, 0.3, 0.01)
output.filename <- paste0("../data/data_DAG_rho_", floor(rho*10), "_type_", floor(type), "_seed_", seed, ".RData")

### Expr 1: Regular shallow DAG
DAG <- DAG.gen(250, type = "reg_layers", 4, nlinks = 3)
n <- length(V(DAG))
V(DAG)$name <- 1:n
H0 <- rep(FALSE, n)
H0.ind <- DAG.strong.null(DAG, 1:16, 50)
H0[H0.ind] <- TRUE
mu <- rep(0, n)
mu[H0.ind] <- 3

DAG.expr1 <- DAG.expr(DAG, mu, H0, rho, alpha.list, repeat.times)
data.DAG <- list(DAG.expr1)
save(file = output.filename, data.DAG)

### Expr 2: Regular deep DAG
DAG <- DAG.gen(100, type = "reg_layers", 10, nlinks = 3)
n <- length(V(DAG))
V(DAG)$name <- 1:n
H0 <- rep(FALSE, n)
H0.ind <- DAG.strong.null(DAG, 1:14, 50)
H0[H0.ind] <- TRUE
mu <- rep(0, n)
mu[H0.ind] <- 3

DAG.expr2 <- DAG.expr(DAG, mu, H0, rho, alpha.list, repeat.times)
data.DAG <- list(DAG.expr1, DAG.expr2)
save(file = output.filename, data.DAG)

### Expr 3: Irregular DAG (triangle)
DAG <- DAG.gen(c(50, 100, 200, 300, 350),
               type = "nreg", nlinks = 3)
n <- length(V(DAG))
V(DAG)$name <- 1:n
H0 <- rep(FALSE, n)
H0.ind <- DAG.strong.null(DAG, 1:9, 50)
H0[H0.ind] <- TRUE
mu <- rep(0, n)
mu[H0.ind] <- 3

DAG.expr3 <- DAG.expr(DAG, mu, H0, rho, alpha.list, repeat.times)
data.DAG <- list(DAG.expr1, DAG.expr2, DAG.expr3)
save(file = output.filename, data.DAG)

### Expr 4: Irregular DAG (inverted triangle)
DAG <- DAG.gen(c(350, 300, 200, 100, 50),
               type = "nreg", nlinks = 3)
n <- length(V(DAG))
V(DAG)$name <- 1:n
H0 <- rep(FALSE, n)
H0.ind <- DAG.strong.null(DAG, 100:120, 50)
H0[H0.ind] <- TRUE
mu <- rep(0, n)
mu[H0.ind] <- 3

DAG.expr4 <- DAG.expr(DAG, mu, H0, rho, alpha.list, repeat.times)
data.DAG <- list(DAG.expr1, DAG.expr2, DAG.expr3, DAG.expr4)
save(file = output.filename, data.DAG)

### Expr 5: Irregular DAG (spindle)
DAG <- DAG.gen(c(100, 200, 400, 200, 100),
               type = "nreg", nlinks = 3)
n <- length(V(DAG))
V(DAG)$name <- 1:n
H0 <- rep(FALSE, n)
H0.ind <- DAG.strong.null(DAG, 1:10, 50)
H0[H0.ind] <- TRUE
mu <- rep(0, n)
mu[H0.ind] <- 3

DAG.expr5 <- DAG.expr(DAG, mu, H0, rho, alpha.list, repeat.times)
data.DAG <- list(DAG.expr1, DAG.expr2, DAG.expr3, DAG.expr4,
                 DAG.expr5)
save(file = output.filename, data.DAG)
