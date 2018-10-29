source('STAR_convex.R')
source('expr_func.R')
library("MASS")
library("mgcv")
library("RCurl")

## AdaPT
source("AdaPT.R")
source("AdaPT_gam.R")
source('summarize_methods.R')
source("expr_func.R")

rho <- as.numeric(Sys.getenv("rho"))
type <- as.numeric(Sys.getenv("type"))
repeat.times <- as.numeric(Sys.getenv("times"))
seed <- as.numeric(Sys.getenv("seed"))

convex.expr <- function(x, mu, H0, rho, type,
                        cov.formula,
                        alpha.list,
                        repeat.times,
                        ...){
    n <- length(mu)
    m <- length(alpha.list)
    summary.FDP <- list()
    summary.power <- list()
    for (j in 1:13){
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
            STARconvex <- try(STAR.convex(pvals, x,
                                          print.quiet = TRUE,
                                          alpha.list = alpha.list,
                                          fun.list = fun,
                                          type = "model-assist",
                                          update.mask.params =
                                              list(dir = "min"),
                                          ...))
            if (class(STARconvex) != "try-error"){
                STAR.convex.result <- summary.STAR(STARconvex, H0)
                summary.FDP[[k]][, i] <- STAR.convex.result[, 2]
                summary.power[[k]][, i] <- STAR.convex.result[, 3]
            }
        }

        fun <- STAR.fun.list.set[[1]]
        STARconvex <- try(STAR.convex(pvals, x,
                                      print.quiet = TRUE,
                                      alpha.list = alpha.list,
                                      fun.list = fun,
                                      type = "naive"))
        if (class(STARconvex) != "try-error"){
            STAR.convex.result <- summary.STAR(STARconvex, H0)
            summary.FDP[[10]][, i] <- STAR.convex.result[, 2]
            summary.power[[10]][, i] <- STAR.convex.result[, 3]
        }
        
        BH.result <- summary.BH(pvals, H0, alpha.list)
        summary.FDP[[11]][, i] <- BH.result[, 2]
        summary.power[[11]][, i] <- BH.result[, 3]

        storey.result <- summary.storey(pvals, H0,
                                        alpha.list = alpha.list)
        summary.FDP[[12]][, i] <- storey.result[, 2]
        summary.power[[12]][, i] <- storey.result[, 3]

        mux.init <- try(STAR.gam.em.censor(x, pvals, rep(TRUE, n),
                                           cov.formula = cov.formula))
        if (class(mux.init) == "try-error"){
            print(i)
            next
        }
        s0 <- ifelse(mux.init > 1.5, 0.45, 0.05)
        mod <- try(AdaPT(x, pvals, s0 = s0,
                         cov.formula = cov.formula,
                         q.list = alpha.list,
                         num.steps = 5,
                         delta.high = 0.3, delta.low = 0.3,
                         R.max = 20,
                         plot.quiet = TRUE,
                         disp.quiet = TRUE))
        if (class(mux.init) != "try-error"){
            AdaPT.result <- summary.AdaPT(mod, H0, pvals)
            summary.FDP[[13]][, i] <- AdaPT.result[, 2]
            summary.power[[13]][, i] <- AdaPT.result[, 3]
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

####### Generate x
set.seed(seed)
n <- 2500
x1 <- seq(-100, 100, length.out = 50)
x2 <- seq(-100, 100, length.out = 50)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")
cov.formula <- "s(x1, x2)"
alpha.list <- seq(0.01, 0.3, 0.01)
num.steps.update.score <- 10
num.steps.gam <- 5
score.params <- list(cov.formula = cov.formula,
                     num.steps = num.steps.gam)
output.filename <- paste0("../data/data_convex_rho_", floor(rho*10), "_type_", floor(type), "_seed_", seed, ".RData")

## Case 1: a circle in the center
H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
mu <- ifelse(H0, 2, 0)

result1.convex <- convex.expr(x, mu, H0, rho, type,
                              cov.formula = cov.formula,
                              alpha.list = alpha.list,
                              repeat.times = repeat.times,
                              num.steps.update.score =
                                  num.steps.update.score,
                              score.params = score.params)
data.convex <- list(result1.convex)
save(file = output.filename, data.convex)

## Case 2: a circle in the corner
H0 <- apply(x, 1, function(coord){sum((coord - 65)^2) < 900})
mu <- ifelse(H0, 2, 0)

result2.convex <- convex.expr(x, mu, H0, rho, type,
                              cov.formula = cov.formula,
                              alpha.list = alpha.list,
                              repeat.times = repeat.times,
                              num.steps.update.score =
                                  num.steps.update.score,
                              score.params = score.params)
data.convex <- list(result1.convex, result2.convex)
save(file = output.filename, data.convex)


## Case 3: a thin ellipsoid
shape.fun <- function(coord){
    transform.coord <- c(coord[1] + coord[2], coord[2] - coord[1])/sqrt(2)
    transform.coord[1]^2 / 100^2 + transform.coord[2]^2 / 15^2 < 1
}
H0 <- apply(x, 1, shape.fun)
mu <- ifelse(H0, 2, 0)

result3.convex <- convex.expr(x, mu, H0, rho, type,
                              cov.formula = cov.formula,
                              alpha.list = alpha.list,
                              repeat.times = repeat.times,
                              num.steps.update.score =
                                  num.steps.update.score,
                              score.params = score.params)

## Save data
data.convex <- list(result1.convex, result2.convex,
                    result3.convex)
save(file = output.filename, data.convex)
