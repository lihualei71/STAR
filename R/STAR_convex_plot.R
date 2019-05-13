source("FDR_power_plot.R")
source("expr_func.R")
source("STAR_convex.R")
source("summarize_methods.R")
source("AdaPT.R")
source("AdaPT_gam.R")

## Plot of truth
pdf("../figs/convex_expr_truth.pdf", width = 5, height = 1.75)
par(mfrow = c(1, 3), mar = c(2, 2, 2, 2))
n <- 2500
x1 <- seq(-100, 100, length.out = 50)
x2 <- seq(-100, 100, length.out = 50)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")
mains <- c("Circle in the middle",
           "Circle in the corner",
           "Thin ellipse")
for (i in 1:3){
    if (i == 1){
        H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
    } else if (i == 2){
        H0 <- apply(x, 1, function(coord){sum((coord - 65)^2) < 900})
    } else {
        shape.fun <- function(coord){
            transform.coord <- c(coord[1] + coord[2], coord[2] - coord[1])/sqrt(2)
            transform.coord[1]^2 / 100^2 + transform.coord[2]^2 / 15^2 < 1
        }
        H0 <- apply(x, 1, shape.fun)
    }
    plot.rej.convex(x, H0, cex = 0.5,
                    col.bg = "#A9A9A9", col.fg = "#000000",
                    main = mains[i], xaxt = "n", yaxt = "n")
    axis(side = 1, at = c(-100, 0, 100))
    axis(side = 2, at = c(-100, 0, 100))
}
dev.off()


#### Plot \hat{\mu} of STAR at t = 0 and t = infty
#### Plot \hat{\mu} with different methods 
n <- 2500
x1 <- seq(-100, 100, length.out = 50)
x2 <- seq(-100, 100, length.out = 50)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")
cov.formula <- "s(x1, x2)"
alpha.list <- seq(0.01, 1, 0.01)
num.steps.update.score <- 10
num.steps.gam <- 5
score.params <- list(cov.formula = cov.formula,
                     num.steps = num.steps.gam)

## Case 1: a circle in the center
H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
mu <- ifelse(H0, 2, 0)
pvals <- pvals.gen(n, mu, 0, 0)

STAR.obj1 <- STAR.convex(pvals, x, 
                         alpha.list = alpha.list,
                         type = "model-assist",
                         update.mask.params = list(dir = "min"),
                         num.steps.update.score = num.steps.update.score,
                         score.params = score.params)
BH.obj1 <- BH(pvals, 0.1)
AdaPT.obj1 <- AdaPT(x, pvals, cov.formula = cov.formula,
                    q.list = 0.1, plot.quiet = TRUE)
AdaPT.obj1 <- as.numeric(pvals <= AdaPT.obj1$s)

oracle.score1 <- STAR.gam.em.censor(x, pvals, rep(FALSE, n),
                                    cov.formula, num.steps = 30)
corr1 <- apply(STAR.obj1$score, 2, function(x){cor(x, oracle.score1)})

## Case 2: a circle in the corner
H0 <- apply(x, 1, function(coord){sum((coord - 65)^2) < 900})
mu <- ifelse(H0, 2, 0)
pvals <- pvals.gen(n, mu, 0, 0)

STAR.obj2 <- STAR.convex(pvals, x, 
                         alpha.list = alpha.list,
                         type = "model-assist",
                         update.mask.params = list(dir = "min"),
                         num.steps.update.score = num.steps.update.score,
                         score.params = score.params)
BH.obj2 <- BH(pvals, 0.1)
AdaPT.obj2 <- AdaPT(x, pvals, cov.formula = cov.formula,
                    q.list = 0.1, plot.quiet = TRUE)
AdaPT.obj2 <- as.numeric(pvals <= AdaPT.obj2$s)

oracle.score2 <- STAR.gam.em.censor(x, pvals, rep(FALSE, n),
                                    cov.formula, num.steps = 30)
corr2 <- apply(STAR.obj2$score, 2, function(x){cor(x, oracle.score2)})

## Case 3: a thin ellipsoid
shape.fun <- function(coord){
    transform.coord <- c(coord[1] + coord[2], coord[2] - coord[1])/sqrt(2)
    transform.coord[1]^2 / 100^2 + transform.coord[2]^2 / 15^2 < 1
}
H0 <- apply(x, 1, shape.fun)
mu <- ifelse(H0, 2, 0)
pvals <- pvals.gen(n, mu, 0, 0)

STAR.obj3 <- STAR.convex(pvals, x, 
                         alpha.list = alpha.list,
                         type = "model-assist",
                         update.mask.params = list(dir = "min"),
                         num.steps.update.score = num.steps.update.score,
                         score.params = score.params)
BH.obj3 <- BH(pvals, 0.1)

AdaPT.obj3 <- AdaPT(x, pvals, cov.formula = cov.formula,
                    q.list = 0.1, plot.quiet = TRUE)
AdaPT.obj3 <- as.numeric(pvals <= AdaPT.obj3$s)

oracle.score3 <- STAR.gam.em.censor(x, pvals, rep(FALSE, n),
                                    cov.formula, num.steps = 30)
corr3 <- apply(STAR.obj3$score, 2, function(x){cor(x, oracle.score3)})

pdf("../figs/convex_expr_score.pdf", width = 5, height = 3)
scores <- list(STAR.obj1$score[, 100], STAR.obj2$score[, 100],
               STAR.obj3$score[, 100],
               oracle.score1, oracle.score2, oracle.score3)
mains <- c("Initial Score (Case 1)", "Initial Score (Case 2)",
           "Initial Score (Case 3)", "Oracle Score (Case 1)",
           "Oracle Score (Case 2)", "Oracle Score (Case 3)")
par(mfrow = c(2, 3), mar = c(1, 2, 2, 2))
for (i in 1:6){
    plot.strength.convex(x, scores[[i]], cex = 1, cex.main = 1, main = mains[i], xaxt = "n", yaxt = "n", col.bg = "#A9A9A9", col.fg = "#000000")
    axis(side = 1, at = c(-100, 0, 100))
    axis(side = 2, at = c(-100, 0, 100))
}
dev.off()

pdf("../figs/convex_expr_rej.pdf", width = 5, height = 4.5)
rejs <- list(STAR.obj1$mask[, 10], STAR.obj2$mask[, 10],
             STAR.obj3$mask[, 10],
             BH.obj1, BH.obj2, BH.obj3,
             AdaPT.obj1, AdaPT.obj2, AdaPT.obj3
             )
mains <- c("STAR (Case 1)", "STAR (Case 2)",
           "STAR (Case 3)", "BH (Case 1)",
           "BH (Case 2)", "BH (Case 3)",
           "AdaPT (Case 1)", "AdaPT (Case 2)",
           "AdaPT (Case 3)")
par(mfrow = c(3, 3), mar = c(1, 2, 2, 2))
for (i in 1:9){
    plot.rej.convex(x, rejs[[i]], cex = 0.3, cex.main = 1,
                    main = mains[i], xaxt = "n", yaxt = "n")
    axis(side = 1, at = c(-100, 0, 100))
    axis(side = 2, at = c(-100, 0, 100))
}
dev.off()

load("../data/data_convex_0.RData")

result.circ <- data.convex[[1]]
result.corner <- data.convex[[2]]
result.ellipse <- data.convex[[3]]

## Compare Methods
## fdp.expr1 <- list(result.circ[[1]][c(1, 10, 11, 13)],
##                   result.corner[[1]][c(1, 10, 11, 13)],
##                   result.ellipse[[1]][c(1, 10, 11, 13)])
## plot.FDR("../figs/rej_convex_methods_FDR.pdf",
##          fdp.expr1, 
##          mains = c("Circle in the middle",
##              "Circle in the corner",
##              "Thin ellipse"),
##          methods = c("STAR-GAM", "STAR", "BH", "AdaPT"),
##          cols = c("black", "red", "blue", "orange")
##          )


## power.expr1 <- list(result.circ[[2]][c(1, 10, 11, 13)],
##                     result.corner[[2]][c(1, 10, 11, 13)],
##                     result.ellipse[[2]][c(1, 10, 11, 13)])
## plot.power("../figs/rej_convex_methods_power.pdf",
##            power.expr1, 
##            mains = c("Circle in the middle",
##                "Circle in the corner",
##                "Thin ellipse"),
##            methods = c("STAR-GAM", "STAR", "BH", "AdaPT"),
##            cols = c("black", "red", "blue", "orange")
##            )

## For Biometrika submission
fdp.expr1 <- list(result.circ[[1]][c(1, 11, 13)],
                  result.corner[[1]][c(1, 11, 13)],
                  result.ellipse[[1]][c(1, 11, 13)])
plot.FDR("../figs/rej_convex_methods_FDR.pdf",
         fdp.expr1, 
         mains = c("Circle in the middle",
             "Circle in the corner",
             "Thin ellipse"),
         methods = c("ours", "method 1", "method 2"),
         cols = c("black", "red", "blue")
         )


power.expr1 <- list(result.circ[[2]][c(1, 11, 13)],
                    result.corner[[2]][c(1, 11, 13)],
                    result.ellipse[[2]][c(1, 11, 13)])
plot.power("../figs/rej_convex_methods_power.pdf",
           power.expr1, 
           mains = c("Circle in the middle",
               "Circle in the corner",
               "Thin ellipse"),
           methods = c("ours", "method 1", "method 2"),
           cols = c("black", "red", "blue")
           )


## Compare parameters of ISS
fdp.expr2 <- list(result.circ[[1]][c(1:5)],
                  result.corner[[1]][c(1:5)],
                  result.ellipse[[1]][c(1:5)])
plot.FDR("../figs/rej_convex_STARSS_FDR.pdf",
         fdp.expr2,
         mains = c("Circle in the middle",
             "Circle in the corner",
             "Thin ellipse"),
         methods = c("p*=0.5", "p*=0.6", "p*=0.7",
               "p*=0.8", "p*=0.9"),
         cols = c("black", "red", "blue", "orange",
             "magenta")
         )


power.expr2 <- list(result.circ[[2]][c(1:5)],
                    result.corner[[2]][c(1:5)],
                    result.ellipse[[2]][c(1:5)])
plot.power("../figs/rej_convex_STARSS_power.pdf",
           power.expr2,
           mains = c("Circle in the middle",
               "Circle in the corner",
               "Thin ellipse"),
           methods = c("p*=0.5", "p*=0.6", "p*=0.7",
               "p*=0.8", "p*=0.9"),
           cols = c("black", "red", "blue", "orange",
               "magenta")
           )


## Compare h-function
fdp.expr3 <- list(result.circ[[1]][c(1, 6:9)],
                  result.corner[[1]][c(1, 6:9)],
                  result.ellipse[[1]][c(1, 6:9)])
plot.FDR("../figs/rej_convex_h_FDR.pdf",
         fdp.expr3, 
         mains = c("Circle in the middle",
             "Circle in the corner",
             "Thin ellipse"),
         methods = c("SS (p*=0.5)", "FS",
             "HE (p*=0.5)",
             "HE (p*=0.7)", "HE (p*=0.9)"),
         cols = c("black", "red", "blue", "orange",
             "magenta"),
         )

power.expr3 <- list(result.circ[[2]][c(1, 6:9)],
                    result.corner[[2]][c(1, 6:9)],
                    result.ellipse[[2]][c(1, 6:9)])

plot.power("../figs/rej_convex_h_power.pdf",
           power.expr3, 
           mains = c("Circle in the middle",
               "Circle in the corner",
               "Thin ellipse"),
           methods = c("SS (p*=0.5)", "FS",
               "HE (p*=0.5)",
               "HE (p*=0.7)", "HE (p*=0.9)"),
           cols = c("black", "red", "blue", "orange",
               "magenta"),
           )

#### Sensitivity Analysis
for (rho in c(5, -5)){
    load(paste0("../data/data_convex_", rho, ".RData"))

    result.circ <- data.convex[[1]]
    result.corner <- data.convex[[2]]
    result.ellipse <- data.convex[[3]]

    ## ## Compare Methods
    ## fdp.expr1 <- list(result.circ[[1]][c(1, 10, 11, 13)],
    ##                   result.corner[[1]][c(1, 10, 11, 13)],
    ##                   result.ellipse[[1]][c(1, 10, 11, 13)])
    ## plot.FDR(paste0("../figs/rej_convex_methods_FDR_",
    ##                 rho, ".pdf"),
    ##          fdp.expr1, 
    ##          mains = c("Circle in the middle",
    ##              "Circle in the corner",
    ##              "Thin ellipse"),
    ##          methods = c("STAR-GAM", "STAR", "BH", "AdaPT"),
    ##          cols = c("black", "red", "blue", "orange")
    ##          )


    ## power.expr1 <- list(result.circ[[2]][c(1, 10, 11, 13)],
    ##                     result.corner[[2]][c(1, 10, 11, 13)],
    ##                     result.ellipse[[2]][c(1, 10, 11, 13)])
    ## plot.power(paste0("../figs/rej_convex_methods_power_",
    ##                   rho, ".pdf"),
    ##            power.expr1, 
    ##            mains = c("Circle in the middle",
    ##                "Circle in the corner",
    ##                "Thin ellipse"),
    ##            methods = c("STAR-GAM", "STAR", "BH", "AdaPT"),
    ##            cols = c("black", "red", "blue", "orange"),
    ##            ylim = c(0, 1.1)
    ##            )

    ## For Biometrika submission    
    fdp.expr1 <- list(result.circ[[1]][c(1, 11, 13)],
                      result.corner[[1]][c(1, 11, 13)],
                      result.ellipse[[1]][c(1, 11, 13)])
    plot.FDR(paste0("../figs/rej_convex_methods_FDR_",
                    rho, ".pdf"),
             fdp.expr1, 
             mains = c("Circle in the middle",
                 "Circle in the corner",
                 "Thin ellipse"),
             methods = c("ours", "method 1", "method 2"),
             cols = c("black", "red", "blue")
             )


    power.expr1 <- list(result.circ[[2]][c(1, 11, 13)],
                        result.corner[[2]][c(1, 11, 13)],
                        result.ellipse[[2]][c(1, 11, 13)])
    plot.power(paste0("../figs/rej_convex_methods_power_",
                      rho, ".pdf"),
               power.expr1, 
               mains = c("Circle in the middle",
                   "Circle in the corner",
                   "Thin ellipse"),
               methods = c("ours", "method 1", "method 2"),
               cols = c("black", "red", "blue"),
               ylim = c(0, 1.1)
               )
}
