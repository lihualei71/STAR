source('STAR_convex.R')
source('expr_func.R')
library("MASS")
library("mgcv")
library("RCurl")

set.seed(2017)
n <- 400
x1 <- rep(1:20, 20) - 10.5
x2 <- rep(1:20, each = 20) - 10.5
x <- data.frame(x1 = x1, x2 = x2)
H0 <- apply(x, 1, function(item){sum(item^2) < 36})
mu.alt <- 3 
mu <- mu.alt * H0

pvals <- pvals.gen(n, mu, 0, 0)
cov.formula <- "s(x1, x2)"
add.plot <- expression(plot.rej.convex(x, H0, main = "Truth",
    cex = 6, pch = 15))

illustrate1 <- STAR.convex(pvals, x,
                           plot.params =
                               list(main.strength = "Canonical Score", pch = 15, cex = 6, add.plot = add.plot, cex.main = 4.5),
                           gif.root.title = "../movies/naiveconvex",
                           type = "naive", delay = 20,
                           width = 1800, height = 600)

illustrate2 <- STAR.convex(pvals, x,
                           plot.params =
                               list(main.strength = "Model-Assisted Score (GAM)", pch = 15, cex = 6, add.plot = add.plot, col.fg.strength = "#000080", col.bg.strength = "#ADD8E6", cex.main = 4.5),
                           score.params =
                               list(cov.formula = cov.formula),
                           update.mask.params = list(dir = "min"),
                           gif.root.title = "../movies/gamconvex",
                           type = "model-assist", delay = 20,
                           width = 1800, height = 600)
