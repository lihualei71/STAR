source('STAR_DAG.R')
library("igraph")
library("MASS")

set.seed(2017)
DAG <- DAG.gen(20, type = "reg_layers", 6, nlinks = 3)
layer <- rep(1:6, each = 20)
n <- length(V(DAG))
V(DAG)$name <- 1:n
strong.H0 <- rep(FALSE, n)
strong.H0.ind <- DAG.strong.null(DAG, 1:10)
strong.H0[strong.H0.ind] <- TRUE

mu <- rep(0, n)
mu[strong.H0.ind] <- 3
pvals <- pvals.gen(n, mu, 0, 0)
add.plot1 <- expression(DAG.plot(DAG, strong.H0, main = "Truth (Strong Heredity)"))

illustrate1 <- STAR.DAG(pvals, DAG,
                        gif.root.title = "../movies/DAG_strong",
                        type = "naive",
                        criterion = "strong",
                        plot.params = list(layer = layer,
                            add.plot = add.plot1,
                            cex.main = 4),
                        delay = 20)

weak.H0 <- rep(FALSE, n)
weak.H0.ind <- c(1, 10,
                 32, 24, 33, 34,
                 54, 55, 47, 56, 57, 58,
                 78, 70, 79, 80, 61,
                 93, 82, 83)
weak.H0[weak.H0.ind] <- TRUE

mu <- rep(0, n)
mu[weak.H0.ind] <- 3
pvals <- pvals.gen(n, mu, 0, 0)
add.plot2 <- expression(DAG.plot(DAG, weak.H0, main = "Truth (Weak Heredity)"))

illustrate2 <- STAR.DAG(pvals, DAG,
                        gif.root.title = "../movies/DAG_weak",
                        type = "naive",
                        criterion = "weak",
                        plot.params = list(layer = layer,
                            add.plot = add.plot2,
                            cex.main = 4),
                        delay = 20)
