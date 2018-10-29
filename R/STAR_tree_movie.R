source('STAR_tree.R')
library("isotone")
library("ape")
library("igraph")
library("MASS")

set.seed(2017)
n <- 100
ns <- 20
tree <- make_tree(n)
V(tree)$name <- 1:n
H0 <- c(rep(TRUE, ns), rep(FALSE, n - ns))
mu <- c(rep(2, ns), rep(0, n - ns))
null.dist <- runif
alt.dist <- function(n){1 - pnorm(rnorm(n, mean = mu))}
pvals <- ifelse(H0, alt.dist(n), null.dist(n))

add.plot <- list()
layout1 <- layout_as_tree(tree, circular = TRUE)
add.plot[[1]] <- expression(plot.rej.tree(tree, H0, layout1,
    main = "Truth (circular layout)"))
layout2 <- layout_as_tree(tree, circular = FALSE)
add.plot[[2]] <- expression(plot.rej.tree(tree, H0, layout2,
    main = "Truth (hierarchical layout)"))

illustrate1 <- STAR.tree(pvals, tree,
                         plot.params = list(
                             add.plot = add.plot,
                             cex.main = 4),
                         gif.root.title = "../movies/naivetree",
                         type = "naive",
                         height = 1440)

illustrate2 <- STAR.tree(pvals, tree,
                         plot.params = list(
                             add.plot = add.plot,
                             cex.main = 4),
                         gif.root.title = "../movies/isotree",
                         type = "model-assist",
                         height = 1440)
