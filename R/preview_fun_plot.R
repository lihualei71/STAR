library("latex2exp")

source("generic_STAR.R")
x <- seq(0, 1, 0.01)

pdf("../figs/preview_fun.pdf", width = 7, height = 2.5)
par(mfrow = c(1, 3), cex = 0.75)
fun1 <- create.fun("SeqStep", 0.5)
plot(x, pmin(x, fun1$g(x)), type = "l", xlab = "p-value", ylab = "g(p)", main = TeX("SeqStep ($p_{*}=0.5$)"))
fun2 <- create.fun("ForwardStop")
plot(x, pmin(x, fun2$g(x)), type = "l", xlab = "p-value", ylab = "g(p)", main = TeX("ForwardStop"))
fun3 <- create.fun("HingeExp", 0.9)
plot(x, pmin(x, fun3$g(x)), type = "l", xlab = "p-value", ylab = "g(p)", main = TeX("HingeExp ($p_{*}=0.9$)"))
dev.off()    
