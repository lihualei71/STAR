plot.results <- function(filename, vals,
                         mains, methods, cols,
                         ylab, ylim,
                         alpha.list = seq(0.01, 0.3, 0.01),
                         width = 7, height = 2,
                         ltys = 1:length(methods),
                         pchs = 1:length(methods)){
    pdf(filename, width = width, height = height)
    num.plots.per.row <- length(mains)
    num.methods <- length(methods)
    par(mfrow = c(1, num.plots.per.row),
        mar = c(4, 4, 2, 2), oma=c(0, 0, 0, 5))

    for (i in 1:(num.plots.per.row)){
        plot(0:1, 0:1, type = 'n',
             xlim = range(alpha.list),
             ylim = ylim,
             xlab = expression(paste('Target FDR level ',alpha)),
             ylab = ylab,
             axes = FALSE, font.main = 1, main = mains[i])
        axis(side = 1,at = c(0, 0.1, 0.2, 0.3))
        axis(side = 2)
        alpha.at <- c(10, 20, 30)
        val.temp <- vals[[i]]
        for (j in 1:num.methods){
            points(alpha.list, val.temp[[j]],
                   type = 'l', col = cols[j],
                   lty = ltys[j])
            points(alpha.list[alpha.at],
                   val.temp[[j]][alpha.at], col = cols[j],
                   pch = pchs[j])
        }
        if(i == num.plots.per.row) {
            legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
                   methods, col = cols, lty = ltys, pch = pchs)
        }
    }

    dev.off()

}

plot.FDR <- function(filename, vals,
                     mains, methods, cols,
                     ylim = c(0, 0.35),
                     ...){
    plot.results(filename, vals, mains, methods, cols,
                 ylab = "FDR", ylim = ylim, ...)
}

plot.power <- function(filename, vals,
                       mains, methods, cols,
                       ylim = c(0, 1),
                     ...){
    plot.results(filename, vals, mains, methods, cols,
                 ylab = "Power", ylim = ylim, ...)
}
