##############################################################################
## STAR for Wavelet Thresholding
##################################################################
############

library("waveslim")
library("igraph")
library("EBImage")
source("STAR_tree.R")
source("Lynch.R")

est.sigma <- function(wd.obj, type=1){
    if (type == 1){
        sigma.HL <- mad(wd.obj$HL1)
        sigma.LH <- mad(wd.obj$LH1)
        sigma.HH <- mad(wd.obj$HH1)
    } else if (type == 2){
        sigma.HL <- sd(wd.obj$HL1)
        sigma.LH <- sd(wd.obj$LH1)
        sigma.HH <- sd(wd.obj$HH1)
    }
    return(list(HL = sigma.HL, LH = sigma.LH, HH = sigma.HH))
}

dwt2mat <- function(wd.obj, nlevels, sigma.type = 1){
    sigma <- est.sigma(wd.obj, sigma.type)
    n <- 2^nlevels
    coef.mat <- matrix(NA, nrow = n, ncol = n)
    pval.mat <- matrix(NA, nrow = n, ncol = n)
    LL.name <- paste0("LL", nlevels)
    coef.mat[1, 1] <- wd.obj[[LL.name]]
    pval.mat[1, 1] <- 0
    for (level in nlevels:1){
        tmp <- nlevels - level
        
        HL.name <- paste0("HL", level)
        HL.coef.mat <- wd.obj[[HL.name]]        
        HL.row.ind <- 1:(2^tmp)
        HL.col.ind <- 1:(2^tmp) + 2^tmp
        coef.mat[HL.row.ind, HL.col.ind] <- HL.coef.mat
        pval.mat[HL.row.ind, HL.col.ind] <-
            2 * (1 - pnorm(abs(HL.coef.mat / sigma$HL)))
        
        LH.name <- paste0("LH", level)
        LH.coef.mat <- wd.obj[[LH.name]]        
        LH.row.ind <- 1:(2^tmp) + 2^tmp
        LH.col.ind <- 1:(2^tmp)
        coef.mat[LH.row.ind, LH.col.ind] <- LH.coef.mat
        pval.mat[LH.row.ind, LH.col.ind] <-
            2 * (1 - pnorm(abs(LH.coef.mat / sigma$LH)))
        
        HH.name <- paste0("HH", level)
        HH.coef.mat <- wd.obj[[HH.name]]
        HH.row.ind <- 1:(2^tmp) + 2^tmp
        HH.col.ind <- 1:(2^tmp) + 2^tmp
        coef.mat[HH.row.ind, HH.col.ind] <- HH.coef.mat
        pval.mat[HH.row.ind, HH.col.ind] <-
            2 * (1 - pnorm(abs(HH.coef.mat / sigma$HH)))
    }
    return(list(coef = coef.mat, pval = pval.mat))
}


mat2dwt <- function(mat, wavelet, boundary){
    n <- nrow(mat)
    nlevels <- floor(log(n, 2))
    obj <- list()
    attributes(obj)$J <- nlevels
    attributes(obj)$wavelet <- wavelet
    attributes(obj)$boundary <- boundary
    class(obj) <- "dwt.2d"

    for (level in 1:nlevels){
        tmp <- nlevels - level

        LH.name <- paste0("LH", level)
        LH.row.ind <- 1:(2^tmp) + 2^tmp
        LH.col.ind <- 1:(2^tmp)
        obj[[LH.name]] <- mat[LH.row.ind, LH.col.ind,
                              drop = FALSE]

        HL.name <- paste0("HL", level)
        HL.row.ind <- 1:(2^tmp)
        HL.col.ind <- 1:(2^tmp) + 2^tmp
        obj[[HL.name]] <- mat[HL.row.ind, HL.col.ind,
                              drop = FALSE]
        
        HH.name <- paste0("HH", level)
        HH.row.ind <- 1:(2^tmp) + 2^tmp
        HH.col.ind <- 1:(2^tmp) + 2^tmp
        obj[[HH.name]] <- mat[HH.row.ind, HH.col.ind,
                              drop = FALSE]
    }
    LL.name <- paste0("LL", nlevels)
    obj[[LL.name]] <- mat[1, 1, drop = FALSE]

    return(obj)
}

gen.edge <- function(mat){
    n <- nrow(mat)
    m <- floor((n - 1) / 2)
    ind.list <- expand.grid(0:m, 0:m)
    edgelist <- sapply(2:nrow(ind.list), function(j){
        ind1 <- ind.list[j, 1]
        ind2 <- ind.list[j, 2]
        parent <- ind1 * n + ind2 + 1
        child.LL <- 2 * ind1 * n + 2 * ind2 + 1
        child.HL <- child.LL + 1
        child.LH <- child.LL + n
        child.HH <- child.LH + 1
        c(parent, child.LL, parent, child.HL,
          parent, child.LH, parent, child.HH)
    })
    edgelist <- cbind(edgelist[1:2,], edgelist[3:4,],
                      edgelist[5:6,], edgelist[7:8,])
    edgelist <- cbind(c(1, 2), c(1, 3), c(1, 4), edgelist)
    edgelist <- data.frame(t(edgelist))
    names(edgelist) <- c("from", "to")
    edgelist <- edgelist[order(edgelist$from), ]
    return(edgelist)
}

STAR.wavelet <- function(image, alpha.list,
                         wavelet = "haar",
                         boundary = "periodic",
                         sigma.type = 1,
                         plot.posit = FALSE,
                         print.quiet = TRUE,
                         ...){
    n <- dim(image)[1]
    m <- floor(log(n, 2))
    nalphas <- length(alpha.list)
    im.dwt <- dwt.2d(image, wavelet, m, boundary)
    wavelet.mat <- dwt2mat(im.dwt, m, sigma.type)
    pvals <- as.numeric(wavelet.mat$pval)
    coefs <- wavelet.mat$coef

    edgelist <- gen.edge(coefs)
    tree <- graph_from_data_frame(edgelist)
    ind <- as.numeric(V(tree)$name)
    pvals <- pvals[ind]

    mod <- STAR.tree(pvals, tree,
                     fun.list = create.fun("SeqStep", 0.5),
                     alpha.list = alpha.list,
                     print.quiet = print.quiet,
                     ...)
    num.rej <- apply(mod$mask, 2, sum)

    rec.images <- lapply(1:nalphas, function(i){
        remove.ind <- V(tree)$name[!mod$mask[, i]]
        mat.copy <- coefs
        mat.copy[as.numeric(remove.ind)] <- 0
        rec.dwt <- mat2dwt(mat.copy, wavelet, boundary)
        img <- idwt.2d(rec.dwt)
        return(img)
    })

    if (plot.posit){
        par(mfrow = c(nalphas, 1))
        for (i in 1:nalphas){
            ind <- V(tree)$name[mod$mask[, i]]
            ind <- as.numeric(ind)
            mat <- matrix(0, nrow = n, ncol = n)
            mat[ind] <- 1
            image(as.Image(mat)[1:32, 32:1])
        }
    }
    return(list(images = rec.images, num.rej = num.rej))
}

BH.wavelet <- function(image, alpha.list,
                       wavelet = "haar",
                       boundary = "periodic",
                       sigma.type = 1,
                       plot.posit = FALSE){
    n <- dim(image)[1]
    m <- floor(log(n, 2))
    nalphas <- length(alpha.list)    
    im.dwt <- dwt.2d(image, wavelet, m, boundary)
    wavelet.mat <- dwt2mat(im.dwt, m, sigma.type)
    pvals <- as.numeric(wavelet.mat$pval)
    coefs <- wavelet.mat$coef

    rec.images <- list()
    num.rej <- rep(NA, nalphas)
    rejs <- matrix(NA, nrow = n^2, ncol = nalphas)
    for (i in 1:nalphas){
        alpha <- alpha.list[i]
        khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
        alpha <- alpha * khat / length(pvals)
        rejs[, i] <- (pvals <= alpha)
        mat.copy <- coefs
        mat.copy[!rejs[, i]] <- 0
        rec.dwt <- mat2dwt(mat.copy, wavelet, boundary)
        rec.images[[i]] <- idwt.2d(rec.dwt)
        num.rej[i] <- sum(rejs[, i])
    }

    if (plot.posit){
        par(mfrow = c(nalphas, 1))
        for (i in 1:nalphas){
            ind <- rejs[, i]
            mat <- matrix(0, nrow = n, ncol = n)
            mat[ind] <- 1
            image(as.Image(mat)[1:32, 32:1])
        }
    }    
    return(list(images = rec.images, num.rej = num.rej))
}

    
## thresh.wavelet <- function(image, 
##                            wavelet = "haar",
##                            boundary = "periodic",
##                            sigma.type = 1){
##     n <- dim(image)[1]
##     m <- floor(log(n, 2))
##     im.dwt <- dwt.2d(image, wavelet, m, boundary)    
##     wavelet.mat <- dwt2mat(im.dwt, m)
##     sigma <- est.sigma(wavelet.mat, sigma.type)    
##     zvals <- abs(as.numeric(wavelet.mat / sigma))

##     ## Hard Thresholding
##     thresh <-  sigma * sqrt(2 * log(n))
##     ind <- (abs(zvals) >= thresh)
##     num.rej <- sum(ind)
##     wavelet.mat[!ind] <- 0
##     rec.dwt <- mat2dwt(wavelet.mat, wavelet, boundary)
##     rec.images.hard <- idwt.2d(rec.dwt)

##     ## Soft Thresholding
##     wavelet.mat[ind] <- sign(wavelet.mat[ind]) *
##         (abs(wavelet.mat[ind]) - thresh)
##     rec.dwt <- mat2dwt(wavelet.mat, wavelet, boundary)
##     rec.images.soft <- idwt.2d(rec.dwt)
    
##     return(list(hard.images = rec.images.hard,
##                 soft.images = rec.images.soft,
##                 num.rej = num.rej))
## }

thresh.wavelet <- function(image, 
                           wavelet = "haar",
                           boundary = "periodic",
                           sigma.type = 1){
    soft <- function(x, delta) sign(x) * pmax(abs(x) - delta, 
        0)
    hard <- function(x, delta) ifelse(abs(x) > delta, x, 0)
    
    n <- dim(image)[1]
    m <- floor(log(n, 2))    
    im.dwt <- dwt.2d(image, wavelet, m, boundary)
    sigma <- est.sigma(im.dwt, sigma.type)
    thresh <- list(HL = sigma$HL * 2 * sqrt(log(n)),
                   LH = sigma$LH * 2 * sqrt(log(n)),
                   HH = sigma$HH * 2 * sqrt(log(n)))

    im.dwt.hard <- im.dwt
    im.dwt.soft <- im.dwt
    num.rej <- 0
    
    for (level in 1:m) {
        HL.name <- paste("HL", level, sep = "")
        im.dwt.hard[[HL.name]] <- hard(im.dwt[[HL.name]], thresh$HL)
        im.dwt.soft[[HL.name]] <- soft(im.dwt[[HL.name]], thresh$HL)
        LH.name <- paste("LH", level, sep = "")
        im.dwt.hard[[LH.name]] <- hard(im.dwt[[LH.name]], thresh$LH)
        im.dwt.soft[[LH.name]] <- soft(im.dwt[[LH.name]], thresh$LH)
        HH.name <- paste("HH", level, sep = "")
        im.dwt.hard[[HH.name]] <- hard(im.dwt[[HH.name]], thresh$HH)
        im.dwt.soft[[HH.name]] <- soft(im.dwt[[HH.name]], thresh$HH)
        num.rej <- num.rej +
            sum(im.dwt.hard[[HL.name]] != 0) +
            sum(im.dwt.hard[[LH.name]] != 0) +
            sum(im.dwt.hard[[HH.name]] != 0)
    }
    
    rec.images.hard <- idwt.2d(im.dwt.hard)
    rec.images.soft <- idwt.2d(im.dwt.soft)
    
    return(list(hard.images = rec.images.hard,
                soft.images = rec.images.soft,
                num.rej = num.rej))
}

gen.noisy.image <- function(image, sigma = 10){
    n <- nrow(image)
    image.noisy <- image + sigma * rnorm(n^2)
    image.noisy[image.noisy < 0] <- 0
    image.noisy[image.noisy > 1] <- 1
    return(image.noisy)
}

compare.wavelet <- function(img.file, output.file = NULL,
                            SNR = 10,
                            STAR.alpha = 0.2,
                            BH.alpha = 0.05,
                            wavelet = "haar",
                            save.plot = FALSE){
    img <- readImage(img.file)
    img <- channel(img, "gray")
    n <- dim(img)[2]
    img <- img[, n:1]

    img.noisy <- gen.noisy.image(img, sd(img) / sqrt(SNR))

    result.STAR <- STAR.wavelet(img.noisy, STAR.alpha, wavelet)
    result.BH <- BH.wavelet(img.noisy, BH.alpha, wavelet)
    result.thresh <- thresh.wavelet(img.noisy, wavelet)

    img.STAR <- result.STAR$images[[1]]
    img.BH <- result.BH$images[[1]]
    img.hard <- result.thresh$hard.images
    img.soft <- result.thresh$soft.images

    CR.STAR <- round(n^2 / result.STAR$num.rej[1], 0)
    CR.BH <- round(n^2 / result.BH$num.rej[1], 0)
    CR.hard <- round(n^2 / result.thresh$num.rej, 0)
    CR.soft <- round(n^2 / result.thresh$num.rej, 0)    
        
    SNR.STAR <- round(20 * log(sd(img)/sd(img.STAR - img), 10)
, 2)
    SNR.BH <- round(20 * log(sd(img)/sd(img.BH - img), 10), 2)
    SNR.hard <- round(20 * log(sd(img)/sd(img.hard - img), 10), 2)
    SNR.soft <- round(20 * log(sd(img)/sd(img.soft - img), 10), 2)

    img.name <- strsplit(img.file, ".png")[[1]][1]
    if (is.null(output.file)){
        output.file <- paste0(img.name, "_", wavelet, "_SNR_",
                              SNR, "_STAR_", STAR.alpha, "_BH_",
                              BH.alpha, ".pdf")
    }

    if (save.plot){
        pdf(output.file, width = 10.5, height = 7)
        par(mfrow = c(2, 3), mar = c(2, 2, 2, 2), cex.main = 2)
        image(img, main = "Original")
        image(as.Image(img.STAR),
              main = paste0("STAR (SNR = ", SNR.STAR,
                  ", CR = ", CR.STAR, ")"))
        image(as.Image(img.BH),
              main = paste0("BH (SNR = ", SNR.BH,
                  ", CR = ", CR.BH, ")"))
        image(as.Image(img.noisy),
              main = paste0("Noisy (SNR = ", SNR, ")"))
        image(as.Image(img.hard),
              main = paste0("H.T. (SNR = ", SNR.hard,
                  ", CR = ", CR.hard, ")"))
        image(as.Image(img.soft),
              main = paste0("S.T. (SNR = ", SNR.soft,
                  ", CR = ", CR.soft, ")"))
        dev.off()
    }
    
    CR <- c(CR.STAR, CR.BH, CR.hard, CR.soft)
    SNR <- c(SNR.STAR, SNR.BH, SNR.hard, SNR.soft)
    return(list(CR = CR, SNR = SNR))
}
