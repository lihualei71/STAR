source("STAR_wavelet.R")

SNR <- 0.5
STAR.alpha <- 0.2
BH.alpha <- 0.05
wavelet <- "la8"
CR.low.list <- list()
SNR.low.list <- list()
for (i in 1:48){
    img.file <- paste0("../images/", i, ".png")
    save.plot <- (i == 8)
    result <- compare.wavelet(img.file, SNR = SNR,
                              wavelet = wavelet,
                              STAR.alpha = STAR.alpha,
                              BH.alpha = BH.alpha,
                              save.plot = save.plot)
    print(paste0(i, "th image finished!"))
    CR.low.list[[i]] <- result$CR
    SNR.low.list[[i]] <- result$SNR
}

SNR <- 10
STAR.alpha <- 0.2
BH.alpha <- 0.05
wavelet <- "la8"
CR.high.list <- list()
SNR.high.list <- list()
for (i in 1:48){
    img.file <- paste0("../images/", i, ".png")
    save.plot <- (i == 8)    
    result <- compare.wavelet(img.file, SNR = SNR,
                              wavelet = wavelet,
                              STAR.alpha = STAR.alpha,
                              BH.alpha = BH.alpha,
                              save.plot = save.plot)
    print(paste0(i, "th image finished!"))
    CR.high.list[[i]] <- result$CR
    SNR.high.list[[i]] <- result$SNR
}

CR.list <- list(CR.low.list, CR.high.list)
SNR.list <- list(SNR.low.list, SNR.high.list)
save(file = "../data/STAR_wavelet_CR.RData", CR.list)
save(file = "../data/STAR_wavelet_SNR.RData", SNR.list)
