source("STAR_wavelet.R")
library("ggplot2")
library("dplyr")
library("reshape2")

load("../data/STAR_wavelet_CR.RData")
load("../data/STAR_wavelet_SNR.RData")

CR.low.list <- CR.list[[1]]
CR.high.list <- CR.list[[2]]
SNR.low.list <- SNR.list[[1]]
SNR.high.list <- SNR.list[[2]]

BH.low.CR.ratio <- sapply(CR.low.list, function(x){
    x[2] / x[1]
})
hard.low.CR.ratio <- sapply(CR.low.list, function(x){
    x[3] / x[1]
})
soft.low.CR.ratio <- sapply(CR.low.list, function(x){
    x[4] / x[1]
})
BH.low.SNR.ratio <- sapply(SNR.low.list, function(x){
    x[2] / x[1]
})
hard.low.SNR.ratio <- sapply(SNR.low.list, function(x){
    x[3] / x[1]
})
soft.low.SNR.ratio <- sapply(SNR.low.list, function(x){
    x[4] / x[1]
})

low.CR.ratio <- data.frame(BH.low.CR.ratio,
                           hard.low.CR.ratio,
                           soft.low.CR.ratio)
names(low.CR.ratio) <- c("BH", "hard", "soft")
low.CR.ratio <- melt(low.CR.ratio)
low.CR.ratio <- data.frame(type = "CR", low.CR.ratio)

low.SNR.ratio <- data.frame(BH.low.SNR.ratio,
                            hard.low.SNR.ratio,
                            soft.low.SNR.ratio)
names(low.SNR.ratio) <- c("BH", "hard", "soft")
low.SNR.ratio <- melt(low.SNR.ratio)
low.SNR.ratio <- data.frame(type = "SNR", low.SNR.ratio)

low.ratio <- rbind(low.CR.ratio, low.SNR.ratio)

low.CR.plot <- ggplot(low.CR.ratio, aes(x = variable, y = value)) +
    geom_boxplot(width = 0.5) + xlab("Methods") + ylab("Ratio of CR") + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20)) + geom_hline(yintercept = 1, color = "red") +
    ggtitle("Comparison of Compression Ratio")
ggsave("../figs/low_CR.pdf", low.CR.plot, width = 7, height = 7)
low.SNR.plot <- ggplot(low.SNR.ratio, aes(x = variable, y = value)) +
    geom_boxplot(width = 0.5) + xlab("Methods") + ylab("Ratio of SNR") + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20)) + geom_hline(yintercept = 1, color = "red") +
    ggtitle("Comparison of Signal-to-Noise Ratio")
ggsave("../figs/low_SNR.pdf", low.SNR.plot, width = 7, height = 7)



BH.high.CR.ratio <- sapply(CR.high.list, function(x){
    x[2] / x[1]
})
hard.high.CR.ratio <- sapply(CR.high.list, function(x){
    x[3] / x[1]
})
soft.high.CR.ratio <- sapply(CR.high.list, function(x){
    x[4] / x[1]
})
BH.high.SNR.ratio <- sapply(SNR.high.list, function(x){
    x[2] / x[1]
})
hard.high.SNR.ratio <- sapply(SNR.high.list, function(x){
    x[3] / x[1]
})
soft.high.SNR.ratio <- sapply(SNR.high.list, function(x){
    x[4] / x[1]
})

high.CR.ratio <- data.frame(BH.high.CR.ratio,
                           hard.high.CR.ratio,
                           soft.high.CR.ratio)
names(high.CR.ratio) <- c("BH", "hard", "soft")
high.CR.ratio <- melt(high.CR.ratio)
high.CR.ratio <- data.frame(type = "CR", high.CR.ratio)

high.SNR.ratio <- data.frame(BH.high.SNR.ratio,
                            hard.high.SNR.ratio,
                            soft.high.SNR.ratio)
names(high.SNR.ratio) <- c("BH", "hard", "soft")
high.SNR.ratio <- melt(high.SNR.ratio)
high.SNR.ratio <- data.frame(type = "SNR", high.SNR.ratio)

high.ratio <- rbind(high.CR.ratio, high.SNR.ratio)

high.CR.plot <- ggplot(high.CR.ratio, aes(x = variable, y = value)) +
    geom_boxplot(width = 0.5) + xlab("Methods") + ylab("Ratio of CR") + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20)) + geom_hline(yintercept = 1, color = "red") +
    ggtitle("Comparison of Compression Ratio")        
ggsave("../figs/high_CR.pdf", high.CR.plot, width = 7, height = 7)
high.SNR.plot <- ggplot(high.SNR.ratio, aes(x = variable, y = value)) +
    geom_boxplot(width = 0.5) + xlab("Methods") + ylab("Ratio of CR") + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20)) + geom_hline(yintercept = 1, color = "red") +
    ggtitle("Comparison of Signal-to-Noise Ratio")        
ggsave("../figs/high_SNR.pdf", high.SNR.plot, width = 7, height = 7)
