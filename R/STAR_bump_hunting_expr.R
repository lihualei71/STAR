library("arules")
library("xtable")
library("ggplot2")
source("STAR_bump_hunting.R")

rej.to.box <- function(rej, x){
    lapply(x, function(covar){
        new.levels <- levels(factor(covar[rej]))
        all.levels <- levels(factor(covar))
        if (length(new.levels) < length(all.levels)){
            return(new.levels)
        } else {
            return()
        }
    })
}

find.most.frequent <- function(items){
    n <- length(items)
    freq <- sapply(1:n, function(i){
        sum(sapply(1:n, function(j){
            identical(items[[i]], items[[j]])
        }))
    })
    if (max(freq) < 0.5 * n) {
        return(NA)
    } else {
        ind <- which.max(freq)[1]
        return(list(item = items[[ind]], freq = max(freq) / n))
    }
}

data(IncomeESL)
data <- IncomeESL[complete.cases(IncomeESL), ]
y <- data$income
levels(y)

prop <- summary(y) / length(y)
n <- length(y)
cum.prop <- cumsum(summary(y)) / length(y)
x <- data[, -1]

box20 <- list()
box10 <- list()
box5 <- list()
for (seed in 1:100){
    set.seed(seed)
    print(seed)
    pvals <- as.numeric(1 - cum.prop[y] + prop[y] * runif(n))
    result <- STAR.bump.hunting(
        pvals, x,
        fun.list = create.fun("SeqStep", 0.7),
        print.quiet = TRUE)
    box20[[seed]] <- rej.to.box(result$mask[, 20], x)
    box10[[seed]] <- rej.to.box(result$mask[, 10], x)
    box5[[seed]] <- rej.to.box(result$mask[, 5], x)
}

box20.result <- lapply(names(box20[[1]]), function(item){
    lapply(box20, function(box){
        box[[item]]
    })
})
box20.box <- lapply(box20.result, function(x){
    find.most.frequent(x)$item
})
names(box20.box) <- names(box20[[1]])
box20.freq <- lapply(box20.result, function(x){
    find.most.frequent(x)$freq
})
names(box20.freq) <- names(box20[[1]])

box10.result <- lapply(names(box10[[1]]), function(item){
    lapply(box10, function(box){
        box[[item]]
    })
})
box10.box <- lapply(box10.result, function(x){
    find.most.frequent(x)$item
})
names(box10.box) <- names(box10[[1]])
box10.freq <- lapply(box10.result, function(x){
    find.most.frequent(x)$freq
})
names(box10.freq) <- names(box10[[1]])

box5.result <- lapply(names(box5[[1]]), function(item){
    lapply(box5, function(box){
        box[[item]]
    })
})
box5.box <- lapply(box5.result, function(x){
    find.most.frequent(x)$item
})
names(box5.box) <- names(box5[[1]])
box5.freq <- lapply(box5.result, function(x){
    find.most.frequent(x)$freq
})
names(box5.freq) <- names(box5[[1]])

df <- data.frame(box20_freq= unlist(box20.freq),
                 box10_freq= unlist(box10.freq),
                 box5_freq= unlist(box5.freq))

box <- c("-", "married/single", "[18, 54]", "higher than high school", "professional/managerial/student", ">10", "not married/yes", "[2,4]", "<=2", "own", "house", "white", "english")
df$box <- box

#### Produce Table 1
## xtable(df)
## capture.output(xtable(df), file = "bump_hunting.txt")

select.index <- sapply(names(box20.box), function(item){
    select.attr <- box20.box[[item]]
    if (!is.null(select.attr)){
        x[[item]] %in% select.attr
    } else {
        rep(TRUE, nrow(x))
    }
})
select.index <- apply(select.index, 1, all)

select.samples <- y[select.index]
other.samples <- y[!select.index]
hist.select <- table(select.samples) / length(select.samples)
hist.others <- table(other.samples) / length(other.samples)

df.select.hist <- data.frame(
    category = "Selected", 
    income = as.integer(select.samples))
df.others.hist <- data.frame(
    category = "Others",
    income = as.integer(other.samples))
df.hist <- rbind(df.select.hist, df.others.hist)
## df.hist$income <- as.factor(df.hist$income)

bump.hunting.plot <- ggplot(df.hist, aes(income, fill = category)) +
    geom_histogram(alpha = 0.7, aes(y = ..density..),
                   position = "identity",
                   bins = 9) + theme_bw() +
    scale_x_continuous(expand = c(0, 0),
                       breaks = 1:9,
                       labels = levels(y)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.grid = element_blank()) +
    xlab("Income Class") + ylab("Frequency") +
    ggtitle("Comparison of income distributions")
ggsave(filename = "../figs/bump_hunting.pdf", bump.hunting.plot, height = 5)
