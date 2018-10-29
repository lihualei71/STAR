##############################################################################
## STAR for Bump Hunting
##############################################################################

source("generic_STAR.R")

#### Find_Candidate_Fun
find.variable <- function(covar, mask, prop,
                          ord.cut = 1){
    temp.var <- covar[mask, ]
    n <- sum(mask)
    candids <- lapply(temp.var, function(x){
        type <- class(x)
        if (length(type) == 1 && type == "factor"){
            sapply(levels(factor(x)), function(val){
                x == val
            })
        } else if (length(type) == 2 && type[1] == "ordered"){
            val.set <- levels(factor(x))
            m <- length(val.set)
            if (length(val.set) <= ord.cut){
                return(rep(TRUE, n))
            } else {
                left.peel <- x %in% val.set[1:ord.cut]
                right.peel <- x %in% val.set[m:(m-ord.cut+1)]
                cbind(left.peel, right.peel)
            }
        }
    })
    candids <- Reduce(cbind, candids)
    ind <- which(mask)
    candids <- lapply(1:ncol(candids), function(i){
        temp <- rep(FALSE, length(mask))
        candid <- candids[, i]
        temp[ind[candid]] <- TRUE
        temp
    })
    return(candids)
}

#### Update_Mask_Fun
## Naive Update
bump.hunting.mask.update <- function(candid, score, mask, prop,
                                     dir = c("min", "max")){
    dir <- dir[1]
    candid.vals <- sapply(candid, function(set){
        mean(score[!set])
    })
    if (dir == "max"){
        reveal.inds <- candid[[which.max(candid.vals)]]
    } else {
        reveal.inds <- candid[[which.min(candid.vals)]]
    }
    mask[reveal.inds] <- FALSE
    return(mask)
}

####
STAR.bump.hunting <- function(
    pvals, x, 
    fun.list = create.fun("SeqStep", pstar = 0.5),
    alpha.list = seq(0.01, 0.3, 0.01),
    type = c("naive", "model-assist"),
    ...){

    type <- type[1]
    find.candid.fun <- find.variable
    update.mask.fun <- bump.hunting.mask.update
    if (type == "naive"){
        score.fun <- NULL
    } else if (type == "model-assist"){
        score.fun <- NULL
    }

    generic.STAR(pvals = pvals, covar = x, type = type,
                fun.list = fun.list, alpha.list = alpha.list,
                score.fun = score.fun,
                find.candid.fun = find.candid.fun,
                update.mask.fun = update.mask.fun,
                ...)
}
