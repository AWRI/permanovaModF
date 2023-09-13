#' Calculates pairwise modified F statistic (Andersson 2017) and permutation p-value for
#' each comparison of factors within a group
#'
#' @description performs a pairwise one way ANOVA and calculates modified F statistic to account for unbalanced designs
#' with heterogeneous multivariate dispersion. Carries out a permutaiton test to calculate significance.
#' Based on the permanova function in Fathom Toolbox.
#'
#' @details Authours: Chris Ward
#'
#'
#' @param mat matrix of ordination data
#' @param metaData data_frame like object with one column with same name as \code{group}, must include
#' a \code{sample_id} column with sample id in the same order as \code{mat}
#' @param group character vector containing colname to test across
#' @param permutations integer number of permutations (including initial observation)
#'
#' @importFrom readr read_csv
#' @import dplyr
#' @import Rcpp
#' @import RcppArmadillo
#' @import RcppEigen
#'
#' @return A \code{data_frame} containing modified F and permutation p values
#'
#' @rdname permanovaModF_PW
#' @export

permanovaModF_PW <- function(mat, metaData, group = c(), permutations = 10000L){

  metaData <- metaData[c("sample_id", group)]
  colnames(metaData)[2] <- "group"
  mdList <- split(metaData, metaData$group)

  pairs <- .generatePairs(mdList)

  all <- lapply(pairs, function(x){

    md_x <- filter(metaData, group %in% x)
    index <- which(colnames(mat) %in% md_x$sample_id)
    mat_x <- mat[index, index]

    df_x <- permanovaModF(mat = mat_x, metaData = md_x, groups = "group", permutations = permutations)

    tibble(P1 = x[1], P2 = x[2], variable = group, F2 = df_x$F2, p = df_x$p)


  }) %>% bind_rows()

  all$t <- sqrt(all$F2)
  all$holm.p <- p.adjust(all$p, method = "holm")
  return(all[c(1:4,6,5,7)])
}


.generatePairs <- function(mdList, ...){
  nam <- names(mdList)
  pairs <- outer(nam, nam, paste, sep = "///")
  pairs <- lapply(pairs,function(z){
    pair <- unlist(strsplit(z, split = "///"))
    if(pair[1] == pair[2]) pair <- NULL
    pair <- sort(pair)
    return(pair)
  })
  pairs <- pairs[!duplicated(pairs)]
  Filter(Negate(is.null), pairs)
}

