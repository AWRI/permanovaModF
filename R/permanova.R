#' Calculates modified F statistic (Andersson 2017) and permutation p-value
#'
#' @description performs a one way ANOVA and calculates modified F statistic to account for unbalanced designs
#' with heterogeneous multivariate dispersion. Carries out a permutaiton test to calculate significance.
#' Based on the permanova function in Fathom Toolbox.
#'
#' @details Authours: Chris Ward
#'
#'
#' @param mat matrix of ordination data
#' @param metaData data_frame like object with columns with same name as \code{groups}
#' @param groups character vector containing colnames to test across
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
#' @rdname permanovaModF
#' @export



permanovaModF <- function(mat, metaData, groups = c(), permutations = 1000L){


  iter <- permutations
  n <- dim(mat)[1]

  A <- -0.5*(mat^2)
  I <- matrix(0, n, n)
  diag(I) <- 1

  uno <- matrix(1,n,n)

  G <- permanovaModF:::e((I-(1/n)*(uno)), A, (I-(1/n)*(uno)))

  all <- lapply(colnames(metaData[groups]), function(var){

    grp <- as.numeric(factor(metaData[[var]]))
    uGrp   <-  unique(grp)
    nGrp   <- length(uGrp)

    idxGrp <- sapply(uGrp, function(x){
      as.integer(grp == x)
    })


    sGrp <- colSums(idxGrp)

    input_matrix <- cbind(uno[1], idxGrp[, 1:(ncol(idxGrp) - 1)])

    qr_result <- qr(input_matrix, LAPACK = TRUE)

    Q1 <- qr.Q(qr_result)
    H <-  Q1%*%t(Q1)


    R_dia <- diag(permanovaModF:::e((I - H), G, (I - H)))

    V <- sapply(1:nGrp, function(y) {
      sum(R_dia[as.logical(idxGrp[,y])]) / (sGrp[y] - 1)
    })


    HGH <- sum(diag(permanovaModF:::e(H, G, H)))


    F_denom <- sum((1 - (sGrp / n)) * V)

    F2 <- HGH / F_denom

    permute_iteration <- function(j) {
      idxS <- sample(1:n)

      R_dia <- diag(permanovaModF:::e((diag(n) - H), G[idxS, idxS], (diag(n) - H)))

      sapply(1:nGrp, function(i){
        sum(R_dia[idxGrp[, i]]) / (sGrp[i] - 1)
      })

      diag_perm <- diag(permanovaModF:::e(H, G[idxS, idxS], H))
      sum(diag_perm) / sum((1 - (sGrp / n)) * V)
    }


    Fperm <- sapply(2:iter, permute_iteration)

    p <- sum(c(F2, Fperm) >= F2) / iter

    tibble(variable = var, F2, p)

  }) %>% bind_rows()

  return(all)

}




#######
