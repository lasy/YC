
#' Compute the Yue and Clayton similarity between the rows of one or two matrices
#'
#' @param mat1 a \code{matrix} or \code{data.frame} where each row contains the data for a sample.
#' @param mat2 (optional) a \code{matrix} or \code{data.frame} where each row contains the data for a sample.
#'
#' @return a similarity matrix.
#' If \code{mat2} is \code{NULL}, \code{YC} returns a square matrix of dimension
#' \code{n}, the number of rows of \code{mat1}.
#' Otherwise, \code{YC} returns a rectangular matrix with \code{n1 = nrow(mat1)} rows
#' and \code{n2 = nrow(mat2)} columns. Each element \code{[i,j]} gives the
#' Yue and Clayton similarity between \code{mat1[i,]} and \code{mat2[j,]}.
#' @export
#'
#' @importFrom dplyr bind_rows
#' @import magrittr
#'
YC <- function(mat1, mat2 = NULL) {

  if (any(mat1 < 0, na.rm = TRUE) | any(mat1 > 1, na.rm = TRUE))
    stop("`mat1` values should stricktly be in [0, 1].\n")

  if (is.null(mat2)) mat2 <- mat1

  # TODO: check that mat1 and mat2 have colnames (and rownames)

  if (any(mat2 < 0, na.rm = TRUE) | any(mat2 > 1, na.rm = TRUE))
    stop("`mat2` values should stricktly be in [0, 1].\n")


  # first, if mat1 and mat2 don't have the same features (colnames), we harmonize them
  mat <-
    bind_rows(
      mat1 %>% as.data.frame(),
      mat2 %>% as.data.frame()
    ) %>%
    as.matrix()
  mat[is.na(mat)] = 0

  m1 <- mat[1:nrow(mat1),] %>% matrix(nrow = nrow(mat1), ncol = ncol(mat))
  m2 <- mat[nrow(mat1) + (1:nrow(mat2)),] %>% matrix(nrow = nrow(mat2), ncol = ncol(mat))

  # Yue and Clayton similarity is
  # theta = \sum_i p_i q_i / (\sum_i (p_i - q_i)^2  + \sum_i p_i q_i)

  # N is \sum_i p_i q_i
  # M is \sum_i (p_i - q_i)^2
  # theta = N / (M + N)

  N <- m1 %*% t(m2)

  dims <- c(nrow(m1), ncol(mat), nrow(m2))
  M1 <- array(m1, dim = dims)
  M2 <- array(rep(t(m2), each = nrow(m1)), dim = dims)
  M <- apply((M1 - M2)^2, c(1,3), sum)

  D <- N / (M + N)
  rownames(D) <- rownames(mat1)
  colnames(D) <- rownames(mat2)
  D
}
