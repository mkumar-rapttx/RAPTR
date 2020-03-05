#' Returns the first or last rows and columns of a subsettable object (matrix, data frame, tibble, data table). Alternative to head or tail.
#'
#' @param x An object
#' @param m a single positive integer. The number of rows to return.
#' @param n a single positive integer. The number of columns to return.
#' @export
corner = function(x, m = 5, n = 5)
{
    if (ncol(x) < n ) {n = ncol(x)}
    if (nrow(x) < m ) {m = nrow(x)}
    print(x[1:m,1:n])
}
