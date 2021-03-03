#' Calculates a robust Z-score that does not have "exploding" values due to invariant/nearly-invariant vectors
#'
#' @param x a matrix/table of values to calculate z-scores from
#' @param min a minimum standard deviation to use
#' @param min.quantile a quantile of the non-zero standard deviations to use as a minimum
#' @param by calculate z-scores across "row" or "col"
#' @return A matrix of z-scores.
#' If code{min} is provided, then that value is used as the minimum standard deviation for calculating z-scores.
#' If code{min} is not provided, but code{min.quantile} is, then that value is used as the determine the minimum standard deviation for calculating z-scores by choosing the SD at that quantile amongs all the non-zero SDs.
#' If neither is provided, then a standard z-score (code{(x - mean(x))/standard deviation(x)}) is calculated
#' @export
#'
#' @examples
#' x <- matrix(c(rnorm(30, mean=10, sd = 1), rnorm(30, mean=10, sd = 0.1)), nrow=10)
#' 
#' # calculate a standard z-score
#' out <- robust_z(x, by="col")
#' apply(out, 2, max) # look at max z-score per column
#'
#' # calculate a robust z-score, using a lower bound of 0.1 for standard deviations
#' out2 <- robust_z(x, min=0.1, by="col")
#' apply(out2, 2, max) # look at max z-score per column
#'
#' # calculate a robust z-score, using a lower bound of 25th percentile for standard deviations
#' out3 <- robust_z(x, min.quantile=0.25, by="col")
#' apply(out3, 2, max) # look at max z-score per column
robust_z <- function(x, min=NULL, min.quantile=NULL, by="row") {
	direction <- if_else(by == "row", 1, 2)
	means <- apply(x, direction, mean, na.rm=TRUE)
	y <- sweep(x, direction, means, "-")
	sds <- apply(y, direction, sd, na.rm=TRUE)
	if (!is.null(min)) { 
		sds <- trunc_range(sds, low=min) 
	} else if (!is.null(min.quantile)) {
		min <- quantile(sds[sds > 0], min.quantile, na.rm=T)
		sds <- trunc_range(sds, low=min)
	}
	z <- sweep(y, direction, sds, "/")
	z
}
#' Truncates a values to lower and/or upper bounds
#'
#' @param x input vector/matrix of values
#' @param low lower limit to truncate to
#' @param high upper limit to truncate to
#' @param filter if code{TRUE} filters out truncated values, otherwise replace them
#' @param low.val an alternate value to replace the lower limit with, rather than the lower limit itself
#' @param high.val an alternate value to replace the upper limit with, rather than the upper limit itself
#' @examples
#' x <- rnorm(500, mean=0, sd=5)
#' low.val/high.val is useful for differentiating truncated values on a plot
#' x1 <- trunc_range(x, low=-5, high=5, low.val=-6, high.val=6)
#' plot(x1 ~ x)
trunc_range <- function(x, low=NA, high=NA, filter=FALSE, low.val=NA, high.val=NA) {
	if (filter) {
		if (!is.na(low)) x <- x[x>=low]
		if (!is.na(high)) x <- x[x<=high]
	} else {
		if (!is.na(low)) x[x<low] <- ifelse(is.na(low.val), low, low.val)
		if (!is.na(high)) x[x>high] <- ifelse(is.na(high.val), high, high.val)
	}
	return(x)
}
