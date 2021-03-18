#' Customized version of scales::col_numeric with extra control of color ranges
#' @usage col_numeric_min_max("heat", min.val=0, max.val=10)
#' @usage col_numeric_min_max("heat", extend.min=0, extend.max=100)
#' @usage col_numeric_min_max(c("pink","white","forestgreen"), symmetric=TRUE)
#'
#' @param palette a valid palette name from `RColorBrewer::brewer.pal.info`, from `viridis` ("viridis", "magma", "inferno", "plasma") or two or more colors for a color ramp
#' @param na.color the color value to use for displaying NAs
#' @param alpha the alpha (transparency) value to use, or FALSE
#' @param reverse a logical. If TRUE, reverse the color direction
#' @param min.val,max.val a double. `NULL` or a minimum/maximum value at which to truncate the displayed data.  Values above/below this value will all be displayed at the min/max color.
#' @param extend.min,extend.max a double. `NULL` or a minimum/maximum value to which the color range extends.  This allows the color range to span outside of the data range.
#' @param symmetric a logical. If TRUE, sets the color range to -/+ the absolute maximum value of the data.  Such that the range spans equal lengths less than and greather than zero.
#'
#' @details
#' This is a modified version of scales::col_numeric with extra arguments.  
#' `min/max.val` allow you to "crush" the color range so that the min/max colors are pinned to these values rather than the true min/max of the data.
#' `extend.min/max` do the opposite and allow you to extend the color range outside of the data range.  This is useful for placing the data range in a greater context than the data itself allows.
#' `symmetric` is useful for tri-color gradients (e.g. red/white/green) for representing negative/zero/positive values.  It dynamically calculates `extend.min/max` such that the minimum/maximum values for the color range are `c(-1,+1) * max(abs(data))`
#' If `extend.min/max` are within the actual range of the data, they have no effect.

col_numeric_min_max <- function (palette, na.color = "white", alpha = FALSE, reverse = FALSE, 
								 min.val=NULL, max.val=NULL, 
								 extend.min=NULL, extend.max=NULL, 
								 symmetric=FALSE) {
	
	require(scales, quietly=TRUE)

		# Note a bunch of functions don't get exported by scales.  Need to access them with `scales:::`
	pf <- scales:::safePaletteFunc(palette, na.color, alpha)

	scales:::withColorAttr("numeric", list(na.color = na.color), function(x) {
		if (length(x) == 0 || all(is.na(x))) return(pf(x))
		# truncate to min.val and max.val
		if (!is.null(min.val)) x[!is.na(x) & x < min.val] <- min.val
		if (!is.null(max.val)) x[!is.na(x) & x > max.val] <- max.val
		
		# extend to extend.min and extend.max
		if (!is.null(extend.min)) {min. <- min(c(extend.min,x),na.rm=T)} else { min. <- min(x,na.rm=T)}
		if (!is.null(extend.max)) {max. <- max(c(extend.max,x),na.rm=T)} else { max. <- max(x,na.rm=T)}
		# mirror min and mix if `symmetric`
		if (symmetric) { max. <- max(abs(c(min.,max.))); min. <- -max. }

		rescaled <- (x - min.)/(max.-min.)
		if (reverse) rescaled <- 1 - rescaled
		pf(rescaled)
	})
}



