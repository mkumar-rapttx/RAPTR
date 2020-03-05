#' theme_custom
#'
#' @param base_size
#' @param x_label_size
#' @param y_label_size
#' @param x_tick_label_size
#' @param y_tick_label_size
#' @param y_label_angle Angle to rotate y axis label.
#' @param title_size
#' @param title_just
#' @param legend_position
#'
#' @return
#' @export
#'
theme_custom <- function(base_size = 16,
                         x_label_size = NA, y_label_size = NA,
                         x_tick_label_size = NA, y_tick_label_size = NA,
                         y_label_angle = 90,
                         title_size = NA, title_just = 0.5,
                         legend_position = "top") {
    # Check size inputs
    x_label_size = ifelse(is.na(x_label_size),base_size,x_label_size)
    y_label_size = ifelse(is.na(y_label_size),base_size,y_label_size)
    x_tick_label_size = ifelse(is.na(x_tick_label_size),0.67*base_size,x_tick_label_size)
    y_tick_label_size = ifelse(is.na(y_tick_label_size),0.67*base_size,y_tick_label_size)
    title_size = ifelse(is.na(title_size),base_size,title_size)

    # Modity theme_bw
    theme_bw( base_size = base_size,
	          base_line_size = base_size / 32,
	          base_rect_size = base_size / 32) %+replace%
        theme(
	     # panel.border = element_rect(fill = NA, colour = "grey95")
	      panel.grid = element_line(colour = "grey95"),
	      panel.grid.minor = element_line(size = rel(0.5)),
	      #strip.background = element_rect(fill = "grey90", colour = "grey80"),
	      legend.position = legend_position,
	      legend.background = element_rect(fill = "transparent",color = NA),
	      legend.key = element_rect(fill = "transparent", color = NA),
	      plot.title = element_text(size= title_size, hjust = title_just),
	      panel.background = element_blank(),
	      axis.text.x = element_text(size = x_tick_label_size),
	      axis.text.y = element_text(size = y_tick_label_size),
	      axis.title.x = element_text(size = x_label_size),
	      axis.title.y = element_text(size = y_label_size, angle = y_label_angle)

	    #  complete = TRUE
	    )
}



#' Title
#'
#' @param color
#' @param range
#' @param space
#'
#' @return
#' @export
#'
heatmap_color = function(colors = c("#0A95A7","white","#E5541B"), range = c(-1,0,1), space = "LAB")
{
    return(circlize::colorRamp2(range, colors))
}


#rapt.extended <- c(gray="#5B666F", blue="#0A95A7", orange="#E5541B", yellow="#FE9F33", green="#4DA167", lilac="#947EB0", `dark blue`="#083D77")
#rapt.darker.extended <- c(gray2="#444D54",blue2="#06707E",orange2="#AC3F12",yellow2="#D07920", green2="#356E46", lilac2="#69597D", `dark blue2`="#052344")

# rapt_colors <- function(...) {
# 	cols <- c(...)
# 	full <- c(rapt.extended, rapt.darker.extended)
# 	if (is.null(cols)) return(full)
# 	return(full[cols])
# }
# rapt_pal <- function(palette = "main", reverse = FALSE, ...) {
# 	if (stringr::str_detect(palette, "dark")) {
# 		pal <- rapt.darker.extended
# 	} else if (stringr::str_detect(palette, "both")) {
# 		pal <-  c(rapt.extended, rapt.darker.extended)
# 	} else {
# 		pal <- rapt.extended
# 	}
# 	if (reverse) pal <- rev(pal)
#
# 	function(n) {
# 		if (n > length(pal)) stop("Palette ",palette," has a total of ",length(pal)," colors, but you need ",n,"!")
# 		return(unname(pal)[seq_len(n)])
# 	}
# }
#
# scale_color_rapt <- function(palette="main", reverse = FALSE, ...) {
# 	pal <- rapt_pal(palette = palette, reverse = reverse)
# 	discrete_scale("colour", paste0("rapt_", palette), palette = pal, ...)
# }
# scale_colour_rapt <- scale_color_rapt
#
#
# scale_fill_rapt <- function(palette = "main", reverse = FALSE, ...) {
# 	pal <- rapt_pal(palette = palette, reverse = reverse)
# 	discrete_scale("fill", paste0("rapt_", palette), palette = pal, ...)
# }
#
# # data(mtcars)
# # ggplot(mtcars, aes(x=wt, y=drat, color=as.character(cyl))) + geom_point(size=3) + facet_wrap(~cyl) +
# #	theme_rapt() + scale_color_rapt()
