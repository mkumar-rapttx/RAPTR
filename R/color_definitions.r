# Some RAPT heatmap color definitions

colors <- list(
	phase = c(
		M=rapt_colors()[["yellow"]],
		C=rapt_colors()[["orange"]],
		X=rapt_colors()[["blue"]]
	),
	phase_txt = c(M="black",C="white",X="black"),
	
	# Teal
	dose = c(
		`25mg`= "#D2EEEA",
		`50mg`= "#88C3C8",
		`75mg`= "#498EA4",
		`100mg`= "#2A5676"
	),
	dose_txt = c(`25mg`="black", `50mg`="black",`75mg`="black",`100mg`="white"),
	
	# BluYl colors
	biopsy = c(
		Screen =  "#F9FFAF",
		Treat =  "#9CDAA1",
		Paired =  "#028190"
	),
	biopsy_txt = c(Screen="black",Treat="black",Paired="white"),
	
	# Oranges
	charged = c(
		Low = "#FEF5EC",
		Med = "#FFBB80",
		High = "#E96200"
	),
	charged_txt = c(Low="black",Med="black",High="black"),
	
	# Terrain colors
	bor = c(
		NE = "#F1F1F1",
		PD = "#FFC59E",
		SD = "#E1BB4E",
		PR = "#9BB306",
		CR = "#26A63A"
	),
	bor_txt = c(NE="black",PD="white",SD="black",PR="black",CR="black")
	
)

require(colorspace)
# use colorRampPalette to get the first color after white
purple_range <- c(colorRampPalette(c("white","#DA95CC"))(10)[c(2,10)])
yellow_range <- c(colorRampPalette(c("white",rapt_colors()[["yellow"]]))(10)[c(2,10)], rapt_colors()[["orange"]])
blue_range <- c(colorRampPalette(c("white",rapt_colors()[["blue"]]))(10)[c(2,10)], rapt_colors()[["dark blue"]])
red_green_range <- c("#C1478A","white","#5CB35D")
heat_range <- rev(sequential_hcl(6,"Heat"))
bugn_range <- rev(sequential_hcl(6,"BuGn"))
pkyl_range <- rev(sequential_hcl(4,"PinkYl"))
lajo_range <- sequential_hcl(6,"Lajolla")
powderblue_range <- c(colorRampPalette(c("white","#9FA2FF"))(10)[c(2,10)])

