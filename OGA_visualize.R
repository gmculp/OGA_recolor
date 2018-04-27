###load dichromat simulation functions###
source('BVM_recolor.R')
    
###load Orange-Gray-Azure Re-coloring algorithm###
source('OGA_recolor.R')

############################
###load required packages###
############################
library("png")
library("jpeg")		
library("grid")
library("gridExtra")
    
###################################################################################	
###function to display original and re-colored images along with CVD simulations### 	###################################################################################

OGA_image_compare <- function(in_image)	{	
	
	###CVD simulation of original image###
    img.deutan <- sim_img(in_image,"deutan")
    img.protan <- sim_img(in_image,"protan")
    img.tritan <- sim_img(in_image,"tritan")

    ###re-color image###
    in_image2 <- OGA_recolor_image(in_image)
    
    ###CVD simulation of re-colored image###
    img.deutan2 <- sim_img(in_image2,"deutan")
    img.protan2 <- sim_img(in_image2,"protan")
    img.tritan2 <- sim_img(in_image2,"tritan")
    
    ###display output###
    grid.arrange(rasterGrob(in_image), rasterGrob(in_image2), 
        rasterGrob(img.deutan), rasterGrob(img.deutan2), 
        rasterGrob(img.protan), rasterGrob(img.protan2), 
        rasterGrob(img.tritan), rasterGrob(img.tritan2), nrow=4)
	

}
	
#####################################################################################
###function to display original and re-colored palettes along with CVD simulations### 	
#####################################################################################

OGA_palette_compare <- function(in_palette)	{

    oga_pal <- OGA_recolor_hex_palette(in_palette)

    pal_list <- list()
    pal_list["normal"] <- ""
    pal_list["normal.original"] <- list(in_palette)
    pal_list["normal.recolored"] <- list(oga_pal)
    pal_list["protan"] <- ""
    pal_list["protan.original"] <- list(CVD_p(in_palette))
    pal_list["protan.recolored"] <- list(CVD_p(oga_pal))
    pal_list["deutan"] <- ""
    pal_list["deutan.original"] <- list(CVD_d(in_palette))
    pal_list["deutan.recolored"] <- list(CVD_d(oga_pal))
    pal_list["tritan"] <- ""
    pal_list["tritan.original"] <- list(CVD_t(in_palette))
    pal_list["tritan.recolored"] <- list(CVD_t(oga_pal))
    pal_list <- rev(pal_list)
    
    nr <- length(pal_list)
    nc <- length(in_palette)

    plot(1, 1, xlim = c(0, nc), ylim = c(0, nr), type = "n", 
    axes = FALSE, bty = "n", xlab = "", ylab = "")

    for (i in 1:nr) {
        this_pal <- unlist(pal_list[i], use.names = FALSE)
        ni <- length(this_pal)
	
          if (ni == 1) next
	
          rect(xleft = 0:(ni - 1), ybottom = i - 1, xright = 1:ni, 
	    ytop = i - 0.2, col = this_pal, border = "light grey")
   }

   text(rep(-0.1, nr), (1:nr) - 0.6, 
    labels = ifelse(grepl("\\.",names(pal_list)),gsub("^.*\\.","",names(pal_list)),""), 
    xpd = TRUE, adj = 1)

   text(rep(nc/2, nr), (1:nr) - 0.6, 
    labels = ifelse(grepl("\\.",names(pal_list)),"",names(pal_list)), 
    xpd = TRUE)

	#return(oga_pal)
	
}	