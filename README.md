# OGA_recolor
This R function re-colors images with an orange-gray-blue palette in an effort to make them more accessible to people with color vision deficiency (CVD) of either the red/green (protan, deutan) or blue/green (tritan) type.  To read more about this algorithm, please refer to my [dissertation](http://academicworks.cuny.edu/gc_etds/1243/).  An additional function is provided which uses the algorithm described in [Brettel et al. (1997)](http://vision.psychol.cam.ac.uk/jdmollon/papers/Dichromat_simulation.pdf) and [Vie´not et al. (1999)](http://vision.psychol.cam.ac.uk/jdmollon/papers/colourmaps.pdf)  to simulate the three types of CVD. 

Using the R function to re-color a palette:

        ```R
        ###load dichromat simulation functions###
        source('BVM_recolor.R')
        
        ###load Orange-Gray-Azure Re-coloring algorithm###
        source('OGA_recolor.R')
        
        my_pal <- c("#009300","#FF7300","#6666FF","#00FF00","#FFBF33","#66D9FF")
        oga_pal <- OGA_recolor_hex_palette(my_pal)
	
        pal_list <- list()
        pal_list["normal"] <- ""
	pal_list["normal.original"] <- list(my_pal)
	pal_list["normal.recolored"] <- list(oga_pal)
	pal_list["protan"] <- ""
	pal_list["protan.original"] <- list(CVD_p(my_pal))
	pal_list["protan.recolored"] <- list(CVD_p(oga_pal))
	pal_list["deutan"] <- ""
	pal_list["deutan.original"] <- list(CVD_d(my_pal))
	pal_list["deutan.recolored"] <- list(CVD_d(oga_pal))
	pal_list["tritan"] <- ""
	pal_list["tritan.original"] <- list(CVD_t(my_pal))
	pal_list["tritan.recolored"] <- list(CVD_t(oga_pal))
	pal_list <- rev(pal_list)
        
        nr <- length(pal_list)
	nc <- length(my_pal)
  
	plot(1, 1, xlim = c(0, nc), ylim = c(0, nr), type = "n", 
	    axes = FALSE, bty = "n", xlab = "", ylab = "")
	
	for (i in 1:nr) {
		
		this_pal <- unlist(pal_list[i], use.names = FALSE)
		
		ni <- length(this_pal)
		
		if (ni == 1) 
			next
		
		rect(xleft = 0:(ni - 1), ybottom = i - 1, xright = 1:ni, 
		    ytop = i - 0.2, col = this_pal, border = "light grey")
	}
	
	text(rep(-0.1, nr), (1:nr) - 0.6, 
	    labels = ifelse(grepl("\\.",names(pal_list)),gsub("^.*\\.","",names(pal_list)),""), 
	    xpd = TRUE, adj = 1)
	
	text(rep(nc/2, nr), (1:nr) - 0.6, 
	    labels = ifelse(grepl("\\.",names(pal_list)),"",names(pal_list)), 
	    xpd = TRUE)
        ```


Using the R function to re-color an image:

        ```R
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
        
        ###load image###
        my.image <- readPNG("candies_color.png")
        
        ###CVD simulation of original image###
        img.deutan <- sim_img(my.image,"deutan")
        img.protan <- sim_img(my.image,"protan")
        img.tritan <- sim_img(my.image,"tritan")

        ###re-color image###
        my.image2 <- OGA_recolor_image(my.image)
        
        ###CVD simulation of re-colored image###
        img.deutan2 <- sim_img(my.image2,"deutan")
        img.protan2 <- sim_img(my.image2,"protan")
        img.tritan2 <- sim_img(my.image2,"tritan")
        
        ###display output###
        grid.arrange(rasterGrob(my.image), rasterGrob(my.image2), 
            rasterGrob(img.deutan), rasterGrob(img.deutan2), 
            rasterGrob(img.protan), rasterGrob(img.protan2), 
            rasterGrob(img.tritan), rasterGrob(img.tritan2), nrow=4)
        ```
