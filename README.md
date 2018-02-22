# OGA_recolor
This R function re-colors images with an orange-gray-blue palette in an effort to make them more accessible to people with color vision deficiency (CVD) of either the red/green (protan, deutan) or blue/green (tritan) type.  To read more about this algorithm, please refer to my [dissertation](http://academicworks.cuny.edu/gc_etds/1243/).  An additional function is provided which uses the algorithm described in [Brettel et al. (1997)](http://vision.psychol.cam.ac.uk/jdmollon/papers/Dichromat_simulation.pdf) and [Vie´not et al. (1999)](http://vision.psychol.cam.ac.uk/jdmollon/papers/colourmaps.pdf)  to simulate the three types of CVD. 

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
        grid.arrange(rasterGrob(my.image), rasterGrob(img.deutan), rasterGrob(img.protan), rasterGrob(img.tritan), rasterGrob(my.image2), rasterGrob(img.deutan2), rasterGrob(img.protan2), rasterGrob(img.tritan2), nrow=4)
        ```
