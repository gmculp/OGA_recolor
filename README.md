# OGA_recolor
This set of R functions re-color palettes and images with an orange-gray-blue palette in an effort to make them more accessible to people with color vision deficiency (CVD) of either the red/green (protan, deutan) or blue/green (tritan) type.  To read more about this algorithm, please refer to my [dissertation](http://academicworks.cuny.edu/gc_etds/1243/).  An additional function is provided which uses the algorithm described in [Brettel et al. (1997)](http://vision.psychol.cam.ac.uk/jdmollon/papers/Dichromat_simulation.pdf) and [Vie´not et al. (1999)](http://vision.psychol.cam.ac.uk/jdmollon/papers/colourmaps.pdf)  to simulate the three types of CVD. 

Using the R function to re-color a palette:

        ```R
       ###load Orange-Gray-Azure Re-coloring algorithm visualization functions###
       source('OGA_visualize.R')

       ###display original and re-colored palettes along with CVD simulations###

       ###good use case for use of re-coloring algorithm###
       ###input palette contains CVD inaccessible color pairs###
       my_palette <- c("#009300","#FF7300","#6666FF","#00FF00","#FFBF33","#66D9FF")
       OGA_palette_compare(my_palette)

       ###get re-colored palette###
       oga_pal <- OGA_recolor_hex_palette(my_palette)

       ###bad use case for use of re-coloring algorithm###
       ###input palette (in this case, the ColorBrewer PRGn palette) is already CVD accessible###
       my_palette <- c("#762A83","#AF8DC3","#E7D4E8","#D9F0D3","#7FBF7B","#1B7837")
       OGA_palette_compare(my_palette)
       ```


Using the R function to re-color an image:

        ```R
        ###load Orange-Gray-Azure Re-coloring algorithm visualization functions###
       source('OGA_visualize.R')
        
       ###good use case for use of re-coloring algorithm###
       ###input image contains CVD inaccessible color pairs###
       my.image <- readPNG("candies_color.png")
       OGA_image_compare(my.image)

       ###get re-colored image###
       oga_image <- OGA_recolor_hex_palette(my_palette)
       writePNG(oga_image,"oga_image.png")

       ###bad use case for use of re-coloring algorithm###
       ###input image is already CVD accessible###
       my.image <- readJPEG("Magenta_flower.JPG")
       OGA_image_compare(my.image)
        ```
