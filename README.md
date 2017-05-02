# OGA_recolor
This function re-colors images with an orange-gray-blue palette in an effort to make them more accessible to people with color vision deficiency (CVD) of either the red/green (protan, deutan) or blue/green (tritan) type.  To read more about this algorithm, please refer to my [dissertation](http://academicworks.cuny.edu/gc_etds/1243/).  

Using the function to recolor an image:

        ```R
        #requires devtools
        library(devtools)
        #make sure Rtools is TRUE
        devtools::find_rtools()
        #install package
        devtools::install_github("gmculp/rGBAT16AB")
        #load package
        library(rGBAT16AB)
        ```
