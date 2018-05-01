library(data.table)
library(pracma)


###gamma value###
gamma <- 2.4

###sRGB to XYZ matrix Y column values for ref. white D65###
Y_r <- 0.2126

Y_g <- 0.7152

Y_b <- 0.0722

###white point###
w_pt <- 1/3

###rounding factor###
rnd_fac <- 0.02

inter_calc1 <- function(n){
  return(ifelse(n > 0.04045, (((n + 0.055)/1.055) ^ gamma), (n/12.92)))
}

inter_calc2 <- function(n){
  return(ifelse(n > 0.0031308, ((1.055 * (n ^ (1/gamma))) - 0.055), (n * 12.92)))
}

calc_pd <- function(xa,ya,xb,yb) {
  return((((xb - xa) ^ 2) + ((yb - ya) ^ 2)) ^ (1/2))
}


#chan_func <- function(rgb,i){
#	return(inter_calc1(rgb[i])/sum(inter_calc1(rgb)))
#}

chan_func <- function(chan,tot) {
  return(ifelse(chan > 0, chan/tot, ifelse(tot==0,1/3,0)))
}

calc_Y <- function(R,G,B) {
  (R * Y_r) + (G * Y_g) + (B * Y_b)
}


#calc_Y <- function(r,g,b) {
#	return((inter_calc1(r) * Y_r) + (inter_calc1(g) * Y_g) + (inter_calc1(b) * Y_b))
#}

prop_chan <- function(chan,dY){
  
  #dY <- Y1/Y2
  
  out_chan <- pmin(pmax(as.numeric(inter_calc2(chan * dY)),0),1)
  
  return(out_chan)
}


####
####
####

func_rr.a <- function(t) inter_calc1(t)/(inter_calc1(t)+inter_calc1((0.6*t)+0.4)+inter_calc1(1))
func_bb.a <- function(t) inter_calc1(1)/(inter_calc1(t)+inter_calc1((0.6*t)+0.4)+inter_calc1(1))
func_Y.a <- function(t) calc_Y(inter_calc1(t), inter_calc1((t*0.6)+0.4), inter_calc1(1))

func_rr.o <- function(t) inter_calc1(1)/(inter_calc1(t)+inter_calc1((0.6*t)+0.4)+inter_calc1(1))
func_bb.o <- function(t) inter_calc1(t)/(inter_calc1(t)+inter_calc1((0.6*t)+0.4)+inter_calc1(1))
func_Y.o <- function(t) calc_Y(inter_calc1(1), inter_calc1((t*0.6)+0.4), inter_calc1(t))

####
####
####







find_root <- function(in.rr,in.bb,in.Y,func_Y,func_rr,func_bb) {
	
	t1 <- try(pracma::brent(function(t) func_Y(t) - in.Y, 0, 1)$root,silent = TRUE)
	t1 <- ifelse(class(t1) == "try-error",0,t1)

	t2.func <- function(t) (((in.rr - func_rr(t))^2) + ((in.bb - func_bb(t))^2))^(1/2)
	t2 <- pracma::fibsearch(t2.func, t1, 1, endp = TRUE)$xmin
	
	return(t2)
	
}


OGA_recolor <- function(test.dt) {

  #in_vec <- RColorBrewer::brewer.pal(12, "Paired"); test.dt <- data.table(t(col2rgb(in_vec))/255)
  
  test.dt[, col_id := .I]
  test.dt[ , alpha := 1]
  test.dt[, in.r2 := inter_calc1(red)]
  test.dt[, in.g2 := inter_calc1(green)]
  test.dt[, in.b2 := inter_calc1(blue)]
  test.dt[, in.rgb := in.r2 + in.g2 + in.b2]
  test.dt[, in.rr := chan_func(in.r2,in.rgb)]
  test.dt[, in.bb := chan_func(in.b2,in.rgb)]
  test.dt[, in.Y := calc_Y(in.r2,in.g2,in.b2)]
  
  test.dt1 <- copy(test.dt)
  test.dt1[, new.t := unlist(lapply(1:nrow(test.dt1), function (x) find_root(test.dt1[x,]$in.rr,test.dt1[x,]$in.bb,test.dt1[x,]$in.Y,func_Y.a,func_rr.a,func_bb.a)))]
  test.dt1[, new.rr := func_rr.a(new.t)]
  test.dt1[, new.bb := func_bb.a(new.t)]
  
  test.dt2 <- copy(test.dt)
  test.dt2[, new.t := unlist(lapply(1:nrow(test.dt2), function (x) find_root(test.dt2[x,]$in.rr,test.dt2[x,]$in.bb,test.dt2[x,]$in.Y,func_Y.o,func_rr.o,func_bb.o)))]
  test.dt2[, new.rr := func_rr.o(new.t)]
  test.dt2[, new.bb := func_bb.o(new.t)]
  
  test.dt <- rbindlist(list(test.dt1,test.dt2))
  test.dt[, pt_d := calc_pd(new.rr,new.bb,in.rr,in.bb)]
  test.dt <- test.dt[test.dt[, .I[which.min(pt_d)], by = col_id]$V1]
  
  test.dt[, new.gg := 1 - new.rr - new.bb] 
  test.dt[, new.Y := calc_Y(new.rr,new.gg,new.bb)]

  test.dt[, new_red := prop_chan(new.rr,in.Y/new.Y)]
  test.dt[, new_green := prop_chan(new.gg,in.Y/new.Y)]
  test.dt[, new_blue := prop_chan(new.bb,in.Y/new.Y)]
  
  return(test.dt[, c("col_id","red","green","blue","alpha","new_red","new_green","new_blue"), with = FALSE])

}

###function to re-color vector of hex values###
OGA_recolor_hex_palette <- function(in_vec) {
  pal.dt <- OGA_recolor(data.table(t(col2rgb(in_vec))/255))
  return(rgb(pal.dt$new_red,pal.dt$new_green,pal.dt$new_blue))
}

########################################
###wrapper function to recolor images###
########################################	

OGA_recolor_image <- function(img) {
  
  #img <- readPNG("candies_color.png")
  
  ###convert matrix to data.table###
  ###round colors by rnd_fac###
  img.dt <- data.table(
    red = round(matrix(img[,,1], ncol=1)/rnd_fac)*rnd_fac,
    green = round(matrix(img[,,2], ncol=1)/rnd_fac)*rnd_fac,
    blue = round(matrix(img[,,3], ncol=1)/rnd_fac)*rnd_fac
  )
  
  setnames(img.dt, gsub(".V1","",colnames(img.dt)))
  
  ###run recolor algorithm on data.table of unique pixel colors###
  img.dt2A <- OGA_recolor(unique(img.dt))
  
  img.dt[,pix.order := .I]
  
  img.dt2 <- merge(img.dt,img.dt2A,by=c("red","green","blue"))
  
  rm(img.dt2A)
  
  setorder(img.dt2, pix.order)
  
  ##################################
  ###return back to matrix format###
  ##################################
  
  R = matrix(as.numeric(img.dt2$new_red), nrow=dim(img)[1])
  
  G = matrix(as.numeric(img.dt2$new_green), nrow=dim(img)[1])
  
  B = matrix(as.numeric(img.dt2$new_blue), nrow=dim(img)[1])
  
  img.new = array(dim=dim(img))
  img.new[,,1] = R
  img.new[,,2] = G
  img.new[,,3] = B
  
  if (dim(img)[3]==4) {
    ###cat("4-dim image")
    img.new[,,4] = matrix(as.numeric(img.dt2$alpha), nrow=dim(img)[1])
  }
  
  rm(img, img.dt, img.dt2, R, G, B)
  
  invisible(gc())
  
  return(img.new)
  
}



