
###DISABLE SCIENTIFIC NOTATION###
options(scipen = 999)

############################
###load required packages###
############################
suppressMessages(library(data.table))


####################
###some constants###
####################

###gamma value###
gamma <- 2.4

###sRGB to XYZ matrix Y column values for ref. white D65###
Y_r <- 0.2126

Y_g <- 0.7152

Y_b <- 0.0722

###white point###
w_pt <- 1/3

###azure###
a_r <- 0; a_g <- 0.4; a_b <- 1

###orange###
o_r <- 1; o_g <- 0.4; o_b <- 0

###rounding factor###
rnd_fac <- 0.02

###segment factor###
seg_fac <- 0.1

####################
###some functions###
####################

calc_Y <- function(R,G,B) {
  (R * Y_r) + (G * Y_g) + (B * Y_b)
}

inter_calc1 <- function(n){
  return(ifelse(n > 0.04045, (((n + 0.055)/1.055) ^ gamma), (n/12.92)))
}

inter_calc2 <- function(n){
  return(ifelse(n > 0.0031308, ((1.055 * (n ^ (1/gamma))) - 0.055), (n * 12.92)))
}

chan_func <- function(chan,tot) {
  return(ifelse(chan > 0, chan/tot, ifelse(tot==0,1/3,0)))
}


calc_u <- function(x1,y1,x2,y2,x3,y3) {
  (((x3 - x1) * (x2 - x1)) + ((y3 - y1) * (y2 - y1))) / (((x2 - x1) ^ 2) + ((y2 - y1) ^ 2))
}


calc_ux <- function(x1,y1,x2,y2,x3,y3) {
  ##xu = x1 + u (x2 - x1)
  u <- calc_u(x1,y1,x2,y2,x3,y3)
  ux <- ifelse(u >= 0 & u <=1 & !is.na(u), x1 + (u * (x2 - x1)),1000)
  return(ux) 
}

calc_uy <- function(x1,y1,x2,y2,x3,y3) {
  ##yu = y1 + u (y2- y1)
  u <- calc_u(x1,y1,x2,y2,x3,y3)
  uy <- ifelse(u >= 0 & u <=1 & !is.na(u), y1 + (u * (y2 - y1)),1000)
  return(uy) 
}

calc_pd <- function(xa,ya,xb,yb) {
  (((xb - xa) ^ 2) + ((yb - ya) ^ 2)) ^ (1/2)
}

prop_chan <- function(chan,Y1,Y2){
  
  out_chan <- pmin(pmax(as.numeric(inter_calc2(chan * (Y1/Y2))),0),1)
  
  return(out_chan)
}

###http://stackoverflow.com/questions/10600060/how-to-do-cross-join-in-r###
cjdt <- function(a,b){
  cj = CJ(1:nrow(a),1:nrow(b))
  cbind(a[cj[[1]],],b[cj[[2]],])
}


##########################################
###assemble data table of line segments###
##########################################

func_rr.a <- function(t) inter_calc1(t)/(inter_calc1(t)+inter_calc1((0.6*t)+0.4)+inter_calc1(1))
func_bb.a <- function(t) inter_calc1(1)/(inter_calc1(t)+inter_calc1((0.6*t)+0.4)+inter_calc1(1))
func_Y.a <- function(t) (Y_r * (inter_calc1(t))) + (Y_g * (inter_calc1((0.6*t)+0.4))) + (Y_b * (inter_calc1(1)))

func_rr.o <- function(t) inter_calc1(1)/(inter_calc1(t)+inter_calc1((0.6*t)+0.4)+inter_calc1(1))
func_bb.o <- function(t) inter_calc1(t)/(inter_calc1(t)+inter_calc1((0.6*t)+0.4)+inter_calc1(1))
func_Y.o <- function(t) (Y_r * (inter_calc1(1))) + (Y_g * (inter_calc1((0.6*t)+0.4))) + (Y_b * (inter_calc1(t)))


make_segments <- function(t,func_rr,func_bb,func_Y){
	dt.seg <- data.table(t.2=t,rr.2=func_rr(t),bb.2=func_bb(t),Y.2=func_Y(t))
	dt.seg[, t.1 := shift(.(t.2), type = "lag")]
	dt.seg[, rr.1 := shift(.(rr.2), type = "lag")]
	dt.seg[, bb.1 := shift(.(bb.2), type = "lag")]	
	dt.seg[, Y.1 := shift(.(Y.2), type = "lag")]	
	dt.seg <- dt.seg[!(is.na(rr.1))]
	return(dt.seg)
}

###select points that are equidistant###
opt_make_segments <- function(seg_num, func_rr, func_bb, func_Y){
	t <- seq(0,1,by=0.001)
	dt.seg <- make_segments(t,func_rr,func_bb,func_Y)
	dt.seg[, dist := (((rr.2-rr.1)^2) + ((bb.2-bb.1)^2) + ((Y.2-Y.1)^2)) ^ (1/2)]
	dt.seg[, dist_cumsum := cumsum(dist)]
	div_pt <- sum(dt.seg$dist)/seg_num
	t_pix <- rbindlist(lapply(1:(seg_num-1), function(i) dt.seg[which.min(abs(dist_cumsum - (div_pt*i)))]))$t.2
	dt.seg <- make_segments(c(0,t_pix,1),func_rr,func_bb,func_Y)
	return(dt.seg)
}

###construct segment look-up table###
nn <- 10
seg.dt <- rbindlist(list(opt_make_segments(nn,func_rr.o,func_bb.o,func_Y.o),opt_make_segments(nn,func_rr.a,func_bb.a,func_Y.a)))
#seg.dt <- rbindlist(list(make_segments(seq(0,1,by=(1/nn)),func_rr.o,func_bb.o,func_Y.o),make_segments(seq(0,1,by=(1/nn)),func_rr.a,func_bb.a,func_Y.a)))

###############################
###main re-coloring function###
###############################

OGA_recolor <- function (test.dt) {
  
  test.dt[ , col_id := 1:nrow(test.dt)]
  
  test.dt[ , alpha := 1]
  
  #setcolorder(test.dt, c("col_id","red","green","blue","alpha"))
  
  test.dt[ , in.r2 := inter_calc1(red)]
  test.dt[ , in.g2 := inter_calc1(green)]
  test.dt[ , in.b2 := inter_calc1(blue)]
  test.dt[ , in.rgb2 := in.r2 + in.g2 + in.b2]
  
  test.dt[ , in.Y := calc_Y(in.r2,in.g2,in.b2)]
  
  test.dt[ , in.rr := chan_func(in.r2,in.rgb2)]
  
  test.dt[ , in.bb := chan_func(in.b2,in.rgb2)]
  
  test.dt[, c("in.r2","in.g2","in.b2","in.rgb2"):=NULL]
  
  ###many-to-many merge with line segments###
  test.dt <- cjdt(test.dt,seg.dt)
  
  ###remove segments that have lower lightness than pixel###
  test.dt <- test.dt[in.Y <= Y.1 | in.Y <= Y.2]
  
  ###sort by pixel id
  setorder(test.dt, col_id)
  
  #########################################################################
  ###find point at which pixel color lightness channel intersects curves###
  #########################################################################
  
  test.dt[ , pt.rr := rr.2 + (((in.Y - Y.2)/(Y.1 - Y.2))*(rr.1-rr.2))]
  test.dt[ , pt.bb := bb.2 + (((in.Y - Y.2)/(Y.1 - Y.2))*(bb.1-bb.2))]
  
  test.dt[ , rr.1 := round(ifelse(in.Y >= Y.1, pt.rr, rr.1),6)]
  test.dt[ , bb.1 := round(ifelse(in.Y >= Y.1, pt.bb, bb.1),6)]
  
  test.dt[, c("pt.rr","pt.bb","Y.1","Y.2"):=NULL] 
  
  ##################################################
  ###closest point on chromaticity lines to color###
  ##################################################
  
  ###data table containing closest point on line... excludes line's start and end point###
  test.dt1 <- copy(test.dt)
  test.dt1[, new.rr := calc_ux(rr.1,bb.1,rr.2,bb.2,in.rr,in.bb)]
  test.dt1[, new.bb := calc_uy(rr.1,bb.1,rr.2,bb.2,in.rr,in.bb)]
  test.dt1[, pt_d := calc_pd(new.rr,new.bb,in.rr,in.bb)]
  
  ###data table for start point distance###
  test.dt2 <- copy(test.dt)
  test.dt2[, new.rr := rr.1]
  test.dt2[, new.bb := bb.1]
  test.dt2[, pt_d := calc_pd(rr.1,bb.1,in.rr,in.bb)]
  
  ###data table for end point distance###
  test.dt3 <- copy(test.dt)
  test.dt3[, new.rr := rr.2]
  test.dt3[, new.bb := bb.2]
  test.dt3[, pt_d := calc_pd(rr.2,bb.2,in.rr,in.bb)]
  
  ###bind all potential points together###
  test.dt <- rbindlist(list(test.dt1,test.dt2,test.dt3))
  
  ###extract rows with minimum point distances###
  test.dt <- test.dt[test.dt[, .I[which.min(pt_d)], by = col_id]$V1]
  
  ###calculate greenness channel###
  test.dt[ , new.gg :=  1 - new.rr - new.bb]
  
  #########################
  ###convert back to RGB###
  #########################
  
  test.dt[ , temp.Y := calc_Y(new.rr, new.gg, new.bb)]
  
  test.dt[ , new_red :=  prop_chan(new.rr,in.Y,temp.Y)]
  test.dt[ , new_green :=  prop_chan(new.gg,in.Y,temp.Y)]
  test.dt[ , new_blue :=  prop_chan(new.bb,in.Y,temp.Y)]
  
  ###remove temporary columns###
  d_c <- colnames(test.dt)[!colnames(test.dt) %in% c("col_id","red","green","blue","alpha","new_red","new_green","new_blue")]
  
  test.dt[, (d_c):=NULL]
  
  ###sort by color id
  setorder(test.dt, col_id)
  
  invisible(gc())
  
  return(test.dt)
}



########################################
###wrapper function to recolor images###
########################################	

OGA_recolor_image <- function(img) {
  
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
  
  img.dt[,pix.order := 1:nrow(img.dt)]
  
  img.dt2 <- merge(img.dt,img.dt2A,by=c("red","green","blue"))
  
  rm(img.dt2A)
  
  img.dt2 <- setorder(img.dt2, pix.order)
  
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

###function to re-color vector of hex values###
OGA_recolor_hex_palette <- function(in_vec) {
  pal.dt <- OGA_recolor(data.table(t(col2rgb(in_vec))/255))
  return(rgb(pal.dt$new_red,pal.dt$new_green,pal.dt$new_blue))
}

###function to convert from hex to rbY###
hex2rbY <- function(hex) {
  
  temp.dt <- data.table(t(col2rgb(hex))/255)
  temp.dt$hex <- hex
  
  temp.dt$Y <- calc_Y(inter_calc1(temp.dt$red),inter_calc1(temp.dt$green),inter_calc1(temp.dt$blue))
  
  temp.dt$rgb2 <- inter_calc1(temp.dt$red) + inter_calc1(temp.dt$green) + inter_calc1(temp.dt$blue)
  
  temp.dt$rr <- chan_func(inter_calc1(temp.dt$red),temp.dt$rgb2)
  
  temp.dt$bb <- chan_func(inter_calc1(temp.dt$blue),temp.dt$rgb2)
  
  return(temp.dt[,c("hex","rr","bb","Y")])
  
}

###function to convert from rbY to hex###
rbY2hex <- function(rr,bb,Y){
	
	Y1 <- calc_Y(rr,(1 - rr - bb),bb) 
	new_r <- as.numeric(inter_calc2(rr * (Y/Y1)))
	new_g <- as.numeric(inter_calc2((1 - rr - bb) * (Y/Y1)))
	new_b <- as.numeric(inter_calc2(bb * (Y/Y1)))
	
	if(max(new_r,new_g,new_b) > 1) {
		return("out of range")
	} else{
		return(rgb(new_r,new_g,new_b))
	}
	
}
