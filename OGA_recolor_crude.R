
###DISABLE SCIENTIFIC NOTATION###
options(scipen = 999)

############################
###load required packages###
############################
suppressMessages(library(data.table))


####################
###some constants###
####################

###rounding factor###
pix_rnd <- 0.01

###gamma value###
gamma <- 2.4

###sRGB to XYZ matrix Y column values for ref. white D65###
Y_r <- 0.2126

Y_g <- 0.7152

Y_b <- 0.0722

###white point###
w_pt <- 1/3

##################
###hue 1: azure###
##################
rr.a <- 0
bb.a <- 1/1.4

###################
###hue 2: orange###
###################
rr.o <- 1/1.4
bb.o <- 0


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
  return(ifelse(chan > 0, chan/tot, ifelse(tot==0,w_pt,0)))
}

prop_chan <- function(chan,Y1,Y2){
  
  out_chan <- pmin(pmax(as.numeric(inter_calc2(inter_calc1(chan) * (Y1/Y2))),0),1)
  
  return(out_chan)
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



###http://stackoverflow.com/questions/10600060/how-to-do-cross-join-in-r###
cjdt <- function(a,b){
  cj = CJ(1:nrow(a),1:nrow(b))
  cbind(a[cj[[1]],],b[cj[[2]],])
}

##########################################
###assemble data table of line segments###
##########################################

make_segments <- function(rr,bb){

	gg <- 1 - rr - bb

	###sequence of 10###
	t <- seq(0, 1, by = 0.1)

	val.v <- c(rr,gg,bb)
	name.v <- c('r','g','b')

	###check for corner colors###
	min.c <- name.v[which( round(val.v,2) == min(round(val.v,2)))]
	max.c <- name.v[which( round(val.v,2) == max(round(val.v,2)))]
	max.c <- max.c[!(max.c %in% min.c)] #for white
	mid.c <- name.v[!(name.v %in% min.c) & !(name.v %in% max.c)]

	rgb.list <- append(setNames(lapply(min.c,function(x) t),min.c), setNames(lapply(max.c,function(x) rep(1,length(t))),max.c))

	if(length(mid.c)==1){
		f1 <- (val.v[name.v == mid.c] - val.v[name.v == max.c])/(val.v[name.v == min.c] - val.v[name.v == max.c])
		f2 <- 1 - f1
		rgb.list <- append(rgb.list,setNames(lapply(mid.c,function(x) (t*f1)+f2),mid.c))
	} 
		
	rgb.dt <- as.data.table(rgb.list)	
	setcolorder(rgb.dt, c("r","g","b"))	
		
	rgb.dt[, rr.2 := r/(r+g+b)]
	rgb.dt[, bb.2 := b/(r+g+b)]
	rgb.dt[, Y.2 := calc_Y(inter_calc1(r),inter_calc1(g),inter_calc1(b))]
	
	rgb.dt[, c('r','g','b') := NULL]
	
	rgb.dt[, rr.1 := shift(.(rr.2), type = "lag")]
	rgb.dt[, bb.1 := shift(.(bb.2), type = "lag")]	
	rgb.dt[, Y.1 := shift(.(Y.2), type = "lag")]	
	
	###add line for below saturated color###
	rgb.dt[, rr.1 := ifelse(is.na(rr.1),rr.2,rr.1)]
	rgb.dt[, bb.1 := ifelse(is.na(bb.1),bb.2,bb.1)]
	rgb.dt[, Y.1 := ifelse(is.na(Y.1),0,Y.1)]

	return(rgb.dt)
}


###############################
###main re-coloring function###
###############################

OGA_recolor <- function (test.dt) {
  
  test.dt[ , col_id := 1:nrow(test.dt)]
  
  test.dt[ , alpha := 1]
  
  #setcolorder(test.dt, c("col_id","red","green","blue","alpha"))
  
  test.dt[ ,in.rr := chan_func(red,(red + green + blue))]
  
  test.dt[ ,in.bb := chan_func(blue,(red + green + blue))]
  
  test.dt[ ,in.Y := calc_Y(inter_calc1(red),inter_calc1(green),inter_calc1(blue))]
  
  ###build look-up of line segments###
  seg.dt <- rbindlist(list(make_segments(rr.a,bb.a),make_segments(rr.o,bb.o)))
  
  ###many-to-many merge with line segments###
  test.dt <- cjdt(test.dt,seg.dt)
  
  ###remove out-of-range segments###	
  test.dt <- test.dt[in.Y >= Y.1 & in.Y <= Y.2]
  
  #########################################################################
  ###find point at which pixel color lightness channel intersects curves###
  #########################################################################
  
  test.dt[ , rr.1 := rr.2 + (((in.Y - Y.2)/(Y.1 - Y.2))*(rr.1-rr.2))]
  test.dt[ , bb.1 := bb.2 + (((in.Y - Y.2)/(Y.1 - Y.2))*(bb.1-bb.2))]
  
  #################################
  ###extend lines to white point###
  #################################
  
  test.dt[, rr.2 := w_pt]
  test.dt[, bb.2 := w_pt]
  
  ###remove columns that are no longer necessary###
  test.dt[, c("Y.1","Y.2"):=NULL] 
  
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
  
  test.dt[ , temp.Y := calc_Y(inter_calc1(new.rr),inter_calc1(new.gg),inter_calc1(new.bb))]
 
  test.dt[ , new_red := prop_chan(new.rr,in.Y,temp.Y)]
  test.dt[ , new_green := prop_chan(new.gg,in.Y,temp.Y)]
  test.dt[ , new_blue := prop_chan(new.bb,in.Y,temp.Y)]
  
  ###remove unnecessary columns###
  d_c <- colnames(test.dt)[!colnames(test.dt) %in% c("col_id","red","green","blue","alpha","new_red","new_green","new_blue")]
  test.dt[, (d_c):=NULL]
  
  ###sort by pixel id###
  setorder(test.dt, col_id)
  
  ###flush memory###
  invisible(gc())
  
  return(test.dt)
}



########################################
###wrapper function to recolor images###
########################################	

OGA_recolor_image <- function(img) {
  
  ###convert matrix to data.table###
  ###round colors by pix_rnd###
  img.dt <- data.table(
    red = round(matrix(img[,,1], ncol=1)/pix_rnd)*pix_rnd,
    green = round(matrix(img[,,2], ncol=1)/pix_rnd)*pix_rnd,
    blue = round(matrix(img[,,3], ncol=1)/pix_rnd)*pix_rnd
  )
  
  colnames(img.dt) <- gsub(".V1","",colnames(img.dt))
  
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

OGA_recolor_hex_palette <- function(in_vec) {
  pal.dt <- OGA_recolor(data.table(t(col2rgb(in_vec))/255))
  return(rgb(pal.dt$new_red,pal.dt$new_green,pal.dt$new_blue))
}

hex2rbY <- function(hex) {
  
  temp.dt <- data.table(t(col2rgb(hex))/255)
  temp.dt$hex <- hex
  
  temp.dt$Y <- calc_Y(inter_calc1(temp.dt$red),inter_calc1(temp.dt$green),inter_calc1(temp.dt$blue))
  
  temp.dt$rgb2 <- temp.dt$red + temp.dt$green + temp.dt$blue
  
  temp.dt$rr <- chan_func(temp.dt$red,temp.dt$rgb2)
  
  temp.dt$bb <- chan_func(temp.dt$blue,temp.dt$rgb2)
  
  return(temp.dt[,c("hex","rr","bb","Y")])
  
}


rbY2hex <- function(rr,bb,Y,warp=FALSE){
	
	Y1 <- calc_Y(inter_calc1(rr),inter_calc1(1 - rr - bb),inter_calc1(bb)) #z1
	
	new_r <- as.numeric(inter_calc2(inter_calc1(rr) * (Y/Y1)))
	new_g <- as.numeric(inter_calc2(inter_calc1(1 - rr - bb) * (Y/Y1)))
	new_b <- as.numeric(inter_calc2(inter_calc1(bb) * (Y/Y1)))
	
	if(max(new_r,new_g,new_b) > 1) {
	
		if(warp==TRUE){
			t.dt <- make_segments(in.rr,in.bb)
			t.dt <- t.dt[Y >= Y.1 & Y <= Y.2]
			t.dt <- t.dt[1,]
			
			rr.pt <- t.dt$rr.2 + (((Y - t.dt$Y.2)/(t.dt$Y.1 - t.dt$Y.2))*(t.dt$rr.1-t.dt$rr.2))
			bb.pt <- t.dt$bb.2 + (((Y - t.dt$Y.2)/(t.dt$Y.1 - t.dt$Y.2))*(t.dt$bb.1-t.dt$bb.2))
			
			return(rbY2hex(rr.pt,bb.pt,Y))
			
		} else {
			return("out of range")
		}
	} else{
		return(rgb(new_r,new_g,new_b))
	}
}
