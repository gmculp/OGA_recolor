
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

##################
###hue 1: azure###
##################
###azure hue's green channel = (b_fac.a*b) + (r_fac.a*r)###
b_fac.a <- 0.4 
r_fac.a <- 0.6
hue.a <- "azure"


###################
###hue 2: orange###
###################
###orange hue's green channel = (b_fac.b*b) + (r_fac.b*r)###
b_fac.b <- 0.6
r_fac.b <- 0.4
hue.b <- "orange"

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

##########################################
###assemble data table of line segments###
##########################################

###sequence of 10###
t <- seq(0, 1, by = 0.1)

seg.dt <- rbindlist(list(
  data.table(
    rr = inter_calc1(t),
    gg = inter_calc1(b_fac.a + (t*r_fac.a)),
    bb = inter_calc1(1),
    hue=rep(hue.a,length(t))),
  data.table(
    rr = inter_calc1(1),
    gg = inter_calc1(r_fac.b + (t*b_fac.b)),
    bb = inter_calc1(t),
    hue=rep(hue.b,length(t)))
))

seg.dt[, x1 := chan_func(rr,(rr + gg + bb))]
seg.dt[, y1 := chan_func(bb,(rr + gg + bb))]
seg.dt[, z1 := calc_Y(rr,gg,bb)]

seg.dt[, x2 := shift(.(x1), type = "lead"), by = hue]
seg.dt[, y2 := shift(.(y1), type = "lead"), by = hue]	
seg.dt[, z2 := shift(.(z1), type = "lead"), by = hue]	


seg.dt <- seg.dt[!(is.na(x2))]	
seg.dt <- seg.dt[order(x1,y1),]	
seg.dt[, seg_id := 	1:nrow(seg.dt)]
seg.dt[, c('rr','gg','bb','hue') := NULL]




#########################
###some more functions###
#########################


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
  
  #out_chan <- pmin(pmax(as.numeric((chan * (Y1/Y2)) ^ (1/gamma)),0),1)
  
  out_chan <- pmin(pmax(as.numeric(inter_calc2(chan * (Y1/Y2))),0),1)
  
  return(out_chan)
}

###http://stackoverflow.com/questions/10600060/how-to-do-cross-join-in-r###
cjdt <- function(a,b){
  cj = CJ(1:nrow(a),1:nrow(b))
  cbind(a[cj[[1]],],b[cj[[2]],])
}


###############################
###main re-coloring function###
###############################

OGA_recolor <- function (test.dt) {
  
  test.dt[ , pix_id := 1:nrow(test.dt)]
  
  test.dt[ , alpha := 1]
  
  setcolorder(test.dt, c("pix_id","red","green","blue","alpha"))
  
  test.dt[ , in_r2 := inter_calc1(red)]
  test.dt[ , in_g2 := inter_calc1(green)]
  test.dt[ , in_b2 := inter_calc1(blue)]
  
  test.dt[ , Y_chan := calc_Y(in_r2,in_g2,in_b2)]
  
  test.dt[ , in_rgb2 := in_r2 + in_g2 + in_b2]
  
  test.dt[ , in_rr2 := chan_func(in_r2,in_rgb2)]
  
  test.dt[ , in_bb2 := chan_func(in_b2,in_rgb2)]
  
  ###many-to-many merge with line segments###
  test.dt <- cjdt(test.dt,seg.dt)
  
  ###remove segments that have lower lightness than pixel###
  test.dt <- test.dt[Y_chan <= z1 | Y_chan <= z2]
  
  ###sort by pixel id and segment id
  setorder(test.dt, pix_id, seg_id)
  
  #########################################################################
  ###find point at which pixel color lightness channel intersects curves###
  #########################################################################
  
  test.dt[ , pt_x := x2 + (((Y_chan - z2)/(z1 - z2))*(x1-x2))]
  test.dt[ , pt_y := y2 + (((Y_chan - z2)/(z1 - z2))*(y1-y2))]
  
  test.dt[ , x1 := round(ifelse(z1 < z2 & Y_chan >= z1 & Y_chan <= z2, pt_x, x1),4)]
  test.dt[ , y1 := round(ifelse(z1 < z2 & Y_chan >= z1 & Y_chan <= z2, pt_y, y1),4)]
  test.dt[ , x2 := round(ifelse(z1 > z2 & Y_chan >= z2 & Y_chan <= z1, pt_x, x2),4)]
  test.dt[ , y2 := round(ifelse(z1 > z2 & Y_chan >= z2 & Y_chan <= z1, pt_y, y2),4)]
  
  test.dt[, c("pt_x","pt_y","z1","z2"):=NULL] 
  
  ##################################################
  ###closest point on chromaticity lines to color###
  ##################################################
  
  test.dt[ , ptu_x := calc_ux(x1,y1,x2,y2,in_rr2,in_bb2)]
  
  test.dt[ , ptu_y := calc_uy(x1,y1,x2,y2,in_rr2,in_bb2)]
  
  test.dt[ , ptu_d := calc_pd(ptu_x,ptu_y,in_rr2,in_bb2)]
  
  test.dt[ , pt1_d := calc_pd(x1,y1,in_rr2,in_bb2)]
  
  test.dt[ , pt2_d := calc_pd(x2,y2,in_rr2,in_bb2)]
  
  test.dt[ , new_rr2 := ifelse(pt1_d <= pt2_d & pt1_d < ptu_d,x1, ifelse(pt2_d < pt1_d & pt2_d < ptu_d,x2,ptu_x))]
  
  test.dt[ , new_rr2 := ifelse(new_rr2 < 0, 0, new_rr2)]
  
  test.dt[ , new_bb2 := ifelse(pt1_d <= pt2_d & pt1_d < ptu_d,y1, ifelse(pt2_d < pt1_d & pt2_d < ptu_d,y2,ptu_y))]
  
  test.dt[ , new_bb2 := ifelse(new_bb2 < 0, 0, new_bb2)]
  
  test.dt[ , new_gg2 :=  1 - new_rr2 - new_bb2]
  
  test.dt[ , pt_d :=  pmin(pt1_d,pt2_d,ptu_d)]
  
  #########################
  ###convert back to RGB###
  #########################
  
  test.dt[ , new_red :=  prop_chan(new_rr2,Y_chan,calc_Y(new_rr2, new_gg2, new_bb2))]
  
  test.dt[ , new_green :=  prop_chan(new_gg2,Y_chan,calc_Y(new_rr2, new_gg2, new_bb2))]
  
  test.dt[ , new_blue :=  prop_chan(new_bb2,Y_chan,calc_Y(new_rr2, new_gg2, new_bb2))]
  
  test.dt <- test.dt[test.dt[, .I[which.min(pt_d)], by = pix_id]$V1]
  
  d_c <- colnames(test.dt)[!colnames(test.dt) %in% c("pix_id","red","green","blue","alpha","new_red","new_green","new_blue")]
  
  test.dt[, (d_c):=NULL]
  
  ###sort by pixel id
  setorder(test.dt, pix_id)
  
  invisible(gc())
  
  return(test.dt)
}



########################################
###wrapper function to recolor images###
########################################	

OGA_recolor_image <- function(img) {
  
  ###convert matrix to data.table###
  ###round colors by 0.02###
  img.dt <- data.table(
    red = round(matrix(img[,,1], ncol=1)/0.02)*0.02,
    green = round(matrix(img[,,2], ncol=1)/0.02)*0.02,
    blue = round(matrix(img[,,3], ncol=1)/0.02)*0.02
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

hex2rbY <- function(in_vec) {
  
  temp.df <- data.frame(t(col2rgb(in_vec))/255)
  temp.df$hex <- in_vec
  temp.df$r2 <- temp.df$red ^ gamma
  temp.df$g2 <- temp.df$green ^ gamma
  temp.df$b2 <- temp.df$blue ^ gamma
  temp.df$Y <- calc_Y(temp.df$r2,temp.df$g2,temp.df$b2)
  
  temp.df$rgb2 <- temp.df$r2 + temp.df$g2 + temp.df$b2
  
  temp.df$rr <- chan_func(temp.df$r2,temp.df$rgb2)
  
  temp.df$bb <- chan_func(temp.df$b2,temp.df$rgb2)
  
  return(temp.df[,c("hex","rr","bb","Y")])
  
}
