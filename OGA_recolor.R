############################
###load required packages###
############################
library("png")	
library("jpeg")		
library("grid")
library("gridExtra")
library(data.table)


####################
###some constants###
####################

gamma <- 2.4

Y_r <- 0.2126

Y_g <- 0.7152

Y_b <- 0.0722

w_pt <- 1/3


###################################
###generate line segment look-up###
###################################

##################################
###old option for opponent pair###
##################################

###b_fac <- 0.6
###r_fac <- 0.4 

###t <- seq(0, 1, by = 0.05)

###x1 <- (t^gamma)/((t^gamma) + (((r_fac*t)+(b_fac*(1-t)))^gamma) + ((1-t)^gamma)) 
###y1 <- ((1-t)^gamma)/((t^gamma) + (((r_fac*t)+(b_fac*(1-t)))^gamma) + ((1-t)^gamma)) 

###seg.df <- data.frame(x1=x1[1:length(x1)-1],y1=y1[1:length(y1)-1],x2=x1[2:length(x1)],y2=y1[2:length(y1)])

###seg.df <- seg.df[order(seg.df$x1,seg.df$y1),]

###seg.df$seg_id <- 1:nrow(seg.df)

###seg.df$m <- (seg.df$y2-seg.df$y1)/(seg.df$x2-seg.df$x1)

###seg.df$b <- seg.df$y2 - (seg.df$m * seg.df$x2)

###seg.df$hue <- ifelse(round(seg.df$x1,4) > 0.3333 | round(seg.df$x2,4) > 0.3333, "orange", "azure")

######################################
###new option for non opponent pair###
######################################

t <- seq(0, 1, by = 0.1)

###################
###hue 1: bluish###
###################
###bluish hue's green channel = (b_fac.a*b) + (r_fac.a*r)###
b_fac.a <- 0.4 
r_fac.a <- 0.6
hue.a <- "azure"
x.a <- (t ^ gamma) / ((t ^ gamma) + ((b_fac.a +(t*r_fac.a)) ^ gamma) + 1)
y.a <- 1/((t ^ gamma) + ((b_fac.a + (t*r_fac.a)) ^ gamma) + 1)
z.a <- Y_b + ((t ^ gamma) * Y_r) + (((b_fac.a +(t*r_fac.a)) ^ gamma) * Y_g) 


####################
###hue 2: reddish###
####################
###reddish hue's green channel = (b_fac.b*b) + (r_fac.b*r)###
b_fac.b <- 0.6
r_fac.b <- 0.4
hue.b <- "orange"
x.b <- rev(1 / ((t ^ gamma) + (((t*b_fac.b)+ r_fac.b) ^ gamma) + 1))
y.b <- rev((t ^ gamma)/((t ^ gamma) + (((t*b_fac.b) + r_fac.b) ^ gamma) + 1))
z.b <- rev(((t ^ gamma) * Y_b) + Y_r + ((((t*b_fac.b)+ r_fac.b) ^ gamma) * Y_g))

#########################
###assemble data frame###
#########################
seg.df <- rbind(data.frame(x1=x.a[1:length(t)-1],
y1=y.a[1:length(t)-1],
z1=z.a[1:length(t)-1],
x2=x.a[2:length(t)],
y2=y.a[2:length(t)],
z2=z.a[2:length(t)],
hue=rep(hue.a,length(t)-1)),
data.frame(x1=x.b[1:length(t)-1],
y1=y.b[1:length(t)-1],
z1=z.b[1:length(t)-1],
x2=x.b[2:length(t)],
y2=y.b[2:length(t)],
z2=z.b[2:length(t)],
hue=rep(hue.b,length(t)-1)))

seg.df <- seg.df[order(seg.df$x1,seg.df$y1),]

seg.df$seg_id <- 1:nrow(seg.df)

seg.df$m <- (seg.df$y2-seg.df$y1)/(seg.df$x2-seg.df$x1)

seg.df$b <- seg.df$y2 - (seg.df$m * seg.df$x2)



####################
###some functions###
####################

calc_Y <- function(R,G,B) {
	(R * Y_r) + (G * Y_g) + (B * Y_b)
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
	out_chan <- as.numeric((chan * (Y1/Y2)) ^ (1/gamma))
	out_chan <- ifelse(out_chan > 1, 1, out_chan)
	out_chan <- ifelse(out_chan < 0, 0, out_chan)
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

	# Start the clock!
	ptm <- proc.time()
	
	#test.dt <- img.dt
	
	test.dt$pix_id <- 1:nrow(test.dt)
	
	test.dt[ , alpha := 1]
	
	setcolorder(test.dt, c("pix_id","red","green","blue","alpha"))
	
	test.dt[ , in_r2 := red ^ gamma]
	test.dt[ , in_g2 := green ^ gamma]
	test.dt[ , in_b2 := blue ^ gamma]
	
	test.dt[ , Y_chan := calc_Y(in_r2,in_g2,in_b2)]
	
	test.dt[ , in_rgb2 := in_r2 + in_g2 + in_b2]
	
	test.dt[ , in_rr2 := chan_func(in_r2,in_rgb2)]
	
	test.dt[ , in_bb2 := chan_func(in_b2,in_rgb2)]
	
	
	
	###many-to-many merge with line segments###
	test.dt <- cjdt(test.dt,data.table(seg.df))
	
	###remove segments that have lower lightness than pixel###
	test.dt <- test.dt[Y_chan <= z1 | Y_chan <= z2]
	
	###sort by pixel id and segment id
	test.dt <- setorder(test.dt, pix_id, seg_id)
	
	################
	###azure line###
	################
	
	test.dt[ , pt_azr_y := ((m * Y_g) + (b * Y_g) - (b * Y_r))/((Y_chan * m) + (m * Y_g) - (m * Y_b) + Y_g - Y_r) ]
	
	test.dt[ , pt_azr_y := ifelse(pt_azr_y > y1, y1, pt_azr_y) ]

	test.dt[ , pt_azr_x := (pt_azr_y - b)/m]

	#################
	###orange line###
	#################
	
	test.dt[ , pt_org_x := ((b * Y_b) + Y_g - (b * Y_g))/ ((m * Y_g) - (m * Y_b) + Y_chan + Y_g - Y_r)]
	
	test.dt[ , pt_org_x := ifelse(pt_org_x > x2, x2, pt_org_x)]

	test.dt[ , pt_org_y := (pt_org_x * m) + b]
	
	######################################
	###assign point values based on hue###
	######################################
	
	test.dt[ , x1 := round(ifelse(hue=="azure",pt_azr_x,x1),4)]
	
	test.dt[ , y1 := round(ifelse(hue=="azure",pt_azr_y,y1),4)]
	
	test.dt[ , x2 := round(ifelse(hue !="azure",pt_org_x,x2),4)]
	
	test.dt[ , y2 := round(ifelse(hue !="azure",pt_org_y,y2),4)]
	
	test.dt[, c("pt_azr_y","pt_org_y","pt_azr_x","pt_org_x"):=NULL] 
	
	###remove records where point is outside range###
	test.dt <- test.dt[x1 <= x2 & y1 >= y2]
	
	
	#####################################
	###closest point on lines to color###
	#####################################
	
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
	
	test.dt[ , new_red :=  prop_chan(new_rr2,Y_chan,calc_Y(new_rr2, new_gg2, new_bb2))]
	
	test.dt[ , new_green :=  prop_chan(new_gg2,Y_chan,calc_Y(new_rr2, new_gg2, new_bb2))]

	test.dt[ , new_blue :=  prop_chan(new_bb2,Y_chan,calc_Y(new_rr2, new_gg2, new_bb2))]
	
	test.dt <- test.dt[test.dt[, .I[which.min(pt_d)], by = pix_id]$V1]
	
	d_c <- colnames(test.dt)[!colnames(test.dt) %in% c("pix_id","red","green","blue","alpha","new_red","new_green","new_blue")]
	
	test.dt <- test.dt[, (d_c):=NULL]
	
	###sort by pixel id
	test.dt <- setorder(test.dt, pix_id)
	
	invisible(gc())
	
	#cat(paste0("\n",round(as.numeric((proc.time() - ptm)[3]),2)," seconds\n"))
	
	return(test.dt)
}



########################################
###wrapper function to recolor images###
########################################	

OGA_recolor_image <- function(img) {

	#img <- readPNG("C:/Users/gmculp/Google Drive/various/gb_budgies.png")
	
	#img <- readPNG("J:/candies_color3.png")
	
	# reshape image into a data table
	#img.dt <- data.table(
	#	red = matrix(img[,,1], ncol=1),
	#	green = matrix(img[,,2], ncol=1),
	#	blue = matrix(img[,,3], ncol=1)
	#)
	
	#img.dt <- data.table(
	#	red = round(matrix(img[,,1], ncol=1),2),
	#	green = round(matrix(img[,,2], ncol=1),2),
	#	blue = round(matrix(img[,,3], ncol=1),2)
	#)
	
	img.dt <- data.table(
		red = round(matrix(img[,,1], ncol=1)/0.02)*0.02,
		green = round(matrix(img[,,2], ncol=1)/0.02)*0.02,
		blue = round(matrix(img[,,3], ncol=1)/0.02)*0.02
	)
	
	colnames(img.dt) <- gsub(".V1","",colnames(img.dt))
	
	#img.dt2 <- OGA_recolor(img.dt)
	img.dt2A <- OGA_recolor(unique(img.dt))
	img.dt$pix.order <- 1:nrow(img.dt)
	
	#ptm <- proc.time()
	###fast merge with data.table###
	img.dt2 <- merge(img.dt,img.dt2A,by=c("red","green","blue"))
	rm(img.dt2A)
	
	img.dt2 <- setorder(img.dt2, pix.order)
	#cat(paste0("\n",round(as.numeric((proc.time() - ptm)[3]),2)," seconds\n"))
	

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
	
	return(img.new)
	invisible(gc())
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
