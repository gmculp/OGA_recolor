
#######################

exponent <- function(a, pow) (abs(a)^pow)*sign(a)

#######################

calc_sc <- function(this_m){
	#get row containing negative value
	m_neg <- which(round(this_m[],7)<0, arr.ind=T)
	m_row <- as.numeric(m_neg[1,"row"])
	m_col <- as.numeric(m_neg[1,"col"])

	sc_a <- (sum(this_m[m_row,])/sum(this_m[m_row,this_m[m_row,]>0])) + (this_m[m_row,m_col]/sum(this_m[m_row,this_m[m_row,]>0]))

	sc_b <- -(this_m[m_row,m_col]/sum(this_m[m_row,this_m[m_row,]>0]))
	
	return(c(sc_a,sc_b))
}



RGB2XYZ <- matrix(c(40.9568, 35.5041, 17.9167, 21.3389, 70.6743, 7.98680, 1.86297, 11.4620, 91.2367), ncol = 3, byrow = TRUE)

XYZ2LMS <- matrix(c(0.15514, 0.54312, -0.03286, -0.15514, 0.45684, 0.03286, 0, 0, 0.01608), ncol = 3, byrow = TRUE)

RGB2LMS <- XYZ2LMS %*% RGB2XYZ 

L_w <- sum(RGB2LMS[1,])
M_w <- sum(RGB2LMS[2,])
S_w <- sum(RGB2LMS[3,])

L_r <- RGB2LMS[1,1]
M_r <- RGB2LMS[2,1]
S_r <- RGB2LMS[3,1]

#L_c <- sum(RGB2LMS[1,2:3])
#M_c <- sum(RGB2LMS[2,2:3])
#S_c <- sum(RGB2LMS[3,2:3])

L_b <- RGB2LMS[1,3]
M_b <- RGB2LMS[2,3]
S_b <- RGB2LMS[3,3]

a_RG <- (M_w * S_b) - (M_b * S_w)
b_RG <- (S_w * L_b) - (S_b * L_w)
g_RG <- (L_w * M_b) - (L_b * M_w)

a_BG <- (M_w * S_r) - (M_r * S_w)
b_BG <- (S_w * L_r) - (S_r * L_w)
g_BG <- (L_w * M_r) - (L_r * M_w)

LMS2pLMS <- matrix(c(0, (-b_RG/a_RG), (-g_RG/a_RG), 0, 1, 0, 0, 0, 1), ncol = 3, byrow = TRUE)

LMS2dLMS <- matrix(c(1, 0, 0, (-a_RG/b_RG), 0, (-g_RG/b_RG), 0, 0, 1), ncol = 3, byrow = TRUE)

LMS2tLMS <- matrix(c(1, 0, 0, 0, 1, 0, (-a_BG/g_BG), (-b_BG/g_BG), 0), ncol = 3, byrow = TRUE)





RGB2pRGB <- solve(RGB2LMS) %*% LMS2pLMS %*% RGB2LMS

RGB2dRGB <- solve(RGB2LMS) %*% LMS2dLMS %*% RGB2LMS

RGB2tRGB <- solve(RGB2LMS) %*% LMS2tLMS %*% RGB2LMS




##########################
###lightness correction###
##########################

RGB2XYZ <- matrix(c(0.4124,0.3576,0.1805,0.2126,0.7152,0.0722,0.0193,0.1192,0.9505), ncol = 3, byrow = TRUE)

hex1 <- c("#FF00FF","#FF0000")

hex2 <- c("#FF0000","#FF00FF")

correct_Y <- function(hex1,hex2) {
	xyz1 <- as.data.frame(t(RGB2XYZ %*% ((col2rgb(hex1)/255) ^ 2.4)))
	xyz2 <- as.data.frame(t(RGB2XYZ %*% ((col2rgb(hex2)/255) ^ 2.4)))

	xyz2_sum <- xyz2$V1 + xyz2$V2 + xyz2$V3

	xx2 <-  ifelse(xyz2_sum>0, xyz2$V1/xyz2_sum, 1/3)
	yy2 <-  ifelse(xyz2_sum>0, xyz2$V2/xyz2_sum, 1/3)

	new.X <- xx2 * ( xyz1$V2 / yy2 )
	new.Y <- xyz1$V2
	new.Z <- (1 - xx2 - yy2) * ( xyz1$V2 / yy2 )
	
	#return(rgb(t(round(exponent((solve(RGB2XYZ) %*% matrix(c(new.X,new.Y,new.Z), nrow = 3, byrow=TRUE)),(1/2.4)),2))))
	
	return(rgb(pmax(pmin(t(round(exponent((solve(RGB2XYZ) %*% matrix(c(new.X,new.Y,new.Z), nrow = 3, byrow=TRUE)),(1/2.4)),2)),1),0)))
}





CVD_p <- function(in_hex) {
	out_hex <- rgb(round(t(exponent((RGB2pRGB %*% ((((col2rgb(in_hex)/255) ^ 2.2)*calc_sc(RGB2pRGB)[1]) + calc_sc(RGB2pRGB)[2])),(1/2.2))),3))
	#return(correct_Y(in_hex,out_hex))
	return(out_hex)
}

CVD_d <- function(in_hex) {
	out_hex <- rgb(round(t(exponent((RGB2dRGB %*% ((((col2rgb(in_hex)/255) ^ 2.2)*calc_sc(RGB2dRGB)[1]) + calc_sc(RGB2dRGB)[2])),(1/2.2))),3))
	#return(correct_Y(in_hex,out_hex))
	return(out_hex)
}

CVD_t <- function(in_hex) {
	out_hex <- rgb(round(t(exponent((RGB2tRGB %*% ((((col2rgb(in_hex)/255) ^ 2.2)*calc_sc(RGB2tRGB)[1]) + calc_sc(RGB2tRGB)[2])),(1/2.2))),3))
	return(correct_Y(in_hex,out_hex))
	#return(out_hex)
}


sim_img <- function(img,CVD_type){

	#cat(CVD_type)
	if (CVD_type == "tritan"){
		CVD_fun <- CVD_t
	} else if (CVD_type == "protan"){
		CVD_fun <- CVD_p
	} else {
		CVD_fun <- CVD_d
	}

	img.dt <- data.table(
			red = round(matrix(img[,,1], ncol=1),2),
			green = round(matrix(img[,,2], ncol=1),2),
			blue = round(matrix(img[,,3], ncol=1),2),
			alpha = 1
	)
		
	colnames(img.dt) <- gsub(".V1","",colnames(img.dt))

	img.dt$pix_id <- 1:nrow(img.dt)

	img.dt[ , in_hex := rgb(red,green,blue, maxColorValue = 1)]
	img.dt[ , out_hex := CVD_fun(in_hex)]
	img.dt[ , new_red := t(col2rgb(out_hex))[,"red"]/255]
	img.dt[ , new_green := t(col2rgb(out_hex))[,"green"]/255]
	img.dt[ , new_blue := t(col2rgb(out_hex))[,"blue"]/255]


	img.new = array(dim=dim(img))
	img.new[,,1] = matrix(as.numeric(img.dt$new_red), nrow=dim(img)[1])
	img.new[,,2] = matrix(as.numeric(img.dt$new_green), nrow=dim(img)[1])
	img.new[,,3] = matrix(as.numeric(img.dt$new_blue), nrow=dim(img)[1])

	if (dim(img)[3]==4) {
		#cat("4-dim image")
		img.new[,,4] = matrix(as.numeric(img.dt$alpha), nrow=dim(img)[1])
	}
	return(img.new)
}
