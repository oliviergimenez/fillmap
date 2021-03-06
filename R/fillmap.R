#' Create maps with SpatialPolygons objects
#'
#' This function offers several options for creating plots with SpatialPolygons objects. This can be a single plot or multiple plots. These multi-plots can have the same legend or individual legends. Additionally, the layout call within this function offers many options for improving the resulting figures.
#' @param map The SpatialPolygon object.
#' @param figtitle The title to be displayed on the map(s). Can be of length 1 or more.
#' @param y The vector or matrix object to be mapped.
#' @param n.col The number of legend levels.
#' @param bk The break specifications defaults to "e" for equal. Other options are "c" to specify your own breaks and "q" for quantile breaks.
#' @param cuts The cuts to be specified when bk="c".
#' @param legendtxt The text to be displayed in the legend (must be NA (which removes the text altogether) or equal length n.col+1). Defaults to "" which creates text based on the cuts. 
#' @param leg.loc The loation for the legend defaults to "bottomright". Can also be given as a list of x,y coordinates
#' @param leg.cex Specifies the font size in the legend. Defaults to 1.5.
#' @param main.cex Specifies the font size for the title. Defaults to 1.5.
#' @param main.line Specifies the line for the title of the plot to appear on. Defaults to -2.
#' @param leg.horiz Indicates if the legend should be horizontal. Defaults to FALSE.
#' @param map.lty Specifies the line type for the plot. Defaults to 1.
#' @keywords disease mapping 
#' @export
#' @examples
#' Sr1 = Polygon(cbind(c(0,0,1,1,0),c(0,1,1,0,0)))
#' Sr2 = Polygon(cbind(c(0,1,1,0,0),c(0,0,-1,-1,0)))
#' Sr3 = Polygon(cbind(c(0,-1,-1,0,0),c(0,0,1,1,0)))
#' Sr4 = Polygon(cbind(c(0,0,-1,-1,0),c(0,-1,-1,0,0)))
#' Sr5 = Polygon(cbind(c(1,1,2,2,1),c(0,1,1,0,0)))
#' Sr6 = Polygon(cbind(c(0,2,2,1,1),c(0,0,-1,-1,0)))
#' Sr7 = Polygon(cbind(c(-1,-1,0,0,-1),c(1,2,2,1,1)))
#' Sr8 = Polygon(cbind(c(-1,-2,-2,-1,-1),c(1,1,2,2,1)))
#' Sr9 = Polygon(cbind(c(0,0,1,1,0),c(1,2,2,1,1)))
#' Sr10 = Polygon(cbind(c(-2,-2,-1,-1,-2),c(-2,-1,-1,-2,-2)))
#' Sr11 = Polygon(cbind(c(-2,-3,-3,-2,-2),c(-2,-2,-1,-1,-2)))
#' Sr12 = Polygon(cbind(c(-1,-1,0,0,-1),c(-2,-1,-1,-2,-2)))
#' Srs1 = Polygons(list(Sr1), "s1")
#' Srs2 = Polygons(list(Sr2), "s2")
#' Srs3 = Polygons(list(Sr3), "s3")
#' Srs4 = Polygons(list(Sr4), "s4")
#' Srs5 = Polygons(list(Sr5), "s5")
#' Srs6 = Polygons(list(Sr6), "s6")
#' Srs7 = Polygons(list(Sr7), "s7")
#' Srs8 = Polygons(list(Sr8), "s8")
#' Srs9 = Polygons(list(Sr9), "s9")
#' Srs10 = Polygons(list(Sr10), "s10")
#' Srs11 = Polygons(list(Sr11), "s11")
#' Srs12 = Polygons(list(Sr12), "s12")
#' SpP = SpatialPolygons(list(Srs1,Srs2,Srs3,Srs4,Srs5,Srs6,Srs7,Srs8,Srs9,
#' 	Srs10,Srs11,Srs12), 1:12)
#' 
#' v=as.matrix(runif(12,0,1))
#' fillmap(SpP,"Example 1",v,n.col=3)
#' fillmap(SpP,"Example 2",v,n.col=3,bk="q",leg.loc="left",
#' 	legendtxt=c("low","middle","high"))
#' 


fillmap<-function(map, figtitle, y , n.col, bk="e", cuts,legendtxt="",
	leg.loc="bottomright",leg.cex=1.5,main.cex=1.5,main.line=-2,leg.horiz=F,map.lty=1){

if(bk=="q"){if (min(y)<min(y) | max(y)>max(y)){
	print("The minimum or maximum	values of y fall outside of those for y")
	} else {p <- seq(0,1, length=n.col+1)
      br <- round(quantile(y, probs=p),2)}}
if(bk=="e"){if (min(y)<min(y) | max(y)>max(y)){
	print("The minimum or maximum values of y fall outside of those for y")
	} else {br <- round(seq(min(y), max(y), length=n.col+1),6)}}
if(bk=="c"){if (length(cuts)!= (n.col+1)) {cat("Cut off and color categories 
	do not match. ", "\n")
      break}  else {br <- cuts}  }

# 0: dark 1: light light Current shading ranges from darkest to light gray white (to distinguish with lakes).
shading<-gray(rev(0:(n.col-1)/(n.col-1)))
#shading<-hsv(.6,alpha=0:(n.col-1)/(n.col-1))
y.grp<-findInterval(y, vec=br, rightmost.closed = TRUE, all.inside = TRUE)
y.shad<-shading[y.grp]

plot(map,col=y.shad,axes=F, lty=map.lty)
title(main=figtitle,cex.main=main.cex,line=main.line) 

br<-round(br, 2)

if (is.na(legendtxt[1])){print("No legend specifed")
} else if (legendtxt[1]==""){
 	leg.txt<-paste("[",br[n.col],",",br[n.col+1],"]",sep="")
	for(j in (n.col-1):1){ 
		leg.txt<-append(leg.txt,paste("[",br[j],",",br[j+1],")",sep="")) }
		leg.txt<-rev(leg.txt)
	legend(leg.loc,legend=leg.txt,fill=shading,cex=leg.cex,ncol=1,bty="n",
		horiz=leg.horiz)
} else if (length(legendtxt) != n.col){cat("Length of lengendtxt must equal 
		n.col", "\n")
            break
} else {leg.txt<-legendtxt
	legend(leg.loc,legend=leg.txt,fill=shading,cex=leg.cex,ncol=1,bty="n",
		horiz=leg.horiz)}
}

