#' Assesses spatial random effects, such as those mapped with fillmap.
#'
#' This function offers several options for secondary assessments of spatial random effects.
#' @param u The random effect to be assessed.
#' @param v The second random effect to be assessed such that u + v.
#' @param x The the risk factor for exploration.
#' @param bk The categorization specification for the random effect. Defaults to "quart." Other options are: "quint", "tert", "med", "trend", "user".
#' @param user The user specified categoriation of the random effect.
#' @param graph The INLA graph file neccessary for the secondary assessment.
#' @keywords disease mapping secondary assessment
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
#' nb.r <- poly2nb(SpP, queen=F)    
#' mat <- nb2mat(nb.r, style="B",zero.policy=TRUE) # mat is the 0/1 adjacency matrix
#' n.site <- dim(mat)[1]          # n.site: number of areas
#' n.edge <- sum(mat)/2           # n.edge: number of unique pairs
#' nbToInlaGraph(nb.r,graphFile="~\\INLAgraph.dat")
#' 
#' v=as.matrix(runif(12,0,1))
#' x=as.matrix(runif(12,0,1))
#' library(INLA)
#' assessmap(u=v,x=x,graph="~\\INLAgraph.dat")
#' 


assessmap=function(u,v=rep(0,length(u)),x,bk="quart",user,graph){
 
  if (bk=="quart"){
  G=ifelse(u+v<=sort(u+v)[round(length(u)/4)],1,
      ifelse(u+v>sort(u+v)[round(length(u)*3/4)],4,
        ifelse((u+v>sort(u+v)[round(length(u)/4)])&(u+v<=sort(u+v)[round(length(u)/2)]),2,3)))
  }
  else if (bk=="quint"){
    G=ifelse(u+v<=sort(u+v)[round(length(u)/5)],1,
        ifelse(u+v>sort(u+v)[round(length(u)*4/5)],5,
          ifelse((u+v>sort(u+v)[round(length(u)/5)])&(u+v<=sort(u+v)[round(length(u)*2/5)]),2,
            ifelse((u+v>sort(u+v)[round(length(u)*2/5)])&(u+v<=sort(u+v)[round(length(u)*3/5)]),3,4))))
  }
  else if (bk=="tert"){
    G=ifelse(u+v<=sort(u+v)[round(length(u)/3)],1,
             ifelse(u+v>sort(u+v)[round(length(u)*2/3)],3,2))
  }
  else if (bk=="med"){
    G=ifelse(u+v<=sort(u+v)[round(length(u)/2)],1,2)
  }
  else if (bk=="trend"){
    G=u+v
  }
  else if (bk=="user"){
    if (length(user)!=length(u)){
      print("Length of user must be the same as length of u. Quartiles used instead.")
      G=ifelse(u+v<=sort(u+v)[round(length(u)/4)],1,
               ifelse(u+v>sort(u+v)[round(length(u)*3/4)],4,
                      ifelse((u+v>sort(u+v)[round(length(u)/4)])&(u+v<=sort(u+v)[round(length(u)/2)]),2,3)))
    }
    else {G=user}
  }
  else {
    print("User must specify bk= 'med', 'tert', 'quart', 'quint', 'trend', or 'user' for comaprison levels. Quartiles used instead.")
    G=ifelse(u+v<=sort(u+v)[round(length(u)/4)],1,
             ifelse(u+v>sort(u+v)[round(length(u)*3/4)],4,
                    ifelse((u+v>sort(u+v)[round(length(u)/4)])&(u+v<=sort(u+v)[round(length(u)/2)]),2,3)))
  }
 
  if (bk=="trend"){
    id=id2=1:length(u)
    formG=y~1+G+f(id,model="iid",param=c(2,1))+f(id2,model="besag",graph=graph,param=c(2,1),scale.model=TRUE)
    resG<-inla(formG,data=list(y=x,G=G),
               control.fixed     = list(prec.intercept = 1, prec = 1),
               control.inla      = list(strategy = "laplace"))
    result=resG$summary.fixed[-1,c(1,3,5)]
    row.names(result)="Testing the linear trend"
    return(result)
  }
  else {
  id=id2=1:length(u)
  formG=y~1+as.factor(G)+f(id,model="iid",param=c(2,1))+f(id2,model="besag",graph=graph,param=c(2,1),scale.model=TRUE)
  resG<-inla(formG,data=list(y=x,G=G),
           control.fixed     = list(prec.intercept = 1, prec = 1),
           control.inla      = list(strategy = "laplace"))
  result=resG$summary.fixed[-1,c(1,3,5)]
  for (i in 1:(length(table(G))-1)){
    row.names(result)[i]=paste(levels(as.factor(G))[1]," vs. ",levels(as.factor(G))[i+1])
  }
  return(result)
  }
}