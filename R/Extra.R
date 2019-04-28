#
#' @useDynLib AbdominalEITextra
#' @export
systembfort <- function(numnp, mband, b, a, mxbnd, mxnode) {
  out <- .Fortran("systembf", numnp = as.integer(numnp), mband = as.integer(mband), 
                  b = as.double(as.vector(b)), lb = as.integer(length(b)), a = as.double(as.vector(as.matrix(a))), 
                  la = as.integer(length(as.vector(a))), mxbnd = as.integer(mxbnd), mxnode = as.integer(mxnode))
  rb <- as.vector(out$b)
  return(rb)
}
#
#' hoge.
#' @useDynLib AbdominalEITextra
#' @export
systemafort <- function(numnp, mband, a, mxbnd, mxnode) {
  nr <- nrow(a)
  nc <- ncol(a)
  out <- .Fortran("systemaf", numnp = as.integer(numnp), mband = as.integer(mband), 
                  a = as.double(as.vector(as.matrix(a))), la = as.integer(length(as.vector(a))), 
                  mabnd = as.integer(mxbnd), mxnode = as.integer(mxnode))
  ra <- matrix(as.vector(out$a), nrow = nr, ncol = nc)
  return(ra)
} 
#
#############################################
################################################################################################ 
#' FUNCTION fieldfort.
#' @useDynLib AbdominalEITextra
#' @export
fieldfort <- function(nex, nez, cex, bb, nnodez, nodex, nodez, mxnode, mze, mxe) {
  exyz <- array(0, c(3, mze, mxe))
  out <- .Fortran("fieldf", nex = as.integer(nex), nez = as.integer(nez), lcex = as.integer(length(cex)), 
                  cex = as.double(as.vector(cex)), lbb = as.integer(length(bb)), bb = as.double(as.vector(bb)), 
                  nnodez = as.integer(nnodez), lnodex = as.integer(length(nodex)), nodex = as.integer(as.vector(nodex)), 
                  lnodez = as.integer(length(nodez)), nodez = as.integer(as.vector(nodez)), 
                  lexyz = as.integer(length(exyz)), exyz = as.double(as.vector(exyz)), mxe = as.integer(mxe), 
                  mze = as.integer(mze), mxnode = as.integer(mxnode))
  exyzarray <- array(as.vector(out$exyz), c(3, mze, mxe))
  return(exyzarray)
}
################################################################################################ 
#' FUNCTION field.
#' @useDynLib AbdominalEITextra
#' @export
field <- function(nex, nez, cex, bb, nnodez, nodex, nodez) {
  exyzl <- array(0, dim = c(3, nez, nex))
  for (n in 1:nex) {
    for (l in 1:nez) {
      for (j in 1:3) {
        exyzl[j, l, n] <- -sum(cex[j, 1:4, l, n] * (bb[nnodez * (nodex[1:4, 
                                                                       l, n] - 1) + nodez[1:4, l, n]]))
      }
    }
  }
  return(exyzl)
}
################################################################################################ FUNCTION invcurfort
invcurfort <- function(ncur, nex, nez, vol, ezfac, exyz, exyz0, sigma, pod, mxe, 
                       mze, mxnode, jcur, kvol) {
  sensemattemp <- array(0, c(mxe, jcur))
  out <- .Fortran("invcur", ncur = as.integer(ncur), nex = as.integer(nex), nez = as.integer(nez), 
                  vol = as.double(as.vector(vol)), ezfac = as.double(as.vector(ezfac)), as.double(as.vector(exyz)), 
                  exyz0 = as.double(as.vector(exyz0)), sigma = as.double(as.vector(sigma)), 
                  pod = as.double(as.vector(pod)), sensemat = as.double(as.vector(sensemattemp)), 
                  mxe = as.integer(mxe), mze = as.integer(mze), mxnode = as.integer(mxnode), 
                  jcur = as.integer(jcur), kvol = as.integer(kvol))
  sensemattemp <- array(as.vector(out$sensemat), c(mxe, jcur))
  return(list(sense = sensemattemp))
}
################################################################################################ FUNCTION invcur
invcur <- function(ncur, nex, nez, vol, ezfac, exyz, exyz0, sigma, pod, mxe, mze, 
                   mxnode, jcur, kvol) {
  sensematl <- array(0, c(mxe, jcur))
  for (j in 1:jcur) {
    for (n in 1:nex) {
      summ <- sum(2 * vol[1:nez, n] * (exyz[1, 1:nez, n] * exyz0[1, 1:nez, 
                                                                 n, j] + exyz[2, 1:nez, n] * exyz0[2, 1:nez, n, j] + ezfac[n] * exyz[3, 
                                                                                                                                     1:nez, n] * exyz0[3, 1:nez, n, j]))
      sensematl[n, j] <- -summ * sigma[n]/pod[ncur, j]
    }
  }
  return(list(sense = sensematl))
} 
###############################################################################################
#' FUNCTION resist.
#' @useDynLib AbdominalEITextra
#' @export
resist <- function(j,exyz0,exyz1,height,xcod,ycod,nodes,element,sigmaxy,sigmaz,face,iter,est64){
  #browser()
  interresl <- array(0,c(12,64,32,32))
  for (i in 1:32){
    #cat("i=",i,"\n")
    #DF
    nvect <- c(-diff(ycod[nodes[5,j,i,]]),diff(xcod[nodes[5,j,i,]]))
    nvect <- nvect/norm(as.matrix(nvect),type="f")
    if (nvect %*% c(xcod[nodes[5,j,i,2]],ycod[nodes[5,j,i,2]]) > 0) nvect <- -nvect
    #if (FALSE){
    #if (j<=6){
    interresl[5,(2*i-1)%%64+1,1:32,1:32] <- (
      drop(height[j]*(0*nvect%*%c(-sin(pi*(i-1)/16),cos(pi*(i-1)/16)) + ifelse(j<=6,2,1)*nvect%*%c(-sin(pi*(i-0.5)/16),cos(pi*(i-0.5)/16)))/6)*(
        (sigmaxy[element[4,i]]-sigmaxy[element[6,i]])*exyz0[1,face[5,i,1,j],element[4,i],1:32,j]%o%exyz1[1,face[5,i,3,j],element[6,i],1:32,j] +
          (sigmaxy[element[4,i]]-sigmaxy[element[6,i]])*exyz0[2,face[5,i,1,j],element[4,i],1:32,j]%o%exyz1[2,face[5,i,3,j],element[6,i],1:32,j] +
          ( sigmaz[element[4,i]]- sigmaz[element[6,i]])*exyz0[3,face[5,i,1,j],element[4,i],1:32,j]%o%exyz1[3,face[5,i,3,j],element[6,i],1:32,j])+
        drop(height[j]*(0*nvect%*%c(-sin(pi*(i-1)/16),cos(pi*(i-1)/16)) + ifelse(j<=6,1,2)*nvect%*%c(-sin(pi*(i-0.5)/16),cos(pi*(i-0.5)/16)))/6)*(
          (sigmaxy[element[4,i]]-sigmaxy[element[6,i]])*exyz0[1,face[5,i,2,j],element[4,i],1:32,j]%o%exyz1[1,face[5,i,4,j],element[6,i],1:32,j] +
            (sigmaxy[element[4,i]]-sigmaxy[element[6,i]])*exyz0[2,face[5,i,2,j],element[4,i],1:32,j]%o%exyz1[2,face[5,i,4,j],element[6,i],1:32,j] +
            ( sigmaz[element[4,i]]- sigmaz[element[6,i]])*exyz0[3,face[5,i,2,j],element[4,i],1:32,j]%o%exyz1[3,face[5,i,4,j],element[6,i],1:32,j]))
    #EG
    nvect <- c(-diff(ycod[nodes[6,j,(i-1)%%32+1,]]),diff(xcod[nodes[6,j,(i-1)%%32+1,]]))
    nvect <- nvect/norm(as.matrix(nvect),type="f")
    if (nvect %*% c(xcod[nodes[6,j,(i-1)%%32+1,1]],ycod[nodes[6,j,(i-1)%%32+1,1]]) > 0) nvect <- -nvect
    #if (TRUE){
    #if (j<=6){
    interresl[6,(2*i-1)%%64+1,1:32,1:32] <- (
      drop(height[j]*(ifelse(j<=6,1,2)*nvect%*%c(-sin(pi*(i-0.5)/16),cos(pi*(i-0.5)/16)) + 0*nvect%*%c(-sin(pi*(i-0)/16),cos(pi*(i-0)/16)))/6)*(
        (sigmaxy[element[5,(i-1)%%32+1]]-sigmaxy[element[7,(i-1)%%32+1]])*exyz0[1,face[6,(i-1)%%32+1,1,j],element[5,(i-1)%%32+1],1:32,j]%o%exyz1[1,face[6,(i-1)%%32+1,3,j],element[7,(i-1)%%32+1],1:32,j] +
          (sigmaxy[element[5,(i-1)%%32+1]]-sigmaxy[element[7,(i-1)%%32+1]])*exyz0[2,face[6,(i-1)%%32+1,1,j],element[5,(i-1)%%32+1],1:32,j]%o%exyz1[2,face[6,(i-1)%%32+1,3,j],element[7,(i-1)%%32+1],1:32,j] +
          ( sigmaz[element[5,(i-1)%%32+1]]- sigmaz[element[7,(i-1)%%32+1]])*exyz0[3,face[6,(i-1)%%32+1,1,j],element[5,(i-1)%%32+1],1:32,j]%o%exyz1[3,face[6,(i-1)%%32+1,3,j],element[7,(i-1)%%32+1],1:32,j])+
        drop(height[j]*(ifelse(j<=6,2,1)*nvect%*%c(-sin(pi*(i-0.5)/16),cos(pi*(i-0.5)/16)) + 0*nvect%*%c(-sin(pi*(i-0)/16),cos(pi*(i-0)/16)))/6)*(
          (sigmaxy[element[5,(i-1)%%32+1]]-sigmaxy[element[7,(i-1)%%32+1]])*exyz0[1,face[6,(i-1)%%32+1,2,j],element[5,(i-1)%%32+1],1:32,j]%o%exyz1[1,face[6,(i-1)%%32+1,4,j],element[7,(i-1)%%32+1],1:32,j] +
            (sigmaxy[element[5,(i-1)%%32+1]]-sigmaxy[element[7,(i-1)%%32+1]])*exyz0[2,face[6,(i-1)%%32+1,2,j],element[5,(i-1)%%32+1],1:32,j]%o%exyz1[2,face[6,(i-1)%%32+1,4,j],element[7,(i-1)%%32+1],1:32,j] +
            ( sigmaz[element[5,(i-1)%%32+1]]- sigmaz[element[7,(i-1)%%32+1]])*exyz0[3,face[6,(i-1)%%32+1,2,j],element[5,(i-1)%%32+1],1:32,j]%o%exyz1[3,face[6,(i-1)%%32+1,4,j],element[7,(i-1)%%32+1],1:32,j]))  
    #DF-
    nvect <- c(-diff(ycod[nodes[5,j,(i)%%32+1,]]),diff(xcod[nodes[5,j,(i)%%32+1,]]))
    nvect <- nvect/norm(as.matrix(nvect),type="f")
    if (nvect %*% c(xcod[nodes[5,j,(i)%%32+1,1]],ycod[nodes[5,j,(i)%%32+1,1]]) > 0) nvect <- -nvect
    #if (FALSE){
    #if (j<=6){
    interresl[5,(2*i-0)%%64+1,1:32,1:32] <- (
      drop(height[j]*(ifelse(j<=6,1,2)*nvect%*%c(-sin(pi*(i)/16),cos(pi*(i)/16)) + 0  *nvect%*%c(-sin(pi*(i+0.5)/16),cos(pi*(i+0.5)/16)))/6)*(
        (sigmaxy[element[4,(i)%%32+1]]-sigmaxy[element[6,(i)%%32+1]])*exyz0[1,face[5,(i)%%32+1,1,j],element[4,(i)%%32+1],1:32,j]%o%exyz1[1,face[5,(i)%%32+1,3,j],element[6,(i)%%32+1],1:32,j] +
          (sigmaxy[element[4,(i)%%32+1]]-sigmaxy[element[6,(i)%%32+1]])*exyz0[2,face[5,(i)%%32+1,1,j],element[4,(i)%%32+1],1:32,j]%o%exyz1[2,face[5,(i)%%32+1,3,j],element[6,(i)%%32+1],1:32,j] +
          ( sigmaz[element[4,(i)%%32+1]]- sigmaz[element[6,(i)%%32+1]])*exyz0[3,face[5,(i)%%32+1,1,j],element[4,(i)%%32+1],1:32,j]%o%exyz1[3,face[5,(i)%%32+1,3,j],element[6,(i)%%32+1],1:32,j])+
        drop(height[j]*(ifelse(j<=6,2,1)*nvect%*%c(-sin(pi*(i)/16),cos(pi*(i)/16)) + 0  *nvect%*%c(-sin(pi*(i+0.5)/16),cos(pi*(i+0.5)/16)))/6)*(
          (sigmaxy[element[4,(i)%%32+1]]-sigmaxy[element[6,(i)%%32+1]])*exyz0[1,face[5,(i)%%32+1,2,j],element[4,(i)%%32+1],1:32,j]%o%exyz1[1,face[5,(i)%%32+1,4,j],element[6,(i)%%32+1],1:32,j] +
            (sigmaxy[element[4,(i)%%32+1]]-sigmaxy[element[6,(i)%%32+1]])*exyz0[2,face[5,(i)%%32+1,2,j],element[4,(i)%%32+1],1:32,j]%o%exyz1[2,face[5,(i)%%32+1,4,j],element[6,(i)%%32+1],1:32,j] +
            ( sigmaz[element[4,(i)%%32+1]]- sigmaz[element[6,(i)%%32+1]])*exyz0[3,face[5,(i)%%32+1,2,j],element[4,(i)%%32+1],1:32,j]%o%exyz1[3,face[5,(i)%%32+1,4,j],element[6,(i)%%32+1],1:32,j]))
    #EG-
    nvect <- c(-diff(ycod[nodes[6,j,i,]]),diff(xcod[nodes[6,j,i,]]))
    nvect <- nvect/norm(as.matrix(nvect),type="f")
    if (nvect %*% c(xcod[nodes[6,j,i,2]],ycod[nodes[6,j,i,2]]) > 0) nvect <- -nvect
    #if (FALSE){
    #if (j<=6){
    interresl[6,(2*i-0)%%64+1,1:32,1:32] <- (
      drop(height[j]*(0  *nvect%*%c(-sin(pi*(i-0.5)/16),cos(pi*(i-0.5)/16)) + ifelse(j<=6,2,1)*nvect%*%c(-sin(pi*(i)/16),cos(pi*(i)/16)))/6)*(
        (sigmaxy[element[5,i]]-sigmaxy[element[7,i]])*exyz0[1,face[6,i,1,j],element[5,i],1:32,j]%o%exyz1[1,face[6,i,3,j],element[7,i],1:32,j] +
          (sigmaxy[element[5,i]]-sigmaxy[element[7,i]])*exyz0[2,face[6,i,1,j],element[5,i],1:32,j]%o%exyz1[2,face[6,i,3,j],element[7,i],1:32,j] +
          ( sigmaz[element[5,i]]- sigmaz[element[7,i]])*exyz0[3,face[6,i,1,j],element[5,i],1:32,j]%o%exyz1[3,face[6,i,3,j],element[7,i],1:32,j])+
        drop(height[j]*(0  *nvect%*%c(-sin(pi*(i-0.5)/16),cos(pi*(i-0.5)/16)) + ifelse(j<=6,1,2)*nvect%*%c(-sin(pi*(i)/16),cos(pi*(i)/16)))/6)*(
          (sigmaxy[element[5,i]]-sigmaxy[element[7,i]])*exyz0[1,face[6,i,2,j],element[5,i],1:32,j]%o%exyz1[1,face[6,i,4,j],element[7,i],1:32,j] +
            (sigmaxy[element[5,i]]-sigmaxy[element[7,i]])*exyz0[2,face[6,i,2,j],element[5,i],1:32,j]%o%exyz1[2,face[6,i,4,j],element[7,i],1:32,j] +
            ( sigmaz[element[5,i]]- sigmaz[element[7,i]])*exyz0[3,face[6,i,2,j],element[5,i],1:32,j]%o%exyz1[3,face[6,i,4,j],element[7,i],1:32,j]))
  }
  return(interres=interresl)
}
#################################################################################################################################################################
