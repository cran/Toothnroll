#' ToothPCA
#'
#' Perform the Principal Component Analysis on a list of tooth.shape objects
#' @param mpShapeList list: tooth.shape objects
#' @param gamMap logical: if TRUE gamMap spline method is applied
#' @param nrow numeric: number of rows if gamMap is TRUE
#' @param ncol numeric: number of columns if gamMap is TRUE
#' @param gdl numeric: degree of freedom (if gamMap is TRUE)
#' @param rem.out logical: if TRUE outliers are removed
#' @param scaleThick logical: if TRUE thickness values are scaled from 0 to 1
#' @param relThick logical: if TRUE the thickness values are scaled by the diameter from the centroid to the external outline
#' @param fac.out numeric: threshold to define an outlier observation
#' @param method character: "equiangular" or "closest" to define the thickness from evenly spaced or closest semilandmarks between the external and internal outline
#' @param scalePCA logical: indicate whether the variables should be scaled to have unit variance
#' @return PCscores matrix of PC scores 
#' @return PCs principal components 
#' @return variance table of the explained Variance by the PCs 
#' @return meanMap mean map
#' @return CorMaps maps of thickness used as input in the PCA
#' @author Antonio Profico; Mathilde Augoyard
#' @examples
#' ### Example on the canine crown
#' data("UCcrown")
#' require(morphomap)
#' shapeList<-UCcrown
#' PCA<-ToothPCA(shapeList,gamMap = FALSE,scaleThick = TRUE,scalePCA = TRUE ,relThick = FALSE)
#'  \donttest{
#' #gamMap set on TRUE
#' PCA<-ToothPCA(shapeList,gamMap = TRUE,scaleThick = TRUE,scalePCA = TRUE,relThick = FALSE)
#' }
#' otu<-substr(names(shapeList),1,2)
#' pchs <- ifelse(otu == "MH", 16, 17)
#' cols <- ifelse(otu == "MH", "orange", "darkblue")
#' 
#' plot(PCA$PCscores,col=cols,cex=1, pch = pchs,
#'      xlab=paste("PC1 (",round(PCA$Variance[1,2],2),"%)"),
#'      ylab=paste("PC2 (",round(PCA$Variance[2,2],2),"%)"),
#'      cex.lab=1,cex.axis=1)
#' title (main="UC (radicular dentine)", font.main= 1,adj = 0, cex.main = 1.2)
#' legend("topright", legend = c("MH", "NE"), col = c("orange", "darkblue"), pch = c(16,17), cex = 0.8)
#' abline(v=0,h=0,col="black",lwd=2,lty=3)
#' 
#' hpts1 <- chull(PCA$PCscores[which(otu=="MH"),1:2])
#' hpts1 <- c(hpts1, hpts1[1])
#' polygon(PCA$PCscores[which(otu=="MH")[hpts1],1:2 ], col = adjustcolor("orange", 0.3), border = NA)
#' 
#' hpts2 <- chull(PCA$PCscores[which(otu=="NE"),c(1:2)])
#' hpts2 <- c(hpts2, hpts2[1])
#' polygon(PCA$PCscores[which(otu=="NE")[hpts2], ], col = adjustcolor("darkblue", 0.3), border = NA)
#' 
#' PC1min<-ToothVariations(PCA,min(PCA$PCscores[,1]),PCA$PCs[,1],asp=0.5,meanmap = FALSE)
#' PC1max<-ToothVariations(PCA,max(PCA$PCscores[,1]),PCA$PCs[,1],asp=0.5,meanmap = FALSE)
#' 
#' ### Example on the canine root
#' data("UCroot")
#' require(morphomap)
#' shapeList<-UCroot
#' PCA<-ToothPCA(shapeList,gamMap = FALSE,scaleThick = TRUE,scalePCA = TRUE,relThick = FALSE)
#' otu<-substr(names(shapeList),1,2)
#' pchs <- ifelse(otu == "MH", 16, 17)
#' cols <- ifelse(otu == "MH", "orange", "darkblue")
#' 
#' plot(PCA$PCscores,col=cols,cex=1, pch = pchs,
#'      xlab=paste("PC1 (",round(PCA$Variance[1,2],2),"%)"),
#'      ylab=paste("PC2 (",round(PCA$Variance[2,2],2),"%)"),
#'      cex.lab=1,cex.axis=1)
#' title (main="UC (radicular dentine)", font.main= 1,adj = 0, cex.main = 1.2)
#' legend("topright", legend = c("MH", "NE"), col = c("orange", "darkblue"), pch = c(16,17), cex = 0.8)
#' abline(v=0,h=0,col="black",lwd=2,lty=3)
#' 
#' hpts1 <- chull(PCA$PCscores[which(otu=="MH"),1:2])
#' hpts1 <- c(hpts1, hpts1[1])
#' polygon(PCA$PCscores[which(otu=="MH")[hpts1],1:2 ], col = adjustcolor("orange", 0.3), border = NA)
#' 
#' hpts2 <- chull(PCA$PCscores[which(otu=="NE"),c(1:2)])
#' hpts2 <- c(hpts2, hpts2[1])
#' polygon(PCA$PCscores[which(otu=="NE")[hpts2], ], col = adjustcolor("darkblue", 0.3), border = NA)
#' 
#' PC1min<-ToothVariations(PCA,min(PCA$PCscores[,1]),PCA$PCs[,1],asp=0.5,meanmap = FALSE)
#' PC1max<-ToothVariations(PCA,max(PCA$PCscores[,1]),PCA$PCs[,1],asp=0.5,meanmap = FALSE)
#'
#' ### Example on the central upper incisor (crown)
#' data("UI1crown")
#' require(morphomap)
#' shapeList<-UI1crown
#' PCA<-ToothPCA(shapeList,gamMap = FALSE,scaleThick = TRUE,scalePCA = TRUE ,relThick = FALSE)
#' otu<-substr(names(shapeList),1,2)
#' pchs <- ifelse(otu == "MH", 16, 17)
#' cols <- ifelse(otu == "MH", "orange", "darkblue")
#' 
#' plot(PCA$PCscores,col=cols,cex=1, pch = pchs,
#'      xlab=paste("PC1 (",round(PCA$Variance[1,2],2),"%)"),
#'      ylab=paste("PC2 (",round(PCA$Variance[2,2],2),"%)"),
#'      cex.lab=1,cex.axis=1)
#' title (main="UC (radicular dentine)", font.main= 1,adj = 0, cex.main = 1.2)
#' legend("topright", legend = c("MH", "NE"), col = c("orange", "darkblue"), pch = c(16,17), cex = 0.8)
#' abline(v=0,h=0,col="black",lwd=2,lty=3)
#' 
#' hpts1 <- chull(PCA$PCscores[which(otu=="MH"),1:2])
#' hpts1 <- c(hpts1, hpts1[1])
#' polygon(PCA$PCscores[which(otu=="MH")[hpts1],1:2 ], col = adjustcolor("orange", 0.3), border = NA)
#' 
#' hpts2 <- chull(PCA$PCscores[which(otu=="NE"),c(1:2)])
#' hpts2 <- c(hpts2, hpts2[1])
#' polygon(PCA$PCscores[which(otu=="NE")[hpts2], ], col = adjustcolor("darkblue", 0.3), border = NA)
#' 
#' PC1min<-ToothVariations(PCA,min(PCA$PCscores[,1]),PCA$PCs[,1],asp=0.5,meanmap = FALSE)
#' PC1max<-ToothVariations(PCA,max(PCA$PCscores[,1]),PCA$PCs[,1],asp=0.5,meanmap = FALSE)
#' 
#' ### Example on the upper central incisor (root)
#' data("UI1root")
#' require(morphomap)
#' shapeList<-UI1root
#' PCA<-ToothPCA(shapeList,gamMap = FALSE,scaleThick = TRUE,scalePCA = TRUE,relThick = FALSE)
#' otu<-substr(names(UI1root),1,2)
#' pchs <- ifelse(otu == "MH", 16, 17)
#' cols <- ifelse(otu == "MH", "orange", "darkblue")
#' 
#' plot(PCA$PCscores,col=cols,cex=1, pch = pchs,
#'      xlab=paste("PC1 (",round(PCA$Variance[1,2],2),"%)"),
#'      ylab=paste("PC2 (",round(PCA$Variance[2,2],2),"%)"),
#'      cex.lab=1,cex.axis=1)
#' title (main="UC (radicular dentine)", font.main= 1,adj = 0, cex.main = 1.2)
#' legend("topright", legend = c("MH", "NE"), col = c("orange", "darkblue"), pch = c(16,17), cex = 0.8)
#' abline(v=0,h=0,col="black",lwd=2,lty=3)
#' 
#' hpts1 <- chull(PCA$PCscores[which(otu=="MH"),1:2])
#' hpts1 <- c(hpts1, hpts1[1])
#' polygon(PCA$PCscores[which(otu=="MH")[hpts1],1:2 ], col = adjustcolor("orange", 0.3), border = NA)
#' 
#' hpts2 <- chull(PCA$PCscores[which(otu=="NE"),c(1:2)])
#' hpts2 <- c(hpts2, hpts2[1])
#' polygon(PCA$PCscores[which(otu=="NE")[hpts2], ], col = adjustcolor("darkblue", 0.3), border = NA)
#' 
#' PC1min<-ToothVariations(PCA,min(PCA$PCscores[,1]),PCA$PCs[,1],asp=0.5,meanmap = FALSE)
#' PC1max<-ToothVariations(PCA,max(PCA$PCscores[,1]),PCA$PCs[,1],asp=0.5,meanmap = FALSE)
#'
#' ### Example on the lateral upper incisor (crown)
#' data("UI2crown")
#' require(morphomap)
#' shapeList<-UI2crown
#' PCA<-ToothPCA(shapeList,gamMap = FALSE,scaleThick = TRUE,scalePCA = TRUE ,relThick = FALSE)
#' otu<-substr(names(shapeList),1,2)
#' pchs <- ifelse(otu == "MH", 16, 17)
#' cols <- ifelse(otu == "MH", "orange", "darkblue")
#' 
#' plot(PCA$PCscores,col=cols,cex=1, pch = pchs,
#'      xlab=paste("PC1 (",round(PCA$Variance[1,2],2),"%)"),
#'      ylab=paste("PC2 (",round(PCA$Variance[2,2],2),"%)"),
#'      cex.lab=1,cex.axis=1)
#' title (main="UC (radicular dentine)", font.main= 1,adj = 0, cex.main = 1.2)
#' legend("topright", legend = c("MH", "NE"), col = c("orange", "darkblue"), pch = c(16,17), cex = 0.8)
#' abline(v=0,h=0,col="black",lwd=2,lty=3)
#' 
#' hpts1 <- chull(PCA$PCscores[which(otu=="MH"),1:2])
#' hpts1 <- c(hpts1, hpts1[1])
#' polygon(PCA$PCscores[which(otu=="MH")[hpts1],1:2 ], col = adjustcolor("orange", 0.3), border = NA)
#' 
#' hpts2 <- chull(PCA$PCscores[which(otu=="NE"),c(1:2)])
#' hpts2 <- c(hpts2, hpts2[1])
#' polygon(PCA$PCscores[which(otu=="NE")[hpts2], ], col = adjustcolor("darkblue", 0.3), border = NA)
#' 
#' PC1min<-ToothVariations(PCA,min(PCA$PCscores[,1]),PCA$PCs[,1],asp=0.5,meanmap = FALSE)
#' PC1max<-ToothVariations(PCA,max(PCA$PCscores[,1]),PCA$PCs[,1],asp=0.5,meanmap = FALSE)
#' 
#' ### Example on the upper lateral incisor (root)
#' data("UI2root")
#' require(morphomap)
#' shapeList<-UI2root
#' PCA<-ToothPCA(shapeList,gamMap = FALSE,scaleThick = TRUE,scalePCA = TRUE,relThick = FALSE)
#' otu<-substr(names(UI2root),1,2)
#' pchs <- ifelse(otu == "MH", 16, 17)
#' cols <- ifelse(otu == "MH", "orange", "darkblue")
#' 
#' plot(PCA$PCscores,col=cols,cex=1, pch = pchs,
#'      xlab=paste("PC1 (",round(PCA$Variance[1,2],2),"%)"),
#'      ylab=paste("PC2 (",round(PCA$Variance[2,2],2),"%)"),
#'      cex.lab=1,cex.axis=1)
#' title (main="UC (radicular dentine)", font.main= 1,adj = 0, cex.main = 1.2)
#' legend("topright", legend = c("MH", "NE"), col = c("orange", "darkblue"), pch = c(16,17), cex = 0.8)
#' abline(v=0,h=0,col="black",lwd=2,lty=3)
#' 
#' hpts1 <- chull(PCA$PCscores[which(otu=="MH"),1:2])
#' hpts1 <- c(hpts1, hpts1[1])
#' polygon(PCA$PCscores[which(otu=="MH")[hpts1],1:2 ], col = adjustcolor("orange", 0.3), border = NA)
#' 
#' hpts2 <- chull(PCA$PCscores[which(otu=="NE"),c(1:2)])
#' hpts2 <- c(hpts2, hpts2[1])
#' polygon(PCA$PCscores[which(otu=="NE")[hpts2], ], col = adjustcolor("darkblue", 0.3), border = NA)
#' 
#' PC1min<-ToothVariations(PCA,min(PCA$PCscores[,1]),PCA$PCs[,1],asp=0.5,meanmap = FALSE)
#' PC1max<-ToothVariations(PCA,max(PCA$PCscores[,1]),PCA$PCs[,1],asp=0.5,meanmap = FALSE)
#' @export

ToothPCA<-function (mpShapeList, gamMap = FALSE, nrow = 120, ncol = 80, gdl=250,
          rem.out = TRUE, scaleThick = FALSE, relThick=FALSE,fac.out = 1.5, method = "equiangular", 
          scalePCA = TRUE) 
  
  
 
{ 
  if (length(unique(unlist(lapply(mpShapeList, function(x) dim(x$"out3D")[1])))) != 
      1) {
    stop("the number of landmarks is different")
  }
  if (length(unique(unlist(lapply(mpShapeList, function(x) dim(x$"out3D")[3])))) != 
      1) {
    stop("the number of cross-sections is different")
  }
  
  
 
  
  if (isTRUE(gamMap)) {
    nrows <- nrow
    ncols <- ncol}else{nrow <-dim(mpShapeList[[1]]$"out3D")[1];
        ncol <-dim(mpShapeList[[1]]$"out3D")[3]}
  
  if (isTRUE(gamMap)) {
    ncol<-nrows
    nrow<-ncols
  }
  
  if (isTRUE(gamMap)) {
  ThickMat <- array(NA, dim = c(ncol, nrow, length(mpShapeList)))}else{
  ThickMat <- array(NA, dim = c(nrow, ncol, length(mpShapeList)))}
  
  for(i in 1:length(mpShapeList)){
  thickarray_i <- ToothThickness(mpShapeList[[i]])
  data_t <- ToothDF(thickarray_i, rem.out = rem.out, fac.out = fac.out, 
                    smooth = FALSE, scale = scaleThick, smooth.iter = NULL, 
                    method = method,relThick=relThick)

  data<-data_t$XYZ
  
  if (gamMap == TRUE) {
    m1 <- gam(Z ~ s(X, Y, k = gdl), data = data)
    minx <- min(data$X, na.rm = T)
    maxx <- max(data$X, na.rm = T)
    miny <- min(data$Y, na.rm = T)
    maxy <- max(data$Y, na.rm = T)
    wx <- maxx - minx
    wy <- maxy - miny
    sx <- wx/ncol
    sy <- wy/nrow
    xx <- seq(minx, maxx, length.out = ncol)
    yy <- seq(miny, maxy, length.out = nrow)
    xp <- vector()
    yp <- vector()
    xc <- vector()
    yc <- vector()
    for (m in seq(1, nrow)) for (j in seq(1, ncol)) {
      xp[(m - 1) * ncol + j] <- xx[j]
      yp[(m - 1) * ncol + j] <- yy[m]
      yc[(m - 1) * ncol + j] <- yy[m]
      xc[(m - 1) * ncol + j] <- xx[j]
    }
    fit <- predict.gam(m1, list(X = xp, Y = yp))  
    levelplot(matrix(fit, nrow = ncol, ncol = nrow),asp=.5)
    ThickMat[, , i] <- matrix(fit, nrow = ncol, ncol = nrow)
  }
  if (gamMap == FALSE){
  mati<-matrix(data[,3], nrow = nrow, ncol = ncol)
  ThickMat[, , i] <- matrix(mati, nrow = nrow, ncol = ncol)
  }
  }

  CortThick_mat <- vecx(ThickMat)
  Thick_PCA <- prcomp(CortThick_mat, scale. = scalePCA)
  values <- 0
  eigv <- Thick_PCA$sdev^2
  values <- eigv[which(eigv > 1e-16)]
  lv <- length(values)
  PCs <- Thick_PCA$rotation[, 1:lv]
  PCscores <- as.matrix(Thick_PCA$x[, 1:lv])
  Variance <- cbind(sqrt(eigv), eigv/sum(eigv), cumsum(eigv)/sum(eigv)) * 
    100
  Variance <- Variance[1:lv, ]
  out <- list(PCscores = PCscores, PCs = PCs, Variance = Variance, 
              meanMap = Thick_PCA$center, CorMaps = ThickMat)
  return(out)
}
