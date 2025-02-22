#' Tooth2Dmap
#'
#' Create 2D morphometric maps of enamel/dentin thickness 
#' @param tooth.shape list: output from the function ToothShape 
#' @param input list: output from the function ToothAlignment 
#' @param rem.out logical: if TRUE outliers will be removed 
#' @param fac.out numeric: parameter to set the threshold in outlier detection
#' @param smooth logical: if TRUE a smoothing filter is applied
#' @param scale logical: if TRUE the thichkness matrix is scaled from 0 to 1 
#' @param smooth.iter numeric: number of smoothing iterations 
#' @param gamMap logical: if TRUE gam smoothing is applied
#' @param nrow numeric: number of rows for gam smoothing matrix
#' @param ncol numeric: number of columns for gam smoothing matrix
#' @param gdl numeric: number of degree of freedom for gam smoothing matrix
#' @param method character: if set on "equiangular" the dentine or enamel thickness is meant as the distance of the segment intersecting the external and internal outline starting from the centroid of the section. If set on "closest" the dentine or enamel thickness is calculated at each point as the closest distance between external and internal outlines
#' @param plot logical: if TRUE the 2D morphometric map is plotted
#' @param pal character vector: colors to be used in the map production
#' @param aspect numeric: axis ratio for 2D morphometric map
#' @param labels character vector: names for x labels in the morphometric map
#' @param ylab character vector: label for y axis in the morphometric map
#' @return dataframe dataframe for colormap production
#' @return 2Dmap thickness color map
#' @return gamoutput output from GAM
#' @return data input used to build the GAM map
#' @author Antonio Profico; Mathilde Augoyard
#' @examples
#' \donttest{
#' data("URI1_tooth")
#' require(morphomap)
#' Enamel<-URI1_tooth$mesh1
#' Dentin<-URI1_tooth$mesh2
#' Pulp<-URI1_tooth$mesh3
#' outline<-URI1_tooth$outline
#' set<-URI1_tooth$set
#' #Map of the crown
#' AlignMeshes<-ToothAlignment(mesh1=Enamel,mesh2=Dentin,mesh3=Pulp,set,outline,analyse = "c")
#' #Virtual sectioning dentine-pulp
#' External<-AlignMeshes$almesh1$mesh
#' Internal<-AlignMeshes$almesh2$mesh
#' #Define 16 cross-sections from the 30% to the 90% along the crown
#' Core<-ToothCore(External,Internal,num.points = 1000,num.sect =16,
#'                     bio.len = AlignMeshes$length,start=0.3,end=0.9)
#' #Extract 25 equiangular semilandmarks from each cross-section (anticlockwise)
#' Shape<-ToothShape(Core,25,sects_vector = NULL,cent.out = "E",direction="a")
#' Tooth2Dmap(Shape,AlignMeshes,rem.out =TRUE,scale=FALSE,smooth = FALSE,aspect = 0.5,gamMap = FALSE,
#'            nrow = 100,ncol = 100,gdl = 250,method="equiangular")
#' 
#' #Map of the root
#' AlignMeshes<-ToothAlignment(mesh1=Enamel,mesh2=Dentin,mesh3=Pulp,set,outline,analyse = "r")
#' #Virtual sectioning dentine-pulp
#' External<-AlignMeshes$almesh2$mesh
#' Internal<-AlignMeshes$almesh3$mesh
#' #Define 16 cross-sections from the 10% to the 50% along the root
#' Core<-ToothCore(External,Internal,num.points = 1000,num.sect =16,
#'                     bio.len = AlignMeshes$length,start=0.1,end=0.5)
#' #Extract 25 equiangular semilandmarks from each cross-section (anticlockwise)
#' Shape<-ToothShape(Core,25,sects_vector = NULL,cent.out = "E",direction="a")
#' Tooth2Dmap(Shape,AlignMeshes,rem.out =FALSE,scale=FALSE,smooth = FALSE,aspect = 0.5,gamMap = FALSE,
#'            nrow = 100,ncol = 100,gdl = 250,method="equiangular")
#' }
#' @export

Tooth2Dmap<-function (tooth.shape, input, rem.out = FALSE, fac.out = 0.5, 
                      smooth = FALSE, scale = TRUE, smooth.iter = 5, gamMap = FALSE, 
                      nrow = 120, ncol = 80, gdl = 250, method = "equiangular", 
                      plot = TRUE, pal = blue2green2red(101), aspect = 0.6,
                      labels=c("Li","Mes", "Bu", "D","Li"),ylab="") 
{
  thickarray <- ToothThickness(tooth.shape)
  data_t <- ToothDF(thickarray, rem.out = rem.out, fac.out = fac.out, 
                    smooth = smooth, scale = scale, smooth.iter = smooth.iter, 
                    method = method,labels=labels)
  data <- data_t$XYZ
  labels <- data_t$labels
  if (gamMap == FALSE) {
    peri <- sum(sqrt(rowSums(diff(input$aloutline)^2)))
    selm <- input$margins_sel
    if (selm[1] != 1) {
      dist1 <- (sum(sqrt(rowSums(diff(input$aloutline[1:selm[1], 
      ])^2))) * 100)/peri
    }
    else {
      dist1 = 1
    }
    dist2 <- (sum(sqrt(rowSums(diff(input$aloutline[1:selm[2], 
    ])^2))) * 100)/peri
    dist3 <- (sum(sqrt(rowSums(diff(input$aloutline[1:selm[3], 
    ])^2))) * 100)/peri
    dist4 <- (sum(sqrt(rowSums(diff(input$aloutline[1:selm[4], 
    ])^2))) * 100)/peri
    minx <- min(data$X, na.rm = T)
    miny <- min(data$Y, na.rm = T)
    maxx <- max(data$X, na.rm = T)
    maxy <- max(data$Y, na.rm = T)
    maxz <- max(data$Z, na.rm = T)
    minz <- min(data$Z, na.rm = T)
    atvalues <- seq(minz, maxz, length.out = 100)
    xat <- seq(1, 100, length.out = 5)
    ylabels <- round(seq(tooth.shape$start, tooth.shape$end, 
                         length.out = 5) * 100, 2)
    yat <- round(seq(miny, maxy, length.out = 5), 2)
    map <- levelplot(data[, 3] ~ data[, 1] + data[, 2], 
                     panel = function(...) {
                       panel.levelplot(...)
                       panel.abline(v = c(dist1, dist2, dist3, dist4), 
                                    lwd = 3)
                     }, xlab = "", aspect = aspect, 
                     main = "2D thickness map", at = atvalues, col.regions = pal, 
                     scales = list(x = list(at = xat, labels = labels, 
                                            rot = 90, alternating = 1), y = list(at = yat, 
                                                                                 labels = ylabels, rot = 90, alternating = 1)))
    if (plot == TRUE) {
      graphics::plot(map)
    }
    mat <- data
  }
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
    for (i in seq(1, nrow)) for (j in seq(1, ncol)) {
      xp[(i - 1) * ncol + j] <- xx[j]
      yp[(i - 1) * ncol + j] <- yy[i]
      yc[(i - 1) * ncol + j] <- yy[i]
      xc[(i - 1) * ncol + j] <- xx[j]
    }
    fit <- predict.gam(m1, list(X = xp, Y = yp))
    data <- data.frame(X = xc, Y = yc, Z = fit)
    if (scale == TRUE) {
      x <- data$Z
      thickvector <- (x - min(x))/(max(x) - min(x))
      data$Z <- thickvector
    }
    minx <- min(data$X, na.rm = T)
    miny <- min(data$Y, na.rm = T)
    maxx <- max(data$X, na.rm = T)
    maxy <- max(data$Y, na.rm = T)
    maxz <- max(data$Z, na.rm = T)
    minz <- min(data$Z, na.rm = T)
    atvalues <- seq(minz, maxz, length.out = 100)
    xat <- seq(1, 100, length.out = 5)
    ylabels <- round(seq(tooth.shape$start, tooth.shape$end, 
                         length.out = 5) * 100, 2)
    yat <- round(seq(miny, maxy, length.out = 5), 2)
    peri <- sum(sqrt(rowSums(diff(input$aloutline)^2)))
    selm <- input$margins_sel
    if (selm[1] != 1) {
      dist1 <- (sum(sqrt(rowSums(diff(input$aloutline[1:selm[1], 
      ])^2))) * 100)/peri
    }
    else {
      dist1 = 1
    }
    dist2 <- (sum(sqrt(rowSums(diff(input$aloutline[1:selm[2], 
    ])^2))) * 100)/peri
    dist3 <- (sum(sqrt(rowSums(diff(input$aloutline[1:selm[3], 
    ])^2))) * 100)/peri
    dist4 <- (sum(sqrt(rowSums(diff(input$aloutline[1:selm[4], 
    ])^2))) * 100)/peri
    map <- levelplot(data[, 3] ~ data[, 1] + data[, 2], 
                     panel = function(...) {
                       panel.levelplot(...)
                       panel.abline(v = c(dist1, dist2, dist3, dist4), 
                                    lwd = 3)
                     }, xlab = "", ylab=ylab,aspect = aspect, 
                     main = "2D thickness map", at = atvalues, col.regions = pal, 
                     scales = list(x = list(at = xat, labels = labels, 
                                            rot = 90, alternating = 1), y = list(at = yat, 
                                                                                 labels = ylabels, rot = 90, alternating = 1)))
    if (plot == TRUE) {
      graphics::plot(map)
    }
    mat <- data
  }
  if (gamMap == TRUE) {
    out <- list(dataframe = mat, "map2D" = map, gamoutput = m1, 
                data = data)
  }
  if (gamMap == FALSE) {
    out <- list(dataframe = mat, "map2D" = map, gamoutput = NULL, 
                data = data)
  }
  return(out)
}