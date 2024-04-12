#' ToothVariations
#' 
#' Calculate and return morphometric map variation
#' @param PCA list: output from ToothPCA
#' @param scores numeric: principal component value
#' @param PC numeric: principal component eigenvalue
#' @param asp numeric: x,y ratio of the morphometric map
#' @param pal character vector: vector of colors 
#' @param meanmap logical: if TRUE the mean map corresponds to the real mean morphometric maps, otherwise the mean map is 0
#' @return XYZ data.frame of morphometric map variation
#' @author Antonio Profico; Mathilde Augoyard
#' @export

ToothVariations<-function (PCA, scores, PC, asp = 2,pal = blue2green2red(101), meanmap=TRUE) 
{
  if(isFALSE(meanmap)){PCA$meanMap<-rep(0,length(PCA$meanMap))}
  mapvar <- matrix(restoreFromPCA(scores, PC, PCA$meanMap), 
                   nrow = dim(PCA$CorMaps)[1], ncol = dim(PCA$CorMaps)[2])
  map <- levelplot(mapvar, asp = asp, col.regions = pal)
  graphics::plot(map)
  return(mapvar)
}