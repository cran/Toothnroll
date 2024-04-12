#' ToothDF
#' 
#' Tool to build a data.frame suitable for morphometric maps
#' @param tooth.thickness list: tooth.Thickness object
#' @param rem.out logical: if TRUE the outlier will be removed 
#' @param fac.out numeric: parameter to set the threshold in outliers detection
#' @param smooth logical: if TRUE the smooth algorithm is applied
#' @param scale logical: if TRUE the thichkness matrix is scaled from 0 to 1 
#' @param smooth.iter numeric: number of smoothing iterations 
#' @param method character: if set on "equiangular" the dental thickness is meant as the distance of the segment intersecting the external and internal outline starting from the centroid of the section. If set on "closest" the dental thickness is calculated at each point as the closest distance between external and internal outlines
#' @param labels character vector: names for x labels in the morphometric map
#' @param relThick logical: if TRUE the dental thickness is converted into relative dental thickness
#' @return XYZ data.frame for morphometric map
#' @return labels character vector for x labels in the morphometric map
#' @author Antonio Profico; Mathilde Augoyard
#' @export

ToothDF<-function (tooth.thickness, rem.out = TRUE, fac.out = 0.5, 
                   smooth = TRUE, scale = TRUE, smooth.iter = 5, method = "equiangular",labels =c("Li","Mes", "Bu", "D","Li"),
                   relThick=FALSE) {
  tooth.shape <- tooth.thickness$tooth.shape
  if (method == "equiangular") {
    thicks <- tooth.thickness$sect_thickness_eq
  if(isTRUE(relThick)) {thicks<-tooth.thickness$sect_Rthickness_eq}
  }
  
  if (method == "closest") {
    thicks <- tooth.thickness$sect_thickness_cp
    if(isTRUE(relThick)) {thicks<-tooth.thickness$sect_Rthickness_cp}
    
  }
  X <- NULL
  Y <- NULL
  Z <- NULL
  summary <- dim(thicks)
  num.sect <- summary[3]
  num.land <- summary[1]
  if (scale == TRUE) {
    x <- as.vector(thicks)
    thickvector <- (x - min(x))/(max(x) - min(x))
    Thick_arr <- array(thickvector, dim = c(num.land, 1, 
                                            num.sect))
    thicks <- Thick_arr
  }
  if (rem.out == TRUE) {
    x <- as.vector(thicks)
    qnt <- stats::quantile(x, probs = c(0.25, 0.75))
    H <- fac.out * IQR(x)
    y <- x
    y[x < (qnt[1] - H)] <- min(y[x > (qnt[1] - H)])
    y[x > (qnt[2] + H)] <- max(y[x < (qnt[2] + H)])
    thick_values <- y
    Thick_arr <- array(thick_values, dim = c(num.land, 1, 
                                             num.sect))
    thicks <- Thick_arr
    if (scale == TRUE) {
      x <- as.vector(thicks)
      thickvector <- (x - min(x))/(max(x) - min(x))
      Thick_arr <- array(thickvector, dim = c(num.land, 
                                              1, num.sect))
      thicks <- Thick_arr
    }
  }
  if (smooth == TRUE) {
    thick_smoo <- apply(thicks, 2, function(x) (matrixSmooth(x, 
                                                             passes = smooth.iter)))
    thick_out <- array(thick_smoo, dim = c(num.land, 1, num.sect))
    thicks <- thick_out
  }
  for (i in 1:num.sect) {
    X_i <- seq(1, 100, length.out = num.land)
    Y_i <- tooth.shape$"out3D"[, 3, i]
    Z_i <- thicks[, , i]
    X <- c(X, X_i)
    Y <- c(Y, Y_i)
    Z <- c(Z, Z_i)
  }

  XYZ <- data.frame(X, Y, Z)
  return(list(XYZ = XYZ, labels = labels))
}
