#' ToothShape
#' 
#' Tool for the extraction of equiangular landmarks on the entire dental region of interest 
#' @param tooth.core list: tooth.core object 
#' @param num.land numeric: number of landmarks defining each section
#' @param sects_vector numeric: number of sections
#' @param cent.out how to define the center of each section. The method allowed are "CCA" (center of cortical area), "E" (barycenter of the external outline) and "I" (barycenter of the internal outline)
#' @param delta pixel size used to calculate the CCA
#' @param direction clockwise if "c", anticlockwise if"a" 
#' @param out.sur mesh: if provided, the external outlines will be projected back on the surface
#' @param inn.sur mesh: if provided, the internal outlines will be projected back on the surface 
#' @return out3D num.pointsx3xnum.sect array in which the external outlines are stored
#' @return inn3D num.pointsx3xnum.sect array in which the internal outlines are stored
#' @return out2D num.pointsx2xnum.sect array in which the external outlines are stored
#' @return inn2D num.pointsx2xnum.sect array in which the interal outlines are stored
#' @return ALPM_inn array with the coordinates of ALPM coordinates on the external outline
#' @return ALPM_out array with the coordinates of ALPM coordinates on the internal outline
#' @return mech_length length of the selected region of interest
#' @return start percentage of the mechanical length from which the first section is defined
#' @return end percentage of the mechanical length from which the last section is defined
#' @return centroids of the cross-sections
#' @author Antonio Profico; Mathilde Augoyard
#' @examples
#' \donttest{
#' data("URI1_tooth")
#' require(morphomap)
#' require(rgl)
#' Enamel<-URI1_tooth$mesh1
#' Dentin<-URI1_tooth$mesh2
#' Pulp<-URI1_tooth$mesh3
#' outline<-URI1_tooth$outline
#' set<-URI1_tooth$set
#' 
#' #Unrolling the crown
#' AlignMeshes<-ToothAlignment(mesh1=Enamel,mesh2=Dentin,mesh3=Pulp,set,outline,analyse = "c")
#' #Virtual sectioning dentine-pulp
#' External<-AlignMeshes$almesh1$mesh
#' Internal<-AlignMeshes$almesh2$mesh
#' #Define 16 cross-sections from the 30% to the 90% along the crown
#' Core<-ToothCore(External,Internal,num.points = 1000,num.sect =16,
#'                 bio.len = AlignMeshes$length,start=0.3,end=0.9)
#' Shape<-ToothShape(Core,num.land = 21,sects_vector = NULL,direction = "a")
#' 
#' plot3d(morphomapArray2matrix(Shape$"out3D"),type="s",radius = 0.1,aspect=FALSE,
#' xlab="x",ylab="y",zlab="z")
#' plot3d(morphomapArray2matrix(Shape$"inn3D"),type="s",radius = 0.1,aspect=FALSE,
#' add=TRUE)
#' wire3d(AlignMeshes$almesh2$mesh,col="white")
#' wire3d(AlignMeshes$almesh1$mesh,col="violet")
#' 
#' #Unrolling the rooth
#' AlignMeshes<-ToothAlignment(mesh1=Enamel,mesh2=Dentin,mesh3=Pulp,set,outline,
#' analyse = "r")
#' #Virtual sectioning dentine-pulp
#' External<-AlignMeshes$almesh2$mesh
#' Internal<-AlignMeshes$almesh3$mesh
#' #Define 16 cross-sections from the 10% to the 50% along the root
#' Core<-ToothCore(External,Internal,num.points = 1000,num.sect =16,
#'                 bio.len = AlignMeshes$length,start=0.1,end=0.5)
#' Shape<-ToothShape(Core,num.land = 21,sects_vector = NULL,direction = "a")
#' 
#' plot3d(morphomapArray2matrix(Shape$"out3D"),type="s",radius = 0.1,aspect=FALSE,
#' xlab="x",ylab="y",zlab="z")
#' plot3d(morphomapArray2matrix(Shape$"inn3D"),type="s",radius = 0.1,aspect=FALSE,
#' add=TRUE)
#' wire3d(AlignMeshes$almesh3$mesh,col="red")
#' wire3d(AlignMeshes$almesh2$mesh,col="lightblue")
#' wire3d(AlignMeshes$almesh1$mesh,col="white")
#' }
#' @export

ToothShape<-function (tooth.core, num.land, sects_vector, cent.out = "E", 
                      delta = 0.1,direction="c",out.sur = NULL, inn.sur = NULL) {
  
 
  mech.len <- tooth.core$mech_length
  start <- tooth.core$start
  end <- tooth.core$end
  if (is.null(sects_vector) == TRUE) {
    sects_vector <- 1:dim(tooth.core$"out3D")[3]
  }
  num.sect <- length(sects_vector)
  sect_poi <- seq(mech.len * start, mech.len * end, length = num.sect)
  out_coo_3D <- array(NA, dim = c(num.land, 3, num.sect))
  inn_coo_3D <- array(NA, dim = c(num.land, 3, num.sect))
  out_coo_2D <- array(NA, dim = c(num.land, 2, num.sect))
  inn_coo_2D <- array(NA, dim = c(num.land, 2, num.sect))
  ALPM_i_tot <- array(NA, dim = c(4, 2, num.sect))
  ALPM_o_tot <- array(NA, dim = c(4, 2, num.sect))
  centroids<- array(NA, dim = c(1, 2, num.sect))
  for (m in sects_vector) {
    # print(m)
    out_outline <- tooth.core$"out2D"[, , m]
    inn_outline <- tooth.core$"inn2D"[, , m]
    if (cent.out == "CCA") {
      centroid <- morphomapCentroid(out_outline, inn_outline, 
                                    delta)
      centroids[,,m]<-centroid
      ids_o <- ToothRegradius(out_outline, centroid, 
                              num.land,direction=direction)
      ALPM_o <- ToothRegradius(out_outline, centroid, 
                               4,direction=direction)
      
      fho <- ids_o[which.min(abs(ids_o - ALPM_o[1])):(which.min(abs(ids_o - ALPM_o[3])) - 1)]
      sho <- ids_o[which.min(abs(ids_o - ALPM_o[3])):length(ids_o)]
      
      ids_i <- ToothRegradius(inn_outline, centroid, num.land,direction=direction)
      ALPM_i <- ToothRegradius(inn_outline, centroid, 4,direction=direction)
      
      
      fhi <- ids_i[which.min(abs(ids_i - ALPM_i[1])):(which.min(abs(ids_i - ALPM_i[3])) - 1)]
      shi <- ids_i[which.min(abs(ids_i - ALPM_i[3])):length(ids_i)]
      
      out_coo_2D[, , m] <- out_outline[c(sho,fho), ]
      out_coo_3D[, , m] <- cbind(out_coo_2D[, , m],sect_poi[m])
      
        out_coo_3D[,,m]<-equidistantCurve(out_coo_3D[, , m],n = num.land,open = FALSE)
        if(!is.null(out.sur)){
          out_p<-t(projRead(out_coo_3D[,,m],out.sur)$vb)[,c(1:3)]
          out_coo_2D[, , m] <- out_p[,c(1,2)]
          out_coo_3D[, , m] <- cbind(out_p[,c(1,2)],sect_poi[m])
        }
      
      
      ALPM_o_tot[, , m] <- out_outline[ALPM_o[c(1, 2, 3, 4)], ]
      
      inn_coo_2D[, , m] <- inn_outline[c(shi, fhi), ]
      inn_coo_3D[, , m] <- cbind(inn_coo_2D[, , m], sect_poi[m])
      
        inn_coo_3D[,,m]<-equidistantCurve(inn_coo_3D[, , m],n = num.land,open = FALSE)
        if(!is.null(inn.sur)){
          inn_p<-t(projRead(inn_coo_3D[,,m],inn.sur)$vb)[,c(1:3)]
          inn_coo_2D[, , m] <- inn_p[,c(1,2)]
          inn_coo_3D[, , m] <- cbind(inn_p[,c(1,2)],sect_poi[m])
        }
      
      
      ALPM_i_tot[, , m] <- inn_outline[ALPM_i[c(1, 2, 3, 4)], ]
    }
    if (cent.out == "E") {
      centroid <- colMeans(out_outline)
      centroids[,,m]<-centroid
      ids_o <- ToothRegradius(out_outline, centroid, num.land,direction=direction)
      ALPM_o <- ToothRegradius(out_outline, centroid, 4,direction=direction)
      fho <- ids_o[which.min(abs(ids_o - ALPM_o[1])):(which.min(abs(ids_o - ALPM_o[3])) - 1)]
      sho <- ids_o[which.min(abs(ids_o - ALPM_o[3])):length(ids_o)]
      out_coo_2D[, , m] <- out_outline[c(sho,fho), ]
      out_coo_3D[, , m] <- cbind(out_coo_2D[, , m],sect_poi[m])
      
        out_coo_3D[,,m]<-equidistantCurve(out_coo_3D[, , m],n = num.land,open = FALSE)
        if(!is.null(out.sur)){
          out_p<-t(projRead(out_coo_3D[,,m],out.sur)$vb)[,c(1:3)]
          out_coo_2D[, , m] <- out_p[,c(1,2)]
          out_coo_3D[, , m] <- cbind(out_p[,c(1,2)],sect_poi[m])
        }
      
      
      ALPM_o_tot[, , m] <- out_outline[ALPM_o[c(1, 2, 3, 4)], ]
      
      
      
      ids_i <- ToothRegradius(inn_outline, centroid, num.land,direction=direction)
      ALPM_i <- ToothRegradius(inn_outline, centroid, 4,direction=direction)
         
      fhi <- ids_i[which.min(abs(ids_i - ALPM_i[1])):(which.min(abs(ids_i - ALPM_i[3])) - 1)]
      shi <- ids_i[which.min(abs(ids_i - ALPM_i[3])):length(ids_i)]
      
      inn_coo_2D[, , m] <- inn_outline[c(shi, fhi), ]
      inn_coo_3D[, , m] <- cbind(inn_coo_2D[, , m], sect_poi[m])
      
        inn_coo_3D[,,m]<-equidistantCurve(inn_coo_3D[, , m],n = num.land,open = FALSE)
        if(!is.null(inn.sur)){
          inn_p<-t(projRead(inn_coo_3D[,,m],inn.sur)$vb)[,c(1:3)]
          inn_coo_2D[, , m] <- inn_p[,c(1,2)]
          inn_coo_3D[, , m] <- cbind(inn_p[,c(1,2)],sect_poi[m])
        }
      
      
      ALPM_i_tot[, , m] <- inn_outline[ALPM_i[c(1, 2, 3, 4)], ]
      
      
    }
    if (cent.out == "I") {
      centroid <- colMeans(inn_outline)
      centroids[,,m]<-centroid
      ids_o <- ToothRegradius(out_outline, centroid, num.land,direction=direction)
      ALPM_o <- ToothRegradius(out_outline, centroid, 4,direction=direction)
      
      fho <- ids_o[which.min(abs(ids_o - ALPM_o[1])):(which.min(abs(ids_o - ALPM_o[3])) - 1)]
      sho <- ids_o[which.min(abs(ids_o - ALPM_o[3])):length(ids_o)]
      
      
      ALPM_o_tot[, , m] <- out_outline[ALPM_o, ]
      out_coo_2D[, , m] <- out_outline[c(sho,fho), ]
      out_coo_3D[, , m] <- cbind(out_coo_2D[, , m],sect_poi[m])
      
        out_coo_3D[,,m]<-equidistantCurve(out_coo_3D[, , m],n = num.land,open = FALSE)
        if(!is.null(out.sur)){
          out_p<-t(projRead(out_coo_3D[,,m],out.sur)$vb)[,c(1:3)]
          out_coo_2D[, , m] <- out_p[,c(1,2)]
          out_coo_3D[, , m] <- cbind(out_p[,c(1,2)],sect_poi[m])
        }
      
      
      ids_i <- ToothRegradius(inn_outline, centroid, 
                              num.land,direction=direction)
      ALPM_i <- ToothRegradius(inn_outline, centroid, 4,direction=direction)
      ALPM_i_tot[, , m] <- inn_outline[ALPM_i, ]
      
      fhi <- ids_i[which.min(abs(ids_i - ALPM_i[1])):(which.min(abs(ids_i - ALPM_i[3])) - 1)]
      shi <- ids_i[which.min(abs(ids_i - ALPM_i[3])):length(ids_i)]
      
      inn_coo_2D[, , m] <- inn_outline[c(shi, fhi), ]
      inn_coo_3D[, , m] <- cbind(inn_coo_2D[, , m], sect_poi[m])
      
        inn_coo_3D[,,m]<-equidistantCurve(inn_coo_3D[, , m],n = num.land,open = FALSE)
        if(!is.null(inn.sur)){
          inn_p<-t(projRead(inn_coo_3D[,,m],inn.sur)$vb)[,c(1:3)]
          inn_coo_2D[, , m] <- inn_p[,c(1,2)]
          inn_coo_3D[, , m] <- cbind(inn_p[,c(1,2)],sect_poi[m])
      }
    }
  }
  
  
  
  
  out <- list("out3D" = out_coo_3D, "inn3D" = inn_coo_3D, 
              "out2D" = out_coo_2D, "inn2D" = inn_coo_2D, 
              ALPM_inn = ALPM_i_tot, ALPM_out = ALPM_o_tot, start = start, 
              end = end, mech.len = mech.len,centroids=centroids)
  return(out)
}
