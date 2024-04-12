#' ToothCore
#'
#' Tool to build 3D and 2D cross sections 
#' @param out.sur object of class mesh3d
#' @param inn.sur object of class mesh3d
#' @param num.sect number of sections
#' @param bio.len length of the selected region of interest along with extracting the digital section
#' @param clean_int_out_O logical if TRUE the outer section will be cleaned by using spherical flipping
#' @param clean_int_out_I logical if TRUE the inner section will be cleaned by using spherical flipping
#' @param param1_out numeric parameter for clean_int_out_O spherical flipping operator (how much the section will be deformed)
#' @param radius.fact_out logical if TRUE the inner section will be cleaned by using spherical flipping
#' @param param1_inn numeric parameter for clean_int_out_I spherical flipping operator (how much the section will be deformed)
#' @param radius.fact_inn logical if TRUE the inner section will be cleaned by using spherical flipping
#' @param npovs numeric: number of points of view defined around the section
#' @param num.points number of equiengular points to be defined on each section
#' @param start percentage of the mechanical length from which the first section is defined
#' @param end percentage of the mechanical length from which the last section is defined
#' @param curved logical: if TRUE the cutting planes will follow the mesh curvature (beta version) 
#' @param print.progress logical: if TRUE a progress bar is printed to the screen
#' @return out3D num.pointsx3xnum.sect array of the external outlines 
#' @return inn3D num.pointsx3xnum.sect array of the internal outlines
#' @return out3D num.pointsx2xnum.sect array of the external outlines
#' @return inn3D num.pointsx2xnum.sect array of the internal outlines
#' @return mech_length mechanical length of the long bone
#' @return start percentage of the mechanical length from which the first section is defined
#' @return end percentage of the mechanical length from which the last section is defined
#' @author Antonio Profico; Mathilde Augoyard
#' @examples
#' \donttest{
#' data("URI1_tooth")
#' require(morphomap)
#' require(Morpho)
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
#' 
#' plot3d(morphomapArray2matrix(Core$"out3D"),aspect=FALSE,xlab="x",ylab="y",zlab="z")
#' plot3d(morphomapArray2matrix(Core$"inn3D"),aspect=FALSE,add=TRUE)
#' wire3d(AlignMeshes$almesh2$mesh,col="white")
#' wire3d(AlignMeshes$almesh1$mesh,col="violet")
#' 
#' #Unrolling the rooth
#' AlignMeshes<-ToothAlignment(mesh1=Enamel,mesh2=Dentin,mesh3=Pulp,set,outline,analyse = "r")
#' #Virtual sectioning dentine-pulp
#' External<-AlignMeshes$almesh2$mesh
#' Internal<-AlignMeshes$almesh3$mesh
#' #Define 16 cross-sections from the 10% to the 50% along the root
#' Core<-ToothCore(External,Internal,num.points = 1000,num.sect =16,
#'                 bio.len = AlignMeshes$length,start=0.1,end=0.5)
#' 
#' plot3d(morphomapArray2matrix(Core$"out3D"),aspect=FALSE,xlab="x",ylab="y",zlab="z")
#' plot3d(morphomapArray2matrix(Core$"inn3D"),aspect=FALSE,add=TRUE)
#' wire3d(AlignMeshes$almesh3$mesh,col="red")
#' wire3d(AlignMeshes$almesh2$mesh,col="lightblue")
#' wire3d(AlignMeshes$almesh1$mesh,col="white")
#' }
#' @export

ToothCore<-function (out.sur = out.sur, inn.sur = inn.sur, num.sect = 61, 
                     bio.len, clean_int_out_O = TRUE,clean_int_out_I = TRUE, param1_out = 0.5, radius.fact_out = 2.5,
                     param1_inn = 0.5, radius.fact_inn = 2.5, 
                     npovs = 100, num.points = 500, start = 0.2, end = 0.8, curved=FALSE,print.progress = TRUE) 
{
  
  ext_raw_sects_2D <- array(NA, dim = c(num.points, 2, num.sect))
  inn_raw_sects_2D <- array(NA, dim = c(num.points, 2, num.sect))
  ext_raw_sects_3D <- array(NA, dim = c(num.points, 3, num.sect))
  inn_raw_sects_3D <- array(NA, dim = c(num.points, 3, num.sect))
  sect_poi <- seq(bio.len * start, bio.len * end, length = num.sect)
  if (print.progress == TRUE) {
    pb <- txtProgressBar(min = 0, max = num.sect - 1, initial = 0, 
                         style = 3)
    step <- 0
  }
  
  if(isTRUE(curved)){
    core_temp<-ToothCore(out.sur,inn.sur,num.sect = num.sect+1,bio.len = bio.len,start=start,end=end)
    shapi<-ToothShape(core_temp,num.land=48,sects_vector = NULL,cent.out = "CCA",delta = 0.05)
    
    matArray<-matrix(NA,nrow=num.sect+1,ncol=3)
    for(j in 1:(num.sect+1)){
      matArray[j,]<-c(morphomapCentroid(shapi$"3D_out"[,,j],
                                        shapi$"3D_inn"[,,j],delta = 0.05),unique(core_temp$"out3D"[,3,j]))
      close3d()
    }
  }
  
  for (m in 1:num.sect) {
    
    if(isTRUE(curved)){
      sect_poi <- seq(bio.len * start, bio.len * end, length = num.sect)[m]
      Diff<-diff(seq(bio.len * start, bio.len * end, length = num.sect))[1]
      p1 <- c(0, 0, sect_poi)
      p2 <- c(10, 0, sect_poi)
      p3 <- c(0, 10, sect_poi)
      yincrement <- Diff
      xincrement <- matArray[m+1,1]-matArray[m,1]
      slope <- yincrement/xincrement
      ANGLE<-atan(slope)+pi/2
      if(ANGLE>pi/2){
        ANGLE<-ANGLE-pi
      }
      baryt<-colMeans(rbind(p1,p2,p3))
      barytt<-matArray[m,]-baryt
      p1<-p1+barytt
      p2<-p2+barytt
      p3<-p3+barytt
      p1b<-rotate3d(p1,ANGLE,0,1,0)
      p2b<-rotate3d(p2,ANGLE,0,1,0)
      p3b<-rotate3d(p3,ANGLE,0,1,0)
      baryt<-colMeans(rbind(p1b,p2b,p3b))
      barytt<-matArray[m,]-baryt
      p1b<-p1b+barytt
      p2b<-p2b+barytt
      p3b<-p3b+barytt
      p1<-p1b
      p2<-p2b
      p3<-p3b
    }else{
      p1 <- c(0, 0, sect_poi[m])
      p2 <- c(10, 0, sect_poi[m])
      p3 <- c(0, 10, sect_poi[m])
    }
    
    inter_out <- NULL
    inter_inn <- NULL
    
    
    normal <- crossProduct(p2-p1 ,p3-p1)
    normal <- crossProduct(p2-p1 ,p3-p1)
    
    mesh<-out.sur
    v1<-p1
    v2 = p2
    v3 = p3
    
    pointcloud <- vert2points(mesh)
    updown <- cutSpace(pointcloud, v1 = v1, v2 = v2, v3 = v3,normal = normal)
    upface <- getFaces(mesh, which(updown))
    downface <- getFaces(mesh, which(!updown))
    nit <- 1:ncol(mesh$it)
    facesInter <- which(as.logical((nit %in% upface) * nit %in% downface))
    mesh$it <- mesh$it[, facesInter]
    mesh <- rmUnrefVertex(mesh, silent = TRUE)
    edges <- as.matrix(vcgGetEdge(mesh)[, 1:2])
    pointcloud <- vert2points(mesh)
    zeroPro <- points2plane(rep(0,3),p1,normal)
    sig <- sign(crossprod(-zeroPro,normal))
    d <- sig*norm(zeroPro,"2")
    inter_out <- pointcloud[, c(1, 2)]
    mesh<-inn.sur
    v1<-p1
    v2 = p2
    v3 = p3
    normal = normal
    pointcloud <- vert2points(mesh)
    updown <- cutSpace(pointcloud, v1 = v1, v2 = v2, v3 = v3,  normal = normal)
    upface <- getFaces(mesh, which(updown))
    downface <- getFaces(mesh, which(!updown))
    nit <- 1:ncol(mesh$it)
    facesInter <- which(as.logical((nit %in% upface) * nit %in%  downface))
    mesh$it <- mesh$it[, facesInter]
    mesh <- rmUnrefVertex(mesh, silent = TRUE)
    edges <- as.matrix(vcgGetEdge(mesh)[, 1:2])
    pointcloud <- vert2points(mesh)
    inter_inn <- pointcloud[, c(1, 2)]
    
    inters <- inter_inn
    if (clean_int_out_O == TRUE) {
      inter_out <- morphomapFlip(inter_out, param1 = param1_out, 
                                 radius.fact = radius.fact_out, npovs = npovs)
    }
    ordered_out_temp <- morphomapSort(inter_out)
    ordered_out_temp <- rbind(ordered_out_temp, ordered_out_temp[1, 
    ])
    ordered_out <- ordered_out_temp
    
    
    if (clean_int_out_I == TRUE) {
      inter_inn <- morphomapFlip(inter_inn, param1 = param1_inn, 
                                 radius.fact = radius.fact_inn, npovs = npovs)
    }
    
    ordered_inn_temp <- morphomapSort(inter_inn)
    ordered_inn_temp <- rbind(ordered_inn_temp, ordered_inn_temp[1, 
    ])
    ordered_inn <- ordered_inn_temp
    
    ev_out <- equidistantCurve(ordered_out, n = num.points, 
                               iterations = 1, increment = 0)
    ev_inn <- equidistantCurve(ordered_inn, n = num.points, 
                               iterations = 1, increment = 0)
    ext_raw_sects_2D[, , m] <- ev_out
    inn_raw_sects_2D[, , m] <- ev_inn
    ext_raw_sects_3D[, , m] <- cbind(ev_out, sect_poi[m])
    inn_raw_sects_3D[, , m] <- cbind(ev_inn, sect_poi[m])
    if (print.progress == TRUE) {
      setTxtProgressBar(pb, step)
      step <- step + 1
    }
  }
  if (print.progress == TRUE) {
    close(pb)
  }
  out <- list("out3D" = ext_raw_sects_3D, "inn3D" = inn_raw_sects_3D, 
              "out2D" = ext_raw_sects_2D, "inn2D" = inn_raw_sects_2D, 
              mech_length = bio.len, start = start, end = end)
  return(out)
}
