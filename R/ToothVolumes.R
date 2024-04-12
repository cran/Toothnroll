#' ToothVolumes
#'
#' Extract volumes from the object ToothShape
#' @param ShapeExt 3D mesh: external mesh
#' @param ShapeInn 3D mesh: internal mesh
#' @param col1 color of the ShapeExt
#' @param col2 color of the ShapeInn
#' @param col3 color of the boolean operation between ShapeExt and ShapeInn
#' @param alpha1 value to set trasparancy of col1
#' @param alpha2 value to set trasparancy of col1
#' @param alpha3 value to set trasparancy of col1
#' @param plot logical: if TRUE the volumes are shown
#' @return meshOut: external selected mesh 
#' @return meshInnT: internal selected mesh 
#' @return meshDiff: differences between selected meshes
#' @return volumeT: volume of the external mesh
#' @return volInn:volume of the internal mesh
#' @return volDiff: difference between volumeT and volInn  
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
#' 
#' #Unrolling the crown
#' AlignMeshes<-ToothAlignment(mesh1=Enamel,mesh2=Dentin,mesh3=Pulp,set,outline,analyse = "c")
#' #Virtual sectioning dentine-pulp
#' External<-AlignMeshes$almesh1$mesh
#' Internal<-AlignMeshes$almesh2$mesh
#' #Define 16 cross-sections from the 30% to the 90% along the crown
#' Core<-ToothCore(External,Internal,num.points = 1000,num.sect =16,
#'                 bio.len = AlignMeshes$length,start=0.3,end=0.9)
#' Shape<-ToothShape(Core,num.land = 100,sects_vector = NULL,direction = "a")
#' volumes<-ToothVolumes(Shape$"out3D",Shape$"inn3D",plot=TRUE)
#' unlist(volumes[4:6])
#' }
#' @export

ToothVolumes<-function(ShapeExt,ShapeInn,col1="gray",col2="red",
                       col3="green",alpha1=1,alpha2=1,alpha3=1,
                       plot=FALSE){
  
  startOut<-morphomapTri2sects(ShapeExt[,,1],repmat(colMeans(ShapeExt[,,1]),dim(ShapeExt)[1],3))
  endOut<-morphomapTri2sects(ShapeExt[,,dim(ShapeExt)[3]],repmat(colMeans(ShapeExt[,,dim(ShapeExt)[3]]),dim(ShapeExt)[1],3))
  meshSOut<-list("vb"=t(cbind(startOut$matrix,1)),"it"=t(startOut$tri))
  meshEOut<-list("vb"=t(cbind(endOut$matrix,1)),"it"=t(endOut$tri))
  meshEOut$it<-meshEOut$it[c(3,2,1),]
  
  class(meshSOut)<-"mesh3d"
  class(meshEOut)<-"mesh3d"
  
  startInn<-morphomapTri2sects(ShapeInn[,,1],repmat(colMeans(ShapeInn[,,1]),dim(ShapeInn)[1],3))
  endInn<-morphomapTri2sects(ShapeInn[,,dim(ShapeInn)[3]],repmat(colMeans(ShapeInn[,,dim(ShapeInn)[3]]),dim(ShapeInn)[1],3))
  meshSInn<-list("vb"=t(cbind(startInn$matrix,1)),"it"=t(startInn$tri))
  meshEInn<-list("vb"=t(cbind(endInn$matrix,1)),"it"=t(endInn$tri))
  meshEInn$it<-meshEInn$it[c(3,2,1),]
  class(meshEInn)<-"mesh3d"
  class(meshSInn)<-"mesh3d"
  
  start<-morphomapTri2sects(ShapeExt[,,1],ShapeInn[,,1])
  end<-morphomapTri2sects(ShapeExt[,,dim(ShapeExt)[3]],ShapeInn[,,dim(ShapeInn)[3]])
  meshS<-list("vb"=t(cbind(start$matrix,1)),"it"=t(start$tri))
  meshE<-list("vb"=t(cbind(end$matrix,1)),"it"=t(end$tri))
  meshE$it<-meshE$it[c(3,2,1),]
  class(meshS)<-"mesh3d"
  class(meshE)<-"mesh3d"
  
  CilExt<-morphomapTriangulate(morphomapArray2matrix(ShapeExt),dim(ShapeExt)[3],close=F)
  CilInt<-morphomapTriangulate(morphomapArray2matrix(ShapeInn),dim(ShapeExt)[3],close=F)
  CilExt$it<-CilExt$it[c(3,2,1),]
  CilInt$it<-CilInt$it[c(3,2,1),]
  
  volumeT<-suppressWarnings(vcgVolume(CilExt))
  volumeI<-suppressWarnings(vcgVolume(CilInt))
  volumeDiff<-volumeT-volumeI
  
  meshOutT<-mergeMeshes(meshSOut,meshEOut,CilExt)
  meshInnT<-mergeMeshes(meshSInn,meshEInn,CilInt)
  CilInt_t<-CilInt
  CilInt_t$it<-CilInt_t$it[c(3,2,1),]
  
  meshDiffT<-mergeMeshes(meshS,meshE,CilExt,CilInt_t)
  hidmat<-as.matrix(meshcube(meshDiffT))
  
  
  if(isTRUE(plot)){
    layout3d(t(c(1:3)),sharedMouse = T)
    spheres3d(hidmat,radius=0)
    triangles3d(t(meshOutT$vb[1:3, meshOutT$it]), col = col2,alpha=alpha2)
    next3d()
    spheres3d(hidmat,radius=0)
    triangles3d(t(meshInnT$vb[1:3, meshInnT$it]), col = col3,alpha=alpha3)
    next3d()
    spheres3d(hidmat,radius=0)
    triangles3d(t(meshS$vb[1:3, meshS$it]), col = col1,alpha=alpha1)
    triangles3d(t(meshE$vb[1:3, meshE$it]), col = col1,alpha=alpha1)
    triangles3d(t(CilExt$vb[1:3, CilExt$it]), col = col2,alpha=alpha2)
    triangles3d(t(CilInt$vb[1:3, CilInt$it]), col = col3,alpha=alpha3)
  }
  
  out<-list("meshOut"=meshOutT,"meshInn"=meshInnT,"meshDiff"=meshDiffT,"volOut"=volumeT,"volInn"=volumeI,"volDiff"=volumeDiff)
  return(out)
}
