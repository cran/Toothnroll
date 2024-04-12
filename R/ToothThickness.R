#' ToothThickness
#' 
#' Tool for the extraction of equiangular landmarks on the selected dental mesh 
#' @param tooth.shape list: tooth.shape object 
#' @return sect_thickness_eq dental thickness (equiangular method)
#' @return sect_thickness_cp dental thickness (closest point method)
#' @return sect_Rthickness_eq relative dental thickness (equiangular method)
#' @return sect_Rthickness_cp relative dental thickness (closest point method)
#' @return ALPM_thickness dental thickness at ALPM quadrants
#' @return tooth.shape tooth.shape object
#' @author Antonio Profico; Mathilde Augoyard
#' @export

ToothThickness<-function (tooth.shape) 
{

  num.land <- dim(tooth.shape$"out3D")[1]
  num.sect <- dim(tooth.shape$"out3D")[3]
  dist_sections_eq <- array(NA, dim = c(num.land, 1, num.sect))
  dist_sections_cp <- array(NA, dim = c(num.land, 1, num.sect))
  dist_sections_Req <- array(NA, dim = c(num.land, 1, num.sect))
  dist_sections_Rcp <- array(NA, dim = c(num.land, 1, num.sect))
  dist_ALPM <- array(NA, dim = c(4, 1, num.sect))
  
  for (m in 1:num.sect) {
    dist_section <- sqrt(rowSums((tooth.shape$"out2D"[, , m] - tooth.shape$"inn2D"[, , m])^2))
    closedPs <- vcgKDtree(tooth.shape$"inn2D"[, , m], tooth.shape$"out2D"[, , m], k = 1)$distance
    disti<-sqrt(rowSums((tooth.shape$"out2D"[, , m]-repmat(tooth.shape$centroids[,,m],dim(tooth.shape$"out2D"[, , m])[1],2))^2))
    dist_sections_cp[, , m] <- closedPs
    dist_sections_eq[, , m] <- dist_section
    dist_sections_Rcp[, , m] <- closedPs/disti
    dist_sections_Req[, , m] <- dist_section/disti
  }
  
  for (m in 1:num.sect) {
    dist_quadr <- sqrt(rowSums((tooth.shape$ALPM_out[, , m] - tooth.shape$ALPM_inn[, , m])^2))
    dist_ALPM[, , m] <- dist_quadr
  }
  out <- list(sect_thickness_eq = dist_sections_eq, sect_thickness_cp = dist_sections_cp,
              sect_Rthickness_eq = dist_sections_Req, sect_Rthickness_cp = dist_sections_Rcp,
              ALPM_thickness = dist_ALPM, tooth.shape = tooth.shape)
  return(out)
}