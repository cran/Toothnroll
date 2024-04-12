#' ToothExport
#'
#' Export the output from ToothAlignement 
#' @param input list: output from ToothAlignement
#' @param id character: label name
#' @param file character: name the output file
#' @return no return value, called for side effects (see description)
#' @author Antonio Profico; Mathilde Augoyard
#' @export
#' 
ToothExport <- function(input,id,file){
  
  
  cat(paste("# Toothnroll ","version ",packageVersion("Toothnroll"), "\n", 
            sep = ""), file = file, append = FALSE, sep = "")
  cat(paste("# ",id, "\n", 
            sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("# ","List of contents", "\r\n", 
            sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("almesh1 ", "{ 3D mesh vb}", " @1", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("almesh1 ", "{ 3D mesh it}", " @2", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("almesh2 ", "{ 3D mesh vb}", " @3", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("almesh2 ", "{ 3D mesh it}", " @4", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("almesh3 ", "{ 3D mesh vb}", " @5", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("almesh3 ", "{ 3D mesh it}", " @6", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("alset ", "{ 3D coordinates }", " @7", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("length ", "{ numeric }", " @8", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("margins ", "{ 3D coordinates }", " @9", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("margins_sel ", "{ integer }", " @10", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("aloutline ", "{ 3D coordinates }", " @11", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("diamBL ", "{ numeric }", " @12", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("diamMD ", "{ numeric }", " @13", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("### Data ###", "\n", sep = ""), file = file, append = TRUE, sep = "")
  if(!is.null(input$almesh1)){
    vb1<-vert2points(input$almesh1$mesh)
    it1<-t(input$almesh1$mesh$it)
  }else {vb1<-NULL;it1<-NULL}
  if(!is.null(input$almesh2)){
    vb2<-vert2points(input$almesh2$mesh)
    it2<-t(input$almesh2$mesh$it)
  }else {vb2<-NULL;it2<-NULL}
  if(!is.null(input$almesh3)){
    vb3<-vert2points(input$almesh3$mesh)
    it3<-t(input$almesh3$mesh$it)
  }else {vb3<-NULL;it3<-NULL}
  
  cat(paste("@1", "\n", sep = ""), file = file, append = TRUE,sep = "")
  write.table(format(vb1, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@2", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(it1, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@3", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(vb2, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@4", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(it2, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@5", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(vb3, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@6", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(it3, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@7", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(input$alset, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@8", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste(input$length, "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@9", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(input$margins, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@10", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(t(input$margins_sel), scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\r\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@11", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(input$aloutline, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@12", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(t(input$diamBL), scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\r\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@13", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(t(input$diamMD), scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\r\n")
}
