#' Load EcoTest Models
#'
#' @param modelno (optional) An integer value representing the model to be loaded. Default is NaN and leads to the loading of all models.
#' @return all package data objects are placed in the global namespace \code{dir}
#' @examples
#' load_mod()
load_mod<-function(modelno=NaN,quiet=F){

  chk <- "package:ETest" %in% search()
  if(chk){
    datadir<- paste(searchpaths()[match("package:ETest",search())],"/data/",sep="")

    if(is.na(modelno)){
      fnams<-list.files(datadir)
      fnams<-fnams[grepl(".keras",fnams)]
      for(i in 1:length(fnams))load_model(fnams[i])
    }else{
      Model = load_model(paste0(datadir,"/Model_",modelno,".keras"))
    }

    if(!quiet)cat(paste0("Model ",modelno, " loaded"))
    cat("\n")
    return(Model)
  }else{
    stop("Can't find ETest in this session: you may have to load the library: library(ETest)")
  }

}
