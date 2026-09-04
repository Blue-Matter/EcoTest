# Indicator 3 functions


prepare_TD = function(sets = NA){
  if(is.na(sets)[1])sets = 1:11
  TDs = list()
  for(i in sets)  TDs[[i]] = get(paste0("TD_",i))
  cat("Training dataset compiled for Indicators 3 and 4 \n")
  as.data.frame(rbindlist(TDs))
}



