# Indicator 3 functions



prepare_TD = function(sets = NA){
  if(is.na(sets)[1])sets = 1:11
  TDs = list()
  for(i in sets)  TDs[[i]] = get(paste0("TD_",i))
  cat("Indicator 3 training dataset compiled \n")
  as.data.frame(rbindlist(TDs))
}

makerawdata_3 = function(cdat, sno=1, isBrel = F,
                         inc_Irel = T, inc_I = T, inc_CR = T, inc_CAL = T, inc_CAA = T,
                         stock_in = NA, fleet_in = NA, Brange = c(0.025,4.5)){
  # sno=1; isBrel = F; clean = T;  inc_Irel = T; inc_I = T;  inc_CR = T; stock_in = NA; fleet_in = NA; inc_CAL = T; inc_CAA = T
  # cdat = as.data.frame(rbindlist(allout))
  dnames = names(cdat)
  Brelcols = grepl("Brel",dnames)
  stocks = sapply(dnames[Brelcols],function(x)strsplit(x,"_")[[1]][2])
  if(is.na(stock_in[1]))stock_in = stocks[stocks!=sno]
  if(is.na(fleet_in[1]))fleet_in = 1:3

  cdat = stock_subsetter3(cdat, sno, stock_in) # keep only those listed in sno and stock_in
  cdat = fleet_subsetter3(cdat, fleet_in)                # keep only those fleets in nfleet
 # cdat = data_subsetter3(cdat,  inc_Irel, inc_I, inc_CR, inc_CAL, inc_CAA) # only include those data that are specified to be available
  dat = resp_subsetter3(cdat, sno, isBrel)                     # gets the right response variable according to sno
  dat = log(dat)                                  # log imperfect fractions
  #dat = dologit(dat,types = "VML")                    # logit proportions (but rescaled 0.05 - 0.95 prior to logit)
  dat = cleandat_3(dat,resprange = log(Brange))       # clean NAs and Infs
  dat = rem_const(dat)                                # remove any independent variables with no variability (constant over simulations)
  dat

}






cleandat_3 = function(dat,resprange=c(-3.689,0.693)){
  keep = dat$Res > resprange[1] & dat$Res < resprange[2]
  cat(paste0("Brange: ",sum(keep)," of ",length(keep)," records kept (",round(sum(keep)/length(keep)*100,2),"%) /n"))
  dats =  dat[keep,]
  # test = (x==0); (1:length(test))[test]
  # goodrow = function(x)(sum(is.na(x)| x == 'Inf')==0)
  goodrow = function(x)(sum(is.na(x)| is.infinite(x))==0)
  # goodrow = function(x)(sum(is.na(x))==0)

  keep = apply(dats,1,goodrow)
  cat(paste0("NAs and Infs: ",sum(keep)," of ",length(keep)," records kept (",round(sum(keep)/length(keep)*100,2),"%) /n"))
  goodrow = function(x)(sum(is.na(x))==0)

  dats[keep,]
}

stock_subsetter3 = function(cdat, sno,stock_in, maxstock = 3){
  if(length(stock_in)<3){
    dnames = names(cdat)
    stock_in_nam = paste0("s",stock_in)
    all_stock=paste0("s",1:maxstock)
    resp_stock = paste0("s",sno)
    sel_stock = unique(c(resp_stock,stock_in_nam))
    rem_stock = all_stock[!(all_stock%in%sel_stock)]
    ns = length(rem_stock)
    nn = length(dnames)
    condgrid = array(F,c(ns,nn))
    for(i in 1:ns)condgrid[i,]= grepl(rem_stock[i],dnames)
    cdat = cdat[,!apply(condgrid,2,any)] # should the input be there based on the stocks?
  }
  cdat
}

fleet_subsetter3 = function(cdat, fleet_in, dnames, maxfleet = 3){
  if(length(fleet_in)<3){
    dnames = names(cdat)
    fleet_in_nam = paste0("f",fleet_in)
    all_fleet=paste0("f",1:maxfleet)
    rem_fleet = all_fleet[!(all_fleet%in%fleet_in_nam)]
    nf = length(rem_fleet)
    nn = length(dnames)
    condgrid = array(F,c(nf,nn))
    for(i in 1:nf)condgrid[i,]= grepl(rem_fleet[i],dnames)
    cdat = cdat[,!apply(condgrid, 2, any)] # should the input be there based on the stocks?
  }
  cdat
}


data_subsetter3 = function(dat, inc_I=T, inc_CR = T, inc_CAL = T, inc_CAA = T){
  if(!inc_I)    dat = dat[,!grepl("I_",names(dat))]
  if(!inc_CR)   dat = dat[,!(grepl("CR_",names(dat)) | grepl("CC",names(dat)))]
  if(!inc_CAL)  dat = dat[,!(grepl("ML_",names(dat))  | grepl("MV",names(dat)) | grepl("FM", names(dat)))]
  if(!inc_CAA)  dat = dat[,!grepl("MA_",names(dat))]
  dat
}

resp_subsetter3 = function(cdat, sno, isBrel){
  dnames2 = names(cdat)
  Brelcols = grepl("Brel",dnames2)
  ns = sum(Brelcols) # os = (1:ns)[!((1:ns)%in%sno)]
  keepcol = as.numeric((1:ncol(cdat))[grepl(paste0("Brel_s",sno),dnames2)])
  Res = cdat[,keepcol]
  if(isBrel){Res = ((ns-1) * Res) /(apply(cdat[,Brelcols],1,sum)-Res)}# level relative to mean level of other stocks
  dat = cbind(Res,cdat[,!Brelcols])
  dat
}

rem_const = function(dat){
  isconst = apply(dat,2,sd)<1E-10
  if(sum(isconst,na.rm=T)>0){
    cat(paste(paste(names(dat)[isconst],collapse=", "), "dropped for sd < 1E-10 \n"))
    dat = dat[,!isconst]
  }
  dat
}

#dolog_3=function(dat){
 # dat=log(dat)
#}



