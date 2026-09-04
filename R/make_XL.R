
makeXL = function(xlfile,         # Blank excel file
                  tofile,
                  spec="A stock",
                  L50 = NA,     # Cruz-Cástan et al., 2019, #
                  Linf = NA,   # da Silva et al. 2024
                  M = NA,      # El-Haweet et al. (2013)
                  K = NA,      # da Silva et al. 2024
                  C_T = NA,
                  CAL_T = NA,
                  CAL_mids = NA,
                  I_T = NA,
                  L5_T = NA,
                  LFS_T = NA,
                  VML = NA,
                  yrs = NA,
                  plot = T){

  pars_sheet = read_xlsx(xlfile, 1)
  ts_sheet = read_xlsx(xlfile, 2)
  np = nrow(pars_sheet)

  # pars sheet
  maxa = -log(0.05)/M
  dats = round(c(maxa, K, Linf, L50/Linf, L5_T/L50, LFS_T/L50, VML),3)
  pars_sheet[,2] = c(spec,dats,rep(NA,np-length(dats)-1))

  # time series sheet
  years = yrs
  ny = length(yrs)
  tsmat = array(NA,c(ny,25))
  tsmat[,1] = years
  tsmat[,2] = C_T
  CALi = CAL_T[1,,]
  mids = CAL_mids
  tsmat[,3] = apply(CALi*rep(mids,each=ny),1,sum)/apply(CALi,1,sum)
  tsmat[,5] = sapply(1:dim(CALi)[1],function(X,CALi,mids){
    sd(rep(mids,CALi[X,]))/mean(rep(mids,CALi[X,]))
    },CALi=CALi,mids=mids)

  tsmat[,6] =  sapply(1:dim(CALi)[1],function(X,CALi,mids,L50){
    indivs = rep(mids,floor(CALi[X,]*1E3)) # we are trying to calc means weighted by nsamp - these can be less than 1 and hence need upscaling by 1E3 to make sure
    mean(indivs > L50,na.rm=T)
    }, CALi=CALi, mids=mids, L50=L50)

  tsmat[,7] = I_T

  tsmat = round(tsmat,3)
  tsmat=as.data.frame(tsmat)
  names(tsmat) = names(ts_sheet)


  writexl::write_xlsx(list(Parameters=pars_sheet, Time_series = tsmat),path=tofile)
  cat(paste0("XLfile written: ",tofile, "\n"))
  if(plot){

    Cs = !all(is.na(tsmat[,2]))
    MLs = !all(is.na(tsmat[,3]))
    MAs = !all(is.na(tsmat[,4]))
    MVs = !all(is.na(tsmat[,5]))
    FMs = !all(is.na(tsmat[,6]))
    Is = !all(is.na(tsmat[,7]))
    npanel = sum(Cs,MLs,MAs, MVs,FMs,Is)
    if(npanel>1){
      yrs = tsmat[,1]
      ncol = ceiling(npanel^0.5)
      nrow = ceiling(npanel/ncol)
      par(mfrow=c(nrow,ncol),mai=c(0.35,0.7,0.05,0.05),omi=c(0.3,0.05,0.3,0.05))
      if(Cs)plot(yrs, tsmat[,2], pch=19, ylab="Catch",ylim=c(0,max(tsmat[,2],na.rm=T))); grid()
      if(MLs)plot(yrs, tsmat[,3], pch=19, ylab="Mean Length",ylim=c(0,max(tsmat[,3],na.rm=T))); grid()
      if(MAs)plot(yrs, tsmat[,4], pch=19, ylab="Mean Age",ylim=c(0,max(tsmat[,4],na.rm=T))); grid()
      if(MVs)plot(yrs, tsmat[,5], pch=19, ylab="CV Lengths",ylim=c(0,max(tsmat[,5],na.rm=T))); grid()
      if(FMs)plot(yrs, tsmat[,6], pch=19, ylab="Fraction Mature",ylim=c(0,max(tsmat[,6],na.rm=T))); grid()
      if(Is)plot(yrs, tsmat[,7], pch=19, ylab="Index of Rel. Abund.",ylim=c(0,max(tsmat[,7],na.rm=T))); grid()
      mtext(spec, 3, line=0.3, font=2, outer=T)
      mtext("Year", 1, line=0.3, outer=T)
    }
  }

}
