#' Process XL input file into an EcoTest data object of class ETData
#'
#' @param xlfile Character string. The location and name of the Excel input file.
#' @param CVs List. A named list of errors than can be applied to the various features of the Excel file.
#' @param nint Positive integer. The number of interpolations for all time series data. Leave this alone as 40 is the default of the training dataset.
#' @param nsim Positive integer. The number of draws of input features. If nsim = 1, features are subject to very small error (essentially deterministic at input values)
#' @param retros Integer. The number of retrospective peels to conduct. Strips away time series data by retros years to obtain status estimates over the most recent years.
#' @param extrapolate Boolean. Should missing data prior to the first year of data or after the last year of data be imputed with the first and last observation? Note that indicator 4 requires complete time series data so extrapolation may be necessary.
#' @param minspan Positive integer. The minimum number of years that time series need to span to be processed.
#' @author T. Carruthers
#' @export
proc_i4=function(xlfile, CVs = list(K = 0.05, maxa = 0.1, L50_Linf = 0.05, Linf = 0.025,
                 Index = 0.15, Catch = 0.1, ML = 0.1, MA = 0.1, MV = 0.05, FM = 0.01,
                 L5_L50 = 0.1, LFS_L50 = 0.1, VML = 0.01),
                 nint = 40, nsim = 100, retros = 6, extrapolate = T, minspan = 10){

  out=list()
  if(nsim==1)for(ll in 1:length(CVs)){CVs[[ll]]=1E-6} # Assign small number
  for(i in 1:(retros+1)){
    out[[i]] <- make_i4(xlfile, CVs, nint, nsim, peel=i-1, extrapolate, minspan)
  }
  names(out) = paste("peel", 0:retros, sep="_")

  # True SSB/SSBMSY
  sheets = excel_sheets(xlfile)
  if("B_BMSY" %in% sheets){
    Brel = as.data.frame(read_xlsx(xlfile,3))
    out[["Brel"]] = Brel
  }
  class(out) = "ETData"
  out
}


make_i4 = function(xlfile, CVs = list(K = 0.1, maxa = 0.1, L50_Linf = 0.1, Linf = 0.025,
                   Index = 0.2, Catch = 0.2, ML = 0.1, MA = 0.1, MV = 0.05, FM = 0.01,
                   L5_L50 = 0.1, LFS_L50 = 0.1, VML = 0.1),
                   nint = 40, nsim = 100, peel = 0, extrapolate = T, minspan = 15){

  # Internal functions
  plu = function(ps, nam) as.numeric(ps[match(nam, ps[,1]), 2])  # parameter look up

  tsresamp = function(ats, CV, nsim, peel, normalz = F){ # time series resample
    nt = length(ats);
    atp = ats[1:(nt-peel)]
    if(normalz)atp=atp/mean(atp,na.rm=T)
    np = length(atp)
    matrix(atp * trlnorm(nsim * np, 1, CV), byrow=T, nrow=nsim)
  }

  getsdsmooth = function(x, enp.mult){
    if(sum(!is.na(x))>5){
      return(sd(smooth2(x, ret = 'resid', enp.mult = enp.mult),na.rm=T))
    }else{
      return(NA)
    }
  }

  ifunc = function(inp, enp.mult, nint = 40, minspan = 15, extrapolate = T, dolog = F, min = NA, max = NA){ # function to interpolate (and potentially extrapolate) time series
    sim1 = inp[1,]
    enoughdata = sum(!is.na(sim1))>=10
    if(enoughdata){
      rng = range((1:length(sim1))[!is.na(sim1)],na.rm=T)
      enoughspan = (rng[2]-rng[1]) >= minspan
      if(enoughspan){
        out = t(apply(inp,1,interpolate, enp.mult=enp.mult, nint=nint, extrapolate=extrapolate))
        if(!is.na(min)) out[out<min] = min
        if(!is.na(max)) out[out>max] = max
        if(dolog) out = log(out)
        return(out)
      }else{
        return(NULL)
      }
    }else{
      return(NULL)
    }
  }

  ps = as.data.frame(read_xlsx(xlfile,1))  # Parameter sheet
  tss = as.data.frame(read_xlsx(xlfile,2)) # Time series sheet

  # Get all params
  K = trlnorm(nsim, plu(ps, "K"),CVs$K)
  maxa = trlnorm(nsim, plu(ps, "max_age"), CVs$maxa)
  M = log(0.05)/-maxa
  M_K = M/K
  L50_Linf = trlnorm(nsim, plu(ps,"L50_Linf"), CVs$L50_Linf)
  Linf = trlnorm(nsim, plu(ps,"Linf"), CVs$Linf)

  CF_f1 = plu(ps, "Fleet_1_CF")
  CF_f2 = plu(ps, "Fleet_2_CF")
  CF_f3 = plu(ps, "Fleet_3_CF")

  ns = data.frame(K, M_K, maxa, L50_Linf, CF_f1, CF_f2, CF_f3)

  # Time series stuff

  fnams = c("Total_","Fleet_1_", "Fleet_2_","Fleet_3_")
  flabs = c("T", "f1", "f2","f3")

  allts = NULL

  for(ff in 1:length(fnams)){

    # Load time series
    Iobs = tsresamp(tss[,paste0(fnams[ff],"Index")], CV = CVs$Index, nsim, peel = peel, normalz = T)    # Mean 1
    Cobs = tsresamp(tss[,paste0(fnams[ff],"catch")], CV = CVs$Catch, nsim, peel = peel, normalz = T)    # Mean 1
    MLobs = tsresamp(tss[,paste0(fnams[ff],"mean_len")], CV = CVs$ML, nsim, peel = peel) / Linf         # Fraction of Linf
    MAobs = tsresamp(tss[,paste0(fnams[ff],"mean_age")], CV = CVs$MA, nsim, peel = peel) / maxa         # Fraction of maxage
    MVobs = tsresamp(tss[,paste0(fnams[ff],"CV_len")], CV = CVs$MV, nsim, peel = peel)                  # Unscaled (CV)
    FMobs = tsresamp(tss[,paste0(fnams[ff],"frac_mat")], CV = CVs$FM, nsim, peel = peel)                # Unscaled frac mature

    # Variability properties
    ind = ncol(Iobs) - (19:0)
    Isd1 = apply(Iobs[,ind], 1, getsdsmooth, enp.mult = 0.4)
    Isd2 = apply(Iobs[,ind], 1, getsdsmooth, enp.mult = 0.2)
    Isd3 = apply(Iobs[,ind], 1, getsdsmooth, enp.mult = 0.1)
    Csd1 = apply(Cobs[,ind],1, getsdsmooth, enp.mult = 0.4)
    Csd2 = apply(Cobs[,ind],1, getsdsmooth, enp.mult = 0.2)
    Csd3 = apply(Cobs[,ind],1, getsdsmooth, enp.mult = 0.1)

    VMLmu = max(min(plu(ps, paste0(fnams[ff],"VML")),0.98),0.02)
    VMLsd = CVs$VML
    VMLr = rep(NA,nsim)
    if(!is.na(VMLmu)){
      VMLr = rbeta(nsim, alphaconv(VMLmu, VMLsd),betaconv(VMLmu, VMLsd))
      VMLr[VMLr<0.01] = 0.025
      VMLr[VMLr>0.99] = 0.975
    }
    df = data.frame(Isd1, Isd2, Isd3, Csd1, Csd2, Csd3,
                    L5_L50 = trlnorm(nsim, plu(ps, paste0(fnams[ff],"L5_L50")), CVs$L5_L50),
                    LFS_L50 = trlnorm(nsim, plu(ps, paste0(fnams[ff],"LFS_L50")), CVs$LFS_L50),
                    VML = VMLr)
    colnames(df) = paste0(c("Isd1","Isd2","Isd3","Csd1","Csd2","Csd3","L5_L50","LFS_L50","VML"),"_",flabs[ff])
    df = df[,!is.na(df[1,])] # remove NA columns (has to be all rows)

    # Interpolated time series
    # if(plotsmooth)  par(mfrow=c(2,3))
    # ifunc extrapolates, interpolates
    Iint = ifunc(Iobs, 0.2, nint, minspan, min = 0.05)
    Cint = ifunc(Cobs, 0.2, nint, minspan, min = 0.01)
    MLint = ifunc(MLobs, 0.125, nint, minspan, min = 0.025, max = 0.975)
    MAint = ifunc(MAobs, 0.125, nint, minspan)
    MVint = ifunc(MVobs, 0.125, nint, minspan)
    FMint = ifunc(FMobs, 0.125, nint, minspan, min = 0.025, max = 0.975)
    dimo = function(x){if(class(x)[1]=="matrix")return(ncol(x)); if(is.null(x))return(0)}
    times=c(dimo(Iint), dimo(Cint), dimo(MAint), dimo(MLint), dimo(MVint), dimo(FMint))

    if(sum(times)>0){
      labs = paste(rep(c("I", "C", "MA", "ML", "MV","FM"), times), 1:40, flabs[ff], sep="_")
      tsdat = as.data.frame(cbind(Iint,Cint,MAint, MLint, MVint, FMint))
      names(tsdat) = labs
      if(ff==1) allts = cbind(df, tsdat)
      if(ff>1) allts = cbind(allts, df, tsdat)
    }else{
      tsdat = NULL
      if(ff==1) allts = df
      if(ff>1) allts = cbind(allts, df)
    }

  }

  out = cbind(ns,allts)
  outna = out[,apply(out,2,function(x)!all(is.na(x)))] # remove any extra NA columns
  log(outna)                          # all variables are logged

}
