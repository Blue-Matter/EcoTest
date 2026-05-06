
what_fleet = function(dir){
  replist= SS_output(dir)
  SSplotCatch(replist)
  SSplotIndices(replist)
  SSplotData(replist)
  
}

ss_names=function(io){
  io$inputs$dat$fleetinfo$fleetname
  #unique(io$inputs$dat$CPUE$index)
  #unique(io$inputs$dat$catch$fleet)
}

ss_fleet_helper = function(io){
  if(!is.null(io)){
    fleetnames =  io$inputs$dat$fleetinfo$fleetname
    nf=length(fleetnames)
    helpy = array("-",c(4,nf))
    
    ct = io$inputs$dat$catch
    allfleet = aggregate(ct$catch,by=list(fleet = ct$fleet),sum)
    helpy[1,allfleet$fleet] = round(allfleet$x/sum(allfleet$x)*100,1)
    
    ind = io$inputs$dat$CPUE
    ny = aggregate(rep(1,nrow(ind)),by=list(ind$index),sum)
    miny = aggregate(ind$year,by=list(ind$index),min)
    maxy = aggregate(ind$year,by=list(ind$index),max)
    helpy[2,ny$Group.1] = ny$x; helpy[3,miny$Group.1]=miny$x; helpy[4,maxy$Group.1]=maxy$x
    colnames(helpy) = fleetnames
    rownames(helpy) = c("Catch%","No yr","first yr","last yr")
    return(helpy)
  }else{
    return("NA")
  }
}

getmucv = function(comps,lbins, ndraw=20){
  
  nrec = nrow(comps)
  nbin = length(lbins)
  suml = t(t(comps)*lbins)
  totl = apply(suml,1,sum)
  mu = totl / apply(comps,1,sum)
  cv = rep(NA,length(mu))
  mlsamp = array(NA,c(nbin,ndraw))
  
  for(rr in 1:nrec){
    nobs = sum(comps[rr,])
    samp=rmultinom(ndraw,nobs,prob=comps[rr,]/sum(comps[rr,]))
    sums = samp*lbins
    tots = apply(sums,2,sum)
    cv[rr] = sd(tots/nobs)/mu[rr]
  }
  
  list(mu=mu,cv=cv)
}

# gets first two fleets selectivity and weighted selectivity of all other fleets (Catch in last year)
getsels = function(io,Find,dat,ploty){
  # selectivities
  selex = io$outputs$sizeselex
  sel = selex[selex$Yr==max(selex$Yr) & selex$Sex==1 & selex$Factor=="Lsel",]
  
  # you were here - need to locate fleet selectivities
  fleets = sel$Fleet
  nf = length(fleets)
  sind = 6:ncol(sel)
  slen =as.numeric(names(sel)[sind])
  sels = array(NA,c(3, length(sind)))
  for(ff in 1:2) sels[ff,] = as.numeric(sel[sel$Fleet==Find[ff],sind,drop=T])
  frest = (1:nf)[!(sel$Fleet%in%Find)]
  clast = dat$catch[dat$catch$year==max(dat$catch$year),]
  swt = rep(0,nf); swt[match(clast$fleet,fleets)]= clast$catch
  srest = apply(sel[frest,sind]*swt[frest],2,sum) 
  sels[3,]=srest # weighted by catch in last year for now
  sels=sels/apply(sels,1,max)
  if(ploty){matplot(slen,t(sels),type="l",lty=1,col=c("red","green","blue"),lwd=2); abline(h=0.05)}
  L5s = LFSs = rep(NA,3)
  for(ff in 1:3){
    ascend = 1:which.max(sels[ff,])
    L5s[ff] = approx(sels[ff,ascend],slen[ascend],0.05)$y
    LFSs[ff] = approx(sels[ff,ascend],slen[ascend],0.95)$y
  }
  VMLs = sels[,ncol(sels)]
  list(L5s = L5s,LFSs=LFSs, VMLs = VMLs)
}

alphaconv = function(m,sd){ 
  m * (((m * (1 - m))/(sd^2)) - 1)
}

betaconv=function(m,sd){
 (1 - m) * (((m * (1 - m))/(sd^2)) - 1)
}


#  dd= 5; io = ios[[dd]]; Fnam = Fnams[[dd]]; Inam = Inams[[dd]]; nsamp = 10; catch_CV = 0.1; L50_CV = 0.05; Linf_CV = 0.025; K_CV = 0.05; M_CV = 0.05; VML_CV = 0.05; peel=0; plotsmooth=T
SS_2_ET = function(io, Fnam = c("F4_JPN_LL","F8_ESP_LL"), Inam = c("S2_JPN_LATE","S4_EU_ESP"),
                   nsamp=20, catch_CV = 0.05, L50_CV = 0.05, Linf_CV = 0.025, 
                   K_CV = 0.05, M_CV = 0.05, VML_CV = 0.05, plotsmooth=F, peel=0){
  
  lastyr = io$inputs$dat$endyr - peel # Retro
  control = io$inputs$control
  Fnames = io$inputs$dat$fleetinfo$fleetname 
  Find = match(Fnam, Fnames)
  Iind = match(Inam, Fnames)
 
  # === Biological info ======================================
  MGp = control$MG_parms
  Linf_mu = MGp[(1:nrow(MGp))[grepl("L_at_Amax",rownames(MGp))][1],3]
  Linf = rep(100, nsamp) #rlnorm(nsamp, Linf_mu ,Linf_CV)
  K_mu = MGp[(1:nrow(MGp))[grepl("VonBert_K",rownames(MGp))][1],3]
  K =  rlnorm(nsamp, log(K_mu), K_CV)
  L50_mu = MGp[grepl("Mat50",rownames(MGp)),3]
  if(L50_mu==0 & "Age_Maturity" %in% names(control)){
    Mata = control$Age_Maturity
    Lata = (Linf_mu*(1-exp(-K_mu*(0:100))))[1:length(Mata)]
    L50_mu = approx(Mata,Lata,0.5)$y
  }
  
  L50 = rlnorm(nsamp, log(L50_mu/Linf_mu*100), L50_CV) # all lengths are % Linf
  L50_Linf = L50/Linf 
  
  if("natM"%in%names(control)){
    Ma = as.numeric(control$natM[1,,drop=T])
    Sa = exp(-cumsum(Ma))
    M_mu = sum(Sa*Ma)/sum(Sa) # survival weighted M
  }else{
    M_mu =  rlnorm(nsamp,log(MGp[(1:nrow(MGp))[grepl("NatM",rownames(MGp))][1],3]),M_CV)
  }
  M = rlnorm(nsamp,log(M_mu),M_CV)
  M_K = M/K
  maxa = -log(0.05)/M # age at 5% cumulative survival
  
  
  # === Fleet specific stuff ================================
  
  # get catches -----------
  dat = io$inputs$dat
  ct = dat$catch
  nc = nrow(ct)
  newind = rep(3,nc) # default to fleet 3
  for(ff in 1:2)newind[ct$fleet == Find[ff]]=ff  
  aggcatch = aggregate(ct$catch,by=list(year = ct$year, fleet = newind),sum)
  aggcatch = aggcatch[aggcatch$year > 0 & aggcatch$year<=lastyr,] # Retro
  allfleet = aggregate(ct$catch,by=list(fleet = ct$fleet),sum)
  fleetcatch = aggregate(ct$catch,by=list(fleet = newind),sum)
  fleetcatch$x = fleetcatch$x/sum(fleetcatch$x)
  #print(allfleet)
  #print(fleetcatch)
  
  # get indices ----------
  ind = dat$CPUE
  ni = nrow(ind)
  newind = rep(3,ni)
  for(ff in 1:2)newind[ind$index == Iind[ff]]=ff  
  aggind = aggregate(ind$obs, by=list(year = ind$year, index = newind),mean)
  aggindcv = aggregate(ind$se_log, by=list(year = ind$year, index = newind),mean)
  aggind=aggind[aggind$year>0 & aggind$year<=lastyr,] # Retro
  aggindcv=aggindcv[aggindcv$year>0 & aggindcv$year<=lastyr,] # Retro
  
  lengthswitch = !is.null(dat$lencomp)
  if(lengthswitch){
    lbinlab = dat$lbin_vector
    lbins = dat$lbin_vector / Linf_mu * 100 # all lengths are scaled to a % of Linf
    lencomp = dat$lencomp[dat$lencomp$year>0,]
    
    labs = names(lencomp)
    bysex = paste0("f",lbinlab[1])%in%labs
    if(bysex){
      comps = lencomp[,match(paste0("f",lbinlab),labs)]+lencomp[,match(paste0("m",lbinlab),labs)]
    }else{
      comps = lencomp[,match(paste0("l",lbinlab),labs)]
    }
    inflate=1E4/max(apply(comps,1,sum))
    comps =ceiling(comps*inflate)
    
    # get mean length ----------
    mucv = getmucv(comps,lbins)  
    nl = nrow(lencomp)
    newind = rep(3,nl)
    for(ff in 1:2)newind[lencomp$fleet == Find[ff]]=ff  
    aggml = aggregate(mucv$mu, by=list(year = lencomp$year, fleet = newind),mean)
    aggmlcv = aggregate(mucv$cv, by=list(year = lencomp$year, fleet = newind),mean) # this isn't ideal as the aggregate fleet 3 should have a CV to reflect combination of many mean lengths
    aggml=aggml[aggml$year>0 & aggml$year<=lastyr,] # Retro
    aggmlcv=aggmlcv[aggmlcv$year>0 & aggmlcv$year<=lastyr,] # Retro
    
    # get cv in lengths
    CVlen = sapply(1:nl,function(X,comps,lbins){
      sd(rep(lbins,comps[X,]),na.rm=T)/mean(rep(lbins,comps[X,]),na.rm=T)
    },comps=comps,lbins=lbins) 
    
    aggsdl = aggregate(CVlen, by=list(year = lencomp$year, fleet = newind),mean)
    aggsdl=aggsdl[!is.na(aggsdl$x),]
    aggsdl=aggsdl[aggsdl$year>0 & aggsdl$year<=lastyr,] # Retro
    
    # get fraction mature
    Fmat = sapply(1:nl,function(X,comps,lbins,L50){ # all these are phrased as % Linf
      indivs = rep(lbins,floor(comps[X,]*1E3)) # we are trying to calc means weighted by nsamp - these can be less than 1 and hence need upscaling by 1E3 to make sure
      mean(indivs > L50,na.rm=T)
    },comps=comps,lbins=lbins,L50=L50)
    
    aggFM = aggregate(Fmat, by=list(year = lencomp$year, fleet = newind),mean)
    aggFM=aggFM[aggFM$year>0 & aggFM$year<=lastyr,] # Retro
  }
  
  selpar = getsels(io,Find,dat,plotsmooth) # selectivity by fleet (MLE, not sampled)
  
  # Blank vectors
    
   
  # summarise fleet specific data
  fleetdat = list()
  
  for(ff in 1:3){
    I_cur<-I_mu<-I_rel<-I_g5<-I_g10<-I_g20<-I_g40<- 
      Isd1 <- Isd2 <- Isd3 <-
      C_cur<- C_mu<-C_rel<-C_g5<-C_g10<- C_g20<- C_g40<-
      Csd <- CF <- 
      ML_cur<-ML_mu<-ML_rel<-ML_g5<-ML_g10<-ML_g20<- ML_g40 <- ML_L50 <- ML_Linf <-
      MA_cur<-MA_mu<-MA_rel<-MA_g5<-MA_g10<-MA_g20<- MA_g40 <-
      MV_cur<- MV_mu<-MV_rel<- MV_g5<-MV_g10<-MV_g20<- MV_g40<-
      FM_cur <- FM_mu <- FM_rel <- FM_g5 <- FM_g10 <-FM_g20 <- FM_g40<-rep(NA,nsamp)
    
    # catch
    fcat = aggcatch$x[aggcatch$fleet==ff]
    nycat = length(fcat)
    
    # cpue
    find = aggind$x[aggind$index==ff]
    findcv = aggindcv$x[aggind$index==ff]
    nyind = length(find)
    
    if(lengthswitch){
      # mean length
      fml = aggml$x[aggml$fleet ==ff]
      fmlcv = aggmlcv$x[aggmlcv$fleet==ff]
      nyml = length(fml)
      
      #sd length
      fsdl = aggsdl$x[aggsdl$fleet ==ff]
      nysdl = length(fsdl)
      
      # frac mature
      fFM = aggFM$x[aggFM$fleet==ff]
      nyFM = length(fFM)
      
    }  
    
    for(i in 1:nsamp){
      # Indices
      if(length(find)>1){
        findsim = rlnorm(nyind, log(find),findcv)
        Is = smooth2(findsim,plot=plotsmooth)
        #Isd[i] = sd(smooth2(findsim,plot=plotsmooth,ret='resid'))
        I_cur[i]<-Is[nyind] 
        #C_mu[i]<-mean(Cobs[i,ny:Iind[i,2]]) # no info
        I_rel[i] <- I_cur[i] / mean(Is[1:nyind],na.rm=T)
        if(nyind>=5)I_g5[i]<-slp3(findsim[nyind-(4:0)])
        if(nyind>=10)I_g10[i]<-slp3(findsim[nyind-(9:0)])
        if(nyind>=20)I_g20[i]<-slp3(findsim[nyind-(19:0)])
        if(nyind>=40)I_g40[i]<-slp3(findsim[nyind-(39:0)])
        ind = nyind - (min(20,nyind):1) +1
        Isd1[i] = sd(smooth2(findsim[ind], ret = 'resid', enp.mult = 0.4, plot=plotsmooth))
        Isd2[i] = sd(smooth2(findsim[ind], ret = 'resid', enp.mult = 0.2, plot=plotsmooth))
        Isd3[i] = sd(smooth2(findsim[ind], ret = 'resid', enp.mult = 0.1, plot=plotsmooth))
        }
      # Catches 
      fcatsim = rlnorm(nycat,log(fcat),catch_CV)
      Cs = smooth2(fcatsim,plot=plotsmooth)
      Csd[i] = sd(smooth2(fcatsim,plot=plotsmooth,ret='resid'))
      C_cur[i]<-Cs[nycat] 
      #C_mu[i]<-mean(Cobs[i,ny:Iind[i,2]]) # no info
      C_rel[i] <- C_cur[i] / mean(Cs[1:nycat],na.rm=T)
      if(nycat>=5)C_g5[i]<-slp3(fcatsim[nycat-(4:0)])
      if(nycat>=10)C_g10[i]<-slp3(fcatsim[nycat-(9:0)])
      if(nycat>=20)C_g20[i]<-slp3(fcatsim[nycat-(19:0)])
      if(nycat>=40)C_g40[i]<-slp3(fcatsim[nycat-(39:0)])
      CF[i] = fleetcatch$x[fleetcatch$fleet==ff]
      
      if(lengthswitch){
        # Mean length
        mlsim = rlnorm(nyml,log(fml),fmlcv)
        Ls = smooth2(mlsim,plot=plotsmooth, enp.mult = 0.125)
        ML_cur[i]=Ls[nyml]
        #ML_mu[i]=mean(mulen[nyears:Iind[i,2]])
        ML_rel[i] <- ML_cur[i] / mean(Ls[1:nyml],na.rm=T)
        if(nyml>=5)ML_g5[i]<-slp3(mlsim[nyml-(4:0)])
        if(nyml>=10) ML_g10[i]<-slp3(mlsim[nyml-(9:0)])
        if(nyml>=20) ML_g20[i]<-slp3(mlsim[nyml-(19:0)])
        if(nyml>=40) ML_g40[i]<-slp3(mlsim[nyml-(39:0)])
        ML_Linf[i]= ML_cur[i]/Linf[i]
        ML_L50[i]= ML_cur[i]/L50[i]
        
        # cv length
        CVs = smooth2(fsdl, plot=plotsmooth,enp.mult = 0.1)
        MV_cur[i]=CVs[nysdl]
        #MV_mu[i] = mean(fsdl[nyears:Iind[i,2]],na.rm=T)
        MV_rel[i] <- MV_cur[i] / mean(CVs[1:nysdl])
        if(nysdl>=5) MV_g5[i]<-slp3(fsdl[nysdl-(4:0)])
        if(nysdl>=10)MV_g10[i]<-slp3(fsdl[nysdl-(9:0)])
        if(nysdl>=20)MV_g20[i]<-slp3(fsdl[nysdl-(19:0)])
        if(nysdl>=40)MV_g40[i]<-slp3(fsdl[nysdl-(39:0)])
        
        # fraction mature
        FMs = smooth2(fFM, plot = plotsmooth, enp.mult = 0.1)
        FM_cur[i]= FMs[nyFM]
        #FM_mu[i] = mean(FMs[nyears:Iind[i,2]])
        FM_rel[i] <- FM_cur[i] / mean(FMs[1:nyFM])
        if(nyFM >=5) FM_g5[i]<-slp3(FMs[nyFM-(4:0)])
        if(nyFM >=10) FM_g10[i]<-slp3(FMs[nyFM-(9:0)])
        if(nyFM >=20) FM_g20[i]<-slp3(FMs[nyFM-(19:0)])
        if(nyFM >=40) FM_g40[i]<-slp3(FMs[nyFM-(39:0)])
        
      }
       
    }
    
    L5_L50 = selpar$L5s[ff]/L50_mu
    LFS_L50 = selpar$LFSs[ff]/L50_mu
    VMLf = selpar$VMLs[ff]
    if(VMLf<0.025)VMLf = 0.025
    if(VMLf>0.975)VMLf = 0.975
    VML = rbeta(nsamp,alphaconv(VMLf,VML_CV),betaconv(VMLf,VML_CV))
    
    df = data.frame(I_rel, I_g5, I_g10, I_g20, I_g40,
                    Isd1, Isd2, Isd3,
                    C_rel, C_g5, C_g10, C_g20, C_g40, Csd, CF, 
                    ML_cur, ML_rel, ML_g5, ML_g10, ML_g20, ML_g40, ML_Linf,
                    MV_cur, MV_rel, MV_g5, MV_g10, MV_g20, MV_g40,
                    FM_cur, FM_rel, FM_g5, FM_g10, FM_g20, FM_g40,
                    MA_cur, MA_rel, MA_g5, MA_g10, MA_g20, MA_g40,
                    L5_L50, LFS_L50, VML)
    
    names(df) = paste0(names(df),"_s1_f",ff)
    fleetdat[[ff]]=df[,!apply(df,2,function(x){all(is.na(x))})]
    
  }
  outdat = fleetdat[[1]]
  for(ff in 2:3)outdat=cbind(outdat,fleetdat[[ff]])
  
  bdf = data.frame(K, M_K, maxa, L50_Linf)
  names(bdf) = paste0(names(bdf),"_s1")
  outdat = cbind(bdf,outdat)
  
  yrs = io$input$dat$styr:io$input$dat$endyr
  SSBlab = paste0("SSB_",yrs)
  SSBest = io$outputs$SSB_est
  Brel = SSBest[match(SSBlab,rownames(SSBest)),2:3] / SSBest[match("SSB_MSY",rownames(SSBest)),2]
 
  list(dat=outdat,Brel=Brel)
  
}


SS_2_ET_Retro=function(io, Fnam , Inam, npeels=8){
  outlist = list()
  for(pp in 0:(npeels-1))    outlist[[pp+1]]= SS_2_ET(io, Fnam, Inam, peel=pp, plotsmooth=F)
  nd = sapply(outlist,function(x)ncol(x[[1]]))
  
  if(length(unique(nd))!=1){
    nams = table(unlist(lapply(outlist,function(x)names(x[[1]]))))
    minlabs = names(nams)[nams==npeels]
    cat(paste0("Pruning from ",max(nd), " to ",length(minlabs)," data inputs \n"))
    #minlabs = names(outlist[[npeels]][[1]])
    for(pp in 1:npeels){
      tab =  outlist[[pp]][[1]] 
      outlist[[pp]][[1]] = tab[,match(minlabs,names(tab))]
    }
  }
  outlist
    
}



all_ss3 = function(dir){
  inputs = get_ss3_inputs(dir)
  SSout = SS_output(dir)
  outputs = extract_SS(SSout)
  list(inputs=inputs,outputs=outputs)
}

get_ss3_inputs = function(dir){

  dat = SS_readdat(paste0(dir,"/data.ss"))
  control = SS_readctl(paste0(dir,"/control.ss"),datlist=dat)
  starter = SS_readstarter(paste0(dir,"/starter.ss"))
  forecast = SS_readforecast(paste0(dir,"/forecast.ss"))
  list(dat=dat, control=control, starter = starter, forecast = forecast)

}


# fit = SS_output(dir = here('assessment/runs/straw_dog'))    
extract_SS<-function(fit, name="Unnamed"){
  if(class(fit) != "list")  fit = readRDS(fit) 
  
  out<-list()
  out$model_name<-name
  # Estimates
  
  out$dep = fit$current_depletion
  
  # dimensions
  ts = fit$timeseries
  out$startyr = fit$startyr
  out$endyr = fit$endyr
  out$ages = fit$agebins
  out$yrs = ts$Yr
  
  # MLE time series
  out$SSB = ts$SpawnBio
  out$B = ts$Bio_all
  out$SSN = ts$mature_num
  out$ind = fit$cpue[,names(fit$cpue) %in% c("Fleet_name","Yr","Vuln_bio","Obs","Exp","Calc_Q","SE","Dev","Like")]
  out$FM = fit$exploitation[,c(1,4:ncol(fit$exploitation))]
  out$sel = fit$ageselex
  out$M = fit$Natural_Mortality
  out$catch = fit$catch[,names(fit$catch) %in% c("Fleet_name","Yr","Obs","Exp","se","F","Like")]
  NAA = fit$natage[,names(fit$natage)%in%c("Sex","Yr","Beg/Mid",0,out$ages)]
  out$NAA = NAA[NAA[,3]=="B",c(1:2,4:ncol(NAA))]
  
  out$npar = fit$N_estimated_parameters
  
  # Estimates with uncertainty
  dq=fit$derived_quants[,1:3]
  out$FM_est = dq[grepl("F_",dq[,1]),]
  out$SSB_est = dq[grepl("SSB_",dq[,1]),]
  out$Rec_est = dq[grepl("Recr_",dq[,1]),]
  out$Dep_est = dq[grepl("Bratio_",dq[,1]),]
  
  
  # Reference points
  out$Refs=list()
  out$Refs$SSB_MSY = dq[dq[,1]=="SSB_MSY",]
  out$Refs$Brel = dq[dq[,1]=="B_MSY/SSB_unfished",]
  out$Refs$R0 = dq[dq[,1] =="Recr_unfished",]
  out$Refs$SSB0 = dq[dq[,1]=="SSB_unfished",]
  out$Refs$B0 = dq[dq[,1]=="Totbio_unfished",]
  out$Refs$OFL =  dq[grepl("OFLCatch",dq[,1]),]
  out$Refs$MSY = dq[dq[,1]=="Ret_Catch_MSY",]
  
  # Objective function
  out$nll = fit$likelihoods_used
  out$nll_ft = fit$likelihoods_by_fleet
  
  out$AIC  <- 2 * (out$nll$values[1] + out$npar) 
  out$CoVar = fit$CoVar
  # Retro
  #if(!is.null(Fit$peels))out$rho<-mohns_rho(Fit)
  out$RecDevs = fit$recruitpars
  #out$conv = fit$opt$convergence
  out$RunTime = fit$RunTime
  out$Nwarnings = fit$Nwarnings
  out$Warnings = fit$warnings
  out$maxgrad = fit$maximum_gradient_component
  out$sigmaR = fit$sigma_R_info
  
  out$par = fit$parameters
  out$par_h = fit$parameters_with_highest_gradients
  
  out$agedbase = fit$agedbase
  out$sizeselex = fit$sizeselex
  
  out$FleetNames = fit$FleetNames
  # object.size(out)/object.size(fit)
  out
  
}

extract_SS_list<-function(mdirs){
  
  out<-list()
  #mods<-list()
  for(i in 1:length(mdirs)){
    Fit<-readRDS(mdirs[[i]])
    out[[i]]<-extract_SS(Fit)
  }
  names(out)<-names(mdirs)
  
  out
  
}


get_ss3_outputs = function(dir){
  fit = SS_output(dir)
  extract_SS(fit)
}