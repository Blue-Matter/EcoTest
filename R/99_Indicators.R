# Functions for extracting indices

slp3<-function(y0){
  y=log(y0)
  x1<-1:length(y0)
  x1[is.na(y)]=NA
  mux<-mean(x1,na.rm=T)
  muy<-mean(y,na.rm=T)
  SS<-sum((x1-mux)^2,na.rm=T)
  (1/SS)*sum((x1-mux)*(y-muy),na.rm=T)
}

slp3_nolog<-function(y){
  x1<-1:length(y)
  x1[is.na(y)]=NA
  mux<-mean(x1,na.rm=T)
  muy<-mean(y,na.rm=T)
  SS<-sum((x1-mux)^2,na.rm=T)
  (1/SS)*sum((x1-mux)*(y-muy),na.rm=T)
}



smooth2<-function(xx,plot=F,enp.mult=0.2,plotname="",ret = "pred"){
  tofill<-!is.na(xx)
  xx[xx==0]<-1E3
  #xx[tofill] = log(xx[tofill])
  predout<-rep(NA,length(xx))
  dat<-data.frame(x=1:length(xx),y=log(xx))
  enp.target<-sum(tofill)*enp.mult
  out<-loess(y~x,dat=dat,enp.target=enp.target)
  #out<-loess(y~x,dat=dat,span=0.5)
  predout[tofill]<-exp(predict(out))
  if(plot){
    plot(xx,type="p",xlab="x",ylab="y",main=plotname)
    lines(predout,col="#ff000090",lwd=2)
  }
  if(ret =="pred")return(predout)
  if(ret =="resid")return(log(xx/predout))
}

smooth3<-function(xx,plot=F,enp.mult=0.3,plotname=""){
  tofill<-!is.na(xx)
  predout<-rep(NA,length(xx))
  dat<-data.frame(x=1:length(xx),y=xx)
  enp.target<-sum(tofill)*enp.mult
  out<-loess(y~x,dat=dat,enp.target=enp.target)
  predout[tofill]<-predict(out)
  if(plot){
    plot(xx,type="p",xlab="x",ylab="y",main=plotname)
    lines(predout,col="#ff000090",lwd=2)
  }
  predout
}

proc_dat<-function(MMSE,Iind=NA,sno = 1, fno=1, plotsmooth=F){
  
  dat = MMSE@PPD[[sno]][[fno]][[1]]
  stock = MMSE@Stocks[[sno]]
  nsim<-MMSE@nsim
  nyears<-MMSE@nyears
  proyears = MMSE@proyears
  allyears = nyears+proyears
 
  datlist<-list()
  
  Iobs<-dat@Ind
  Cobs<-dat@Cat
  mids = dat@CAL_mids
 
  #prob<-seq(1,0.5,length.out=length(allyears))
   
  if(length(Iind)==1){
    Iyr<-sample((nyears+1):(allyears-2),nsim,replace=T)#,prob=prob)
    Byr = Iyr - nyears
    Iind<-cbind(1:nsim,Iyr)
  }else{
    Iyr = Iind[,2]
    Byr = Iyr - nyears
  }
 
  Bt = MMSE@SB_SBMSY
  I_cur<-I_mu<-I_rel<-I_s5<-I_s10<-I_s20<- 
    Isd1 <- Isd2 <- Isd3 <-
    C_cur<- C_mu<-C_rel<-C_s5<-C_s10<- C_s20<- 
    Csd <-
    ML_cur<-ML_mu<-ML_rel<-ML_s5<-ML_s10<-ML_s20<- ML_L50 <- ML_Linf <-
    MV_cur<- MV_mu<-MV_rel<- MV_s5<-MV_s10<-MV_s20<-
    FM_cur <- FM_mu <- FM_rel <- FM_s5 <- FM_s10 <-FM_s20 <- rep(NA,nsim)
  
  matage = MMSE@multiHist[[sno]][[fno]]@AtAge$Maturity[,,1]
  lenage = MMSE@multiHist[[sno]][[fno]]@AtAge$Length[,,1]
  
  OM = MMSE@OM[[sno]][[fno]]
 
  L50 = OM$L50
  Linf = OM$Linf
  M = OM$M
  K = OM$K
  M_K = M/K
  maxa = -log(0.05)/M # age at 5% cumulative survival
  
  L5 = MMSE@multiHist[[sno]][[fno]]@SampPars$Fleet$L5_y[,1]
  LFS = MMSE@multiHist[[sno]][[fno]]@SampPars$Fleet$LFS_y[,1]
  VML = MMSE@multiHist[[sno]][[fno]]@SampPars$Fleet$Vmaxlen_y[,1]
  L5_L50 = L5/L50
  LFS_L50 = LFS/L50

  for(i in 1:nsim){
    # Index
    Is<-smooth2(Iobs[i,],plot=plotsmooth)
    I_cur[i]<-Is[Iind[i,2]]
    I_mu[i] = mean(Iobs[i,nyears:Iind[i,2]])
    I_rel[i]<-I_cur[i] / mean(Is[1:Iind[i,2]])
    I_s5[i]<-slp3(Iobs[i,Iind[i,2]-(5:1)])
    I_s10[i]<-slp3(Iobs[i,Iind[i,2]-(10:1)])
    I_s20[i]<-slp3(Iobs[i,Iind[i,2]-(20:1)])
    
    ind = Iind[i,2] - (20:1)
    Isd1[i] = sd(smooth2(Iobs[i,ind], ret = 'resid', enp.mult = 0.4, plot=plotsmooth))
    Isd2[i] = sd(smooth2(Iobs[i,ind], ret = 'resid', enp.mult = 0.2, plot=plotsmooth))
    Isd3[i] = sd(smooth2(Iobs[i,ind], ret = 'resid', enp.mult = 0.1, plot=plotsmooth))
    
    # Catch
    Cs = smooth2(Cobs[i,],plot=plotsmooth)
    Csd[i] = sd(smooth2(Cobs[i,],plot=plotsmooth,ret='resid'))
    C_cur[i]<-Cs[Iind[i,2]]
    C_mu[i]<-mean(Cobs[i,nyears:Iind[i,2]])
    C_rel[i] <- C_cur[i] / mean(Cs[1:Iind[i,2]],na.rm=T)
    C_s5[i]<-slp3(Cobs[i,Iyr[i]-(5:1)])
    C_s10[i]<-slp3(Cobs[i,Iyr[i]-(10:1)])
    C_s20[i]<-slp3(Cobs[i,Iyr[i]-(20:1)])
   
    # Length
    CALi = dat@CAL[i,,]
    totlen = CALi*t(array(mids,dim(t(CALi))))
    mulen = apply(totlen,1,sum)/apply(CALi,1,sum)
    if(all(is.na(mulen))){
      Ls = rep(NA,length(mulen))
    }else{
      Ls = smooth2(mulen,plot=plotsmooth, enp.mult = 0.125)
    }
    
    ML_cur[i]=Ls[Iind[i,2]]
    ML_mu[i]=mean(mulen[nyears:Iind[i,2]])
    ML_rel[i] <- ML_cur[i] / mean(Ls[1:Iind[i,2]],na.rm=T)
    ML_s5[i]<-slp3(mulen[Iyr[i]-(5:1)])
    ML_s10[i]<-slp3(mulen[Iyr[i]-(10:1)])
    ML_s20[i]<-slp3(mulen[Iyr[i]-(20:1)])
    ML_Linf[i]= ML_cur[i]/Linf[i]
    ML_L50[i]= ML_cur[i]/L50[i]
    
    # standard deviation in length
    CVlen = sapply(1:dim(CALi)[1],function(X,CALi,mids){
      sd(rep(mids,CALi[X,]))/mean(rep(mids,CALi[X,]))
     },CALi=CALi,mids=mids)
    
    if(all(is.na(CVlen))){
      CVs = rep(NA,length(CVlen))
    }else{
      CVs = smooth2(CVlen, plot=plotsmooth,enp.mult = 0.1)
    }
   
    MV_cur[i]=CVs[Iind[i,2]]
    MV_mu[i] = mean(CVlen[nyears:Iind[i,2]],na.rm=T)
    MV_rel[i] <- MV_cur[i] / mean(CVs[1:Iind[i,2]])
    MV_s5[i]<-slp3(CVlen[Iyr[i]-(5:1)])
    MV_s10[i]<-slp3(CVlen[Iyr[i]-(10:1)])
    MV_s20[i]<-slp3(CVlen[Iyr[i]-(20:1)])
    
    L50i = L50[i]
    Fmat = sapply(1:dim(CALi)[1],function(X,CALi,mids,L50i){
      indivs = rep(mids,floor(CALi[X,]*1E3)) # we are trying to calc means weighted by nsamp - these can be less than 1 and hence need upscaling by 1E3 to make sure
      mean(indivs > L50i,na.rm=T)
    },CALi=CALi,mids=mids,L50i=L50i)
    
    if(all(is.na(Fmat))){
      FMs = rep(NA,length(Fmat))
    }else{
      FMs = smooth2(Fmat, plot=plotsmooth, enp.mult=0.1)
    }
    
    FM_cur[i]= FMs[Iind[i,2]]
    FM_mu[i] = mean(FMs[nyears:Iind[i,2]])
    FM_rel[i] <- FM_cur[i] / mean(FMs[1:Iind[i,2]])
    FM_s5[i]<-slp3(Fmat[Iyr[i]-(5:1)])
    FM_s10[i]<-slp3(Fmat[Iyr[i]-(10:1)])
    FM_s20[i]<-slp3(Fmat[Iyr[i]-(20:1)])

  }
  
  Bind<-cbind(1:nsim,rep(sno,nsim),rep(1,nsim),Byr)
  Brel = MMSE@SB_SBMSY[Bind]
 
  data.frame(Brel, I_rel, I_s5, I_s10, I_s20, 
             Isd1, Isd2, Isd3,
             C_rel, C_s5, C_s10, C_s20, Csd, 
             ML_cur, ML_rel, ML_s5, ML_s10, ML_s20, ML_L50, ML_Linf,
             MV_cur, MV_rel, MV_s5, MV_s10, MV_s20,
             FM_cur, FM_rel, FM_s5, FM_s10, FM_s20,
             M_K, maxa, L5_L50, LFS_L50, VML)
  
}

cat_ratios = function(MMSE, Iind, fno=1, plotsmooth=F, refyrs = 20){
  
  ns = MMSE@nstocks
  nsim = MMSE@nsim
  allyears = MMSE@proyears+MMSE@nyears
  nyears = MMSE@nyears
  
  CRtab = NULL
  
  for(i in 1:nsim){ 
    
    catsmth = catsmth2 = array(NA,c(ns,allyears-1))
    
    for(sno in 1:ns){
      catsmth[sno,] = smooth2(MMSE@PPD[[sno]][[fno]][[1]]@Cat[i,],plot=plotsmooth)
      catsmth2[sno,] = smooth2(MMSE@PPD[[sno]][[fno]][[1]]@Cat[i,],enp.mult = 0.075,ret='resid',plot=plotsmooth) # log residuals after broad detrending
    }
    
    ncomp = ((ns-1)*ns)/2
    
    CR_rel =  CR_s5 = CR_s10 = CR_mu = CR = namy = rep(NA,ncomp)
    
    sampyr = Iind[i,2]
    j = 0
    for(ss in 1:(ns-1)){
      for(s2 in (ss+1):ns){
        
        j=j+1
        namy[j] = paste(ss,s2,sep="_")
        rat = catsmth[ss,]/catsmth[s2,]
        CR[j] = rat[sampyr]
        invalidrat = rat<(0.01*mean(rat))&((1:length(rat))<refyrs)
        rat2 = rat; rat2[invalidrat] = NA
        CR_mu[j] = mean(rat2[1:refyrs],na.rm=T)
        CR_rel[j] = rat2[sampyr] / mean(rat2[1:refyrs],na.rm=T)
        CR_s5[j] = slp3(rat2[sampyr-(5:1)])
        CR_s10[j] = slp3(rat2[sampyr-(10:1)])
        
      }
    }  
    
    CRvec = c(CR,CR_mu,CR_rel,CR_s5,CR_s10)
    names(CRvec) = paste(rep(c("CR","CR_mu","CR_rel","CR_s5","CR_s10"),each=ncomp),
                       rep(namy,5),sep="_")
    
    
    CC20 = CC40 = CC60 = rep(NA, ncomp)
    
    j = 0
    ind20 = Iind[i,2] - (20:1)
    ind40 = Iind[i,2] - (40:21)
    ind60 = Iind[i,2] - (60:41)
    
    for(ss in 1:(ns-1)){
      for(s2 in (ss+1):ns){
        
        j=j+1
        CC20[j] = cor(catsmth2[ss,ind20],catsmth2[s2,ind20])
        CC40[j] = cor(catsmth2[ss,ind40],catsmth2[s2,ind40]) 
        CC60[j] = cor(catsmth2[ss,ind60],catsmth2[s2,ind60]) 
        
      }
    }  
    
    CCvec = c(CC20,CC40,CC60)
    names(CCvec) = paste(rep(c("CC20","CC40","CC60"),each=ncomp),
                         rep(namy,3),sep="_")
    
    
                   
    CRtab = rbind(CRtab, c(CRvec,CCvec))  
    #cat(".")
    
  } # simulation
  #cat("\n")
  CRtab
  
}


filt=function(x){
  x[x%in%c("Inf","-Inf")] = NA
  x
} 


corplot = function(out,quanty=0.5,ptcex=0.8,maxn=20, lasinv=F,lab=""){
  
  nams = names(out)
  ndat = dim(out)[2]
  nobs = dim(out[1])
  ni <- min(ndat,maxn)
  
  ncol = 200
  cols <- rainbow(ncol,start=0.1,end=0.3,alpha=0.7)
  cols <- c(cols,rep(cols[ncol],1000-ncol))
  Brel = out$Brel
  Bcol = cols[ceiling(Brel*100)]; plot(Brel[order(Brel)],col=Bcol[order(Brel)],pch=19)
  par(mfrow=c(ni-1,ni-1),mai=rep(0,4),omi=c(0.55,0.75,0.05,0.05))
  cutoff= c(quanty/100,(100-quanty)/100)
  # labys = c(letters,paste0(letters,letters),paste0(letters,letters,letters),paste0(letters,letters,letters,letters))
  labys = paste0(rep(letters,9),rep(1:9,each=length(letters)))
  labys = paste0("(",labys,")")
  k = 0
  
  for(i in 2:ni){ # i<-i+1
    
    for(j in 1:(ni-1)){
      
      if(j==i|j>i){
        
        plot(1,1,col='white',axes=F)
        
      }else{
        xs = filt(out[,j])
        ys = filt(out[,i])
        xlim<-range(quantile(xs,cutoff,na.rm=T))
        ylim<-range(quantile(ys,cutoff,na.rm=T))
        plot(xs,ys,pch=19,xlim=xlim,ylim=ylim,cex=ptcex,col=Bcol,axes=F)
        k=k+1; mtext(labys[k],3,adj=0.02,line=-1,cex=0.6)
        axis(1,c(-1000,1000))
        axis(2,c(-1000,1000))
        axis(3,c(-1000,1000))
        axis(4,c(-1000,1000))
        #points(obs[j],obs[i],pch=19,cex=1.2,col=cols[2])
        if(nams[j]=="Brel")abline(v=c(0.5,1),lty=c(2,1),col="#99999999")
        
      }
      if(i==2&j==(ni-1)){
        legend('center',legend=lab,text.col="blue",bty='n',text.font=2,cex=1.2)
      }
      
      if(!lasinv){
        if(j==1)mtext(nams[i],2,line=2,cex=0.6,las=2)
        if(i==ni)mtext(nams[j],1,line=1,cex=0.6,las=2)
      }else{
        if(j==1)mtext(nams[i],2,line=2,cex=0.6,srt=90)
        if(i==ni)mtext(nams[j],1,line=1,cex=0.6,las=1)
      }
      #if(j==1)mtext(i,2,line=2,cex=0.5,las=2)
      #if(i==nplotted)mtext(j,1,line=1,cex=0.5,las=2)
      
    }
    
  }
  
#}


}




simulate_tc_lm = function(object, Brel){
  ftd <- predict(object,newdata = data.frame(x=log(Brel)))   # == napredict(*, object$fitted)
  vars <- deviance(object)/ df.residual(object)
  ftd + rnorm(length(Brel), sd = sqrt(vars))
}

sim_spatial_dat = function(MMSE,Iind, spat_mods){
  
  Iyr = Iind[,2]
  nyears= MMSE@nyears
  ns = MMSE@nstocks
  nsim = MMSE@nsim
  ny = 15
  
  spatstat = array(NA,c(nsim,ns))
  
  for(sno in 1:ns){
    histSSB = apply(MMSE@multiHist[[sno]][[1]]@TSdata$SBiomass,1:2,sum)
    futureSSB = MMSE@SSB[,sno,1,]
    SSB = cbind(histSSB,futureSSB)
    SSBMSY = (futureSSB / MMSE@SB_SBMSY[,sno,1,])[1,]
    SBrel = SSB / SSBMSY
    object = spat_mods[[sno]]
    for(i in 1:nsim){
      Bind<-cbind(rep(i,ny),Iyr[i]-(1:ny))
      Brel = SBrel[Bind]
      obs = simulate_tc_lm(object, Brel)
      spatstat[i,sno] = smooth3(obs,plot=F)[length(obs)]  
    }
  }
  
  stat= as.data.frame(spatstat)
  names(stat) = paste0("spat_",1:sno)
  stat
  
}

get_sim_data = function(ff,filelocs, spat_mods){
  
  set.seed(ff)
  MMSE = readRDS(filelocs[ff])
  nsim<-MMSE@nsim
  nyears<-MMSE@nyears
  proyears = MMSE@proyears
  allyears = nyears+proyears
  Iyr<-sample((nyears+1):(allyears-2),nsim,replace=T)
  Byr = Iyr - nyears
  Iind<-cbind(1:nsim,Iyr)
  
  ns=MMSE@nstocks
  outs = list()
  for(sno in 1:ns)outs[[sno]] = proc_dat(MMSE,Iind=Iind,sno=sno)
  outs[[ns+1]] = cat_ratios(MMSE, Iind=Iind)
  outs[[ns+2]] = sim_spatial_dat(MMSE,Iind=Iind,spat_mods)

  for(i in 1:(ns+2)){
    out = outs[[i]]
    if(i <=ns) names(out) = paste0(names(out),"_",i)
    if(i==1) one_tab=out
    if(i>1) one_tab=cbind(one_tab,out)
  }
  row.names(one_tab) = paste0("sim_",1:nrow(one_tab))
  
  cat(".")
  one_tab
  
}

process_sim_data = function(MSEdir, spat_mods, parallel=T, cores = NA){
  
  files = list.files(MSEdir)
  keep = grepl("MMSE",files)
  filelocs = list.files(MSEdir,full.names=T)[keep]
  nfile = length(filelocs)
  
  if(parallel){
    library(snowfall)
    library(parallel)
     if(is.na(cores))cores = detectCores()/2
     sfInit(parallel=T,cpus = cores)
     sfExport("proc_dat"); sfExport("smooth2"); sfExport('slp3'); sfExport('cat_ratios'); 
     sfExport("smooth3"); sfExport("sim_spatial_dat"); sfExport("simulate_tc_lm")
     allout = sfLapply(1:nfile,get_sim_data,filelocs=filelocs, spat_mods=spat_mods)
  }else{
    allout = lapply(1:2, get_sim_data, filelocs=filelocs, spat_mods=spat_mods)
    allout=list()
    for(i in 1:nfile)allout[[i]]=get_sim_data(i,filelocs=filelocs, spat_mods=spat_mods)
    
  }
  
  cat("\n")
  allout
  
}



