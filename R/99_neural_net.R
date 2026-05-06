

# exp(c(-3.689,0.693))
cleandat = function(dat,resprange=c(-3.689,0.693)){
 # test = (x==0); (1:length(test))[test]
  goodrow = function(x)(sum(is.na(x)| x==0 | x == 'Inf')==0)&x[1]>resprange[1] & x[1]<resprange[2]
  keep = apply(dat,1,goodrow)
  dat[keep,]
}

reformat = function(x){
  
  renam = grepl("-",x)
  charmatch("-",x[162])
  
}

dithfunc = function(dat,type = "VML",ub = 0.99, lb = 0.01,sd=0.005){
  dithcols = grep(type,names(dat))
  for(cc in dithcols){
    ind = (dat[,cc] > ub)
    if(sum(ind)>0)  dat[ind,cc] = rbeta(sum(ind),MSEtool::alphaconv(ub,sd), MSEtool::betaconv(ub,sd))
    ind = (dat[,cc] < lb)
    if(sum(ind)>0)  dat[ind,cc] = rbeta(sum(ind),MSEtool::alphaconv(lb,sd), MSEtool::betaconv(lb,sd))
  }
  dat
}

dolog = function(dat, types){
  boolarray = sapply(types,function(x)grepl(x,names(dat)))
  tolog = apply(boolarray,1,any) 
  dat[,tolog] = log(dat[,tolog])
  dat
}

dologit = function(dat, types = "VML", bnds = c(0.05,0.95)){
  if(types[1] %in% names(dat)){
    boolarray = sapply(types,function(x)grepl(x,names(dat)))
    tologit = apply(boolarray,1,any) 
    rng = bnds[2] - bnds[1]
    dat[,tologit] = bnds[1]+dat[,tologit]*rng
    dat[,tologit] = MSEtool:::logit(dat[,tologit])
  }
  dat
}

makerawdata = function(allout, sno=1, isBrel = F, clean = T, 
                       inc_spat = T, inc_Irel = T, inc_I = T,
                       inc_CR = T, stock_in = NA){
  
  
  cdat = as.data.frame(rbindlist(allout))
  dnames = names(cdat)
  Brelcols = grepl("Brel",dnames)
  
  # get stocks to be included (default to all)
  stocks = sapply(dnames[Brelcols],function(x)strsplit(x,"_")[[1]][2])
  ns = length(stocks)
  if(is.na(stock_in[1]))stock_in = stocks[stocks!=sno]
  nsi = length(stock_in)
  
  # keep only those listed in sno and stock_in
  unid = dnames[1:(grep("Brel_2",dnames)-1)]
  ulabs = c(sapply(unid,function(x)strsplit(x,"_1")[[1]][1]),"spat")
  baselabs = paste0(rep(ulabs,1+nsi),"_",rep(c(sno,stock_in),each=length(ulabs)))
  corlab = expand.grid(c(sno,stock_in[1:(nsi-1)]),stock_in,stringsAsFactors = F)
  corlab = corlab[corlab[,1]!=corlab[,2],]
  ncomb = nrow(corlab)
  cortypes = c("CR","CR_mu","CR_rel","CR_s5","CR_s10","CC20","CC40","CC60")
  ncor = length(cortypes)
  CRlabs = paste0(rep(cortypes,each=ncomb),"_",corlab[,1],"_",corlab[,2])
  keepcol =  dnames %in% c(baselabs,CRlabs)
  cdat = cdat[,keepcol]
  
  dnames = names(cdat)
  Brelcols = grepl("Brel",dnames)
  ns = sum(Brelcols) # os = (1:ns)[!((1:ns)%in%sno)]
  keepcol = grep(paste0("Brel_",sno),dnames)
  Res = cdat[,keepcol]
 
  if(isBrel){Res = ((ns-1) * Res) /(apply(cdat[,Brelcols],1,sum)-Res)}# level relative to mean level of other stocks
  Res = log(Res)
  dat = cbind(Res,cdat[,!Brelcols])
  
  if(!inc_spat) dat = dat[,!grepl('spat',names(dat))]
  if(!inc_Irel) dat = dat[,!grepl('I_rel',names(dat))]
  if(!inc_I)    dat = dat[,!grepl("I",names(dat))]
  if(!inc_CR)   dat = dat[,!(grepl("CR",names(dat))|grepl("CC",names(dat)))]
 
  # dither asymptotic selectivity (when exactly 0 or 1)
  dat = dithfunc(dat,type = "VML")
  
  # log fractions 
  dat = dolog(dat, types = c("I_rel","C_rel","CR_mu","FM_cur","FM_rel","ML_cur","ML_rel","MV_cur","MV_rel","ML_Linf",
                             "ML_L50","CR_rel","maxa","M_K","L5_L50", "LFS_L50","Csd","Isd"))
  
  # logit proportions (but rescaled 0.05 - 0.95 prior to logit)
  dat = dologit(dat,types = "VML")
  
  if(clean) dat=cleandat(dat)
  isconst = apply(dat,2,sd)<1E-10 # sum(isconst,na.rm=T); names(dat)[is.na(isconst)]
  if(sum(isconst,na.rm=T)>0){
    cat(paste(paste(names(dat)[isconst],collapse=", "), "dropped for sd < 1E-10 \n"))
    dat = dat[,!isconst]
  }
  
  dat
  
}




power_tab = function(sim, pred, lev = c(0.5,1),asprob = T){
  
  nlev = length(lev)
  ncat = nlev+1
  labs=rep(NA,ncat)
  labs[1] = paste("Brel <",lev[1])
  labs[ncat] = paste(lev[nlev], "< Brel")
  if(nlev>1){ for(i in 2:nlev){
    labs[i] = paste(lev[i-1],"< Brel <",lev[i])
  }}
  
  tab = array(NA,c(ncat,ncat))
  row.names(tab) = paste("(Pred)",labs)
  colnames(tab) = paste("(Sim)",labs)
  
  alllev = c(-Inf,lev,Inf)
  for(i in 1:ncat){ #pred
    for(j in 1:ncat){ #sim
      cond = sim > alllev[j] & sim < alllev[j+1] & pred > alllev[i] & pred < alllev[i+1]
      tab[i,j] = sum(cond)
    }  
  }
  
  if(asprob)  tab = tab / rep(colSums(tab),each=ncat)
  tab
}


NN_fit = function(sim, pred, history,lev, addpow=T, nepoch,plot=T){
  tab = power_tab(sim, pred, lev)
  if(plot) pred_plot(inlist=list(sim=sim,pred=pred,grid=grid,lev =lev,tab=tab,r2 = cor(sim,pred)^2,MAE=mean(abs(sim-pred))))
  tab
}


pred_plot = function(inlist,axlim=c(0,2),newplot=T,lab=NA, adj=0.25){
  
  if(newplot)par(mfrow=c(1,1),mai=c(0.7,0.7,0.05,0.05),omi=c(0,0,0,0))
  sim = inlist$sim
  pred=inlist$pred
  grid=inlist$grid
  lev=inlist$lev
  tab=inlist$tab
  r2 = inlist$r2
  MAE = inlist$MAE
  
  nlev = length(lev)
  alllev = c(-Inf,lev,Inf)
  ncat = nlev+1
  colt = c("red","orange","green") #rainbow(ncat,start=0,end=0.35)
  cols = rep(colt[1],length(sim))
  for(i in 2:ncat){
    cols[pred > alllev[i] & pred < alllev[i+1]] = colt[i]
  }

 # mins_pred = c(min(pred),lev)
#  maxs_pred = c(lev,max(pred))
#  difs_pred = maxs_pred-mins_pred
#  muy = lev + adj
  
 # mins_sim = c(min(sim),lev)
#  maxs_sim = c(lev,max(sim))
 # difs_sim = maxs_sim-mins_sim
  muy = mux = lev[c(1,1:length(lev))] + c(-adj,rep(adj,length(lev))) #maxs_sim-difs_sim/4
  
  grid = expand.grid(muy,mux)
  plot(sim, pred, xlab="",ylab="",pch=19,cex=1.2,col="white",ylim=axlim, xlim=axlim)
  mtext("SSB/SSBMSY (simulated)",1,line=2.3)
  mtext("SSB/SSBMSY (predicted)",2,line=2.3)
  lines(c(0,1E10),c(0,1E10),col='black',lwd=1,lty=2)
  abline(v=0,h=0,lty=2)
  points(sim, pred, pch=19,cex=1.2,col=cols)
  abline(h=lev,v=lev,lty=2)
  text(grid[,2],grid[,1],round(as.vector(tab)*100,1),font=2,cex=1.2)
  legend('topleft',legend = c(paste("MAE =",round(inlist$MAE,3)),
                              paste("R-squared =",round(inlist$r2,3))))
  if(!is.na(lab)) mtext(lab,line=0.5,cex=0.9)
  
}

calc_importance = function(model, testy, adj = 0.25, barno = 30){
  nin = ncol(testy)
  impdf = rbind(testy,testy,testy,testy,testy)[1:(nin+1),]
  impdf[]=0.0
  ind = as.matrix(cbind(2:(nin+1),1:nin))
  impdf[ind] = adj
  sense = exp((model %>% predict(impdf))[,1])
  dif = abs(sense[2:nin] - sense[1])
  ord = order(dif,decreasing=T)
  out = data.frame(input=colnames(testy)[ord],output_diff = dif[ord])
  y = out[1:barno,2]
  x=barplot(y,xaxt='n');grid(); barplot(y,xaxt='n',add=T)
  yadj = max(y)/8
  text(x=x[,1]-0.25, y=rep(-yadj,barno), labels=out[1:barno,1], xpd=TRUE, srt=90,cex=0.85)
  mtext("Relative weight of input",2,line=2.3)
  out
}

