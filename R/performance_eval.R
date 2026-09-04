
perf_plot = function(inlist,axlim=c(0,2),newplot=T,lab=NA, adj=0.25){

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
