
retro_ind = function(pred, newplot=F){

  if(class(pred)=="ETRetro"){

    Prel = t(apply(pred$pred,2,quantile, p=c(0.05,0.5,0.95)))
    ylim = c(0,max(Prel))
    yind = yind2 = 1:nrow(Prel)
    if("Brel" %in% names(pred)){
      Bin = pred$Brel
      Brel = t(sapply(1:nrow(Bin),function(x,Bin)qlnorm(c(0.05, 0.5, 0.95),log(Bin[x,2]),Bin[x,3]/Bin[x,2]),Bin=Bin))
      ylim[2] = max(ylim[2], max(Brel))
      yind = 1:nrow(Brel)
      yind2 = nrow(Brel) - (nrow(Prel)-1):0

    }
    if(newplot)par(mai=c(1.2,0.8,0.2,0.05))
    matplot(yind, Brel, col="red",lty=1,type="l",ylim=ylim,axes=F,xlab="",ylab="SSB/SSBMSY"); grid()
    axis(2)
    axis(1,at=yind,labels=pred$Brel$Year,las=2)
    legend('topright',c("Assessed","EcoTest"),text.col=c("red","blue"),bty="n")
    matplot(yind2, Prel, col="blue",lty=2,add=T,type="l",axes=F)

  }else{
    stop("The EcoTest prediction object does not have multiple peels of data")
  }

}
