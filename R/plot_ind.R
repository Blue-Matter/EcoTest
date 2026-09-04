
plot_ind = function(pred, newplot=F, curY = 2024){

    Prel = t(apply(pred$pred,2,quantile, p=c(0.05,0.5,0.95)))
    ylim = c(0,max(Prel))
    yind2 = curY -nrow(Prel):1 +1
    if(newplot)par(mai=c(1.2,0.8,0.2,0.05))
    matplot(yind2, Prel, col="blue",lty=2,type="l",xlab="",ylab="SSB/SSBMSY",ylim=ylim); grid()
    legend('topright',c("EcoTest"),text.col=c("blue"),bty="n")

}
