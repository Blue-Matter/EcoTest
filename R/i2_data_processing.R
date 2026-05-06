
getdat=function(PPD,lab,simno=1,ay = 70){
  ns = length(PPD);  nf = length(PPD[[1]]);  mpno = 1
  egdat = slot(PPD[[1]][[1]][[mpno]],lab); yrs = PPD[[1]][[1]][[mpno]]@Year[1:ay]
  ny = length(egdat)
  out=array(NA,c(ns,nf,ay)); dimnames(out)[[1]] = paste0("Stock",1:ns);  dimnames(out)[[2]] = paste0("Fleet",1:ns); dimnames(out)[[3]] = yrs
  for(ss in 1:ns)for(ff in 1:nf)out[ss,ff,]=slot(PPD[[ss]][[ff]][[mpno]],lab)[simno,1:ay]
  out
}

getsdat = function(PPD,lab,simno){
  ns = length(PPD)
  out=sapply(PPD,function(x,lab,simno){slot(x[[1]][[1]],lab)[simno]},lab=lab,simno=simno)
  names(out) = paste0("Stock",1:ns)
  out
}

getfdat = function(PPD,lab,simno){
  ns = length(PPD);  nf = length(PPD[[1]]);  mpno = 1
  out=array(NA,c(ns,nf)); dimnames(out)[[1]] = paste0("Stock",1:ns);  dimnames(out)[[2]] = paste0("Fleet",1:ns)
  for(ss in 1:ns)for(ff in 1:nf)out[ss,ff] = slot(PPD[[ss]][[ff]][[mpno]],lab)[simno]
  out
}


# xmu = x10; ymu = y10; logslp = s10; range = 10
plotline=function(xmu,ymu,logslp,range=5,col="#00ff0099"){
  slp = exp(logslp)
  xs = xmu-c(range,1)+range/2+0.5
  ys = exp(c(1,range)*logslp)
  yscale = ys * mean(ymu)/mean(ys)
  lines(xs,yscale,col=col,lwd=3)
}


#' Extract time-series features for Indicator 2 data inputs
#'
#' @param ts A vector of positive real numbers, a data time series (e.g. 20 years of catch data) 
#' @param lab Character string. The label (code) of the time series (e.g., C for catches, I for CPUE index, ML for mean length, MV for variability in mean length, FM for fraction mature).
#' @param enp.mult Positive fraction. Effective number of parameters multiplier. A value of 0.5 means a loess smoother will be fitted that has half as many parameters as data points. Larger numbers are more flexible models (less smoothing).
#' @param plot Boolean. Should a simulated vs predicted plot / confusion matrix be presented for the independent testing dataset?
#' @param rnd Integer. Number of significant digits to plot. 
#' @param plot Boolean. Should a plot be presented? 
#' @examples
#' ts_features(1.5+sin(seq(1,12,length.out=40)),"C")
#' @author T. Carruthers
#' @export
ts_features=function(ts, lab="", enp.mult=0.2, rnd=3,plot=T){
  namy = deparse(substitute(ts))
  ny = length(ts)
  sts = smooth2(ts, enp.mult=enp.mult, plot = plot, plotname = namy)
  rel = sts[ny] / mean(sts)
  yr5 = ny-(4:0); yr10 = ny-(9:0); yr20 = ny-(19:0); yr40 = ny-(39:0)
  ts5 = ts[yr5]; ts10 = ts[yr10]; ts20 = ts[yr20]; ts40 = ts[yr40]
  g5<-slp3(ts5);   g10<-slp3(ts10); g20<-slp3(ts20); g40 = slp3(ts40)
  
  if(plot){
    abline(h=c(mean(sts),sts[ny]),lty=c(2,1),lwd=c(1,2),col="#ff000099")
    x5 = mean(yr5); x10 = mean(yr10); x20 = mean(yr20); x40 = mean(yr40)
    y5 = mean(ts5); y10 = mean(ts10); y20 = mean(ts20); y40 = mean(ts40)
 
    plotline(x5,y5,g5,5,'#00ff0080')
    plotline(x10,y10,g10,10,'#0000ff80')
    plotline(x20,y20,g20,20,'#99999980')
    plotline(x40,y40,g40,40,'#ff00ff90')
    legend('topright',legend=c(paste0(lab,"_rel = ",round(rel,rnd)),
                             paste0(lab,"_g5 = ",round(g5,rnd)),
                             paste0(lab,"_g10 = ",round(g10,rnd)),
                             paste0(lab,"_g20 = ",round(g20,rnd)),
                             paste0(lab,"_g40 = ",round(g40,rnd))
                             ), text.col=c("red","green","blue","darkgrey","#ff00ff"),bty="n")
  }
  
  outlist=list(rel,g5,g10,g20,g40)
  names(outlist)=  c(paste0(lab,"_rel"),paste0(lab,"_g5"), paste0(lab,"_g10"), paste0(lab,"_g20"), paste0(lab,"_g40"))
  outlist     
}

