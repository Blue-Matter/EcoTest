


ICCATtoGEO<-function(dato){   # Convert ICCAT format (corner closest to GMT/equator) to South West corner

  names(dato)[names(dato)=="GeoStrata"]<-"SquareTypeCode"   # Standardize nomenclature among datasets
  ICCATdat<-subset(dato,dato$SquareTypeCode!="none"&dato$SquareTypeCode!="ICCAT")

  Latadj<-as.numeric(array(unlist(strsplit(as.character(ICCATdat$SquareTypeCode),"x")),dim=c(2,nrow(ICCATdat)))[1,])
  Lonadj<-as.numeric(array(unlist(strsplit(as.character(ICCATdat$SquareTypeCode),"x")),dim=c(2,nrow(ICCATdat)))[2,])

  ICCATdat$Lon[ICCATdat$QuadID==1]<-ICCATdat$Lon[ICCATdat$QuadID==1]+(Lonadj[ICCATdat$QuadID==1]/2)
  ICCATdat$Lat[ICCATdat$QuadID==1]<-ICCATdat$Lat[ICCATdat$QuadID==1]+(Latadj[ICCATdat$QuadID==1]/2)
  ICCATdat$Lon[ICCATdat$QuadID==2]<-ICCATdat$Lon[ICCATdat$QuadID==2]+(Lonadj[ICCATdat$QuadID==2]/2)
  ICCATdat$Lat[ICCATdat$QuadID==2]<--ICCATdat$Lat[ICCATdat$QuadID==2]-(Latadj[ICCATdat$QuadID==2]/2)
  ICCATdat$Lon[ICCATdat$QuadID==3]<--ICCATdat$Lon[ICCATdat$QuadID==3]-(Lonadj[ICCATdat$QuadID==3]/2)
  ICCATdat$Lat[ICCATdat$QuadID==3]<--ICCATdat$Lat[ICCATdat$QuadID==3]-(Latadj[ICCATdat$QuadID==3]/2)
  ICCATdat$Lon[ICCATdat$QuadID==4]<--ICCATdat$Lon[ICCATdat$QuadID==4]-(Lonadj[ICCATdat$QuadID==4]/2)
  ICCATdat$Lat[ICCATdat$QuadID==4]<-ICCATdat$Lat[ICCATdat$QuadID==4]+(Latadj[ICCATdat$QuadID==4]/2)

  ICCATdat

}

ICCATtoGEO2<-function(ICCATdat){   # Convert ICCAT format (corner closest to GMT/equator) to South West corner

  ICCATdat<-subset(ICCATdat,!ICCATdat$GeoStrata%in%c("ICCAT","LatLon"))

  Latadj<-as.numeric(array(unlist(strsplit(as.character(ICCATdat$GeoStrata),"x")),dim=c(2,nrow(ICCATdat)))[1,])
  Lonadj<-as.numeric(array(unlist(strsplit(as.character(ICCATdat$GeoStrata),"x")),dim=c(2,nrow(ICCATdat)))[2,])

  ICCATdat$Lon[ICCATdat$QuadID==1]<-ICCATdat$Lon[ICCATdat$QuadID==1]+(Lonadj[ICCATdat$QuadID==1]/2)
  ICCATdat$Lat[ICCATdat$QuadID==1]<-ICCATdat$Lat[ICCATdat$QuadID==1]+(Latadj[ICCATdat$QuadID==1]/2)
  ICCATdat$Lon[ICCATdat$QuadID==2]<-ICCATdat$Lon[ICCATdat$QuadID==2]+(Lonadj[ICCATdat$QuadID==2]/2)
  ICCATdat$Lat[ICCATdat$QuadID==2]<--ICCATdat$Lat[ICCATdat$QuadID==2]-(Latadj[ICCATdat$QuadID==2]/2)
  ICCATdat$Lon[ICCATdat$QuadID==3]<--ICCATdat$Lon[ICCATdat$QuadID==3]-(Lonadj[ICCATdat$QuadID==3]/2)
  ICCATdat$Lat[ICCATdat$QuadID==3]<--ICCATdat$Lat[ICCATdat$QuadID==3]-(Latadj[ICCATdat$QuadID==3]/2)
  ICCATdat$Lon[ICCATdat$QuadID==4]<--ICCATdat$Lon[ICCATdat$QuadID==4]-(Lonadj[ICCATdat$QuadID==4]/2)
  ICCATdat$Lat[ICCATdat$QuadID==4]<-ICCATdat$Lat[ICCATdat$QuadID==4]+(Latadj[ICCATdat$QuadID==4]/2)

  ICCATdat

}


assign_quarter<-function(dat){

  quarter<-c(rep(1:4,each=3),1:4)
  Subyear<-quarter[dat$TimePeriodID]
  if(!"TimePeriodID"%in%names(dat)&"TimeCatch"%in%names(dat)) Subyear<-quarter[dat$TimeCatch]
  cond<-!is.na(Subyear)
  Subyear<-Subyear[cond]
  cbind(dat[cond,],Subyear)

}

assign_year<-function(dat,years){

  dat<-dat[dat$YearC>(years[1]-1)&dat$YearC<(years[2]+1),]
  Year<-as.numeric(dat$YearC)-years[1]+1
  cbind(dat,Year)

}

assign_spatial_strata <- function(dat,res){

  dat$Lat = floor(dat$Lat/res)*res
  dat$Lon = floor(dat$Lon/res)*res
  uni = paste(dat$Lat,dat$Lon,sep="_")
  Area = match(uni,unique(uni))
  cbind(dat,Area)

}

extract_CE = function(dat, SpecCode, model = 'log(CPUE)~Y + Q + A + F'){

  if(!(SpecCode %in% names(dat))) stop(paste0("Species code ", SpecCode, " not a named column in the dataset"))
  Effort = dat$Eff1
  EType = dat$Eff1Type
  CUnit = dat$CatchUnit
  Area = dat$Area
  Fleet = dat$FleetCode
  Subyear = dat$Subyear
  Year = dat$Year
  Catch = dat[[SpecCode]]
  CPUE = data.frame(Year = Year, Subyear = Subyear, Area = Area, Fleet = Fleet, CUnit = CUnit, EType = EType, Catch = Catch, Effort = Effort)
  Check = apply(CPUE[,1:6],2,unique)

  CPUE = CPUE[CPUE$CUnit != "..",]
  CPUE = CPUE[!is.na(CPUE$EType),]
  CPUE = CPUE[CPUE$EType == "NO.HOOKS",]
  CPUE = CPUE[CPUE$Effort>0,]
  bylist = list(Year = CPUE$Year, SubYear = CPUE$Subyear, Area = CPUE$Area, Fleet = CPUE$Fleet, CUnit = CPUE$CUnit)

  Cagg = aggregate(CPUE$Catch, by=bylist, sum)
  names(Cagg)[ncol(Cagg)] = "Catch"
  Eagg = aggregate(CPUE$Effort, by=bylist, sum)
  names(Eagg)[ncol(Eagg)] = "Effort"
  keep = Cagg$Catch > 0

  CPUEagg = cbind(Cagg[keep,], Eagg$Effort[keep])
  names(CPUEagg)[ncol(CPUEagg)] = "Effort"

  if(nrow(CPUEagg) < 100){
    warning("Fewer than 100 aggregated CPUE records, standardization not attempted")
    return(NULL)
  }else{
    cpue = CPUEagg$Catch / CPUEagg$Effort
    CPUEagg = cbind(CPUEagg,cpue)
    for(cc in 1:5) CPUEagg[,cc] = as.factor(CPUEagg[,cc])
    CPUEagg$cpue = as.numeric(CPUEagg$cpue)
    names(CPUEagg) = c("Y","Q","A","F","U","C","E","CPUE")

    # model = 'log(CPUE)~Y + Q + F'
    mod = lm(model, data = CPUEagg)
    preddat = data.frame(Y = unique(CPUEagg$Y), Q = CPUEagg$Q[1], A = CPUEagg$Q[1], F = CPUEagg$F[1], U = CPUEagg$U[1])
    pred=predict(mod,newdata = preddat,se.fit=T)

    #ymarg = margins::margins_summary(mod,variables = "Y")
    scpue = data.frame(Year = preddat$Y, Index = exp(pred$fit)/mean(exp(pred$fit)), SE = pred$se.fit)

    if(any(is.na(scpue$SE))){
      print(scpue)
      warning("One or more year effects could not be quantified with standard error, standardization halted, simplify the model")
      return(NULL)
    }else{
      #yrs = sapply(ymarg$factor,function(x)strsplit(x,"Y")[[1]][2])
      #scpue = data.frame(Year = as.numeric(yrs), Index = exp(ymarg$AME), CV = ymarg$SE)
      plot(scpue$Year, scpue$Index, pch=19, ylab="Index",xlab="");grid()
      mtext(model,3,cex=0.5,line=0.3)
      plot(scpue$Year, scpue$SE, pch=19, ylab="Index CV",xlab="");grid()
      mtext(model,3,cex=0.5,line=0.3)
      return(scpue)
    }
  }
}

