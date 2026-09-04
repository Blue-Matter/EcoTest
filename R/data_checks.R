# Checks


train_data_comp = function(Ind){
  # transformed data
  strain = Ind$strain
  mintrain = apply(strain,2,min)
  maxtrain = apply(strain,2,max)
  input = Ind$input
  datamu = apply(input[[1]],2,mean)

  data.frame(
    mu_trainorg = round(apply(strain,2,mean),3),
    min_trainorg = round(mintrain,3),
    max_trainorg = round(maxtrain,3),
    mu_data = round(datamu,3),
    inrange = datamu>mintrain & datamu<maxtrain)
}


check_train_data = function(data, TD){
  # transformed data
  lookupnam = names(data[[1]])
  dind = match(lookupnam,names(TD))
  if(any(is.na(dind))) stop(paste0("At least one of your input dataset names does not match the training dataset: ",paste0(lookupnam[is.na(dind)],collapse=", ")," \n"))
  strain = TD[,dind] # training data is reduced to the size of the real data and ordered
  mintrain = apply(strain,2,min)
  maxtrain = apply(strain,2,max)
  datamu = apply(data[[1]],2,mean)

  all = data.frame(
    mu_trainorg = round(apply(strain,2,mean),3),
    min_trainorg = round(mintrain,3),
    max_trainorg = round(maxtrain,3),
    mu_data = round(datamu,3),
    inrange = datamu>mintrain & datamu<maxtrain)

  fail = all[!all$inrange,]
  nf = sum(!all$inrange)
  if(nf>0)cat(paste0("WARNING: ",nf, " features of the dataset are outside the range of the training dataset. Check output list item 'fail' for info"))
  if(nf==0)cat("SUCCESS! All features (mean across samples) within range of training dataset \n")

  retslots = grep("peel",names(data))
  nfeats = sapply(retslots,function(x,data)ncol(data[[x]]),data=data)
  if(length(nfeats)>1)if(length(unique(nfeats))!=1)cat("WARNING: your retro peels do not have the same number of features \n")

  list(all = all, fail = fail, nfeats = nfeats)
}

data_comp=function(td, data){

  mutd = apply(td[,2:ncol(td)],2,mean)
  mintd = apply(td[,2:ncol(td)],2,quantile,p=0.025)
  maxtd = apply(td[,2:ncol(td)],2,quantile,p=0.975)
  mudata = apply(data[[1]],2,mean)

  data.frame(mutd=round(mutd,3),
             mintd=round(mintd,3),
             maxtd=round(maxtd,3),
             mudata=round(mudata,3),
             inside = mudata<maxtd & mudata>mintd)

}
