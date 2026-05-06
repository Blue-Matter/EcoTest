


#' Training a neural network for features specified by user data.
#'
#' A function that examines the user specified data, subsets the training dataset for those features (data types), standardizes the user data and trains a bespoke neural network indicator for the user specified data types.
#'
#' @param TD A training dataset of nsimulation rows and ncolumns features (first column is the response variable, e.g. stock status)
#' @param nodes A vector of 2 positions. Positive integers. There are few general rules for neural network design. However one rule of thumb is that the width of the first and second layer of the neural network. In general, the first layer should be smaller than the number of features (data types) in the user specified dataset. The second layer should be smaller than the first layer. Check for overparameterization in the fitting plots (e.g. a validation MAE that is substantially higher than the training MAE).
#' @param nepoch Positive integer. The number of iterations (passes of the back propagation algorithm) used in the training of the neural network. Should be sufficient that loss (or MAE) has stabilised.
#' @param plot Boolean. Should a simulated vs predicted plot / confusion matrix be presented for the independent testing dataset?
#' @param model A keras model. Optional. You may wish to keep training a partially trained model (needs to be compatible with shape of dataset).
#' @param model_savefile Character vector. Optional. A file location and name for saving the fitted keras model. Must have a .keras type. E.g. C:/temp/mymodel.keras.
#' @param validation_split Fraction. The proportion of the dataset that will be used for cross-validation as the neural network is training.
#' @param test_split Fraction. The proportion of the training dataset reserved as an independent testing dataset (those data of the simulated - predicted plot and confusion matrix).
#' @param lev A vector of 2 positions. The location of the breaks of the confusion matrix.
#' @param seed Numeric value. Random seed for stochastic draws of validation and testing datasets.
#' @examples
#' train_NN(TD)
#' @author T. Carruthers
#' @export
train_NN = function(TD, nodes = c(4, 2), nepoch = 20, plot = T, model=NULL, model_savefile=NA, 
                    validation_split=0.1, test_split = 0.1, lev=c(0.5,1), seed=1){
  set.seed(seed)
  nr<-nrow(TD)
  nc<-ncol(TD)
  p_nontrain = test_split
  train_switch = sample(c(TRUE,FALSE),nr,replace=T,prob=c(1-p_nontrain,p_nontrain))
  
  # Training dataset
  train<-TD[train_switch, 2:nc]
  mu = apply(train, 2, mean)
  sd = apply(train, 2, sd)
  train = scale(train, center=mu, scale=sd)
  train_target<-TD[train_switch, 1]
  
  # Test
  #nnt = sum(!train_switch)
  #breakpt = ceiling(nnt*validation_split/p_nontrain)
  #valnos = allnos[1:breakpt]
  #testnos = allnos[(breakpt+1):length(allnos)]
  ind = (1:nr)[!train_switch]
  testy<-TD[ind,2:nc]
  testy=scale(testy,center=mu,scale=sd)
  testy_target<-TD[ind,1]
  
  if(is.null(model))model = keras_model_sequential()
  
  model %>%
    layer_dense(units = nodes[1], activation = "relu") %>%
    layer_dense(units = nodes[2], activation = "relu") %>%
    layer_dense(units = 1)
  
  model %>%
    compile(
      loss = 'MSE', #mean_squared_error',#"mse",#mae",#"mse",
      optimizer =  optimizer_adam(), #'adam',#optimizer_rmsprop(),'rmsprop',
      metrics = "MAE"
    )
  
  history <- model %>% fit(train, train_target,
                           epochs = nepoch,
                           validation_split = validation_split,
                           verbose = 2
  )
  
  
  pred = exp((model %>% predict(testy))[,1])
  sim = exp(testy_target)
  tab = NN_fit(sim, pred, history, lev=lev, addpow=T,nepoch=nepoch,plot=F)

  
  MAE = history$metrics$val_MAE[nepoch]
  
  if(!is.na(model_savefile))save_model(model,filepath = model_savefile, overwrite=T)
  
  outlist=list(MAE = MAE, sim=sim, pred=pred, grid=grid, tab = tab, model=model, train=train, 
       train_target=train_target, nepoch = nepoch, r2 = cor(sim,pred)^2, history = history,
       lev=lev, mu = mu, sd = sd)
  
  if(plot)(pred_plot(outlist))
  
  outlist
  
}



cleandat_2 = function(dat,resprange=c(-3.689,0.693)){
  # test = (x==0); (1:length(test))[test]
  goodrow = function(x)(sum(is.na(x)| x == 'Inf')==0)&x[1]>resprange[1] & x[1]<resprange[2]
  keep = apply(dat,1,goodrow)
  dat[keep,]
}

stock_subsetter = function(cdat, sno,stock_in, maxstock = 3){
  if(length(stock_in)<3){
    dnames = names(cdat)
    stock_in_nam = paste0("s",stock_in)
    all_stock=paste0("s",1:maxstock)
    resp_stock = paste0("s",sno)
    sel_stock = unique(c(resp_stock,stock_in_nam))
    rem_stock = all_stock[!(all_stock%in%sel_stock)]
    ns = length(rem_stock)
    nn = length(dnames)
    condgrid = array(F,c(ns,nn))
    for(i in 1:ns)condgrid[i,]= grepl(rem_stock[i],dnames)
    cdat = cdat[,!apply(condgrid,2,any)] # should the input be there based on the stocks?
  }
  cdat
}

fleet_subsetter = function(cdat, fleet_in, dnames, maxfleet = 3){
  if(length(fleet_in)<3){
    dnames = names(cdat)
    fleet_in_nam = paste0("f",fleet_in)
    all_fleet=paste0("f",1:maxfleet)
    rem_fleet = all_fleet[!(all_fleet%in%fleet_in_nam)]
    nf = length(rem_fleet)
    nn = length(dnames)
    condgrid = array(F,c(nf,nn))
    for(i in 1:nf)condgrid[i,]= grepl(rem_fleet[i],dnames)
    cdat = cdat[,!apply(condgrid, 2, any)] # should the input be there based on the stocks?
  }
  cdat
}


data_subsetter = function(dat, inc_Irel=T, inc_I = T, inc_CR = T, inc_CAL = T, inc_CAA = T){
  if(!inc_Irel) dat = dat[,!grepl('I_rel',names(dat))]
  if(!inc_I)    dat = dat[,!grepl("I",names(dat))]
  if(!inc_CR)   dat = dat[,!(grepl("CR",names(dat)) | grepl("CC",names(dat)))]
  if(!inc_CAL)  dat = dat[,!(grepl("ML",names(dat))  | grepl("MV",names(dat)) | grepl("FM", names(dat)))]
  if(!inc_CAA)  dat = dat[,!grepl("MA",names(dat))]
  dat
}

resp_subsetter = function(cdat, sno, isBrel){
  dnames2 = names(cdat)
  Brelcols = grepl("Brel",dnames2)
  ns = sum(Brelcols) # os = (1:ns)[!((1:ns)%in%sno)]
  keepcol = grepl(paste0("Brel_s",sno),dnames2)
  Res = cdat[,keepcol]
  if(isBrel){Res = ((ns-1) * Res) /(apply(cdat[,Brelcols],1,sum)-Res)}# level relative to mean level of other stocks
  Res = log(Res)
  dat = cbind(Res,cdat[,!Brelcols])
  dat
} 

rem_const = function(dat){
  isconst = apply(dat,2,sd)<1E-10 
  if(sum(isconst,na.rm=T)>0){
    cat(paste(paste(names(dat)[isconst],collapse=", "), "dropped for sd < 1E-10 \n"))
    dat = dat[,!isconst]
  }
  dat
}

dolog_2=function(dat){
  dolog(dat, types = c("I_rel","C_rel","FM_cur","FM_rel","ML_cur",
                             "ML_rel","MV_cur","MV_rel","ML_Linf", "MA_cur", "MA_rel",
                             "ML_L50","CR_rel","maxa","K","M_K","L50_Linf","L5_L50", 
                             "LFS_L50","Csd","CF","Isd"))
}
                          
                         
makerawdata_2 = function(allout, sno=1, isBrel = F,  
                       inc_Irel = T, inc_I = T, inc_CR = T, inc_CAL = T, inc_CAA = T,
                       stock_in = NA, fleet_in = NA, Brange = c(0.025,4.5)){
  # sno=1; isBrel = F; clean = T;  inc_Irel = T; inc_I = T;  inc_CR = T; stock_in = NA; inc_CAL = T; inc_CAA = T
  cdat = as.data.frame(rbindlist(allout))
  dnames = names(cdat)
  Brelcols = grepl("Brel",dnames)
  stocks = sapply(dnames[Brelcols],function(x)strsplit(x,"_")[[1]][2])
  if(is.na(stock_in[1]))stock_in = stocks[stocks!=sno]
  if(is.na(fleet_in[1]))fleet_in = 1:3
  
  cdat = stock_subsetter(cdat, sno, stock_in) # keep only those listed in sno and stock_in
  cdat = fleet_subsetter(cdat, fleet_in)                # keep only those fleets in nfleet
  cdat = data_subsetter(cdat,  inc_Irel, inc_I, inc_CR, inc_CAL, inc_CAA) # only include those data that are specified to be available
  dat = resp_subsetter(cdat, sno, isBrel)                     # gets the right response variable according to sno
  dat = dolog_2(dat)                                  # log imperfect fractions 
  dat = dologit(dat,types = "VML")                    # logit proportions (but rescaled 0.05 - 0.95 prior to logit)
  dat = cleandat_2(dat,resprange = log(Brange))       # clean NAs and Infs 
  dat = rem_const(dat)                                # remove any independent variables with no variability (constant over simulations)
  dat
  
}


