
doscale = function(mat, mu, sd)  t((t(mat)-mu)/sd)

# input = Shark_1_data; nodes = c(6, 3); nepoch = 20; plot = T; model=NULL; model_savefile=NA; validation_split=0.1; test_split = 0.1; lev=c(0.5,1); seed=1

#' Training a neural network (indicator 3) for features specified by user data.
#'
#' A function that examines the user specified data, subsets the training dataset for those features (data types), standardizes the user data and trains a bespoke neural network indicator for the user specified data types.
#'
#' @param input Either a list object with a single dataset (e.g. Shark_1_data) or a list object npeels long with a dataset for retrospective analysis (e.g. Shark_1_retro).
#' @param TD The training dataset (requires EcoTestData) and created by prepare_TD()
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
#' train_ind_3(Shark_1_data)
#' @author T. Carruthers
#' @export
train_ind_3=function(input, TD, nodes = c(10, 5), nepoch = 20, plot = T, model = NULL, model_savefile = NA,
                    validation_split = 0.1, test_split = 0.1, lev = c(0.5, 1), seed = 1){

  if(class(input)=="ETData"){
    data = input[[1]]
    cat(paste0("Dataset has the following list items: ", paste(names(input),collapse=", ")), " \n")
    npeel = sum(grepl("peel",names(input)))-1
    cat(paste0("Dataset provided with an additional ",npeel, " retrospective datasets \n"))
    if("Brel" %in% names(input)) cat("Dataset provided with a time series of estimated B relative to BMSY (Brel) \n")
  }else{
    stop("You must provide an input object of class 'ETData'")
  }

  cat(paste0("Dataset provided with ",ncol(data), " features: ", paste0(colnames(data),collapse=", ")," \n"))

  lookupnam = c("Brel",names(data))
  dind = match(lookupnam,names(TD))
  if(any(is.na(dind))) stop(paste0("At least one of your input dataset names does not match the training dataset: ",paste0(lookupnam[is.na(dind)],collapse=", ")," \n"))
  td = TD[,dind] # training data is reduced to the size of the real data and ordered

  # Automatic network adjustment to prevent overparameterizatoin
  #if(autoadjust){
  #  nfeat = ncol(td)-1
  #  if(nfeat<(nodes[1]*1.5)){
  #    orgspec = nodes[1]
  #    nodes[1] = floor(nfeat/1.5)
  #    cat(paste0("Adjusting neural network design to prevent overparameterization. Dataset has ",nfeat, " features. The argument nodes[1] determines the width of first layer and was set to ", orgspec," nodes. This is more than 2/3 the number of features and has been automatically reduced to ",nodes[1], " nodes to prevent overparameterization. \n"))
  #    if(nodes[1]<nodes[2]){
  #      nodes[2] = nodes[1]-1
  #      cat(paste0("Adjusting neural network design to prevent overparameterization. Width of second layer reduced to one less than first layer. \n"))
  #    }
  #  }
  #}

  nr<-nrow(td);   nc<-ncol(td)
  p_nontrain = test_split
  train_switch = sample(c(TRUE,FALSE),nr,replace=T,prob=c(1-p_nontrain,p_nontrain))

  # Training dataset
  cat("Standardizing non time-series training dataset features to mean 0, stdev 1. \n")

  train_target<-td[train_switch, 1]
  train<-td[train_switch, 2:nc]
  mu = apply(train, 2, mean)
  sd = apply(train, 2, sd)
  namT = names(train)
  nonts = !(grepl("C_",namT) | grepl("I_",namT) | grepl("ML_",namT) | grepl("MV_",namT) |  grepl("FM_",namT) | grepl("MA_",namT))
  strain = as.matrix(train)
  strain[,nonts]= doscale(train[,nonts], mu[nonts], sd[nonts])

  # Empirical dataset
  cat("Using training dataset mean and stdev to standardize submitted data. \n")
  sdata = input[[1]]
  sdata[,nonts] = doscale(sdata[,nonts], mu[nonts], sd[nonts])

  if(npeel > 0){
    for(pp in 0:npeel) input[[paste0("peel_",pp)]][,nonts] = doscale(input[[paste0("peel_",pp)]][,nonts], mu[nonts], sd[nonts])
  }
  #train_data_comp(train, strain, tdata, sdata)

  # Testing dataset
  cat(paste0(test_split*100, "% of training dataset randomly selected for testing of neural network after fitting (test_split). \n"))
  cat(paste0("Of the remaining training dataset (", (1-test_split)*100, "%), ", validation_split*100, "% are randomly selected for validation during fitting (validation_split). \n"))

  testy_target<-td[!train_switch,1]
  testy<-td[!train_switch,2:nc]
  testy[,nonts] = doscale(testy[,nonts], mu[nonts], sd[nonts])

  cat("Initializing and compiling neural network. \n")
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

  cat("Training neural network. \n")
  history <- model %>% fit(strain, train_target,
                           epochs = nepoch,
                           validation_split = validation_split,
                           verbose = 2)

  # predicted vs simulated

  pred = exp((model %>% predict(as.matrix(testy)))[,1])
  sim = exp(testy_target)

  tab = NN_fit(sim, pred, history, lev=lev, addpow=T,nepoch=nepoch,plot=F)

  MAE = history$metrics$val_MAE[nepoch]

  if(!is.na(model_savefile))save_model(model,filepath = model_savefile, overwrite=T)

  if(!all(colnames(train)==names(sdata))) stop("Internal Error: something weird happened, training dataset does not match size of submitted data object")

  outlist=list(MAE = MAE, sim=sim, pred=pred, grid=grid, tab = tab, model=model, strain=strain,
               train_target=train_target, nepoch = nepoch, r2 = cor(sim,pred)^2, history = history,
               lev=lev, mu = mu, sd = sd, input = input)

  if(plot){
    cat("Plotting simulated stock status (from the testing dataset) against that predicted by the neural network (given the features of the testing dataset). \n")
    perf_plot(outlist)
  }

  outlist

}


