# ------ General functions ------

# function that create all possible predictor combinations and puts them in a list. 
# The target variable is always the first predictor

predictor_combinator <- function(predictors, temperature.included = T, which=NULL, rmTfrompred=T){
  # predictors <- predictors[!is.element(predictors, c(target,"temperature"))]
  if(rmTfrompred) predictors <- predictors[!is.element(predictors, c("temperature"))]
  if(is.null(which)){
    num.predictors <- 1:length(predictors)
  } else {
    num.predictors <- which
  }
  l <- sapply(num.predictors, function(x) {
    comb <- combn(predictors, x, simplify = F)
    if(length(comb)>200) comb <- sample(comb,200,replace = F)
    comb
  })
  
  if(length(which)==1) l <- list(l)
  
  if(temperature.included==T) {
    l2 <- lapply(l, function(lis){lapply(lis, append, values="temperature")})
    l <- c(l,l2)
    l <- do.call(list, unlist(l, recursive=FALSE))
    l <- append(l, "temperature", after=0)
  } else {
    l <- do.call(list, unlist(l, recursive=FALSE))
  }
  # names(l) <- paste0(l)
  return(l)
}

### create libraries in correct format
create_segment <- function(vec){
  vec <- sort(vec)
  diffs <- c(0,diff(vec))
  vec.segs <- matrix(ncol = 2, nrow = length(vec))
  row <- 1
  for(i in 1:length(vec)) {
    if(i==1) {
      vec.segs[row,1] <- vec[i]
    } else if(diffs[i]!=1){
      row <- row + 1
      vec.segs[row,1] <- vec[i]
      vec.segs[row-1,2] <- vec[i-1] + 1 
    } 
    
    if(i==length(vec)) {vec.segs[row,2] <- vec[i] + 1}
  }
  
  vec.segs <- vec.segs %>% na.omit()
  return(vec.segs)
}


### create libraries: what data is used to train and what to predict. 
### use list[[1]] as training, and list[[4]] as pred library
create_libraries <- function(nrows, percent_predict=0.25, seed){
  set.seed(seed)
  rows <- 2:(nrows-1)
  pred <- sample(rows, size = floor(percent_predict*length(rows)))
  lib <- rows[!(rows %in% pred)]
  
  pred.mat <- create_segment(pred)
  lib.mat <- create_segment(lib)
  
  return(list(lib.mat,pred.mat,sort(lib),sort(pred)))
}

# Function that handles the data by creating a dataframe with all information wanted.

mv_housekeeping <- function(predictor_combinations, output_list, data, 
                            target, E){
  t_excluded <- sapply(predictor_combinations, function(y){
    !(target %in% y)
  })
  # if(isTRUE(t_excluded)){
  #   predictor_combinations <- lapply(predictor_combinations, function(x) x[-1])
  # }
  num_predictors <- sapply(predictor_combinations, function(x) length(x))
  names <- unname(sapply(predictor_combinations, function(y) paste(y, collapse = " ")))
  k <- sapply(output_list, "[[", 2)
  RMSE <- sapply(output_list, "[[", 1)
  df <- data.frame(Target=target, 
                   predictor_combination = names,
                   num_pred = num_predictors, 
                   target_excluded = t_excluded,
                   RMSE = RMSE,
                   k = k,
                   E=E)
  return(df)
}

# ------ K-fold cross-validation ------

# wrapper function for the Multiview function

mv_wrapper <- function(data, k=0, columns, target, E=3, max_lag=3, num_neighbors=0, 
                       excludeTarget=F, lib=NULL, pred=NULL, lib.prop=0.75, Tp=1){
  data <- cbind(t=1:nrow(data),data)
  if(is.null(lib)){
    lib <- c(1, floor(nrow(data)*lib.prop))
  }
  if(is.null(pred)){
    pred_actual <- (floor(nrow(data)*lib.prop) + 1):NROW(data)
    pred <- c(1,nrow(data))
  } else {
    pred_actual <- pred[1]:pred[2]
    pred <- c(1,nrow(data))
  }
  
  m <- Multiview(dataFrame = data, lib = lib, pred = pred, E = max_lag, D = E, 
                 knn = num_neighbors, multiview = k, columns = columns, Tp=Tp,
                 target = target, excludeTarget=excludeTarget, numThreads = 1)
  m$Predictions <- m$Predictions[m$Predictions$t %in% pred_actual,]
  return(m)
}

# function that calculates the rmse of all models in list and then only keeps the
# model with the best (smallest) rmse value.
# Returns best model, best RMSE and corresponding k value (=number of views used)

model_selector <- function(mv_wrapper_list,k){
  prediction_list <- lapply(mv_wrapper_list, "[[", 2)
  rmse_list <- lapply(prediction_list, 
                      function(x) sqrt(mean((x$Observations-x$Predictions)^2, 
                                            na.rm = T)))
  idx <- which.min(rmse_list)
  return(list(rmse_list[[idx]], k[idx]))
}

# Main function
# function that for a given dataset, predictor list and target variable calculates
# the range of k values to be tested and then. It calls the mv_wrapper function for 
# each value of k. Is saves the return of the mv_wrapper function as entries in a 
# list and then returns the list. It calls the mv_wrapper function twice: once for 
# when the target is a predictor and once for when it is not.

model_fitter <- function(data, target, predictor_combinations, max_lag=3, E=3,
                         num_neighbors=0, k=0, lib=NULL, pred=NULL, Tp=1, ...){
  
  output_list <- unname(lapply(predictor_combinations, function(y){
    
    y <- unlist(y)
    t_excluded = !(target %in% y)
    kmax <- choose(length(y)*max_lag, E) 
    k_input <- k[k<=kmax]
    input_data <- data %>% ungroup() %>% dplyr::select(all_of(y),target)
    if(t_excluded==T) y <- append(y, target)
    columns <- paste(y, collapse = " ")
    mv_wrapper_list <- lapply(k_input, function(x)
      mv_wrapper(data = input_data, k = x, columns = columns, target = target, E = E,
                 max_lag = max_lag, num_neighbors = num_neighbors, Tp=Tp,
                 excludeTarget = t_excluded, lib=lib, pred=pred))
    model_selector(mv_wrapper_list,k_input)}))
  
  df <- mv_housekeeping(predictor_combinations, output_list, data, target, E)
  return(df)
}



# function that prepares the data that goes into the model_fitter function. It calls the
# model_fitter function for each bottle ID in parallel and then handles the return of the 
# model_fitter function by creating a list with two entries. The first entry is the merged
# dataframe across all bottles for the target species and includes the information about 
# predictors combinations and rmse values. The second entry # is a list with 18 entries
# (one for each bottle ID), and each of these 18 entries contains the best models for each 
# predictor combination tested.

model_fitter_wrapper <- function(whole.data, targets, num.clusters, max_lag=3, 
                                 E=3, num_neighbors=0, k=0, predictor_combinations, 
                                 lib=NULL, pred=NULL, predictors, Tp=1, 
                                 OptimalEsDf=NULL, ...){
  
  model_fitter_output_list <- mclapply(targets, function(target){
    if(!is.null(OptimalEsDf)){
      E <- OptimalEsDf[OptimalEsDf$Target==target,"E"]
    }
    write(paste0("Now running target: ", target, "\n"), write.file, append = T)
    old <- Sys.time()
    out <- model_fitter(data = whole.data, target = target, 
                        predictor_combinations = predictor_combinations, max_lag = max_lag,
                        E = E, num_neighbors = num_neighbors, k = k, 
                        lib=lib, pred=pred, Tp=Tp, ...)
    write(paste0("Time needed for target: ", Sys.time() - old, "\n"), write.file, append = T)
    out
  }, mc.cores = num.clusters, mc.cleanup = T, mc.allow.recursive = F)
  
  whole_df <- do.call("rbind", model_fitter_output_list)
  
  whole_df$predictor_combination <- as.character(whole_df$predictor_combination)
  temp <- matrix(F, nrow = nrow(whole_df), ncol = length(predictors))
  colnames(temp) <- predictors
  whole_df <- cbind(whole_df,temp)
  
  whole_df <- as.data.frame(t(apply(whole_df, 1, function(row) { 
    str <- unlist(strsplit(row["predictor_combination"], " "))
    row[str] <- T
    row
  })))
  
  whole_df$num_pred <- as.numeric(whole_df$num_pred)
  whole_df$RMSE <- as.numeric(whole_df$RMSE)
  whole_df$k <- as.numeric(whole_df$k)
  whole_df$E <- as.numeric(whole_df$E)
  for(i in predictors){  whole_df[,i] <- as.logical(whole_df[,i])}
  
  return(whole_df)
}






K_fold_libraries <- function(data, step=12){
  pred <- list()
  lib <- list()
  
  for(i in seq(1,nrow(data),step)){
    if(i==1){
      p <- c(i, i+step-1)
      l <- c(i+step, nrow(data))
    } else if((nrow(data)-i)<step) {
      p <- c(i, nrow(data))
      l <- c(1,i-1)
    } else {
      p <- c(i,i+step-1)
      l <- matrix(rbind(c(1,i-1),c(i+step,nrow(data))), nrow = 2)
    }
    pred <- append(pred, list(p))
    lib <- append(lib, list(l))
  }
  return(list(pred=pred,lib=lib))
}




####### Univariate Simplex forecasting with optimal E

UniSimplex_determineE <- function(data, Target, maxE, lib, pred) {
  
  mat <- EmbedDimension(dataFrame = data, lib = lib, pred = pred,
                        target = Target, columns = Target, maxE = maxE, showPlot = F)
  best <- mat %>%
    mutate(rho = round(rho,3)) %>%
    filter(rho == max(rho))
  return(best$E)
}

simplexForecasting <- function(data, target, num.clusters, 
                               lib=NULL, pred=NULL, maxE){
  
  data.list <- split(x = data, f = data$ID, drop = T)
  
  model_fitter_output_list <- mclapply(data.list, function(df){
    if(is.null(lib)){
      lib <- c(1, floor(NROW(df)*2/3))
    }
    if(is.null(pred)){
      pred <- c(floor(NROW(df)*2/3) + 1, NROW(df))
    }
    
    ID <- unique(df$ID)
    treat <- unique(df$treat)
    treat2 <- unique(df$treat2)
    treat3 <- unique(df$treat3)
    df <- df[,c("day",target)]
    
    E <- UniSimplex_determineE(df, target, maxE, lib, pred)
    forecast <- Simplex(dataFrame = df, lib = lib, pred = pred, target = target,
                        columns = target, E = E)
    rmse <- sqrt(mean((forecast$Observations-forecast$Predictions)^2, na.rm = T))
    data.frame(Target=target, RMSE=rmse, E=E, ID=ID,
               treat = treat, treat2 = treat2, 
               treat3 = treat3)
  }, mc.cores = num.clusters, mc.cleanup = T, mc.allow.recursive = F)
  
  whole_df <- do.call("rbind", model_fitter_output_list)
  return(whole_df)
}

# ------ LOOCV cross-validation ------

loocv_libraries <- function(data){
  
  SmplIntrvl <- mean(diff(data$time))
  forecastSteps <- 12/SmplIntrvl
  pred <- list()
  lib <- list()
  
  for(i in 2:nrow(data)){
    
    if(SmplIntrvl==1){
      if(i<71){
        p <- c(i, i+1)
        if(i<3){
          l <- matrix(rbind(c(i+1,i+forecastSteps-1),
                            c(i+forecastSteps+1,nrow(data))), ncol =  2)
        } else if(i<69) {
          l <- matrix(rbind(c(1,i-1),
                            c(i+1,i+forecastSteps-1),
                            c(i+forecastSteps+1,nrow(data))), ncol =  2)
        } else {
          l <- matrix(rbind(c(1,i-1),
                            c(i+1,i+forecastSteps-1)), ncol =  2)        
        }
      } else break
    }
    
    if(SmplIntrvl==3){
      if(i<79){
        p <- c(i, i+1)
        if(i<3){
          l <- matrix(rbind(c(i+1,i+forecastSteps-1),
                            c(i+forecastSteps+1,nrow(data))), ncol =  2)
        } else if(i<77) {
          l <- matrix(rbind(c(1,i-1),
                            c(i+1,i+forecastSteps-1),
                            c(i+forecastSteps+1,nrow(data))), ncol =  2)
        } else {
          l <- matrix(rbind(c(1,i-1),
                            c(i+1,i+forecastSteps-1)), ncol =  2)        
        }
      } else break
    }
    
    if(SmplIntrvl==6){
      if(i<81){
        p <- c(i, i+1)
        if(i<3){
          l <- matrix(rbind(c(i+3,nrow(data))), ncol =  2)
        } else if(i<79) {
          l <- matrix(rbind(c(1,i-1),
                            c(i+3,nrow(data))), ncol =  2)
        } else {
          l <- matrix(rbind(c(1,i-1)), ncol =  2)        
        }
      } else break
    }
    
    if(SmplIntrvl==12){
      if(i<82){
        p <- c(i, i+1)
        if(i<3){
          l <- matrix(rbind(c(i+2,nrow(data))), ncol =  2)
        } else if(i<80) {
          l <- matrix(rbind(c(1,i-1),
                            c(i+2,nrow(data))), ncol =  2)
        } else {
          l <- matrix(rbind(c(1,i-1)), ncol =  2)        
        }
      } else break
    }
    
    pred <- append(pred, list(p))
    lib <- append(lib, list(l))
    
  }
  return(list(pred=pred,lib=lib))
}







# wrapper function for the Multiview function

mv_wrapper_loocv <- function(data, k=0, columns, target, E=3, max_lag=3, num_neighbors=0, 
                             excludeTarget=F, lib, pred, Tp=1){
  
  data <- cbind(t=1:nrow(data),data)
  m <- lapply(1:length(pred), function(i){
    L <- Multiview(dataFrame = data, lib = lib[[i]], pred = pred[[i]], E = max_lag, D = E, 
                  knn = num_neighbors, multiview = k, columns = columns, Tp=Tp,
                  target = target, excludeTarget=excludeTarget, numThreads = 1)
    predictions <- L$Predictions
    predictions[nrow(predictions)-1,]
  })
  
  m <- do.call("rbind",m)
  return(m)
}

# function that calculates the rmse of all models in list and then only keeps the
# model with the best (smallest) rmse value.
# Returns best model, best RMSE and corresponding k value (=number of views used)

model_selector_loocv <- function(mv_wrapper_list,k){
  rmse_list <- lapply(mv_wrapper_list, 
                      function(x) sqrt(mean((x$Observations-x$Predictions)^2, 
                                            na.rm = T)))
  idx <- which.min(rmse_list)
  return(list(rmse_list[[idx]], k[idx]))
}

# Main function
# function that for a given dataset, predictor list and target variable calculates
# the range of k values to be tested and then. It calls the mv_wrapper function for 
# each value of k. Is saves the return of the mv_wrapper function as entries in a 
# list and then returns the list. It calls the mv_wrapper function twice: once for 
# when the target is a predictor and once for when it is not.

model_fitter_loocv <- function(data, target, predictor_combinations, max_lag=3, E=3,
                               num_neighbors=0, k=0, lib=NULL, pred=NULL, Tp=1, ...){
  
  output_list <- unname(lapply(predictor_combinations, function(y){
    
    y <- unlist(y)
    t_excluded = !(target %in% y)
    kmax <- choose(length(y)*max_lag, E) 
    k_input <- k[k<=kmax]
    input_data <- data %>% ungroup() %>% dplyr::select(all_of(y),target)
    if(t_excluded==T) y <- append(y, target)
    columns <- paste(y, collapse = " ")
    mv_wrapper_list <- lapply(k_input, function(x)
      mv_wrapper_loocv(data = input_data, k = x, columns = columns, target = target, E = E,
                       max_lag = max_lag, num_neighbors = num_neighbors, Tp=Tp,
                       excludeTarget = t_excluded, lib=lib, pred=pred))
    model_selector_loocv(mv_wrapper_list,k_input)}))
  
  df <- mv_housekeeping(predictor_combinations, output_list, data, target, E)
  return(df)
}


# function that prepares the data that goes into the model_fitter function. It calls the
# model_fitter function for each bottle ID in parallel and then handles the return of the 
# model_fitter function by creating a list with two entries. The first entry is the merged
# dataframe across all bottles for the target species and includes the information about 
# predictors combinations and rmse values. The second entry # is a list with 18 entries
# (one for each bottle ID), and each of these 18 entries contains the best models for each 
# predictor combination tested.

model_fitter_wrapper_loocv <- function(whole.data, targets, num.clusters, max_lag=3, 
                                 E=3, num_neighbors=0, k=0, predictor_combinations, 
                                 lib=NULL, pred=NULL, predictors, Tp=1, 
                                 OptimalEsDf=NULL, ...){
  
  model_fitter_output_list <- mclapply(targets, function(target){
    if(!is.null(OptimalEsDf)){
      E <- OptimalEsDf[OptimalEsDf$Target==target,"E"]
    }
    write(paste0("Now running target: ", target, "\n"), write.file, append = T)
    old <- Sys.time()
    out <- model_fitter_loocv(data = whole.data, target = target, 
                              predictor_combinations = predictor_combinations, max_lag = max_lag,
                              E = E, num_neighbors = num_neighbors, k = k, 
                              lib=lib, pred=pred, Tp=Tp, ...)
    write(paste0("Time needed for target: ", Sys.time() - old, "\n"), write.file, append = T)
    out
  }, mc.cores = num.clusters, mc.cleanup = T, mc.allow.recursive = F)
  
  whole_df <- do.call("rbind", model_fitter_output_list)
  
  whole_df$predictor_combination <- as.character(whole_df$predictor_combination)
  temp <- matrix(F, nrow = nrow(whole_df), ncol = length(predictors))
  colnames(temp) <- predictors
  whole_df <- cbind(whole_df,temp)
  
  whole_df <- as.data.frame(t(apply(whole_df, 1, function(row) { 
    str <- unlist(strsplit(row["predictor_combination"], " "))
    row[str] <- T
    row
  })))
  
  whole_df$num_pred <- as.numeric(whole_df$num_pred)
  whole_df$RMSE <- as.numeric(whole_df$RMSE)
  whole_df$k <- as.numeric(whole_df$k)
  whole_df$E <- as.numeric(whole_df$E)
  for(i in predictors){  whole_df[,i] <- as.logical(whole_df[,i])}
  
  return(whole_df)
}


# ------ Reduced time points ------

forecasting <- function(data, steps_ahead=1, lib_step=14, E=2, OptimalEsDf=NULL, length.out=80,
                        predictors, pred_combs, targets, max_lag=3, write.file="temp.txt", removeLast=F){
  
  ks3 <- unique(round(exp(seq(log(0.8),log(floor(.75*nrow(data))),length.out = length.out))))
  
  if(!is.null(OptimalEsDf)) {
    ks3 <- max(ks3) # should prob change to be square root of max
  }
  
  libraries <- K_fold_libraries(data = data, step = lib_step)
  if(removeLast){
    lib <- libraries$lib[-length(libraries$pred)]
    pred <- libraries$pred[-length(libraries$pred)]
  } else {
    lib <- libraries$lib
    pred <- libraries$pred
  }
  
  forecasts <- lapply(1:length(pred), function(i){
    write(paste0("Fold number: ", i, " out of ", length(pred), "\n"), write.file, append = T)
    df <- model_fitter_wrapper(whole.data = data,
                               predictor_combinations = pred_combs,
                               targets = targets,
                               lib = lib[[i]], pred = pred[[i]],
                               num.clusters = detectCores() - 1,
                               E = E, max_lag = max_lag, k=ks3, 
                               Tp = steps_ahead,
                               predictors = predictors,
                               OptimalEsDf = OptimalEsDf,
                               write.file=write.file)
    df$forecast_subrun <- i
    df
  })
  
  forecasts <- do.call("rbind", forecasts)
  return(forecasts)
}

# ------ Reduced sampling freq ------

forecasting_loocv <- function(data, steps_ahead=1, lib_step=14, E=2, OptimalEsDf=NULL, length.out=80,
                              predictors, pred_combs, targets, max_lag=3, write.file="temp.txt", removeLast=F){
  
  ks3 <- unique(round(exp(seq(log(0.8),log(floor(.75*nrow(data))),length.out = length.out))))
  
  if(!is.null(OptimalEsDf)) {
    ks3 <- max(ks3) # should prob change to be square root of max
  }
  
  libraries <- loocv_libraries(data = data)
  lib <- libraries$lib
  pred <- libraries$pred
  
  forecasts <- model_fitter_wrapper_loocv(whole.data = data,
                                          predictor_combinations = pred_combs,
                                          targets = targets,
                                          lib = lib, pred = pred,
                                          num.clusters = detectCores() - 1,
                                          E = E, max_lag = max_lag, k=ks3, 
                                          Tp = steps_ahead,
                                          predictors = predictors,
                                          OptimalEsDf = OptimalEsDf,
                                          write.file=write.file)
  return(forecasts)
}