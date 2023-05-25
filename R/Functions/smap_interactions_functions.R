# Main function (wrapper) ------------------------------------------------------

SMapInteractions <- function(targets, data, theta=NULL, num.clusters.target=1, 
                             num.clusters.CV, ccm_df, method, long_format=F, 
                             grid=NULL, RegressionType=NULL, Regression.Kernel=NULL){
  
  species_interactions <-mclapply(targets, function(Target){
    
    Embedding <- createEmbedding(ccm_df, Target)
    ll <- housekeeping_smap(data, Target, Embedding)
    block<-ll$block
    Edim<-ll$Edim
    coeff_names<-ll$coeff_names
    
    lib <- 1:dim(block)[1]
    pred <- 1:dim(block)[1]
    
    ## interactions estimation
    if(method=="SMap"){
      coeff <- ClassicSMap(Target, block, Edim, theta, coeff_names, lib, pred)
      parameters <- list()
    }
    if(method=="RegSMap"){
      fit <- RegularizedSMap(num.clusters.CV, block, Target, Embedding, grid, 
                             RegressionType, Regression.Kernel, Edim, coeff_names,
                             lib, pred)
      coeff <- fit$Coefficients
      parameters <- fit$Parameters
    }
    
    coeff$target <- Target
    
    if(long_format) {
      coeff <- coeffs_long(coeff, Target, Embedding)
      if(method=="RegSMap"){
        coeff$alpha <- parameters$ALPHA
        coeff$lambda <- parameters$LM
        coeff$theta <- parameters$TH
      }
    }
    list(coeff,parameters)
  }, mc.cores = num.clusters.target, mc.cleanup = T, mc.allow.recursive = F)
  
  coeff <- map(species_interactions, 1)
  parameters <- map(species_interactions, 2)
  
  if(long_format) return(do.call("rbind", coeff))
  return(list(Coefficients = coeff, Parameters = parameters))
}


# General functions ---------------------------------------------------------

### predict target based on coefficients
predict_smap <- function(data, coefs, method, interactors=NULL, target=NULL, rename_rEDM=T){
  if(method=="rEDM"){
    if(rename_rEDM){
      coefs <- rename_rEDM_SMap_coefs(coefs, interactors, target)[[1]]}
    names <- colnames(coefs)[-c(1,2)]
    names(names) <- gsub(".*/d","",names)
    prediction <- diag(as.matrix(data[,names(names)]) %*% t(coefs[-1,names])) + coefs[-1,"C0"]
  }
  
  if(method %in% c("SMap","RegSMap")){
    names <- colnames(coefs)
    names <- names[!(names %in% c("Time","C0","target"))]
    names(names) <- gsub(".*/d","",names)
    prediction <- diag(as.matrix(data[,names(names)]) %*% t(coefs[,names])) + coefs[,"C0"]
    prediction <- unname(unlist(prediction))
  }
  return(prediction)
}



coeffs_long <- function(coeff, Target, Embedding){
  Embedding <- c("C0",Target,Embedding)
  Edim <- length(Embedding)
  coeff_long <- coeff %>% pivot_longer(cols = c(C0, starts_with("d")),
                                       names_to = "Interaction",
                                       values_to = "Interaction_strength")
  coeff_long$interactor <- rep(Embedding,nrow(coeff_long)/Edim)
  return(coeff_long)
}



## function that puts data in right shape, create coef names, etc.

housekeeping_smap <- function(data, Target, Embedding){
  Embedding <- c(Target,Embedding[Embedding!=Target])
  Edim <- length(Embedding)
  data <- data[,Embedding]
  if(is.null(dim(data))) {
    data <- as.data.frame(data)
    colnames(data) <- Target
  }
  targ_col <- which(colnames(data) == Target)
  block <- cbind(data[2:dim(data)[1],targ_col],data[1:(dim(data)[1]-1),])
  block <- scale(block)
  idx <- apply(block,2,is.nan)
  block[idx] <- 0
  if(ncol(block)==2) {
    colnames(block)[2] <- Target
  }
  coeff_names <- sapply(colnames(data), function(x) paste("d", Target, "/d", x, sep = ""))
  return(list(block=block, Edim=Edim, coeff_names=coeff_names))
}




## compare coefficients of different SMap methods

coefficient_comparison <- function(coeff1, coeff2, target, rEDM=F, rename_rEDM=F, name1, name2){
  if(rEDM){
    if(rename_rEDM){
      coeff1 <- rename_rEDM_SMap_coefs(coeff1, interactors = cols, target = target)[[1]]}
    colnames(coeff1)[1] <- "Time"
  }
  
  cols <- colnames(coeff2)[!(colnames(coeff2) %in% c("Time","target"))]
  
  plot.list <- lapply(cols, function(coef){
    if(rEDM){
      tempdd <- data.frame(Method1=coeff1[,coef][-c(1,nrow(coeff1))],
                           Method2=unlist(coeff2[,coef]))
    } else {
      tempdd <- data.frame(Method1=unlist(coeff1[,coef]),
                           Method2=unlist(coeff2[,coef]))
    }
    
    plot <- tempdd %>%
      ggplot(aes(Method1, Method2))+
      geom_point() + 
      geom_abline(slope = 1, intercept = 0)+
      labs(subtitle = coef, x=name1, y=name2)
  })
  
  plot <- Reduce("+",plot.list) + plot_layout(ncol = 3)
  return(plot)
}


# create embedding

createEmbedding <- function(ccm_df, Target){
  emb <- ccm_df %>%
    ungroup() %>%
    dplyr::filter(target==Target, interactor!=Target) %>%
    dplyr::select(interactor) %>%
    unlist() %>%
    unname()
  return(emb)
}

logspace <- function(d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n)) 


# SMap from Deyle 2016 ---------------------------------------------------------

# SVD solver from Deyle et al. 2016

lm_svdsolve <- function(y, x, block, ws,  subset = seq_along(y)){
  if(is.null(dim(x))){
    x <- x[subset]
    dim_x = 1
  } else {
    x <- x[subset,]
    dim_x = dim(x)[2]
  }
  
  y <- y[subset]
  ws <- ws[subset]
  # prepended column of 1s for constant term in linear model
  A <- cbind(1, x) * ws
  A_svd <- svd(A)
  # >>REMOVE SMALL SINGULAR VALUES<<
  s <- A_svd$d
  s_inv <- matrix(0, nrow = dim_x+1, ncol = dim_x+1)
  for(i in seq_along(s)){
    if(s[i] >= max(s) * 1e-5)
      s_inv[i,i] <- 1/s[i]
  }
  coeff <- A_svd$v %*% s_inv %*% t(A_svd$u) %*% (ws * y)
  coeff <- t(coeff)
  colnames(coeff) <- c("C0",colnames(block)[2:dim(block)[2]])
  return(coeff)
}

ClassicSMap <- function(Target, block, Edim, theta, coeff_names,
                        lib, pred){
  
  coeff <- lapply(1:length(pred), function(ipred){
    #target point is excluded from the fitting procedure
    libs = lib[-pred[ipred]]
    # >>CALCULATE WEIGHTS<<
    q <- matrix(as.numeric(block[pred[ipred],2:dim(block)[2]]),
                ncol=Edim, nrow=length(libs), byrow = T)
    distances <- sqrt(rowSums((block[libs,2:dim(block)[2]] - q)^2))
    dbar <- mean(distances)
    Ws <- exp(-theta*distances/dbar)
    # >>REGRESS<<
    svd_fit <- lm_svdsolve(block[libs,1],block[libs,2:dim(block)[2]],block,Ws)
    svd_fit
  })
  coeff <- do.call("rbind", coeff)
  coeff <- cbind(pred,coeff)
  colnames(coeff) <- c("Time","C0",coeff_names)
  coeff <- as.data.frame(coeff)
  
  return(coeff)
}


DetermineNonlinearity <- function(Target, data, thetas, ccm_df, num.clusters=1){
  
  Embedding <- createEmbedding(ccm_df, Target)
  ll <- housekeeping_smap(data, Target, Embedding)
  block<-ll$block
  Edim<-ll$Edim
  coeff_names<-ll$coeff_names
  
  lib <- 1:dim(block)[1]
  pred <- 1:dim(block)[1]
  
  theta_rho <- mclapply(thetas, function(theta){
    # estimate coefficients for a certain theta value
    coeff <- ClassicSMap(Target = Target, block = block, Edim = Edim, theta = theta,
                         coeff_names = coeff_names, pred = pred, lib = lib)
    # predict target abundance based on estimate coefficients
    predicted <- predict_smap(data, coeff, "SMap")
    # calculate correlation rho between observed and  predicted target abundances
    rho <- cor(predicted, data[-1,Target], use = "complete.obs")
    # save theta and rho
    data.frame(theta=theta, rho=rho)
  }, mc.cores=num.clusters)
  
  theta_rho <- do.call("rbind", theta_rho)
  return(theta_rho)
}




# SMap from rEDM ---------------------------------------------------------

## rename smap coefs (to be checked whether correct...)

rename_rEDM_SMap_coefs <- function(coef_redm, interactors, target){
  new.names <- c()
  ordered.cols <- c()
  for(name in interactors){
    if(name==target){
      split <- strsplit(colnames(coef_redm), "/")
      coef.col <- which(grepl(name, map(split, 1)) & grepl(name, map(split, 2)))
    } else {
      coef.col <- which(grepl(name, colnames(coef_redm)))
    }
    new.name <- paste("d",target, "/d", name, sep = "")
    new.names <- c(new.names, new.name)  
    ordered.cols <- c(ordered.cols, coef.col)
  }
  ordered.cols <- ordered.cols - min(ordered.cols) + 1
  new.names <- c(c("Time","C0"), new.names[ordered.cols]) 
  colnames(coef_redm) <- new.names
  names(new.names) <- c("Time","C0",interactors[ordered.cols])
  return(list(coef_redm, new.names))
}





# Regularized SMap (Cenci et al.) ----------------------------------------------


## Kernel functions  (choice of weights) ---------------------------------------


Exponential.Kernel <- function(dst, theta){
  dbar <- mean(dst)
  krnl <- exp(-theta*dst/dbar)
  return(krnl)
}
Epanechnikov.Kernel <- function(dst, theta){
  tt <- dst/theta * ifelse(abs(dst/theta) > 1, 0, 1)
  krnl <- 0.75*(1 - tt^2)
  return(krnl)
}
TriCubic.Kernel <- function(dst,theta){
  tt <- dst/theta * ifelse(abs(dst/theta) > 1, 0, 1)
  krnl <- (1 - tt^3)^3
  return(krnl)
}
Matern.Kernel <- function(dst,theta){
  krnl = (1 + sqrt(3)*dst/theta)*exp(-(sqrt(3)*dst)/theta)
  return(krnl)
}



## Regularization functions -------------------------------------------------------------------------

ELNET_fit <- function(block, Target, Embedding, theta, lambda, alp, Regression.Kernel,
                      Edim, coeff_names, lib, pred){
  
  nrow <- length(pred)
  
  Krnl = match.fun(Regression.Kernel)
  
  coeff <- lapply(1:nrow, function(ipred){
    libs = lib[-pred[ipred]]
    q <- matrix(as.numeric(block[pred[ipred],2:dim(block)[2]]),
                ncol=Edim, nrow=length(libs), byrow = T)
    distances <- sqrt(rowSums((block[libs,2:dim(block)[2]] - q)^2))
    ### weights
    Ws = Krnl(distances, theta)
    ############ Fit function
    x = as.matrix(block[libs,2:dim(block)[2]])[-nrow, ]
    y = as.matrix(block[libs,1])[-nrow]
    Ws = Ws[-nrow]
    x = Ws * cbind(1, x)
    y = Ws * y
    fit <- enet(x, y, lambda = lambda, normalize = TRUE, intercept = FALSE)
    predict(fit, s = alp, type="coefficients", mode="fraction")$coefficients 
  })
  
  coeff <- do.call("rbind",coeff)
  coeff <- cbind(pred,coeff)
  colnames(coeff) <- c("Time","C0",coeff_names)
  return(as.data.frame(coeff))
}

ridge_fit <- function(block, Target, Embedding, theta, lambda, alp, Regression.Kernel,
                      Edim, coeff_names, lib, pred){
  
  nrow <- length(pred)
  
  Krnl = match.fun(Regression.Kernel)
  
  coeff <- lapply(1:nrow, function(ipred){
    #target point is excluded from the fitting procedure
    libs = lib[-pred[ipred]]
    # q is a (N-1)xE matrix with the ipred-th entry of the rescaled time series (yes it is a N-1 repetation of one value)
    q <- matrix(as.numeric(block[pred[ipred],2:dim(block)[2]]),
                ncol=Edim, nrow=length(libs), byrow = T)
    ######### Here compute the weigths. Wx is going to be a (N-1) vector of weights
    distances <- sqrt(rowSums((block[libs,2:dim(block)[2]] - q)^2))
    Ws = Krnl(distances, theta)
    ####### svd_fit gives me the Jacobian element at each time step ipred
    lm_regularized(block[libs,1], block[libs,2:dim(block)[2]], Ws, lambda, Edim)
  })
  
  coeff <- do.call("rbind",coeff)
  coeff <- cbind(pred,coeff)
  colnames(coeff) <- c("Time","C0",coeff_names)
  return(as.data.frame(coeff))
}
lm_regularized <- function(y, x, ws, lambda, dimension, subset = seq_along(y)){
  if(is.null(dim(x))){
    x <- x[subset]
    dim_x = 1
  } else {
    x <- x[subset,]
    dim_x = dim(x)[2]
  }
  #### y is the target colum (N-1) vector.
  #### x are all the others (N-1) x E vector
  #### Ws is a (N-1) vector
  # x <- x[subset,]
  y <- y[subset]
  ws <- ws[subset]
  WWs = diag(ws)
  Xx = as.matrix(x)
  #### For constant term in linear fit
  Xx = cbind(1, Xx)
  coeff <- solve(t(Xx) %*% WWs %*% Xx + lambda*nrow(Xx)*diag(1,dimension + 1)) %*% t(Xx) %*%(ws * y)
  coeff <- t(coeff)
  
  return(coeff)
}


## Main functions -------------------------------------------------------------------------

## this functions finds best theta and lambda
LOOCV_smap_parameters <- function(num.clusters, block, Target, Embedding, 
                                  grid, RegressionType, Regression.Kernel,
                                  Edim, coeff_names, lib, pred){
  
  S_target <- mclapply(1:nrow(grid), function(row){
    CV(row, block, Target, Embedding, grid, RegressionType, Regression.Kernel,
       Edim, coeff_names, lib, pred)
  }, mc.cores = num.clusters)
  
  error.mat <- do.call("rbind",S_target)
  best <- error.mat %>%
    dplyr::mutate(mse = round(mse,4)) %>%
    dplyr::filter(mse == min(mse)) %>%
    dplyr::filter(bndwth == min(bndwth)) %>%
    dplyr::filter(lmb == min(lmb)) %>%
    dplyr::filter(alpha == min(alpha))
  return(list(BestTH = best$bndwth,
              BestLM = best$lmb,
              BestALPHA = best$alpha,
              min.var.err = best$mse,
              val.err = error.mat))
}

## cross validation: run a parameter combination
CV <- function(row, block, Target, Embedding, grid,
               RegressionType, Regression.Kernel,
               Edim, coeff_names, lib, pred){
  FUN <- match.fun(RegressionType)
  #### Fit the model
  coefficients <- FUN(block, Target, Embedding, grid[row,1], grid[row,2], grid[row,3], 
                      Regression.Kernel, Edim, coeff_names, lib, pred)
  #### Take the forecast
  Data <- predict_smap(block, coefficients, "RegSMap")
  #### Compute the mean square error
  MSE = mean((block[-1,Target]-Data)^2)
  #return(MSE)
  return(data.frame(mse = MSE, bndwth = grid[row,1], lmb = grid[row,2], alpha = grid[row,3]))
}

## main function: regularized Smap
RegularizedSMap <- function(num.clusters, block, Target, Embedding, grid,
                            RegressionType, Regression.Kernel,
                            Edim, coeff_names, lib, pred){
  FUN = match.fun(RegressionType)
  RegularizedParameters <- LOOCV_smap_parameters(num.clusters, block, Target, Embedding, 
                                                 grid, RegressionType, Regression.Kernel,
                                                 Edim, coeff_names, lib, pred)
  ## Now compute the optimum regularized coefficients
  J = FUN(block, Target, Embedding, RegularizedParameters$BestTH, RegularizedParameters$BestLM,
          RegularizedParameters$BestALPHA, Regression.Kernel,
          Edim, coeff_names, lib, pred)
  
  return(list(Coefficients = J, Parameters = list(TH = RegularizedParameters$BestTH,
                                                  LM = RegularizedParameters$BestLM,
                                                  ALPHA = RegularizedParameters$BestALPHA)))
}






# Multiview Regularized SMap ---------------------------------------------------

# to be done 




# Wrapper functions for the monthly lake data  ---------------------------------

LakeInteractions <- function(data, targets, num.clusters.target, num.clusters.CV,
                             long_format = T, method, ccm_df,
                             Regression.Kernel, thetas, lambdas,
                             ELNET0.9=T, ELNET0.5=T, ELNET0.1=T, Ridge=T){
  
  returnlist <- list()
  
  ### ELNET with alpha=0.9
  if(ELNET0.9){
    grid <- expand.grid(tht=thetas, lambda=lambdas, alpha=0.9)
    
    dd_Regsmap_ELNET0.9 <- SMapInteractions(targets = targets, data = data, 
                                            num.clusters.target = num.clusters.target, 
                                            num.clusters.CV = num.clusters.CV, 
                                            ccm_df = ccm_df, method = method,
                                            long_format = long_format, grid = grid, 
                                            RegressionType = "ELNET_fit", 
                                            Regression.Kernel = Regression.Kernel)
    
    
    dd_Regsmap_ELNET0.9 <- dd_Regsmap_ELNET0.9 %>%
      dplyr::filter(Interaction!="C0")
    
    dd_Regsmap_ELNET0.9$method <- "ELNET0.9"
    list_ELNET0.9 <- list(ELNET0.9 = dd_Regsmap_ELNET0.9)
    returnlist <- append(returnlist, list_ELNET0.9)
  }
  
  
  ### ELNET with alpha=0.5
  if(ELNET0.5){
    grid <- expand.grid(tht=thetas, lambda=lambdas, alpha=0.5)
    
    dd_Regsmap_ELNET0.5 <- SMapInteractions(targets = targets, data = data, 
                                            num.clusters.target = num.clusters.target, 
                                            num.clusters.CV = num.clusters.CV, 
                                            ccm_df = ccm_df, method = method,
                                            long_format = long_format, grid = grid, 
                                            RegressionType = "ELNET_fit", 
                                            Regression.Kernel = Regression.Kernel)
    
    
    dd_Regsmap_ELNET0.5 <- dd_Regsmap_ELNET0.5 %>%
      dplyr::filter(Interaction!="C0")
    
    dd_Regsmap_ELNET0.5$method <- "ELNET0.5"
    list_ELNET0.5 <- list(ELNET0.5 = dd_Regsmap_ELNET0.5)
    returnlist <- append(returnlist, list_ELNET0.5)
  }
  
  
  ### ELNET with alpha=0.1
  if(ELNET0.1){
    grid <- expand.grid(tht=thetas, lambda=lambdas, alpha=0.1)
    
    dd_Regsmap_ELNET0.1 <- SMapInteractions(targets = targets, data = data, 
                                            num.clusters.target = num.clusters.target, 
                                            num.clusters.CV = num.clusters.CV, 
                                            ccm_df = ccm_df, method = method,
                                            long_format = long_format, grid = grid, 
                                            RegressionType = "ELNET_fit", 
                                            Regression.Kernel = Regression.Kernel)
    
    
    dd_Regsmap_ELNET0.1 <- dd_Regsmap_ELNET0.1 %>%
      dplyr::filter(Interaction!="C0")
    
    dd_Regsmap_ELNET0.1$method <- "ELNET0.1"
    list_ELNET0.1 <- list(ELNET0.1 = dd_Regsmap_ELNET0.1)
    returnlist <- append(returnlist, list_ELNET0.1)
  }
  
  
  ### Ridge regression (no alpha)
  if(Ridge){
    grid <- expand.grid(tht=thetas, lambda=lambdas, alpha=1)
    
    dd_Regsmap_ridge <- SMapInteractions(targets = targets, data = data, 
                                         num.clusters.target = num.clusters.target, 
                                         num.clusters.CV = num.clusters.CV, 
                                         ccm_df = ccm_df, method = method,
                                         long_format = long_format, grid = grid, 
                                         RegressionType = "ridge_fit", 
                                         Regression.Kernel = Regression.Kernel)
    
    dd_Regsmap_ridge <- dd_Regsmap_ridge %>%
      dplyr::filter(Interaction!="C0")
    
    dd_Regsmap_ridge$method <- "ridge"
    list_ridge <- list(ridge = dd_Regsmap_ridge)
    returnlist <- append(returnlist, list_ridge)
  }
  
  return(returnlist)
}



