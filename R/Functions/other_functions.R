lm_autoplot <- function(model, leg="collect", group=NULL){
  data <- model@frame %>%
    ungroup() %>%
    dplyr::mutate(resid=resid(model),
                  fitted = fitted(model),
                  sqrt_abs_resid = sqrt(abs(resid(model))),
                  std_resid = rstudent(model),
                  leverage = hatvalues(model),
                  qqx = qqnorm(rstudent(model), plot.it = F)$x,
                  qqy = qqnorm(rstudent(model), plot.it = F)$y) 
  
  y <- quantile(data$std_resid, c(0.25,0.75))
  x <- qnorm(c(0.25,0.75))
  
  diag1 <- data %>%
    ggplot(aes(fitted,resid,col=group))+
    geom_point(pch=1) +
    labs(title = "Residuals vs Fitted", x="Fitted values", y="Residuals")
  
  diag3 <- data %>%
    ggplot(aes(fitted,sqrt_abs_resid,col=group))+
    geom_point(pch=1) +
    labs(title = "Scale-Location", x="Fitted values", y=expression(sqrt("Standardized Residuals")))
  
  diag2 <- data %>%
    ggplot(aes(x=qqx,y=qqy,col=group))+
    geom_point(pch=1) +
    geom_abline(slope = diff(y)/diff(x), intercept = y[1L] - diff(y)/diff(x) * x[1L])+
    labs(title = "Normal Q-Q", x="Theoretical Quantiles", y="Studentized Residuals")
  
  diag4 <- data %>%
    ggplot(aes(leverage,std_resid,col=group))+
    geom_point(pch=1) +
    labs(title = "Residuals vs Leverage", x="Leverage", y="Studentized Residuals")
  
  if(leg=="keep"){
    plot <- (diag1 + diag2)/(diag3 + diag4)+ plot_layout(guides = leg, heights = c(4,4)) &
      theme_bw() & theme(legend.position = "bottom") 
  } else {
    plot <- (diag1 + diag2)/(diag3 + diag4)/guide_area() + plot_layout(guides = leg, heights = c(4,4,2)) &
      theme_bw() &  theme(legend.position = "bottom") & guides(col=guide_legend(nrow=3,byrow=TRUE))    
  }
  return(plot)
}


# Standard theme function
standard_theme <- function(legend="none"){ 
  theme_bw() %+replace%   
    theme(
      legend.position = legend,
      axis.title = element_text(size=13),
      axis.text = element_text(size=10),
      legend.text = element_text(size=11),
      strip.text = element_text(size=11, hjust = 0, vjust = 1),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    )
}

standard_theme2 <- function(legend="none"){ 
  theme_bw() %+replace%   
    theme(
      legend.position = legend,
      axis.title = element_text(size=13),
      axis.text = element_text(size=10),
      legend.text = element_text(size=11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    )
}



data_transformation <- function(data, detrending="temporal", time_column, interpolation=TRUE){
  ## Interpolation
  if(interpolation==T){
    data$time <- unlist(data[,time_column])
    time.interp <- time.interp.fun(timestep = data$time)
    data <- data %>% na.omit()
    data <- setDT(data)[, list(time = time.interp,
                               value = chi(time, value, time.interp)), by = .(name)]
  }

  ## fourth root transformation
  data$value <- (data$value)^(1/4)
  ## remove linear trend
  if(detrending=="temporal"){
    data <- setDT(data)[, list(time = time,
                               value = residuals(lm(value~time))), by = .(name)]
  } 
  
  if(detrending=="standardization"){
    data <- setDT(data)[, list(time = time,
                               value = (value-mean(value))/sd(value)), by = .(name)]
  }

  data <- data %>% na.omit()
  ## Standardization
  data <- data %>%
    group_by(name) %>%
    mutate(value = (value-mean(value, na.rm=T))/sd(value, na.rm = T))
  return(data)
}


figure_RMSEvsInteractions <- function(data, group, col, se=FALSE, median.numint=F) {

  data$group <- unlist(data[,group])
  data$col <- unlist(data[,col])
  num.levels <- length(unique(data$col))
  
  if(median.numint){
    a <- data %>%
      ggplot(aes(medianNumberOfInteractions, medianRMSE, group=group, col=col, fill=col))
  } else {
    a <- data %>%
      ggplot(aes(NumberOfInteractions, medianRMSE, group=group, col=col, fill=col))
  }
  
  b <- data %>%
    ggplot(aes(mean_mean_abs, medianRMSE, group=group, col=col, fill=col)) +
    labs(x="Mean interaction strength")

  c <- data %>%
    ggplot(aes(sum_mean_abs, medianRMSE, group=group, col=col, fill=col)) +
    labs(x="Sum of interaction strengths")

  d <- data %>%
    ggplot(aes(NumberOfInteractions,mean_mean_abs, group=group, col=col, fill=col)) +
    labs(y="Mean interaction strength")

  plot <- a + b + c + d & 
    standard_theme2(legend = "right") & 
    geom_smooth(method = "lm", alpha=0.2, se = se) &
    geom_point(size=2, col="black", shape=21)
  plot <- plot + plot_layout(guides = "collect")
  if(num.levels<=8) {
    plot <- plot & scale_color_brewer(palette = "Dark2", aesthetics = c("col","fill"))
  }

  return(plot)
}
