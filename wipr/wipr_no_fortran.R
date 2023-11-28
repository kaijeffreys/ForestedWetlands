# Load packages
library(terra)
library(MultiscaleDTM)
library(randomForest)

# Attempts to do whole WIP tool without fortran
wipr <- function(DEM, len, buffer = TRUE,
                 metrics = c("grad", "plan", "prof", "dev", "twi"),
                 train) {
  
  # Sets up the resolution
  k <- round(len/res(DEM)[1])
  if (k %% 2 == 0) {
    k <- k + 1
  }
  
  # Initialize the inputs for the model
  print("Creating input metrics")
  in_rast <- list()
  
  if("grad" %in% metrics) {
    j <- k/2 - 0.5

    xl.end <- matrix(c(1, rep(NA_real_, times=k-1)), ncol=k, nrow=1)
    xr.end <- matrix(c(rep(NA_real_, times=k-1), 1), ncol=k, nrow=1)

    x.mids <- matrix(NA_real_, ncol=k, nrow=j-1)

    xl.mid <- matrix(c(2, rep(NA_real_, times=k-1)), ncol=k, nrow=1)
    xr.mid <- matrix(c(rep(NA_real_, times=k-1), 2), ncol=k, nrow=1)

    xl.mat <- rbind(xl.end, x.mids, xl.mid, x.mids, xl.end)
    xr.mat <- rbind(xr.end, x.mids, xr.mid, x.mids, xr.end)
    
    yt.end <- matrix(c(1, rep(NA_real_, times=k-1)), ncol=1, nrow=k)
    yb.end <- matrix(c(rep(NA_real_, times=k-1), 1), ncol=1, nrow=k)
    
    y.mids <- matrix(NA_real_, ncol=j-1, nrow=k)
    
    yt.mid <- matrix(c(2, rep(NA_real_, times=k-1)), ncol=1, nrow=k)
    yb.mid <- matrix(c(rep(NA_real_, times=k-1), 2), ncol=1, nrow=k)
    
    yt.mat <- cbind(yt.end, y.mids, yt.mid, y.mids, yt.end)
    yb.mat <- cbind(yb.end, y.mids, yb.mid, y.mids, yb.end)

    dz.dx.l <- terra::focal(r, xl.mat, fun=sum, na.rm=F)
    dz.dx.r <- terra::focal(r, xr.mat, fun=sum, na.rm=F)
    dz.dy.t <- terra::focal(r, yt.mat, fun=sum, na.rm=F)
    dz.dy.b <- terra::focal(r, yb.mat, fun=sum, na.rm=F)

    dz.dx <- (dz.dx.r-dz.dx.l)/(8*j*terra::res(DEM)[1])
    dz.dy <- (dz.dy.t-dz.dy.b)/(8*j*terra::res(DEM)[2])
    
    grad <- sqrt(dz.dx^2 + dz.dy^2)
    in_rast <- c(in_rast, grad = grad)
    
    # grad <- MultiscaleDTM::SlpAsp(DEM, w = k, metrics = "slope")
    # in_rast <- c(in_rast, grad = grad)
  }
  
  if("plan" %in% metrics) {
    if ("prof" %in% metrics) {
      both <- MultiscaleDTM::Qfit(DEM, metrics = c("planc", "profc"), w = k)
      in_rast <- c(in_rast, plan = both[[1]], prof = both[[2]])
    } else {
      plan <- MultiscaleDTM::Qfit(DEM, metrics = "planc", w = k)
      in_rast <- c(in_rast, plan = plan)
    }
  } else if("prof" %in% metrics) {
    prof <- MultiscaleDTM::Qfit(DEM, metrics = "profc", w = k)
    in_rast <- c(in_rast, prof = prof)
  }
  
  if("dev" %in% metrics) {
    dev <- (DEM - focal(DEM, w = k, fun = "mean")) / focal(DEM, w = k,
                                                           fun = "sd") 
    in_rast <- c(in_rast, dev = dev)
  }
  
  if("twi" %in% metrics) {
    topidx <- topmodel::topidx(terra::as.matrix(DEM), res = res(DEM)[1])
    a <- terra::setValues(r, topidx$area)
    twi <- a / tan(terra::terrain(r, unit = "radians"))
    values(twi) <- ifelse(values(twi) < 0, 0, values(twi))
    twi <- terra::focal(twi, w = k, mean)
    
    in_rast <- c(in_rast, twi = twi)
  }
  
  # Set up training data
  print("Building model")
  df_train <- data.frame(class = factor(train$Class))
  for(i in 1:length(in_rast)) {
    vals <- terra::extract(in_rast[[i]], train, ID = F)
    df_train <- cbind(df_train, vals)
  }
  df_train <- na.omit(df_train)
  colnames(df_train) <- c("class", names(in_rast))
  
  # Build model
  mod <- randomForest::randomForest(class ~ ., data = df_train, ntree = 200)
  
  # Run model
  input_raster <- in_rast[[1]]
  if(length(in_rast) > 1) {
    for(i in 2:length(in_rast)) {
      input_raster <- c(input_raster, in_rast[[i]])
    }
  }
  names(input_raster) <- names(in_rast)
  
  print("Running model")
  prob_rast <- terra::predict(input_raster, mod, na.rm = T, type = "prob")
  wet_prob <- prob_rast["WET"]
  return(wet_prob)
}
