# Load packages
library(terra)
library(MultiscaleDTM)
library(randomForest)

# Attempts to do whole WIP tool without fortran
wipr <- function(DEM, len, poly_inputs = list(),
                 elev_dev = c("grad", "plan", "prof", "dev", "twi"),
                 train, export_prob = FALSE, prob_rast_name = NA) {
  
  # Checks if inputs are file names and loads them in
  if(is.character(DEM)) {
    if(!file.exists(DEM)) {
      stop("Cannot find DEM file")
    }
    DEM <- terra::rast(DEM)
  }
  
  if(is.character(train)) {
    if(!file.exists(train)) {
      stop("Cannot find training points file")
    }
    train <- terra::vect(train)
  }
  
  if(length(poly_inputs) > 0) {
    for(i in 1:length(poly_inputs)) {
      if(is.character(poly_inputs[[i]])) {
        if(!file.exists(poly_inputs[[i]])) {
          stop(paste0("Cannot find poly input file:", poly_inputs[i]))
        }
        poly_inputs[[i]] <- terra::vect(poly_inputs[[i]]) 
      }
    }
  }
  
  # Sets up the resolution
  k <- round(len/res(DEM)[1])
  if (k %% 2 == 0) {
    k <- k + 1
  }
  
  # Initialize the inputs for the model
  print("Creating input metrics")
  in_rast <- list()
  
  if("grad" %in% elev_dev) {
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

    dz.dx.l <- terra::focal(DEM, xl.mat, fun=sum, na.rm=T, na.policy = "omit")
    dz.dx.r <- terra::focal(DEM, xr.mat, fun=sum, na.rm=T, na.policy = "omit")
    dz.dy.t <- terra::focal(DEM, yt.mat, fun=sum, na.rm=T, na.policy = "omit")
    dz.dy.b <- terra::focal(DEM, yb.mat, fun=sum, na.rm=T, na.policy = "omit")

    wts.l <- terra::focal(!is.na(DEM), w=xl.mat, fun=sum, na.rm=TRUE,
                            na.policy = "omit")
    wts.r <- terra::focal(!is.na(DEM), w=xr.mat, fun=sum, na.rm=TRUE,
                            na.policy = "omit")
    wts.t <- terra::focal(!is.na(DEM), w=yt.mat, fun=sum, na.rm=TRUE,
                            na.policy = "omit")
    wts.b <- terra::focal(!is.na(DEM), w=yb.mat, fun=sum, na.rm=TRUE,
                            na.policy = "omit")
    dz.dx <- ((dz.dx.r/wts.r) - (dz.dx.l/wts.l))/(2*j*terra::xres(DEM))
    dz.dy <- ((dz.dy.t/wts.t) - (dz.dy.b/wts.b))/(2*j*terra::yres(DEM))
    
    grad <- sqrt(dz.dx^2 + dz.dy^2)
    in_rast <- c(in_rast, grad = grad)
  }
  
  if("plan" %in% elev_dev) {
    if ("prof" %in% elev_dev) {
      both <- MultiscaleDTM::Qfit(DEM, metrics = c("planc", "profc"),
                                  w = k, na.rm = T)
      in_rast <- c(in_rast, plan = both[[1]], prof = both[[2]])
    } else {
      plan <- MultiscaleDTM::Qfit(DEM, metrics = "planc", w = k, na.rm = T)
      in_rast <- c(in_rast, plan = plan)
    }
  } else if("prof" %in% elev_dev) {
    prof <- MultiscaleDTM::Qfit(DEM, metrics = "profc", w = k, na.rm = T)
    in_rast <- c(in_rast, prof = prof)
  }
  
  if("dev" %in% elev_dev) {
    dev <- (DEM - focal(DEM, w = k, fun = "mean", na.rm = T, na.policy = "omit")) / focal(DEM, w = k, fun = "sd", na.rm = T, na.policy = "omit") 
    in_rast <- c(in_rast, dev = dev)
  }
  
  if("twi" %in% elev_dev) {
    topidx <- topmodel::topidx(terra::as.matrix(DEM), res = res(DEM)[1])
    a <- terra::setValues(DEM, topidx$area)
    twi <- a / tan(terra::terrain(DEM, unit = "radians"))
    values(twi) <- ifelse(values(twi) < 0, 0, values(twi))
    twi <- terra::focal(twi, w = k, mean, na.rm = T, na.policy = "omit")
    
    in_rast <- c(in_rast, twi = twi)
  }
  
  if(length(poly_inputs) > 0) {
    for(i in 1:length(poly_inputs)) {
      vr_name <- names(poly_inputs)[i]
      temp_rast <- terra::rasterize(vr_name, DEM, field = vr_name)
      in_rast <- c(in_rast, temp_rast)
    }
  }
  
  # Set up training data
  df_train <- data.frame(class = factor(train$Class))
  for(i in 1:length(in_rast)) {
    vals <- terra::extract(in_rast[[i]], train, ID = F)
    df_train <- cbind(df_train, vals)
  }
  df_train <- na.omit(df_train)
  colnames(df_train) <- c("class", names(in_rast))
  
  # Build the model
  print("Building and running model")
  
  mod <- randomForest::randomForest(class ~ ., data = df_train, ntree = 200)
  
  # Move the rasters into one multilayered raster
  input_raster <- in_rast[[1]]
  if(length(in_rast) > 1) {
    for(i in 2:length(in_rast)) {
      input_raster <- c(input_raster, in_rast[[i]])
    }
  }
  names(input_raster) <- names(in_rast)
  
  # Run the model
  prob_rast <- terra::predict(input_raster, mod, na.rm = T, type = "prob")
  wet_prob <- prob_rast["WET"]
  
  # Return the wetlands probability raster as the output
  if(export_prob) {
    if(!is.na(prob_rast_name)) {
      terra::writeRaster(wet_prob, filename = paste0(prob_rast_name), ".tif")
    } else {
      terra::writeRaster(wet_prob, filename = "wetProb.tif")
    }
  }
  
  return(wet_prob)
}
