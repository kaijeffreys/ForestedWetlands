# Load in packages
library(terra)
library(randomForest)

# Surface metrics function using makeGrids
surface_met1 <- function(len, metrics = c("grad", "plan", "prof", "dev"),
                         dem_dir, exec_dir, out_dir=getwd(), re_sample = NA) {
  
  # Checking to see if directories exist
  if(!file.exists(dem_dir)) {
    stop("DEM directory does not exist!")
  }
  
  if(!dir.exists(exec_dir)) {
    stop("Executable Files directory does not exist!")
  }
  
  # if(!endsWith(dem_dir, ".flt")) {
  #   r <- terra::rast(dem_dir)
  #   dem_dir <- paste0(out_dir, "/temp.flt")
  #   writeRaster(r, dem_dir, overwrite = T)
  #   print("Created .flt file")
  # }
  
  # Prepare inputs
  dem_dir <- normalizePath(dem_dir)
  # dem_dir <- substr(dem_dir, 1, nchar(dem_dir)-4)
  out_dir <- normalizePath(out_dir)
  if(!endsWith(out_dir, "\\")) {
    out_dir <- paste0(out_dir, "\\")
  }
  exec_dir <- normalizePath(exec_dir)
  
  # Write input file
  file_name <- paste0(out_dir, "input_makeGrids.txt")
  file.create(file_name)
  
  writeLines(c("# Input file for makeGrids",
               "",
               paste0("DEM: ", dem_dir),
               paste0("SCRATCH DIRECTORY: ", out_dir),
               paste0("LENGTH SCALE: ", len)), con = file_name)
  
  if("grad" %in% metrics) {
    write(paste0("GRID: GRADIENT, OUTPUT FILE = ", out_dir, "grad", len, ".flt"),
          file = file_name, append = T) 
  }
  
  if("plan" %in% metrics) {
    write(paste0("GRID: PLAN CURVATURE, OUTPUT FILE = ", out_dir,
                 "plan", len), file = file_name, append = T)
  }
  
  if("prof" %in% metrics) {
    write(paste0("GRID: PROFILE CURVATURE, OUTPUT FILE = ", out_dir,
                 "prof", len), file = file_name, append = T)
  }
  
  # Run surface metrics sans DEV
  system(paste0(exec_dir, "\\makeGrids"), input = file_name)
  
  # Writing input file for DEV
  if ("dev" %in% metrics) {
    if(is.na(re_sample)) {
      stop("Set re_sample level")
    }
    
    # Prepare inputs
    file_name <- paste0(out_dir, "input_localRelief.txt")
    rad <- len / 2
    
    # Create and write input file
    file.create(file_name)
    writeLines(c("# Input file for LocalRelief",
                 "# Creating by surfaceMetrics.R",
                 paste0("# On ", Sys.time()),
                 paste0("DEM: ", dem_dir),
                 paste0("SCRATCH DIRECTORY: ", out_dir),
                 paste0("RADIUS: ", rad),
                 paste0("DOWN SAMPLE: ", re_sample),
                 paste0("SAMPLE INTERVAL: ", re_sample),
                 paste0("OUTPUT LOCAL RASTER: ", out_dir, "local", len)),
               con = file_name)
    
    # Run DEV in console
    system(paste0(exec_dir, "\\localRelief"), input = file_name)
  }
}


# Build training data
# train_data <- function(region_file, wet_file, wet_rate, not_wet_rate) {
#   
#   # Sampling points helper function
#   samp_pts <- function(polys, rate) {
#     coords <- NA
#     for(i in 1:length(polys)) {
#       pol_area <- expanse(polys[i], unit = "km")
#       
#       if(pol_area != 0) {
#         samp_size <- ceiling(pol_area * rate)
#         samp_coords <- crs(spatSample(poly, samp=samp_size))
#         coords <- rbind(coords, pol_area)
#       }
#     }
#   }
#   
#   # Load in polygons
#   reg_poly <- vect(region_file)
#   wet_poly <- vect(wet_file)
#   
#   wet_poly <- project(wet_poly, reg_poly)
#   wet_poly <- crop(wet_poly, reg_poly)
#   dry_poly <- erase(reg_poly, wet_poly)
#   
#   # Get the coordinates
#   wet_crds <- samp_pts(wet_poly, wet_rate)
#   dry_crds <- samp_pts(dry_poly, not_wet_rate)
#   
#   crds <- rbind(wet_crds, dry_crds)
#   atts <- data.frame(class = factor(rep(c(1, 0), c(length(wet_crds),
#                                                    length(non_wet_crds)))))
#   train_pts <- vect(crds, atts = atts, crs = crs(reg_poly))
#   
#   return(train_pts)
# }

# Build forest model
build_forest <- function(ref_raster_file, input_raster_files,
                         input_poly_files, train_data) {
  # Load in reference raster
  ref_raster <- rast(ref_raster_file)
  
  # Check to see if raster files exist and load them in 
  raster_list <- list()
  for(i in 1:length(input_raster_files)) {
    file <- input_raster_files[i]
    if(!file.exists(file)) {
      stop(paste0("Could not find raster file:", file))
    }
    raster_list <- c(raster_list, rast(file))
  }
  
  # raster_list <- terra::align(ref_raster, raster_list)
  
  # Check to see if poly files exist, load them in, and then 
  # poly_list <- list()
  # lapply(input_poly_files, function(file) {
  #   if(!file.exists(file)) {
  #     stop(paste0("Could not find polygon file:", file))
  #   }
  #   poly_list[[length(poly_list) + 1]] <- vect(file)
  # })
  # 
  # poly_raster <- list()
  # 
  # for(pol in poly_list) {
  #   pol <- project(pol, ref_raster)
  #   
  #   for(i in 1:length(names(pol))) {
  #     poly_raster[[i]] <- rasterize(pol, ref_raster, field = names(pol)[i])
  #   }
  # }
  # 
  # raster_list <- c(raster_list, poly_raster)
  
  train_df <- data.frame(class = factor(train_data$Class))
  for(i in 1:length(raster_list)) {
    vals <- terra::extract(raster_list[[i]][[1]], train_data, ID = F)
    train_df <- cbind(train_df, vals)
  }
  train_df <- na.omit(train_df)
  
  # Build model
  model <- randomForest::randomForest(class ~ ., data = train_df, ntree = 200)
  return(model)
}
