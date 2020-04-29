# Import a raster brick for a specified time bin and environmental variable.
# Temperature raster files have 19 layers. If using mixed layer depth or BVF,
# then modify code for 1 layer.
getBrik <- function(bin, envNm, mods){
  modRow <- mods$age_1000ka == bin
  id <- mods$id[modRow]
  
  # Load the rasters for only the desired env variables and time step
  allFls <- list.files('Data/', recursive = TRUE)
  txt <- paste0(id,'.*tif')
  modFls <- grep(txt, allFls)
  flNms <- paste0('Data/', allFls[modFls])
  envFlPos <- sapply(envNm, grep, flNms)
  envFlNms <- flNms[envFlPos]
  
  r <- lapply(envFlNms, brick)
  names(r) <- envNm
  r
}
