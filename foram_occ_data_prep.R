# Rasterize to the resolution of GCM data

rEmpt <- raster(ncols=288, nrows= 144, xmn=0, xmx=360, ymn=-90, ymx=90)
  #eck.proj <- "+proj=eck6 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"
  #rEmpt <- projectRaster(rEmpt, crs = eck.proj)
  #rEmpt <- setValues(rEmpt, c(1:ncell(rEmpt))) 
llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
palCoords <- occ[,c('pal.long','pal.lat')]
pts <- SpatialPointsDataFrame(palCoords, data = occ, proj4string = CRS(llPrj))
  #pts_eck6 <- spTransform(pts, eck.proj)

