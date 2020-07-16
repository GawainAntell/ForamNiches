# This function is adapted from the github folder https://github.com/macroecology/mapast 
# so that it doesn't only work on PBDB data. 
# (Though given I wrote it most of a year ago it's possible they've updated 
# the functionality of the original function)
# - IF

# A few lines are commented out because they are unnecessary. In addition, 
# 1. column names have been renamed to match foram data structure
# ($lat and $long were originally $latitude and $longitude)
# 2. the function originally looped through ages in the data,
# but in the context of foram data prep this function is used for
# data already subset to a given age, so the function is modified accordingly.
# - GSA

# possible models: "SETON2012" (default?), "MULLER2016", "GOLONKA", "PALEOMAP" or "MATTHEWS2016".

library(rjson)

pal.coord <- function(data, model, Ma) {
  ## this is from palaeocoords, which is a function in mapast
  # it assumes the age column is 'bin'
  paleolng <- c()
  paleolat <- c()
#  data$order <- seq(1:nrow(data))
#  data <- data[order(data$age), ]
  # subset the data for unique ages
  pts <- ""
  if (nrow(data) > 200) {
    # if there is lots of data in that age category, then run it in subsets
    num <- base::ceiling(nrow(data)/200)
    round <- 1
    while (round <= num) {
      pts <- ""
      if (round < num) {
        pts <- ""
        part2 <- data[((round - 1) * 200 + 1):(round * 200), ]
        for (j in 1:nrow(part2)) {
          pts <- base::paste0(pts, ",", part2$long[j], ",", part2$lat[j])
        }
        pts <- base::substring(pts, 2)
        url <- base::paste0("http://gws.gplates.org/reconstruct/reconstruct_points/?points=", 
                            pts, "&time=", Ma, "&model=", model, 
                            "&return_null_points")
        paleopts <- rjson::fromJSON(file = url)
      for (k in 1:base::length(paleopts$coordinates)) {
          if (base::is.null(paleopts$coordinates[[k]])) {
            paleolng <- c(paleolng, NA)
            paleolat <- c(paleolat, NA)
          } else {
            paleolng <- c(paleolng, paleopts$coordinates[[k]][1])
            paleolat <- c(paleolat, paleopts$coordinates[[k]][2])
          }
        }
      } else {
        pts <- ""
        part2 <- data[((round - 1) * 200 + 1):nrow(data), 
                      ]
        for (j in 1:nrow(part2)) {
          pts <- base::paste0(pts, ",", part2$long[j], 
                              ",", part2$lat[j])
        }
        pts <- base::substring(pts, 2)
        url <- base::paste0("http://gws.gplates.org/reconstruct/reconstruct_points/?points=", 
                            pts, "&time=", Ma, "&model=", model, 
                            "&return_null_points")
        paleopts <- rjson::fromJSON(file = url)
        for (k in 1:base::length(paleopts$coordinates)) {
          if (base::is.null(paleopts$coordinates[[k]])) {
            paleolng <- c(paleolng, NA)
            paleolat <- c(paleolat, NA)
          } else {
            paleolng <- c(paleolng, paleopts$coordinates[[k]][1])
            paleolat <- c(paleolat, paleopts$coordinates[[k]][2])
          }
        }
      }
      round <- round + 1
    }
  } else {
    # if there isn't >200 data points, then run the whole thing
    for (j in 1:nrow(data)) {
      pts <- base::paste0(pts, ",", data$long[j], 
                          ",", data$lat[j])
    }
    pts <- base::substring(pts, 2)
    url <- base::paste0("http://gws.gplates.org/reconstruct/reconstruct_points/?points=", 
                        pts, "&time=", Ma, "&model=", model, 
                        "&return_null_points")
    paleopts <- rjson::fromJSON(file = url)
    for (k in 1:base::length(paleopts$coordinates)) {
      if (base::is.null(paleopts$coordinates[[k]])) {
        paleolng <- c(paleolng, NA)
        paleolat <- c(paleolat, NA)
      } else {
        paleolng <- c(paleolng, paleopts$coordinates[[k]][1])
        paleolat <- c(paleolat, paleopts$coordinates[[k]][2])
      }
    }
  }
 
#  paleolat <- paleolat[order(data$order)]
#  paleolng <- paleolng[order(data$order)]
  return(list(paleolat = paleolat, paleolng = paleolng))
}
