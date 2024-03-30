#' #' @import raster
#' #' @import sf
#' boundary_to_linestrings <- function(boundary_slice) {
#'   # Convert the boundary slice to a raster object
#'   boundary_raster <- raster(boundary_slice)
#'
#'   # Convert the raster object to SpatialLinesDataFrame
#'   boundary_lines <- rasterToContour(boundary_raster, levels = c(0.5))
#'
#'   # Convert the SpatialLinesDataFrame to an sf object
#'   boundary_sf <- sf::st_as_sf(boundary_lines)
#'
#'   return(boundary_sf)
#' }
#'
#' boundary_to_polygons <- function(boundary_slice) {
#'   # Convert the boundary slice to a raster object
#'   boundary_raster <- raster(boundary_slice)
#'
#'   # Convert the raster object to SpatialPolygonsDataFrame
#'   boundary_polygons <- rasterToPolygons(boundary_raster, fun = function(x) x > 0)
#'
#'   # Convert the SpatialPolygonsDataFrame to an sf object
#'   boundary_sf <- st_as_sf(boundary_polygons)
#'
#'   return(boundary_sf)
#' }
#'
#' cluster_to_polgygons <- function(slice) {
#'   # Convert the boundary slice to a raster object
#'   boundary_raster <- raster(slice)
#'
#'   # Convert the raster object to SpatialPolygonsDataFrame
#'   boundary_polygons <- rasterToPolygons(boundary_raster, fun = function(x) x > 0)
#'
#'   # Convert the SpatialPolygonsDataFrame to an sf object
#'   boundary_sf <- st_as_sf(boundary_polygons)
#'
#'   return(boundary_sf)
#'
#' }
#'
#' boundary_linestrings_list <- lapply(1:dim(boundary_volume)[3], function(slice_idx) {
#'   boundary_slice <- boundary_volume[, , slice_idx]
#'   linestrings <- boundary_to_linestrings(boundary_slice)
#'   return(linestrings)
#' })
#'
#' ras <- raster::raster(slice)
#' g <- st_as_sf(st_as_stars(ras), merge=TRUE, connect8=FALSE)
#' g2 = rmapshaper::ms_simplify(g, keep = .2, keep_shapes = FALSE, weighting=1.5)
#' sg2 = smoothr::smooth(g2, method="ksmooth", smoothness=.5,n=15)
#'
#' g <- st_as_sf(st_as_stars(ras), merge=TRUE, connect8=TRUE)
#'
#' boundary_sf <- boundary_to_linestrings(boundary_slice)
#'
#' # Plot the linestrings using ggplot2
#' ggplot() +
#'   geom_sf(data = sg2, aes(fill=factor(id), alpha=.2), colour="black",lwd = 1) +
#'   theme_minimal() +
#'   scale_fill_discrete()+
#'   theme(panel.grid.major = element_blank(),
#'         panel.grid.minor = element_blank(),
#'         axis.text = element_blank(),
#'         axis.ticks = element_blank(),
#'         axis.title = element_blank())
#'
#' smp = smoothr::smooth(g, method="ksmooth", smoothness=1)
#' ggplot() +
#'   geom_sf(data = smp, aes(color=factor(layer)),size = 1) +
#'   theme_minimal() +
#'   theme(panel.grid.major = element_blank(),
#'         panel.grid.minor = element_blank(),
#'         axis.text = element_blank(),
#'         axis.ticks = element_blank(),
#'         axis.title = element_blank())
