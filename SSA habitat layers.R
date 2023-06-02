### Importing habitat layers

library(tidyverse)
library(terra)


ndvi_projected <- rast("mapping/cropped rasters/ndvi_GEE_projected_watermask20230207.tif")

# scaling from the buffalo data
scaling_mean <- as.numeric(attributes(buffalo_CLR_year_harmonics$ndvi_scaled)[1]) # 0.3125089
scaling_sd <- as.numeric(attributes(buffalo_CLR_year_harmonics$ndvi_scaled)[2]) # 0.143179 for July 2018
ndvi_stack_scaled <- (ndvi_projected - scaling_mean) / scaling_sd
# plot(ndvi_stack_scaled)
# time_vector <- terra::time(ndvi_projected)
# terra::time(ndvi_stack_scaled) <- NULL
# terra::time(ndvi_stack_scaled) <- lubridate::date(as.POSIXlt(as.POSIXlt("1970/01/01 00:00:00") + time_vector))
# ndvi_above0_stack_scaled <- rast("mapping/cropped rasters/ndvi_above0_stack_scaled_20230210.tif")

# elev <- rast("mapping/cropped rasters/DEM_H_raster.tif")
# slope <- rast("mapping/cropped rasters/slope_raster.tif")
# veg_herby <- rast("mapping/cropped rasters/veg_herby.tif")
# canopy_cover <- rast("mapping/cropped rasters/canopy_cover.tif")

# names(slope) <- "Slope"
# names(veg_herby) <- "Herbaceous vegetation"
# names(canopy_cover) <- "Canopy cover"

# plot(elev)
# plot(slope)
# plot(veg_herby)
# plot(canopy_cover)
# plot(ndvi_july2018)

# hist(buffalo_CLR_year$slope_scaled)
# slope_scaled <- (slope - attr(buffalo_CLR_year$slope_scaled, "scaled:center")) / attr(buffalo_CLR_year$slope_scaled, "scaled:scale")
# plot(slope_scaled)
# writeRaster(slope_scaled, "mapping/cropped rasters/slope_scaled_by_buffalo_data_20230210.tif")
# slope_scaled <- rast("mapping/cropped rasters/slope_scaled_by_buffalo_data_20230210.tif")

# hist(buffalo_CLR_year$canopy_scaled)
# canopy_scaled <- ((canopy_cover / 100) - attr(buffalo_CLR_year$canopy_scaled, "scaled:center")) / attr(buffalo_CLR_year$canopy_scaled, "scaled:scale")
# plot(canopy_scaled)
# writeRaster(canopy_scaled, "mapping/cropped rasters/canopy_scaled_by_buffalo_data_20230210.tif")
canopy_scaled <- rast("mapping/cropped rasters/canopy_scaled_by_buffalo_data_20230210.tif")

# hist(buffalo_CLR_year$herby_scaled)
# herby_scaled <- ((veg_herby) - attr(buffalo_CLR_year$herby_scaled, "scaled:center")) / attr(buffalo_CLR_year$herby_scaled, "scaled:scale")
# plot(herby_scaled)
# writeRaster(herby_scaled, "mapping/cropped rasters/herby_scaled_by_buffalo_data_20230210.tif")
herby_scaled <- rast("mapping/cropped rasters/herby_scaled_by_buffalo_data_20230210.tif")


# Subsetting the spatial extent to predict over

# finding the minimum and maximum extent (using the pseudo-absences as a buffer around the points)

# buffalo_CLR_year_harmonics %>% dplyr::summarise(xmin = min(x2), xmax = max(x2),
#                                                 ymin = min(y2), ymax = max(y2))
# 
# buffalo_CLR_year_harmonics %>% dplyr::summarise(xmin = min(x1), xmax = max(x1),
#                                                 ymin = min(y1), ymax = max(y1))

# subset_extent <- terra::ext(round(min(buffalo_CLR_year_harmonics$x2), digits = - 2), round(max(buffalo_CLR_year_harmonics$x2), digits = - 2),
#                             round(min(buffalo_CLR_year_harmonics$y2), digits = - 2), round(max(buffalo_CLR_year_harmonics$y2), digits = - 2))

subset_extent <- terra::ext(4200, 54300, -1458800, -1409700)

ndvi_stack_cropped <- terra::crop(ndvi_stack_scaled, subset_extent)
canopy_scaled_cropped <- terra::crop(canopy_scaled, subset_extent)
herby_scaled_cropped <- terra::crop(herby_scaled, subset_extent)

# plot(resources_cropped)
# plot(resources_cropped[[1]])
# points(as.numeric(buffalo_CLR_year %>% 
#                     filter(y == 1) %>% 
#                     dplyr::select(x1) %>% unlist()), 
#        as.numeric(buffalo_CLR_year %>% 
#                     filter(y == 1) %>% 
#                     dplyr::select(y1) %>% unlist()))

