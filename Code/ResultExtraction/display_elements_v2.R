
library(tidyverse)
library(ggplot2)
library(ggmap)
library(tmap)
library(tmaptools)
library(viridis)
library(sp)
library(sf)
library(stars)

#######################
# GET DATA
#######################
CH_canton <- st_read("C:/Users/folly/Desktop/Homeoffice/Geo_Political/gg2019/ggg_2019-LV03/shp/k4l19.shp")
CH_lakes <- st_read("C:/Users/folly/Desktop/Homeoffice/Geo_Political/gg2019/ggg_2019-LV03/shp/k4s19.shp")
# diff models
grid_lm <- readRDS("../_15_noflight/Model_LM_inla/Output/grid_results.rds")
grid_me <- readRDS("../_15_noflight/Model_standard/Output/grid_results.rds")
grid_ex <- readRDS("../_15_noflight/Model_extended_priors/Main/Output/grid_results.rds")
grid_sp <- readRDS("../_15_noflight/Model_SPDE/Output/grid_results.rds")
# combine data from diff models
grid_ex$mean_lm <- grid_lm$mean
grid_ex$sd_lm <- grid_lm$sd
grid_ex$mean_sp <- grid_sp$mean
grid_ex$sd_sp <- grid_sp$sd
grid_ex$mean_me <- grid_me$mean
grid_ex$sd_me <- grid_me$sd
# prepare to go spatial
grid <- grid_ex
grid$ID <- 1:nrow(grid)
grid$X <- 1000 * grid$X
grid$Y <- 1000 * grid$Y
# go spatial
coordinates(grid) <- ~X+Y
proj4string(grid) <- proj4string(as(CH_canton,"Spatial"))
gridded(grid) <- TRUE
grid <- as(grid, "SpatialPolygonsDataFrame")
grid <- st_as_sf(grid)
# tidy data
#grid <- gather(grid, "key","value", -c("geometry"))
#grid$key <- as.factor(grid$key)
#grid$type <- ifelse(grid$key %in% c("mean","mean_lm","mean_me","mean_sp"), "mean", "ID")
#grid$type[grid$key %in% c("sd","sd_lm","sd_me","sd_sp")] <- "sd"
#grid$type[grid$key %in% c("C_ID","GENESE","GESTEINKL","LC09_6","LEG_TEK_1")] <- "covar"

# plot multiples???? takes forever
#map_mean <- ggplot() + geom_sf(data = grid[grid$type %in% c("mean"),],  aes(fill = value),col=NA) + scale_fill_viridis_c() + theme_bw() + labs(fill = " ") + facet_wrap(.~key)
#map_sd <- ggplot() + geom_sf(data = grid[grid$type %in% c("sd"),],  aes(fill = value),col=NA) + scale_fill_viridis_c() + theme_bw() + labs(fill = " ") + facet_wrap(.~key)

#What scale is needed for colors? min = 2.497, max = 5.765
max(grid$mean) # 5.552
max(grid$mean_lm) # 5.765
max(grid$mean_me) # 5.593
max(grid$mean_sp) # 5.497
min(grid$mean) # 2.604
min(grid$mean_lm) # 3.422
min(grid$mean_me) # 2.497
min(grid$mean_sp) # 2.662

#what scale needed for sd? min = 0, max = 0.48
max(grid$sd) # 0.396
max(grid$sd_lm) # 0.136
max(grid$sd_me) # 0.42
max(grid$sd_sp) # 0.478
min(grid$sd) # 0.021
min(grid$sd_lm) # 0.003
min(grid$sd_me) # 0.005
min(grid$sd_sp) # 0.004



# plot single
ggplot() + geom_sf(data = grid,  aes(fill = mean_lm),col=NA) + scale_fill_viridis_c(limits = c(2.497,5.765)) + theme_bw() + labs(fill = " ")

# make and store plots
# mean
p1 <- ggplot() + geom_sf(data = grid,  aes(fill = mean), col = NA) + scale_fill_viridis_c(limits = c(2.497,5.765)) + theme_bw() + labs(fill = " ")
png("extended_mean.png", res = 300, width = 8, height = 6, units = "cm")
p1
dev.off()
p2 <- ggplot() + geom_sf(data = grid,  aes(fill = mean_me), col = NA) + scale_fill_viridis_c(limits = c(2.497,5.765)) + theme_bw() + labs(fill = " ")
png("mixed_mean.png", res = 300, width = 8, height = 6, units = "cm")
p2
dev.off()
p3 <- ggplot() + geom_sf(data = grid,  aes(fill = mean_sp), col = NA) + scale_fill_viridis_c(limits = c(2.497,5.765)) + theme_bw() + labs(fill = " ")
png("spatial_mean.png", res = 300, width = 8, height = 6, units = "cm")
p3
dev.off()
p4 <- ggplot() + geom_sf(data = grid,  aes(fill = mean_lm), col = NA) + scale_fill_viridis_c(limits = c(2.497,5.765)) + theme_bw() + labs(fill = " ")
png("linear_mean.png", res = 300, width = 8, height = 6, units = "cm")
p4
dev.off()

#sd
p5 <- ggplot() + geom_sf(data = grid,  aes(fill = sd), col = NA) + scale_fill_viridis_c(limits = c(0,0.48), option = "A") + theme_bw() + labs(fill = " ")
png("extended_sd.png", res = 300, width = 8, height = 6, units = "cm")
p5
dev.off()
p6 <- ggplot() + geom_sf(data = grid,  aes(fill = sd_me), col = NA) + scale_fill_viridis_c(limits = c(0,0.48), option = "A") + theme_bw() + labs(fill = " ")
png("mixed_sd.png", res = 300, width = 8, height = 6, units = "cm")
p6
dev.off()
p7 <- ggplot() + geom_sf(data = grid,  aes(fill = sd_sp), col = NA) + scale_fill_viridis_c(limits = c(0,0.48), option = "A") + theme_bw() + labs(fill = " ")
png("spatial_sd.png", res = 300, width = 8, height = 6, units = "cm")
p7
dev.off()
p8 <- ggplot() + geom_sf(data = grid,  aes(fill = sd_lm), col = NA) + scale_fill_viridis_c(limits = c(0,0.48), option = "A") + theme_bw() + labs(fill = " ")
png("linear_sd.png", res = 300, width = 8, height = 6, units = "cm")
p8
dev.off()






# map of predictors
p9 <- ggplot() + geom_sf(data = grid,  aes(fill = GESTEINKL), col = NA) + scale_fill_viridis_d() + theme_bw() + labs(fill = "Lithology")
png("gesteinkl.png", res = 300, width = 12, height = 8, units = "cm")
p9
dev.off()
# tectonic - solve issue with legend.....
p10 <- ggplot() + geom_sf(data = grid,  aes(fill = LEG_TEK_1), col = NA) + scale_fill_viridis_d() + theme_bw() + labs(fill = "Tectonic")
p10legend <- cowplot::get_legend(p10)
library(grid)
library(gridExtra)
grid.newpage()
grid.draw(p10legend)
png("tectonic_map.png", res = 300, width = 12, height = 8, units = "cm")
p10 + theme(legend.position = "none")
dev.off()
# land cover
p11 <- ggplot() + geom_sf(data = grid,  aes(fill = LC09_6), col = NA) + scale_fill_viridis_d() + theme_bw() + labs(fill = "Land coverage")
png("landcover.png", res = 300, width = 12, height = 8, units = "cm")
p11
dev.off()
#rainfall
p12 <- ggplot() + geom_sf(data = grid,  aes(fill = rain*1000), col = NA) + scale_fill_viridis_c() + theme_bw() + labs(fill = "Cumulative rainfall [mm]")
png("rainfall.png", res = 300, width = 12, height = 8, units = "cm")
p12
dev.off()
























###################
#Predictors
grid_ex$Preds <- exp(grid_ex$mean)
quilt.plot(grid_ex$X,grid_ex$Y, ifelse(grid_ex$Preds > 40,grid_ex$Preds,40),nx=300,ny=200,asp=1, nlevel=80, main = "Mean of Prediction",axes=FALSE) #screen values lower than 40nSv/h to enhance contrast?


grid_ex$Y <- grid_ex$Y*1000
grid_ex$X <- grid_ex$X*1000
coordinates(grid_ex) <- ~X+Y
proj4string(grid_ex) <- proj4string(CH_canton)
gridded(grid_ex) <-TRUE

tm_shape(CH_canton) + tm_borders(col = "slategrey") +
  tm_shape(CH_lakes) + tm_fill(col = "cadetblue1",alpha=0.4)+
  tm_shape(grid_ex) + tm_raster(col = "Preds", palette = "-plasma" ,style="cont",title = "Terrestrial Radiation [nSv/h]")+
  tm_layout(frame=FALSE, legend.text.size = 1.2, legend.title.size = 2)

terrmap <- tm_shape(CH_canton) + tm_borders(col = "slategrey") +
  tm_shape(CH_lakes) + tm_fill(col = "cadetblue1",alpha=0.4)+
  tm_shape(grid_ex) + tm_raster(col = "Preds", palette = "-plasma" ,style="cont",title = "Terrestrial Radiation [nSv/h]")+
  tm_layout(frame=FALSE, legend.text.size = 1.2, legend.title.size = 2 ,bg.color = "transparent") 
tmap_save(terrmap, "map_terrRad.jpg")

# Lithology
tm_shape(grid_ex) + tm_raster(col = "GESTEINKL", title = "Lithology") +
  tm_layout(legend.outside = TRUE, frame = FALSE, legend.text.size = 0.86, legend.title.size = 1.3)
Litho <- tm_shape(grid_ex) + tm_raster(col = "GESTEINKL", title = "Lithology") +
  tm_layout(legend.outside = TRUE, frame = FALSE, legend.text.size = 0.86, legend.title.size = 1.3)
tmap_save(Litho, "lithology.jpg")

# Tectonics
tm_shape(grid_ex) + tm_raster(col = "LEG_TEK_1", title = "Tectonic units") +
  tm_layout(legend.outside = TRUE, frame = FALSE, legend.title.size = 1.3)
Tecto <- tm_shape(grid_ex) + tm_raster(col = "LEG_TEK_1", title = "Tectonic units") +
  tm_layout(legend.outside = TRUE, frame = FALSE, legend.title.size = 1.3)
tmap_save(Tecto, "tectonics.jpg")

# Land cover
tm_shape(grid_ex) + tm_raster(col = "LC09_6", title = "Land coverage") +
  tm_layout(legend.outside = TRUE, frame = FALSE, legend.title.size = 1.3, legend.text.size = 1.2)
LC <- tm_shape(grid_ex) + tm_raster(col = "LC09_6", title = "Land coverage") +
  tm_layout(legend.outside = TRUE, frame = FALSE, legend.title.size = 1.3, legend.text.size = 1.2)
tmap_save(LC, "landcover.jpg")

# rain fall
ScaleRain <- readRDS("../_15_noflight/Input/SCALE_rain.rds")
grid_ex$rain <- grid_ex$rain*ScaleRain
tm_shape(grid_ex) + tm_raster(col = "rain", title = "Cumulative rainfall [mm]", style = "cont") +
  tm_layout(legend.outside = TRUE, frame = FALSE, legend.title.size = 1.5, legend.text.size = 1.1)
rainfall <- tm_shape(grid_ex) + tm_raster(col = "rain", title = "Cumulative rainfall [mm]", style = "cont") +
  tm_layout(legend.outside = TRUE, frame = FALSE, legend.title.size = 1.5, legend.text.size = 1.1)
tmap_save(rainfall, "rainfall.jpg")
