#######################################
#
# Prepare the data to run the model
#
#######################################



# == load needed packages == 
library(tidyverse)
library(Metrics)
library(blockCV)
library(sf)
library(tmap)
library(tmaptools)
library(sp)
library(gstat)
library(variosig)
library(raster)
library(ncdf4)
library(viridis)
library(fields)
library(INLA)


# == Load Data ==
#aeroradiometrie
aero <- readRDS("../Data/processed/aero_data_clean.rds")
#grid
aero.grid <- readRDS("../Data/processed/aero.grid.rds")
#map of switzerland
#CH_canton <- as(read_shape("Y:/ENVEPI/RESTRICTED/temp/p_background_radiation/origdata/swissBOUNDERIES3D/BOUNDARIES_2018/DATEN/swissBOUNDARIES3D/SHAPEFILE_LV03_LN02/swissBOUNDARIES3D_1_3_TLM_KANTONSGEBIET.shp"),"Spatial")
#swissborder <- CH_canton %>% st_as_sf() %>% st_union() %>% as("Spatial") 






##############################################
# == Some Cleaning & Preparation ==

# subset grid for values on land
idgenese <- (aero.grid$GENESE %in% c("GewÃ¤sser, Seen"))#~17'000, "Gletscher, Firn"
idlanduse <- (aero.grid$LU09_46 %in% c(401)) #401 - Seen; 402 - Flüsse,Bäche; ~32'000 (~17'000 over lakes, Flüsse/Bäche seems to have little influence, but could contai also smaller lakes..)   
idlandcover <- (aero.grid$LC09_27 %in% c(61,63,64)) #61 - Wasser; 62 - Gletscher/Firn,62; 63 - Nassstandorte; 64 - Schilfbestände
idwater <- (idgenese | idlanduse | idlandcover) 
sum(idwater) #3342 measurements
grid.water <- aero.grid[idwater,]
aero.grid <- aero.grid[!idwater,] 
#tm_shape(grid.water)+tm_dots(col = "GESTEINKL")+tm_layout(title = "excluded grid points") 
rm(idgenese, idlanduse, idlandcover, idwater)

# negative measurements are physical nonsense, same for 0. 
# working with log -> some lower limit to guarantee symmetry of distribution can make sense
# sum(aero$TERR<7) only 325 measurements are lower than 7nSv/h.. could well be influenced by water bodies
#tmap_mode("view")
#tm_shape(aero[(aero$TERR < 7),])+tm_dots() #visual inspection: most water or glacier related...
aero <- aero[(aero$TERR > 7),] #325 measurements excluded
plot(density(log(aero$TERR)))

# compare dist with standard normal
x <- seq(0, 8, length=100)
hx <- dnorm(x, mean = mean(log(aero$TERR)), sd = sd(log(aero$TERR)))
plot(density(log(aero$TERR)))
lines(x,hx,col = "red") #TERR slightly more narrow,but larger tails

# factorize
aero$C_ID <- as.factor(aero$C_ID)
aero.grid$C_ID <- as.factor(aero.grid$C_ID)
aero.grid$GESTEINKL <- as.factor(as.character(aero.grid$GESTEINKL))
aero$GESTEINKL <- as.factor(as.character(aero$GESTEINKL))
aero.grid$GENESE <- as.factor(as.character(aero.grid$GENESE))
aero$GENESE <- as.factor(as.character(aero$GENESE))
aero.grid$LU09_46 <- as.factor(as.character(aero.grid$LU09_46))
aero$LU09_46 <- as.factor(as.character(aero$LU09_46))
aero.grid$LU09_10 <- as.factor(as.character(aero.grid$LU09_10))
aero$LU09_10 <- as.factor(as.character(aero$LU09_10))
aero.grid$LU09_4 <- as.factor(as.character(aero.grid$LU09_4))
aero$LU09_4 <- as.factor(as.character(aero$LU09_4))
aero.grid$LC09_27 <- as.factor(as.character(aero.grid$LC09_27))
aero$LC09_27 <- as.factor(as.character(aero$LC09_27))
aero.grid$LC09_6 <- as.factor(as.character(aero.grid$LC09_6))
aero$LC09_6 <- as.factor(as.character(aero$LC09_6))





##################################################################
# == Sample every 15th measurement ==

# indexes to be included
idx <- 20 * (1:as.integer(nrow(aero)/20)) - as.integer(runif(1,min=0, max=20))

# not included in modeling
temp <- aero[-idx,]

# save not included and remove
saveRDS(temp, "Input/excluded_measurements20.rds")
rm(temp)

# subset aero to move on with random sample of 100'000 measurements  
aero <- aero[idx,]







##################################################################
# == Folds for spatial cross validation ==

# folds of data for spatial CV
# not sure how well that actually works....
set.seed(123)
spCV_idx <- spatialBlock(aero, #that one takes a while to compute... ~3-4 minutes
                         selection = "random",
                         theRange = 15000,
                         iteration = 30,k = 4)
aero$fold <- as.factor(spCV_idx$foldID)
# visualize
#tmap_mode("plot")
tm_shape(aero[sample(1:nrow(aero), size = 30000),])+
  tm_dots(col = "fold") + 
  tm_layout(bg.color = "transparent",frame = F,
            panel.labels = "Folds for spatial cross-validation",
            panel.show = T)




#######################################################
# == Data preparation 2 == Data format

# the sf datatype cannot be used for INLA
aero <- cbind(st_set_geometry(aero,NULL),st_coordinates(aero))
aero.grid <- st_as_sf(aero.grid)
aero.grid <- cbind(st_set_geometry(aero.grid,NULL),st_coordinates(aero.grid))

# rescale
# .....of coordinates          
# borders
#swissboundary <- spTransform(swissborder, CRSobj = "+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +units=km +no_defs")
#swissboundary <- swissboundary@polygons[[1]]@Polygons[[1]]@coords
# points
aero[,c("X","Y")] <- aero[,c("X","Y")]/1000
aero.grid[,c("X","Y")] <- aero.grid[,c("X","Y")]/1000
# ......of numeric covariates


# rybach: cs + terr || not needed
#aero$Ryb_TerrCaes <- aero$Ryb_Terrestrial + aero$Chernobyl_Fall_Out
#aero.grid$Ryb_TerrCaes <- aero.grid$Ryb_Terrestrial + aero.grid$Chernobyl_Fall_Out
#rybMEAN <- mean(aero.grid$Ryb_TerrCaes)
#rybSD <- sd(aero.grid$Ryb_TerrCaes)
#aero$Ryb_TerrCaes <- (aero$Ryb_TerrCaes - rybMEAN)/rybSD
#aero.grid$Ryb_TerrCaes <- (aero.grid$Ryb_TerrCaes - rybMEAN)/rybSD
#rybach terr
#rybMEAN_terr <- mean(aero.grid$Ryb_Terrestrial)
#rybSD_terr <- mean(aero.grid$Ryb_Terrestrial)
#aero$Ryb_Terrestrial <- (aero$Ryb_Terrestrial - rybMEAN_terr)/rybSD_terr
#aero.grid$Ryb_Terrestrial <- (aero.grid$Ryb_Terrestrial - rybMEAN_terr)/rybSD_terr

# Redefine precipitation variables
# include 30.april (eastern switzerland) until 5th of may 
# threshold ~10^-1 Bq/m3
# https://www.youtube.com/watch?v=mqu_l29WioM : cloud arrived on 30.april, goes away 5th of may
aero$rain <- aero$Rain_april30 + aero$Rain_may1 + aero$Rain_may2 + aero$Rain_may3 + aero$Rain_may4 + aero$Rain_may5
# needed 1km grid - already done 100m grid
aero.grid$rain <- aero.grid$Rain_april30 + aero.grid$Rain_may1 + aero.grid$Rain_may2 + aero.grid$Rain_may3 + aero.grid$Rain_may4 + aero.grid$Rain_may5

#plot(density(aero.grid$rain))

#quilt.plot(aero.grid$X, aero.grid$Y, aero.grid$rain, nlevel = 20, nx = 100, ny = 50, FUN = mean)

# rescale rain
SCALE_rain <- 1000
aero$rain <- aero$rain/SCALE_rain
aero.grid$rain <- aero.grid$rain/SCALE_rain


# subset for needed variables
aero <- aero[,c("X","Y","fold","LEG_TEK_1","C_ID","GENESE","GESTEINKL","rain","LC09_6","TERR","Flight")]
aero.grid <- aero.grid[,c("X","Y","LEG_TEK_1","C_ID","GENESE","GESTEINKL","rain","LC09_6")]

# add log of terrestrial radiation
aero$LogTerr <- log(aero$TERR)

# save data for use on cluster
saveRDS(aero,"Input/aero_cluster20.rds")
saveRDS(aero.grid, "Input/grid_cluster20.rds")
saveRDS(SCALE_rain, "Input/SCALE_rain20.rds")



