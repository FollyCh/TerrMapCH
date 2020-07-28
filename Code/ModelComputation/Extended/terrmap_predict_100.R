###################################3
#
# Compute predictions on 100x100m grid
#
###################################

# load packages
# load INLA
library(INLA)

# inla load dynamic
INLA:::inla.dynload.workaround() #needed on cluster




# load data
aero <- readRDS("Output/aero_results.rds")
aero.grid <- readRDS("../Input/grid_cluster.rds")
res <- readRDS("Output/results.rds")
grid <- readRDS("../Input/aero.grid_100m_df.rds")
msh1 <- readRDS("../Output/mesh_u1.rds")
msh2 <- readRDS("../Output/mesh_u2.rds")
#library(sf)
#grid <- rbind(st_coordinates(grid)/1000,st_set_geometry(grid,NULL))
#grid <- grid[,names(aero.grid)]
#grid <- grid[!(grid$GESTEINKL %in% "Gewässer, Seen"),]
#grid$GESTEINKL <- as.factor(as.character(grid$GESTEINKL))
#levels(grid$GESTEINKL)
#saveRDS(grid,"../Input/aero.grid_100m_df.rds")



# prepare stack pred100
loc.grid100 <- as.matrix(grid[,c("X","Y")])
loc.grid <- as.matrix(aero.grid[,c("X","Y")])
loc.aero <- as.matrix(aero[,c("X","Y")])
A.grid100_1 <- inla.spde.make.A(mesh = msh1,
                            loc = loc.grid100)

A.grid100_2 <- inla.spde.make.A(mesh = msh2,
                            loc = loc.grid100)
 

# reprepare original data stack
A.train1 <- inla.spde.make.A(mesh = msh1,
                             loc = loc.aero)
A.grid1 <- inla.spde.make.A(mesh = msh1,
                            loc = loc.grid)
A.train2 <- inla.spde.make.A(mesh = msh2,
                             loc = loc.aero)
A.grid2 <- inla.spde.make.A(mesh = msh2,
                            loc = loc.grid)

# == model definition ==
# spde model
rf1 <- inla.spde2.pcmatern(mesh = msh1,
                           prior.range = c(15,0.6),
                           prior.sigma = c(10,0.01))
rf2 <- inla.spde2.pcmatern(mesh = msh2,
                           prior.range = c(2,0.02),
                           prior.sigma = c(10,0.01))
# formula 
formula <- LogTerr ~ 0 + gestein + tectonic + landcover + rain + f(field1, model = rf1) + f(field2, model = rf2)
# stack train
stk.train <- inla.stack(data = list(LogTerr = aero$LogTerr),
                        A = list(A.train1,A.train2, 1),
                        effects = list(field1 = 1:rf1$n.spde,
                                       field2 = 1:rf2$n.spde,
                                       data.frame(gestein = aero$GESTEINKL,
                                                  tectonic = aero$LEG_TEK_1,
                                                  landcover = aero$LC09_6,
                                                  rain = aero$rain)),
                        tag = "train")
# stk grid
stk.grid <- inla.stack(data = list(LogTerr = NA),
                       A = list(A.grid1, A.grid2, 1),
                       effects = list(field1 = 1:rf1$n.spde,
                                      field2 = 1:rf2$n.spde,
                                      data.frame(gestein = aero.grid$GESTEINKL,
                                                 tectonic = aero.grid$LEG_TEK_1,
                                                 landcover = aero.grid$LC09_6,
                                                 rain = aero.grid$rain)),
                       tag = "grid")

####
# stack 100m
stk.grid100 <- inla.stack(data = list(LogTerr = NA),
                        A = list(A.grid100_1,A.grid100_2, 1),
                        effects = list(field1 = 1:rf1$n.spde,
                                       field2 = 1:rf2$n.spde,
                                       data.frame(gestein = grid$GESTEINKL,
                                                  tectonic = grid$LEG_TEK_1,
                                                  landcover = grid$LC09_6,
                                                  rain = grid$rain)),
                        tag = "grid100")
#combine 
fullstk <- inla.stack(stk.train, stk.grid,stk.grid100)


# compute predictions
set.seed(74)
res <- inla(formula,
            family = "Gaussian",
            data = inla.stack.data(fullstk),
            control.predictor = list(A = inla.stack.A(fullstk), link = 1,compute = TRUE),
            control.mode = list(restart = TRUE, theta = res$mode$theta),
            control.inla = list(int.strategy = "eb"),
            quantiles = c(0.05,0.5,0.95),
            control.compute = list(dic = T, waic = T, config = TRUE),
            verbose = T)


# save res
# get indexes
idx.train <- inla.stack.index(fullstk,"train")$data
idx.grid <- inla.stack.index(fullstk, "grid")$data
idx.grid100 <- inla.stack.index(fullstk,"grid100")$data
# extract predicted
aero.grid <- cbind(aero.grid, "mean" = res$summary.fitted.values[idx.grid,c("mean")])
aero.grid <- cbind(aero.grid, "sd" = res$summary.fitted.values[idx.grid,c("sd")])
# extract grid100
grid <- cbind(grid, "mean" = res$summary.fitted.values[idx.grid100,c("mean")])
grid <- cbind(grid, "sd" = res$summary.fitted.values[idx.grid100,c("sd")])
# extract fitted
aero <- cbind(aero, "mean" = res$summary.fitted.values[idx.train,c("mean")])
aero <- cbind(aero, "sd" = res$summary.fitted.values[idx.train,c("sd")])

#save res
saveRDS(res, "Output100/results.rds")
saveRDS(aero, "Output100/aero_results.rds")
saveRDS(aero.grid,"Output100/grid_results.rds")
saveRDS(grid, "Output100/grid100_res.rds")

















