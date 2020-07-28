#################################################
#
#  FIT Model only spatial process
#
##################################################

# load INLA
library(INLA)

# inla load dynamic
INLA:::inla.dynload.workaround() #needed on cluster


# readData
aero <- readRDS("Input/aero_cluster.rds")
aero.grid <- readRDS("Input/grid_cluster.rds")

# == locations ==
loc.grid <- cbind(aero.grid$X, aero.grid$Y)
loc.aero <- cbind(aero$X,aero$Y)

# == mesh, projector matrix == CREATED FOR ITERATIVE FITTING; USE SAME MESH HERE
# create mesh
set.seed(1)
msh1 <- inla.mesh.2d(loc = loc.aero,
                     max.edge = c(5,10),
                     cutoff = 3,
                     min.angle = c(31,30),
                     offset = c(3, 10))
set.seed(1)
msh2 <- inla.mesh.2d(loc = loc.aero,
                     max.edge = c(0.6,5),
                     cutoff = 0.45,
                     min.angle = c(31,30),
                     offset = c(3, 10))
# projector matrix
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
# combine stack
fullstk <- inla.stack(stk.train, stk.grid)

# == compute model and predictions ==
# compute
set.seed(74)
res <- inla(formula,
            family = "Gaussian",
            data = inla.stack.data(fullstk),
            control.predictor = list(A = inla.stack.A(fullstk), compute = TRUE),
            control.inla = list(int.strategy = "eb"),
            quantiles = c(0.05,0.5,0.95),
            control.compute = list(dic = T, waic = T, config = TRUE),
            verbose = T)

# == Extract results == 
# get indexes
idx.train <- inla.stack.index(fullstk,"train")$data
idx.grid <- inla.stack.index(fullstk, "grid")$data
# extract predicted
aero.grid <- cbind(aero.grid, "mean" = res$summary.fitted.values[idx.grid,c("mean")])
aero.grid <- cbind(aero.grid, "sd" = res$summary.fitted.values[idx.grid,c("sd")])
# extract fitted
aero <- cbind(aero, "mean" = res$summary.fitted.values[idx.train,c("mean")])
aero <- cbind(aero, "sd" = res$summary.fitted.values[idx.train,c("sd")])



# = = save = = 
saveRDS(res, "Output/results.rds")
saveRDS(aero, "Output/aero_results.rds")
saveRDS(aero.grid,"Output/grid_results.rds")


