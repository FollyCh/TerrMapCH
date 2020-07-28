#################################################
#
#  FIT Model only spatial process
#
##################################################

# load INLA
library(INLA)

# inla load dynamic
#INLA:::inla.dynload.workaround() #needed on cluster


# readData
aero <- readRDS("../Input/aero_cluster.rds")
aero.grid <- readRDS("../Input/grid_cluster.rds")

# == locations ==
loc.grid <- cbind(aero.grid$X, aero.grid$Y)
loc.aero <- cbind(aero$X,aero$Y)

# == mesh, projector matrix ==
# create mesh 
set.seed(2)
msh <- inla.mesh.2d(loc = loc.aero,
                    max.edge = c(5,15),
                    cutoff = 3.5,
                    min.angle = c(31,29),
                    offset = c(5,15))
# projector matrix
A.train <- inla.spde.make.A(mesh = msh,
                            loc = loc.aero)
A.grid <- inla.spde.make.A(mesh = msh,
                           loc = loc.grid)

# == model definition ==
# spde model
rf <- inla.spde2.pcmatern(mesh = msh,
                          prior.range = c(15,0.5),
                          prior.sigma = c(10,0.01))
# formula 
formula <- LogTerr ~ 0 + gestein + tectonic + landcover + rain + f(field, model = rf)
# stack train
stk.train <- inla.stack(data = list(LogTerr = aero$LogTerr),
                        A = list(A.train, 1),
                        effects = list(field = 1:rf$n.spde,
                                       data.frame(gestein = aero$GESTEINKL,
                                                  tectonic = aero$LEG_TEK_1,
                                                  landcover = aero$LC09_6,
                                                  rain = aero$rain)),
                        tag = "train")
# stk grid
stk.grid <- inla.stack(data = list(LogTerr = NA),
                       A = list(A.grid, 1),
                       effects = list(field = 1:rf$n.spde,
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
            control.compute = list(dic = T, waic = T),
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
saveRDS(msh, "Output/mesh.rds")
saveRDS(res, "Output/results.rds")
saveRDS(aero, "Output/aero_results.rds")
saveRDS(aero.grid,"Output/grid_results.rds")


