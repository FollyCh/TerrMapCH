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
aero <- readRDS("../Input/aero_cluster.rds")
aero$idx <- readRDS("randomCV_index.rds")
# subset
valid <- aero[aero$idx %in% "valid",]
aero <- aero[aero$idx %in% "train",]

# == locations ==
loc.grid <- cbind(valid$X, valid$Y)
loc.aero <- cbind(aero$X,aero$Y)

# == mesh, projector matrix ==
# create mesh
msh1 <- readRDS("../Output/mesh_u1.rds")
msh2 <- readRDS("../Output/mesh_u2.rds")
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
                                      data.frame(gestein = valid$GESTEINKL,
                                                 tectonic = valid$LEG_TEK_1,
                                                 landcover = valid$LC09_6,
                                                 rain = valid$rain)),
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
valid <- cbind(valid, "mean" = res$summary.fitted.values[idx.grid,c("mean")])
valid <- cbind(valid, "sd" = res$summary.fitted.values[idx.grid,c("sd")])
# extract fitted
aero <- cbind(aero, "mean" = res$summary.fitted.values[idx.train,c("mean")])
aero <- cbind(aero, "sd" = res$summary.fitted.values[idx.train,c("sd")])




# = = save = = 
saveRDS(res, "Output/results_cvr.rds")
saveRDS(aero, "Output/aero_results_cvr.rds")
saveRDS(valid,"Output/valid_results_cvr.rds")


