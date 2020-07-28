###############################3
#
# Compute posterior samples
#
#################################

#packages
# load INLA
library(INLA)

# inla load dynamic
INLA:::inla.dynload.workaround() #needed on cluster

##############################
#create stacks used for model
# readData
aero <- readRDS("../Input/aero_cluster.rds")
aero.grid <- readRDS("../Input/grid_cluster.rds")

# == locations ==
loc.grid <- cbind(aero.grid$X, aero.grid$Y)
loc.aero <- cbind(aero$X,aero$Y)

# == mesh, projector matrix == CREATED FOR ITERATIVE FITTING; USE SAME MESH HERE
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
                                      data.frame(gestein = aero.grid$GESTEINKL,
                                                 tectonic = aero.grid$LEG_TEK_1,
                                                 landcover = aero.grid$LC09_6,
                                                 rain = aero.grid$rain)),
                       tag = "grid")
# combine stack
fullstk <- inla.stack(stk.train, stk.grid)

#################

# load data
res <- readRDS("Output/results.rds")


# sample from res
N <- 5
posteriors <- inla.posterior.sample(n = N, res)


# merge with grid
nms <- names(aero.grid)
for(i in 1:N){
  aero.grid <- cbind(aero.grid, posteriors[[i]]$latent[inla.stack.index(fullstk,"grid")$data,])
}
names(aero.grid) <- c(nms, paste0("sample_",1:N))


# save
saveRDS(aero.grid,"Output/posterior_samples_1.rds")
saveRDS(posteriors, "Output/posteriors_1til_5.rds")















