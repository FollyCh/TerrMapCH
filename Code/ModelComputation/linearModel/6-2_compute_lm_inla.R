##############################
# Compute linear model
#############################

# packages
library(INLA)

# work around for issues on cluster
INLA:::inla.dynload.workaround()

# load data
aero <- readRDS("Input/aero_cluster.rds")
aero.grid <- readRDS("Input/grid_cluster.rds")

# define model
formula <- LogTerr ~ -1 + GESTEINKL + LEG_TEK_1+ LC09_6 + rain

# compute model
set.seed(3)
res <- inla(formula = formula, family = "gaussian",
            data = rbind(aero[,c("GESTEINKL","LEG_TEK_1","LC09_6","rain", "LogTerr")],
                         cbind(aero.grid[,c("GESTEINKL","LEG_TEK_1","LC09_6","rain")], "LogTerr" = NA)),
            control.predictor=list(compute=TRUE),
            verbose=T) 

# extract fittet and predicted values
aero$mean <- res$summary.fitted.values[1:nrow(aero),"mean"]
aero$sd <- res$summary.fitted.values[1:nrow(aero),"sd"]
#
aero.grid$mean <- res$summary.fitted.values[(nrow(aero)+1):nrow(res$summary.fitted.values),"mean"]
aero.grid$sd <- res$summary.fitted.values[(nrow(aero)+1):nrow(res$summary.fitted.values),"sd"]

# save results
saveRDS(res, "Output/results.rds")
saveRDS(aero, "Output/aero_results.rds")
saveRDS(aero.grid, "Output/grid_results.rds")