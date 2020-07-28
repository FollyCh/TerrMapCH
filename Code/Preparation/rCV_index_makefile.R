aero <- readRDS("Input/aero_cluster.rds")

idx <- rep("train",nrow(aero))
set.seed(1234567)
idx[sample(1:nrow(aero), size = length(idx)*0.3)] <- "valid"
saveRDS(idx, "randomCV_index.rds")


