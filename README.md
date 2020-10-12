# TerrMapCH

This repository contains:
 - A predicted map for terrestrial radiation in Switzerland as shapefile. 
    - Variable Pred_Terr: posterior mean of predictions, unit is nSv/h
    - Variable Pred_logTerr: posterior mean, on log-scale
    - Variable Unc_logTerr: posterior standard deviation
 - 100 realisations from the joint posterior, delivered in 4 seperate shapefiles.
    - Variables sample_i: realisation i, with i in 1:100, unit is nSv/h
 - the R Code used to generate the map.
 
Reference: [arXiv:2010.00534](arXiv:2010.00534) 


