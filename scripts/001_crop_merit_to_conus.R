## This script crops the downloaded MERIT datasets down to the spatial extent of 
## CONUS. We do this step first since this dataset is one of the major bottlenecks
## due to its size and the number of different computations we need to do on it
## We do need to compute for both the catchments ("cat") and flowine ("riv") files, 
## the catchment files are used to verify we're matching Maavara, and ID any 
## problematic catchments. The flowline file is why we actually need this dataset

