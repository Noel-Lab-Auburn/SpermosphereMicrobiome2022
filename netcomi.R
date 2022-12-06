
# read in rds files. 
assoMat1 <- readRDS("/scratch/aubzxn001/netcomi/soybean_association_network.rds")
assoMat2 <- readRDS("/scratch/aubzxn001/netcomi/cotton_association_network.rds")

# Network construction (pass association matrices to netConstruct)
# - sparsMethod must be set to "none" because sparsification is already included in SpiecEasi
sperm.crossdomain.cottonvsoybean <- netConstruct(data = assoMat1, data2 = assoMat2, 
                                 dataType = "condDependence", 
                                 sparsMethod = "none")

# Network analysis
netprops_sperm.crossdomain.cottonvsoybean <- netAnalyze(sperm.crossdomain.cottonvsoybean, hubPar = "eigenvector")

