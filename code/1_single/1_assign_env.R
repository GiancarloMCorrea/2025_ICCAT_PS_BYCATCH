rm(list = ls())

# Define type of school to be analyzed:
this_type = 'FOB' # or FSC
source('code/1_single/load_libs.R')

# -------------------------------------------------------------------------
source("C:/Use/GitHub/extractOceanVariables/code/copernicus/multiple/matchCOPERNICUS.R")
source("C:/Use/GitHub/extractOceanVariables/code/auxFunctions.R")

# -------------------------------------------------------------------------
# Read data in:
effPoints = readRDS(file.path(data_folder, 'effPoints.rds'))
obsPoints = readRDS(file.path(data_folder, 'obsPoints.rds'))

# Assign environmental information:
# -------------------------------------------------------------------------
# Potential temperature surface:

# Define Lan/Lot and Date column names in your dataset:
lonlat_cols = c("lon", "lat")
date_col = "date"
fields = "thetao"
savedir = "C:/Use/OneDrive - AZTI/Data/ICCAT/Env_Data_ATL_2010-onwards/thetao_depth0-100_P1M-m/"
obsPoints = matchCOPERNICUS(data           = obsPoints, 
                            lonlat_cols    = lonlat_cols,
                            date_col       = date_col,
                            var_label      = fields, 
                            var_path     = savedir,
                            depth_range  = c(0, 5),
                            depth_FUN    = 'mean')
effPoints = matchCOPERNICUS(data           = effPoints, 
                            lonlat_cols    = lonlat_cols,
                            date_col       = date_col,
                            var_label      = fields, 
                            var_path     = savedir,
                            depth_range  = c(0, 5),
                            depth_FUN    = 'mean')

# # -------------------------------------------------------------------------
# # Mixing layer depth:
# 
# # Define Lan/Lot and Date column names in your dataset:
# lonlat_cols = c("longitude", "latitude")
# date_col = "date"
# fields = "mlotst"
# saveEnvDir = "C:/Use/OneDrive - AZTI/Data/ICCAT/Env_Data_ATL_2010-onwards/mlotst_depth0-100_P1M-m"
# obsPoints = matchCOPERNICUS(data           = obsPoints, 
#                             lonlat_cols    = lonlat_cols,
#                             date_col       = date_col,
#                             var_label      = fields, 
#                             varPath        = saveEnvDir)
# effPoints = matchCOPERNICUS(data           = effPoints, 
#                             lonlat_cols    = lonlat_cols,
#                             date_col       = date_col,
#                             var_label      = fields, 
#                             varPath        = saveEnvDir)
# 
# # -------------------------------------------------------------------------
# # NPP:
# 
# # Define Lan/Lot and Date column names in your dataset:
# lonlat_cols = c("longitude", "latitude")
# date_col = "date"
# fields = "nppv"
# saveEnvDir = "C:/Use/OneDrive - AZTI/Data/ICCAT/Env_Data_ATL_2010-onwards/nppv_depth0-100_P1M-m"
# obsPoints = matchCOPERNICUS(data           = obsPoints, 
#                             lonlat_cols    = lonlat_cols,
#                             date_col       = date_col,
#                             summ_fun       = 'sum',
#                             var_label      = fields, 
#                             varPath        = saveEnvDir)
# effPoints = matchCOPERNICUS(data           = effPoints, 
#                             lonlat_cols    = lonlat_cols,
#                             date_col       = date_col,
#                             summ_fun       = 'sum',
#                             var_label      = fields, 
#                             varPath        = saveEnvDir)

# -------------------------------------------------------------------------
# Check env data:
summary(obsPoints)
summary(effPoints)

# Rename:
obsPoints = obsPoints %>% rename(sst = thetao)
effPoints = effPoints %>% rename(sst = thetao)

# -------------------------------------------------------------------------
# Save data with non standardized env information:
saveRDS(obsPoints, file = file.path(data_folder, 'obsPoints_nostd.rds'))
saveRDS(effPoints, file = file.path(data_folder, 'effPoints_nostd.rds'))

# -------------------------------------------------------------------------
# Standardized covariates to explore effect better:

# Select variables to standardize:
std_cov = c('sst', 'trop_catch')

# Std covariates in observations
for(j in seq_along(std_cov)) {
  mean_var = mean(obsPoints[,std_cov[j]], na.rm = T)
  sd_var = sd(obsPoints[,std_cov[j]], na.rm = T)
  obsPoints[,std_cov[j]] = (obsPoints[,std_cov[j]] - mean_var)/sd_var
  effPoints[,std_cov[j]] = (effPoints[,std_cov[j]] - mean_var)/sd_var
}

# -------------------------------------------------------------------------
# Save data with env information:
saveRDS(obsPoints, file = file.path(data_folder, 'obsPoints.rds'))
saveRDS(effPoints, file = file.path(data_folder, 'effPoints.rds'))
