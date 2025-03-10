rm(list = ls())
# -------------------------------------------------------------------------
library(VAST)
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
library(boot)
library(ggeffects)
library(pdp)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Select species and school type:
n_sp = 5 # first N species (only applies to weight or numbers datasets)
n_fac = 2 # Number of factors
this_type = 'FOB'

# Define folder to save results
plot_folder = file.path('figures', this_type, 'VAST-multi')
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Define model folder:
model_folder = file.path('models', this_type, 'VAST-multi')
dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path('data', this_type, 'weight_data.rds'))
extraRegion_VAST = readRDS(file = file.path('data', this_type, 'extraRegion_VAST.rds'))
extraRegion_tinyVAST = readRDS(file = file.path('data', this_type, 'extraRegion_tinyVAST.rds'))

# -------------------------------------------------------------------------
# Filter first N species:
mod_data = weight_data
cumsp_data = mod_data %>% group_by(sp_name) %>% summarise(Catch = sum(Catch))
cumsp_data = arrange(cumsp_data, desc(Catch))
cumsp_data = cumsp_data[1:n_sp, ]
sp_data = mod_data %>% dplyr::filter(sp_name %in% cumsp_data$sp_name) %>%
  mutate(sp_number = factor(sp_name, levels = cumsp_data$sp_name, labels = 1:n_sp))
# Define sp ID:
sp_data = sp_data %>% mutate(sp_number = as.integer(as.character(sp_number)))
sp_data = as.data.frame(sp_data)

# Make settings:
settings <- make_settings(n_x = 100, Region='User',
                          purpose = "ordination", bias.correct = FALSE,
                          use_anisotropy = FALSE, 
                          fine_scale = TRUE,
                          knot_method = 'grid',
                          ObsModel = c(4,1), # Delta-gamma/poisson link delta model
                          n_categories = n_fac)
settings$FieldConfig[c('Omega','Epsilon','Beta'),'Component_1'] = c(0,0,'IID')
settings$FieldConfig[c('Omega','Epsilon','Beta'),'Component_2'] = c(n_fac,0,'IID')
settings$RhoConfig = c("Beta1" = 3, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0)

# Fit VAST model
VModel_j <- fit_model(settings=settings,
                 Lat_i=sp_data$Lat, 
                 Lon_i=sp_data$Lon,
                 t_i=sp_data$Year, 
                 b_i=sp_data$Catch,
                 a_i=sp_data$AreaSwept_km2,
                 c_i=sp_data$sp_number-1,
                 newtonsteps = 0,
                 getsd = FALSE,
                 input_grid=extraRegion_VAST)

VModel_j$parameter_estimates$diagnostics
dim(VModel_j$Report$Omega2_sc)
VModel_j$Report$Omega2_sc[1:10, ]
VModel_j$Report$L_omega2_cf
VModel_j$Report$lowercov_uppercor_omega2
VModel_j$Report$Omegainput2_sf
VModel_j$Report$D_gct

L_list = VModel_j$Report$L_omega2_cf
rownames(L_list) = cumsp_data$sp_name
Psi_sjt = VModel_j$Report$Omegainput2_sf
Psi_gjt = VModel_j$Report$Omegainput2_gf
logkappa = unlist(VModel_j$ParHat[c('logkappa1','logkappa2')])[2]
tau = 1 / (exp(logkappa) * sqrt(4*pi))

Var_rot = rotate_factors( L_pj = L_list,
                          Psi_sjt = Psi_gjt/tau,
                          RotationMethod = "PCA",
                          testcutoff = 1e-4,
                          quiet = TRUE )
Var_rot$Psi_rot


save(VModel_j, file = file.path(model_folder, 'fit.RData'))

# Plot results
plot_results( VModel_j, plot_set = c(3,17), category_names = cumsp_data$sp_name,
              working_dir = file.path(getwd(), plot_folder))


# -------------------------------------------------------------------------
# Plot correlations (Omega2)
Cov_omega2 = fit$Report$L_omega2_cf %*% t(fit$Report$L_omega2_cf)
cor_mat = cov2cor(Cov_omega2)
colnames(cor_mat) = cumsp_data$sp_code
rownames(cor_mat) = cumsp_data$sp_code
corrplot(cor_mat ,method="pie", type="lower")
corrplot.mixed( cor_mat )


# -------------------------------------------------------------------------
# Plot cluster (Omega 2)
Cov_omega2 = fit$Report$L_omega2_cf %*% t(fit$Report$L_omega2_cf)
sel_mat = fit$Report$L_omega2_cf
rownames(sel_mat) = cumsp_data$sp_code
Dist = dist(sel_mat, diag=TRUE, upper=TRUE)
plot(hclust(Dist))