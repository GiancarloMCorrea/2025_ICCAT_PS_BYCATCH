rm(list = ls())
# -------------------------------------------------------------------------
library(tinyVAST)
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
library(boot)
library(ggeffects)
library(pdp)
library(fmesher)
library(sdmTMB)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Select species and school type:
n_sp = 5 # first N species (only applies to weight or numbers datasets)
n_fac = 2 # Number of factors
this_type = 'FOB'

# Define folder to save results
plot_folder = file.path('figures', this_type, 'tinyVAST-multi')
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Define model folder:
model_folder = file.path('models', this_type, 'tinyVAST-multi')
dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path('data', this_type, 'weight_data.rds'))
extraRegion_tinyVAST = readRDS(file = file.path('data', this_type, 'extraRegion_tinyVAST.rds'))

# -------------------------------------------------------------------------
# Filter species:
mod_data = weight_data
cumsp_data = mod_data %>% group_by(sp_name) %>% summarise(Catch = sum(Catch))
cumsp_data = arrange(cumsp_data, desc(Catch))
cumsp_data = cumsp_data[1:n_sp, ]
sp_data = mod_data %>% dplyr::filter(sp_name %in% cumsp_data$sp_name) %>%
            mutate(sp_number = factor(sp_name, levels = cumsp_data$sp_name, labels = 1:n_sp))
# Define sp ID:
sp_data = sp_data %>% mutate(sp_number = as.integer(as.character(sp_number)))
sp_data = as.data.frame(sp_data)

# Make mesh
my_mesh = fm_mesh_2d( sp_data[,c('Lon', 'Lat')], cutoff = 3 )
p1 = ggplot() + geom_fm(data = my_mesh)
ggsave(paste0('map_mesh', img_type), path = plot_folder, 
       plot = p1, width = img_width*0.5, height = 80, units = 'mm', dpi = img_res)

# Plot mesh with observations:
mesh_df = data.frame(Lon = my_mesh$loc[,1], Lat = my_mesh$loc[,2]) %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
plot_dat = sp_data %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
p1 = ggplot() + 
  geom_sf(data = plot_dat, color = 'gray80', fill = 'gray80', shape = 21, size = 0.5) +
  geom_sf(data = mesh_df, color = 'black', fill = 'black', shape = 21, size = 0.5) 
p1 = add_sf_map(p1)
ggsave(paste0('map_mesh_obs', img_type), path = plot_folder, 
       plot = p1, width = img_width*0.75, height = 120, units = 'mm', dpi = img_res)

# SEM model:
sem_mod = "
  f1 -> 1, o1
  f1 -> 2, o2
  f1 -> 3, o3
  f1 -> 4, o4
  f1 -> 5, o5
  f2 -> 2, o6
  f2 -> 3, o7
  f2 -> 4, o8
  f2 -> 5, o9
  f1 <-> f1, sd_o_f1
  f2 <-> f2, sd_o_f2
  1 <-> 1, NA, 0
  2 <-> 2, NA, 0
  3 <-> 3, NA, 0
  4 <-> 4, NA, 0
  5 <-> 5, NA, 0
"

# DSEM model
dsem_mod = "
  f1 -> 1, 0, e1
  f1 -> 2, 0, e2
  f1 -> 3, 0, e3
  f1 -> 4, 0, e4
  f1 -> 5, 0, e5
  f2 -> 2, 0, e6
  f2 -> 3, 0, e7
  f2 -> 4, 0, e8
  f2 -> 5, 0, e9
  f1 -> f1, 1, NA, 0
  f2 -> f2, 1, NA, 0
  1 -> 1, 1, NA, 0
  2 -> 2, 1, NA, 0
  3 -> 3, 1, NA, 0
  4 -> 4, 1, NA, 0
  5 -> 5, 1, NA, 0
  f1 <-> f1, 0, sd_e_f1
  f2 <-> f2, 0, sd_e_f2
  1 <-> 1, 0, NA, 0
  2 <-> 2, 0, NA, 0
  3 <-> 3, 0, NA, 0
  4 <-> 4, 0, NA, 0
  5 <-> 5, 0, NA, 0
"

# Fit tinyVAST model
tinyVModel_j = tinyVAST( data = sp_data,
                    formula = Catch ~ 0 + factor(Year):factor(sp_number),
                    delta_options = list(
                         delta_formula = ~ 0 + factor(Year):factor(sp_number),
                         delta_sem = sem_mod,
                         delta_dsem = dsem_mod),
                    variables = c( "f1", "f2", 1:n_sp ),
                    variable_column = 'sp_number',
                    space_columns = c("Lon", "Lat"),
                    time_column = 'Year',
                    family = delta_lognormal(type="poisson-link"),
                    spatial_graph = my_mesh,
                    control = tinyVASTcontrol(gmrf="proj"))

tinyVModel_j$run_time
tinyVModel_j$sdrep
tinyVModel_j$internal$delta_sem_ram$output$model
tinyVModel_j$internal$delta_sem_ram$output$ram
tinyVModel_j$internal$delta_dsem_ram$output$model
tinyVModel_j$internal$delta_dsem_ram$output$ram

dim(tinyVModel_j$internal$parlist$omega2_sc)
tinyVModel_j$internal$parlist$theta2_z
tinyVModel_j$internal$parlist$omega2_sc[1:10, ]

dim(tinyVModel_j$internal$parlist$epsilon2_stc)
tinyVModel_j$internal$parlist$beta2_z
tinyVModel_j$internal$parlist$epsilon2_stc[1:10, , ]


save(tinyVModel_j, file = file.path(model_folder, 'fit.RData'))


# Check residuals
y_ir = replicate( n = 100, expr = tinyVModel_j$obj$simulate()$y_i )
res = DHARMa::createDHARMa( simulatedResponse = y_ir, 
                            observedResponse = sp_data$Catch, 
                            fittedPredictedResponse = fitted(tinyVModel_j) )
png(file = file.path(plot_folder, 'resid_check.png'), width = img_width, 
    height = 100, res = img_res, units = 'mm')
plot(res)
dev.off()

# -------------------------------------------------------------------------
# Get index:
Index = NULL
# To save annual estimates by species:
plot_dir_sp = file.path(plot_folder, 'Bycatch_est_ts_sp')
dir.create(plot_dir_sp, showWarnings = FALSE)
for(i in 1:n_sp) {
  this_sp = cumsp_data$sp_name[i]
  tmp_pred = extraRegion_tinyVAST
  tmp_pred$sp_number = i
  Est = sapply( sort(unique(sp_data$Year)), FUN=\(t) {
    pred_data = subset(tmp_pred, Year == t)
    integrate_output(tinyVModel_j, newdata = pred_data, area = pred_data$n_sets) 
  })
  tmp = data.frame( Year = sort(unique(sp_data$Year)), t(Est))
  colnames(tmp) = c('Year', 'est', 'se', 'est_corr', 'se_corr')
  # Calculate CI:
  tmp$est = tmp$est_corr # replace uncorrected values by corrected?
  tmp$lwr = tmp[,'est'] - 1.96*tmp[,'se']
  tmp$upr = tmp[,'est'] + 1.96*tmp[,'se']
  tmp$species = this_sp
  Index = rbind(Index, tmp)
  
  # Plot annual estimates with observed estimate:
  obs_df = read.csv(file = file.path('data', this_type, 'bycatch_est_obs', paste0(this_sp, '.csv')))
  p1 = plot_time_predictions(tmp, obs_df, Year, est, lwr, upr, yLab = 'Estimated bycatch (tons)')
  ggsave(filename = paste0(this_sp, img_type), path = plot_dir_sp, plot = p1, 
         width = img_width*0.5, height = 70, units = 'mm', dpi = img_res)
  
  cat('Species:', this_sp, 'done \n')
}
Index$model = 'tinyVAST'
Index$type = 'joint'
write.csv(Index, file = file.path(plot_folder, 'Bycatch_est_ts.csv'), row.names = FALSE)

# -------------------------------------------------------------------------
# Plot Omega 2nd component (Factors):
plot_dat = NULL
for(i in 1:n_fac) {
  tmp = data.frame(Lon = tinyVModel_j$spatial_graph$loc[,1], 
                   Lat = tinyVModel_j$spatial_graph$loc[,2], 
                   omega2 = tinyVModel_j$internal$parlist$omega2_sc[,i],
                   type_fac = paste0('factor_', i))
  plot_dat = rbind(plot_dat, tmp)
}
plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
p1 = ggplot(plot_dat) + geom_sf(aes(color = omega2), size = 2) + scale_colour_gradient2() + labs(color = 'Omega2')
p1 = add_sf_map(p1)
p1 = p1 + facet_wrap(~ type_fac)
ggsave(filename = paste0('omega2_fac', img_type), path = plot_folder, plot = p1, 
       width = img_width, height = 90, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Plot Epsilon 2nd component (Factors):
all_years = sort(unique(sp_data$Year))
for(j in 1:n_fac) {
  plot_dat = NULL
  for(i in seq_along(all_years)) {
    tmp_df = data.frame(year = all_years[i],
                        Lon = tinyVModel_j$spatial_graph$loc[,1], 
                        Lat = tinyVModel_j$spatial_graph$loc[,2], 
                        epsilon2 = tinyVModel_j$internal$parlist$epsilon2_stc[,i,j])
    plot_dat = rbind(plot_dat, tmp_df)
  }
  plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  p1 = ggplot(plot_dat) + geom_sf(aes(color = epsilon2), size = 1) + scale_colour_gradient2() + labs(color = 'Epsilon2')
  p1 = add_sf_map(p1)
  p1 = p1 + facet_wrap(~ year)
  ggsave(filename = paste0('epsilon2_Factor', j, img_type), path = plot_folder, plot = p1, 
         width = img_width, height = 130, units = 'mm', dpi = img_res)
}


# -------------------------------------------------------------------------
# Plot Omega 2nd component (Species):
# plot_dat = NULL
# for(i in 1:n_sp) {
#   tmp = data.frame(Lon = tinyVModel_j$spatial_graph$loc[,1], 
#                    Lat = tinyVModel_j$spatial_graph$loc[,2], 
#                    omega2 = tinyVModel_j$internal$parlist$omega2_sc[,n_fac + i],
#                    type_fac = cumsp_data$sp_name[i])
#   plot_dat = rbind(plot_dat, tmp)
# }
# plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
# p1 = ggplot(plot_dat) + geom_sf(aes(color = omega2), size = 2) + scale_colour_gradient2() + labs(color = 'Omega2')
# p1 = add_sf_map(p1)
# p1 = p1 + facet_wrap(~ type_fac)
# ggsave(filename = paste0('omega2_sp', img_type), path = plot_folder, plot = p1, 
#        width = img_width, height = 90, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Plot predicted values 
pred_df = as.data.frame(extraRegion_tinyVAST)
for(i in 1:n_sp) {
  tmp_pred = pred_df
  this_sp = cumsp_data$sp_name[i]
  tmp_pred$sp_number = i
  tmp_pred$mu_g = predict(tinyVModel_j, newdata = tmp_pred, what = "mu_g")

  plot_dat = tmp_pred %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  p1 = ggplot(plot_dat) + geom_sf(aes(color = mu_g), size = 0.5) + scale_colour_viridis() + labs(color = 'mu_g')
  p1 = add_sf_map(p1)
  p1 = p1 + facet_wrap(~ Year)
  ggsave(filename = paste0('mu_g_', this_sp, img_type), path = plot_folder, plot = p1, 
         width = img_width, height = 130, units = 'mm', dpi = img_res)
  
}

# -------------------------------------------------------------------------
# Association sp to factors (Omega2):
Lhat_cf = matrix(0, nrow = n_sp, ncol = n_fac)
theta2_vec = as.list(tinyVModel_j$sdrep, what="Estimate")$theta2_z[1:(n_fac*n_sp - sum(0:(n_fac-1)))]
Lhat_cf[lower.tri(Lhat_cf, diag=TRUE)] = theta2_vec
Lhat_cf = rotate_pca( L_tf = Lhat_cf, order = "decreasing" )$L_tf
dimnames(Lhat_cf) = list( cumsp_data$sp_name[1:n_sp],
                          paste0("Factor ", 1:ncol(Lhat_cf)) )
load_df = as.data.frame(Lhat_cf)
load_df$species = rownames(load_df)
load_df = pivot_longer(load_df, cols = starts_with("Factor"), names_to = 'Factor', values_to = 'loading')

p1 = ggplot(data = load_df) + geom_col(aes(x = species, y = loading), position = 'dodge') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_rect(fill="white")) +
  labs(y = 'Factor loading', x = NULL) +
  facet_wrap(~ Factor, scales = 'free_y')
ggsave(filename = paste0('Loading_Omega2', img_type), path = plot_folder, plot = p1, 
       width = img_width, height = 130, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Association sp to factors (Epsilon2):
Lhat_cf = matrix(0, nrow = n_sp, ncol = n_fac)
beta2_vec = as.list(tinyVModel_j$sdrep, what="Estimate")$beta2_z[1:(n_fac*n_sp - sum(0:(n_fac-1)))]
Lhat_cf[lower.tri(Lhat_cf, diag=TRUE)] = beta2_vec
Lhat_cf = rotate_pca( L_tf = Lhat_cf, order = "decreasing" )$L_tf
dimnames(Lhat_cf) = list( cumsp_data$sp_name[1:n_sp],
                          paste0("Factor ", 1:ncol(Lhat_cf)) )
load_df = as.data.frame(Lhat_cf)
load_df$species = rownames(load_df)
load_df = pivot_longer(load_df, cols = starts_with("Factor"), names_to = 'Factor', values_to = 'loading')

p1 = ggplot(data = load_df) + geom_col(aes(x = species, y = loading), position = 'dodge') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_rect(fill="white")) +
  labs(y = 'Factor loading', x = NULL) +
  facet_wrap(~ Factor, scales = 'free_y')
ggsave(filename = paste0('Loading_Epsilon2', img_type), path = plot_folder, plot = p1, 
       width = img_width, height = 130, units = 'mm', dpi = img_res)
