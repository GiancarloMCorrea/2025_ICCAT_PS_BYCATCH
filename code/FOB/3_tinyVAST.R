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
this_type = 'FOB'

# Define folder to save results
plot_folder = file.path('figures', this_type, 'tinyVAST')
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Define model folder:
model_folder = file.path('models', this_type, 'tinyVAST')
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
this_sp = cumsp_data$sp_name[36] # change this index if want to make a loop
sp_data = mod_data %>% dplyr::filter(sp_name %in% this_sp) 
# Check sum catch by year:
sp_data %>% group_by(Year) %>% summarise(Catch = sum(Catch))
# Define sp ID:
sp_data = sp_data %>% mutate(sp_id = 'sp')
sp_data = as.data.frame(sp_data)

# Create dir to save plots and model outputs:
dir.create(file.path(model_folder, this_sp), showWarnings = FALSE)
dir.create(file.path(plot_folder, this_sp), showWarnings = FALSE)

# Make mesh
my_mesh = fm_mesh_2d( sp_data[,c('Lon', 'Lat')], cutoff = 2 )
plot(sp_data$Lon, sp_data$Lat, col='red', pch = '.')
points(my_mesh$loc[,1], my_mesh$loc[,2], pch = 19, cex = 0.5)

# Fit tinyVAST model
tVModel = tinyVAST( data = sp_data,
                    formula = Catch ~ 0 + factor(Year),
                    delta_options = list(
                         delta_formula = ~ 0 + factor(Year),
                         delta_sem = "sp <-> sp, spatial_sd",
                         delta_dsem = "sp -> sp, 1, rho"),
                    variable_column = 'sp_id',
                    space_columns = c("Lon", "Lat"),
                    time_column = 'Year',
                    family = delta_lognormal(type="poisson-link"),
                    spatial_graph = my_mesh )
save(tVModel, file = file.path(model_folder, this_sp, 'fit.RData'))


# Check residuals
y_ir = replicate( n = 100, expr = tVModel$obj$simulate()$y_i )
res = DHARMa::createDHARMa( simulatedResponse = y_ir, 
                            observedResponse = sp_data$Catch, 
                            fittedPredictedResponse = fitted(tVModel) )
plot(res)

# Get index:
Est = sapply( sort(unique(sp_data$Year)), FUN=\(t) {
  pred_data = subset(as.data.frame(extraRegion_tinyVAST), Year == t)
  pred_data$sp_id = 'sp'
  integrate_output(tVModel, newdata = pred_data, area = pred_data$n_sets) 
  })
Index = data.frame( Year = sort(unique(sp_data$Year)), t(Est))
colnames(Index) = c('Year', 'est', 'se', 'est_corr', 'se_corr')
Index$lwr = Index[,'est_corr'] - 1.96*Index[,'se']
Index$upr = Index[,'est_corr'] + 1.96*Index[,'se']
Index$type = 'tinyVAST'
write.csv(Index, file = file.path(model_folder, this_sp, 'Bycatch_est_ts.csv'), row.names = FALSE)

# Ratio estimator: (calculate mean catch per set per year, and then multiply by total effort per year)
eff_yr = extraRegion_tinyVAST %>% group_by(Year, ID) %>% summarise(n_sets = sum(n_sets))
obs_yr = sp_data %>% group_by(Year, ID) %>% summarise(avgCatch = mean(Catch))
obs_df = left_join(obs_yr, eff_yr)
# There may be some NA due to ID present in obs data (different flags) but absent in effData (only SPA)
obs_df = obs_df %>% mutate(est = avgCatch*n_sets) %>% group_by(Year) %>% 
        summarise(est = sum(est, na.rm = TRUE))
p1 = plot_time_predictions(Index, obs_df, Year, est, lwr, upr, yLab = paste0('Est (', this_sp, ')'))
ggsave(filename = paste0('Bycatch_est_ts', img_type), path = file.path(plot_folder, this_sp), plot = p1, 
       width = img_width*0.5, height = 70, units = 'mm', dpi = img_res)
