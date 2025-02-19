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
sp_data = sp_data %>% mutate(sp_id = paste0('sp', sp_number),
                             sp_number = as.integer(as.character(sp_number)))
sp_data = as.data.frame(sp_data)

# Make mesh
my_mesh = fm_mesh_2d( sp_data[,c('Lon', 'Lat')], cutoff = 2 )
plot(sp_data$Lon, sp_data$Lat, col='red', pch = '.')
points(my_mesh$loc[,1], my_mesh$loc[,2], pch = 19, cex = 0.5)

# SEM model:
sem = "
  f1 -> 1, l1
  f1 -> 2, l2
  f1 -> 3, l3
  f1 -> 4, l4
  f1 -> 5, l5
  f2 -> 2, l6
  f2 -> 3, l7
  f2 -> 4, l8
  f2 -> 5, l9
  f1 <-> f1, NA, 1
  f2 <-> f2, NA, 1
  1 <-> 1, NA, 0
  2 <-> 2, NA, 0
  3 <-> 3, NA, 0
  4 <-> 4, NA, 0
  5 <-> 5, NA, 0
"

# Fit tinyVAST model
tVModel = tinyVAST( data = sp_data,
                    formula = Catch ~ 0 + factor(Year),
                    delta_options = list(
                         delta_formula = ~ 0 + factor(Year),
                         delta_sem = sem),
                    variables = c( "f1", "f2", 1:n_sp ),
                    variable_column = 'sp_number',
                    space_columns = c("Lon", "Lat"),
                    time_column = 'Year',
                    family = delta_lognormal(type="poisson-link"),
                    spatial_graph = my_mesh,
                    control = tinyVASTcontrol(gmrf="proj"))
save(tVModel, file = file.path(model_folder, 'fit.RData'))


# Check residuals
y_ir = replicate( n = 100, expr = tVModel$obj$simulate()$y_i )
res = DHARMa::createDHARMa( simulatedResponse = y_ir, 
                            observedResponse = sp_data$Catch, 
                            fittedPredictedResponse = fitted(tVModel) )
plot(res)

# Get index:
Index = NULL
for(i in 1:n_sp) {
  tmp_pred = extraRegion_tinyVAST
  tmp_pred$sp_number = i
  Est = sapply( sort(unique(sp_data$Year)), FUN=\(t) {
    pred_data = subset(tmp_pred, Year == t)
    integrate_output(tVModel, newdata = pred_data, area = pred_data$n_sets) 
  })
  tmp = data.frame( Year = sort(unique(sp_data$Year)), t(Est))
  colnames(tmp) = c('Year', 'est', 'se', 'est_corr', 'se_corr')
  tmp$sp_name = cumsp_data$sp_name[i]
  tmp$sp_number = i
  Index = rbind(Index, tmp)
  cat('Species:', cumsp_data$sp_name[i], 'done \n')
}
# Calculate CI:
Index$lwr = Index[,'est_corr'] - 1.96*Index[,'se']
Index$upr = Index[,'est_corr'] + 1.96*Index[,'se']
Index$type = 'tinyVAST'
write.csv(Index, file = file.path(model_folder, 'Bycatch_est_ts.csv'), row.names = FALSE)

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
