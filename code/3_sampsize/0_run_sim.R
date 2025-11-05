rm(list = ls())

# -------------------------------------------------------------------------
# IMPORTANT:
# Note that you need to run single models first (part 1 of this repo)

# Define type of school to be analyzed:
source('code/3_sampsize/load_libs.R')
source('code/1_single/sdmTMB-config.R') # for model-based estimates

# Read effort data:
effPoints = readRDS(file = file.path('data/1_single', this_type, 'effPoints.rds')) # for prediction
effPoints$id_obs = 1:nrow(effPoints) # to do matching later

# Specify number of simulations:
nSims = 100

# Specify sampling fraction vector:
frac_vector = c(0.05, seq(from = 0.1, to = 0.5, by = 0.1), 0.7, 0.9)

# Mesh cutoff:
mesh_cutoff = 1.5 # same as in part 1

# N cores for parallel:
n_cores = 4

# -------------------------------------------------------------------------

# Find species with fitted model:
all_sp = list.files(file.path("model/1_single", this_type))
these_sp = numeric(0)
for(k in seq_along(all_sp)) {
  this_file = file.path("model/1_single", this_type, all_sp[k], 'mod.RData')
  if(file.exists(this_file)) these_sp = c(these_sp, all_sp[k])
}

# Make seeds for simulation in sdmTMB:
set.seed(8675309)
seeds_sim = sample(x = (-1e9):(1e9), size = nSims, replace = FALSE)


# -------------------------------------------------------------------------
# Simulate values for each species using all data:
save_sim_matrix = list()
for(isp in seq_along(these_sp)) {
  rm(mod_upd)
  tot_matrix = matrix(0, nrow = nrow(effPoints), ncol = nSims)
  load(file.path("model/1_single", this_type, these_sp[isp], 'mod.RData'))
  # Select data to simulate:
  mod_dat = mod_upd$data
  yr_keep = unique(mod_dat$year)
  qt_keep = unique(mod_dat$quarter)
  # Filter observed data for simulation
  sim_data = effPoints %>% filter(year %in% yr_keep, quarter %in% qt_keep)
  sim_data = sim_data %>% mutate(fyear = factor(year, levels = sort(unique(sim_data$year))))
  # Simulate:
  sim_matrix = simulate(mod_upd, nsim = nSims, seed = seeds_sim, newdata = sim_data)
  tot_matrix[sim_data$id_obs, ] = sim_matrix 
  save_sim_matrix[[isp]] = tot_matrix
  cat("Simulation for species", isp, "done!\n")
}
# Save simulated values:
save(save_sim_matrix, file = file.path(model_folder, 'sim_sp_matrix.RData'))

# -------------------------------------------------------------------------
# Start running simulation in parallel:

snowfall::sfInit(parallel=TRUE, cpus=n_cores)
snowfall::sfExportAll()
trash = snowfall::sfLapply(1:nSims, function(j) {
  
  # Load required libraries:
  require(dplyr)
  require(ggplot2)
  require(sdmTMB)
  require(sf)
  require(fmesher)
  
  for(i in 1:length(frac_vector)) {
  
    obs_marea = effPoints %>% group_by(year, marea_id) %>% summarise() %>% group_by(year) %>% slice_sample(prop = frac_vector[i]) 
    obs_data = effPoints %>% filter(marea_id %in% obs_marea$marea_id)
    # Obs_data will be used for all species
    
    # Create folder to save sim estimates 
    if(j == 1) { 
      dir.create(file.path(model_folder, 'sim_est', frac_vector[i]), showWarnings = FALSE, recursive = TRUE)
      dir.create(file.path(model_folder, 'obs_sim_data', frac_vector[i]), showWarnings = FALSE, recursive = TRUE)
    }
    
    # Loop over species
    save_sim = list() # to save estimates
    for(isp in seq_along(these_sp)) {
      
      # Read simulated values:
      sim_matrix = save_sim_matrix[[isp]]
      
      # Estimate true values:
      true_data = data.frame(year = effPoints$year, true = sim_matrix[, j])
      true_data = true_data %>% group_by(year) %>% summarise(true = sum(true), .groups = 'drop')
      true_data$sim = paste0("sim_", j)
      true_data$samp_frac = frac_vector[i]
      
      # Now calculate bycatch per year for each sim:
      obs_data$bycatch = sim_matrix[obs_data$id_obs, j]
      obs_data$sp_name = these_sp[isp]
      if(j == 1) {
        obs_data$samp_frac = frac_vector[i]
        saveRDS(obs_data, file = file.path(model_folder, 'obs_sim_data', frac_vector[i], 
                                           paste0("sp_", isp, ".rds")))
      }
      
      # ----------------------------------------
      # Ratio estimator:
      ratio_df = calculate_ratio_bycatch(obs_df = obs_data, eff_df = effPoints, type = 'production')
      ratio_df = ratio_df %>% group_by(year, sp_name) %>% summarise(est_ratio = sum(est), .groups = 'drop')
      estimate_df = left_join(true_data, ratio_df, by = c("year"))
      
      # ----------------------------------------
      # Model-based estimator:
      this_formula = list(as.formula(bycatch ~ 0 + fyear + quarter + trop_catch))
      # Remove years with no presence (if any):
      yr_pa = obs_data %>% group_by(year) %>% summarise(bycatch = sum(bycatch))
      yr_keep = yr_pa$year[which(yr_pa$bycatch > 0)]
      sp_data = obs_data %>% filter(year %in% yr_keep)
      if(length(yr_keep) == 1) {
        this_formula = list(update(this_formula[[1]], ~ . - fyear))
      }
      # Remove quarters with no presence (if any):
      qt_pa = obs_data %>% group_by(quarter) %>% summarise(bycatch = sum(bycatch))
      qt_keep = qt_pa$quarter[which(qt_pa$bycatch > 0)]
      sp_data = sp_data %>% filter(quarter %in% qt_keep)
      if(length(qt_keep) == 1) {
        this_formula = list(update(this_formula[[1]], ~ . - quarter))
      }
      # Add fyear:
      sp_data = sp_data %>% mutate(fyear = factor(year, levels = sort(unique(sp_data$year))))
      # Make mesh:
      sp_mesh = sdmTMB::make_mesh(data = sp_data, xy_cols = c('lon', 'lat'),
                                  mesh = fmesher::fm_mesh_2d( sp_data[,c('lon', 'lat')], 
                                                              cutoff = mesh_cutoff ))
      # Run model:
      sdmtmb_df = run_sdmTMB_model(sp_data = sp_data, this_formula = this_formula, 
                                  sp_mesh = sp_mesh, this_sp = these_sp[isp], 
                                  effPoints = effPoints, yr_keep = yr_keep, qt_keep = qt_keep,
                                  save_model = FALSE, make_plots = FALSE, save_results = FALSE,
                                  check_residuals = FALSE)
      # Merge with results df:
      if(!is.null(sdmtmb_df)) { 
        sdmtmb_df = sdmtmb_df %>% rename(est_model = est)
        estimate_df = left_join(estimate_df, 
                                sdmtmb_df %>% select(year, est_model, category), 
                                by = c("year"))
      } else {
        estimate_df$est_model = NA
        estimate_df$category = NA
      }
      
      # Save results
      save_sim[[isp]] = estimate_df
      
    } # Loop over species
    
    # Save estimates:
    save_sim = bind_rows(save_sim)
    saveRDS(save_sim, file = file.path(model_folder, 'sim_est', frac_vector[i], 
                                        paste0("sim_", j, ".rds")))
      
  } # Loop sampling vector
  
} )# parallel loop over sims

# Stop cluster
snowfall::sfStop()
