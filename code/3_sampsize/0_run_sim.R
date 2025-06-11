rm(list = ls())

# IMPORTANT:
# Note that you need to run single models first (part 1 of this repo)

# Define type of school to be analyzed:
source('code/3_sampsize/load_libs.R')

# Read effort data:
effPoints = readRDS(file = file.path('data/1_single', this_type, 'effPoints.rds')) # for prediction
effPoints$id_obs = 1:nrow(effPoints) # to do matching later

# Specify number of simulations:
nSims = 100

# Specify sampling fraction vector:
frac_vector = seq(from = 0.1, to = 1, by = 0.1)

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
# Loop over sampling effort vector:
for(j in 1:nSims) {
  
  for(i in 1:length(frac_vector)) {
  
    obs_data = effPoints %>% group_by(year, quarter) %>% slice_sample(prop = frac_vector[i]) 
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
      
      # Now calculate bycatch per year for each sim:
      obs_data$bycatch = sim_matrix[obs_data$id_obs, j]
      obs_data$sp_name = these_sp[isp]
      if(j == 1) {
        obs_data$samp_frac = frac_vector[i]
        saveRDS(obs_data, file = file.path(model_folder, 'obs_sim_data', frac_vector[i], 
                                           paste0("sp_", isp, ".rds")))
      }
      # Continue estimating bycatch:
      est_df = calculate_ratio_bycatch(obs_df = obs_data, eff_df = effPoints, type = 'production')
      est_df = est_df %>% group_by(year, sp_name) %>% summarise(est = sum(est), .groups = 'drop')
      est_df$sim = paste0("sim_", j)
      est_df$samp_frac = frac_vector[i]
      est_df = left_join(est_df, true_data, by = c("year"))
      save_sim[[isp]] = est_df
      
    } # Loop over species
    
    # Save estimates:
    save_sim = bind_rows(save_sim)
    saveRDS(save_sim, file = file.path(model_folder, 'sim_est', frac_vector[i], 
                                        paste0("sim_", j, ".rds")))
      
  } # Loop sampling vector
  
  cat("Simulation", j, "done!\n")

} # Loop over sims
