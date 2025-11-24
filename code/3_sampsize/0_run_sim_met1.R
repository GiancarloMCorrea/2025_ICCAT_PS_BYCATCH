rm(list = ls())

# -------------------------------------------------------------------------

# Define type of school to be analyzed:
source('code/3_sampsize/load_libs.R')
source('code/1_single/sdmTMB-config.R') # for model-based estimates

# Read 'effort' data:
# This is the observers data that will be treated as 100% 
effPoints = readRDS(file = file.path('data/1_single', this_type, 'weight_data.rds'))
effPoints = effPoints %>% mutate(marea_id = paste(vessel_code, trip_start_date, sep = '_'))

# Species to analyze:
these_sp = unique(effPoints$sp_name)

# Specify number of simulations:
nSims = 100

# Mesh cutoff:
if(this_type == 'FOB') mesh_cutoff = 1.5
if(this_type == 'FSC') mesh_cutoff = 1

# N cores for parallel:
n_cores = 4

# -------------------------------------------------------------------------
# Estimate true values:
true_data = effPoints %>% group_by(year, sp_name) %>% summarise(true = sum(bycatch), .groups = 'drop')

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
      
      # Filter sp:
      sp_data = obs_data %>% filter(sp_name == these_sp[isp])
      
      # Effort sp:
      effSp = effPoints %>% filter(sp_name == these_sp[isp])
      
      # Save data:
      if(j == 1) {
        sp_data$samp_frac = frac_vector[i]
        saveRDS(sp_data, file = file.path(model_folder, 'obs_sim_data', frac_vector[i], 
                                          paste0("sp_", isp, ".rds")))
      }
      
      # ----------------------------------------
      # Ratio estimator:
      ratio_df = calculate_ratio_bycatch(obs_df = sp_data, eff_df = effSp, type = 'production')
      ratio_df = ratio_df %>% group_by(year, sp_name) %>% summarise(est_ratio = sum(est), .groups = 'drop')
      estimate_df = left_join(true_data %>% filter(sp_name == these_sp[isp]), ratio_df, by = c("year", "sp_name"))
      # Add information:
      estimate_df$sim = paste0("sim_", j)
      estimate_df$samp_frac = frac_vector[i]
      
      # ----------------------------------------
      # Model-based estimator:
      this_formula = list(as.formula(bycatch ~ 0 + fyear + quarter + trop_catch))
      # Remove years with no presence (if any):
      yr_pa = sp_data %>% group_by(year) %>% summarise(bycatch = sum(bycatch))
      yr_keep = yr_pa$year[which(yr_pa$bycatch > 0)]
      mod_data = sp_data %>% filter(year %in% yr_keep)
      if(length(yr_keep) == 1) {
        this_formula = list(update(this_formula[[1]], ~ . - fyear))
      }
      # Remove quarters with no presence (if any):
      qt_pa = sp_data %>% group_by(quarter) %>% summarise(bycatch = sum(bycatch))
      qt_keep = qt_pa$quarter[which(qt_pa$bycatch > 0)]
      mod_data = mod_data %>% filter(quarter %in% qt_keep)
      if(length(qt_keep) == 1) {
        this_formula = list(update(this_formula[[1]], ~ . - quarter))
      }
      if(nrow(mod_data) > 0) { # if there is still data, then run model:
        
        # Add fyear:
        mod_data = mod_data %>% mutate(fyear = factor(year, levels = sort(unique(sp_data$year))))
        # Make mesh:
        sp_mesh = sdmTMB::make_mesh(data = mod_data, xy_cols = c('lon', 'lat'),
                                    mesh = fmesher::fm_mesh_2d( sp_data[,c('lon', 'lat')], 
                                                                cutoff = mesh_cutoff ))
        # Run model:
        sdmtmb_df = run_sdmTMB_model(sp_data = mod_data, this_formula = this_formula, 
                                     sp_mesh = sp_mesh, this_sp = these_sp[isp], 
                                     effPoints = effSp, yr_keep = yr_keep, qt_keep = qt_keep,
                                     save_model = FALSE, make_plots = FALSE, save_results = FALSE,
                                     check_residuals = FALSE)
        # Merge with results df:
        if(!is.null(sdmtmb_df)) { 
          sdmtmb_df = sdmtmb_df %>% rename(est_model = est)
          estimate_df = left_join(estimate_df, 
                                  sdmtmb_df %>% select(year, est_model, category), 
                                  by = c("year"))
          estimate_df$est_model[is.na(estimate_df$est_model)] = 0 # fill NA with zeros
          estimate_df$category = mean(estimate_df$category, na.rm = TRUE) # to fill NA when missing years
        } else { # if model failed, then NA:
          estimate_df$est_model = NA
          estimate_df$category = 3 # model failed identifier
        }
        
      } else { # If missing data, then estimate = 0:
        
        estimate_df$est_model = 0
        estimate_df$category = 4 # no data identifier
        
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
