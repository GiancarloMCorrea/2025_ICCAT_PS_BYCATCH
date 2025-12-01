rm(list = ls())

# -------------------------------------------------------------------------
# Define type of school to be analyzed:
source('code/3_sampsize/load_libs.R')
source('code/1_single/sdmTMB-config.R')

# Read data:
effPoints = readRDS(file = file.path('data/1_single', this_type, 'weight_data.rds'))
effPoints = effPoints %>% mutate(marea_id = paste(vessel_code, trip_start_date, sep = '_'))

# Species to analyze:
these_sp = unique(effPoints$sp_name)

# Specify number of replicates:
nSims = 4

# Mesh cutoff:
if(this_type == 'FOB') mesh_cutoff = 1.5
if(this_type == 'FSC') mesh_cutoff = 1

# N cores for parallel:
n_cores = 4

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
    obs_data = effPoints %>% mutate(cluster_id = if_else(marea_id %in% obs_marea$marea_id, 1, 2))
    # Obs_data will be used for all species
    
    # Create folder to save results:
    if(j == 1) dir.create(file.path(model_folder, 'crossval', frac_vector[i]), showWarnings = FALSE, recursive = TRUE)
    
    # Loop over species
    for(isp in seq_along(these_sp)) {
      
      # Filter sp:
      sp_data = obs_data %>% filter(sp_name == these_sp[isp], cluster_id == 1)
      
      # Effort sp:
      effSp = obs_data %>% filter(sp_name == these_sp[isp]) # predict in the entire dataset
      
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
        mod_data = mod_data %>% mutate(fyear = factor(year, levels = sort(unique(mod_data$year))))
        # Make mesh:
        sp_mesh = sdmTMB::make_mesh(data = mod_data, xy_cols = c('lon', 'lat'),
                                    mesh = fmesher::fm_mesh_2d( mod_data[,c('lon', 'lat')], 
                                                                cutoff = mesh_cutoff ))
        # Run model:
        sdmtmb_df = run_sdmTMB_model(sp_data = mod_data, this_formula = this_formula, 
                                     sp_mesh = sp_mesh, this_sp = these_sp[isp], 
                                     effPoints = effSp, yr_keep = yr_keep, qt_keep = qt_keep,
                                     save_model = FALSE, make_plots = FALSE, save_results = FALSE,
                                     check_residuals = FALSE)
        
        if(!is.null(sdmtmb_df[[1]])) { 
          # Merge predictions:
          pred_data = left_join(effSp, sdmtmb_df[[2]] %>% ungroup() %>% select(id_set, bycatch_est), by = "id_set")
          pred_data$bycatch_est[is.na(pred_data$bycatch_est)] = 0 # zeros for no predictions
          # Calculate RMSE and MAE:
          out_data = data.frame(rmse = sqrt(mean((pred_data$bycatch - pred_data$bycatch_est)^2)),
                                mae = mean(abs(pred_data$bycatch - pred_data$bycatch_est)) )
          out_data$category = mean(sdmtmb_df[[1]]$category)
        } else { # if model failed, then NA:
          pred_data = effSp
          pred_data$bycatch_est = NA # NA values when model failed
          # Calculate RMSE and MAE:
          out_data = data.frame(rmse = NA, mae = NA )
          out_data$category = 3
        }
        
      } else { # If missing data, then predictions = 0:
        
        pred_data = effSp
        pred_data$bycatch_est = 0
        # Calculate RMSE and MAE:
        out_data = data.frame(rmse = sqrt(mean((pred_data$bycatch - pred_data$bycatch_est)^2)),
                              mae = mean(abs(pred_data$bycatch - pred_data$bycatch_est)) )
        out_data$category = 4
        
      }
      
      # Add sim and sample_frac scenarios info:
      out_data$sim = paste0("sim_", j)
      out_data$samp_frac = frac_vector[i]
      out_data$sp_name = these_sp[isp]
      
      # Save results
      saveRDS(out_data, file = file.path(model_folder, 'crossval', frac_vector[i], 
                                         paste0("sp-", isp, "_sim-", j, ".rds")))

    } # Loop over species
    
  } # Loop sampling vector
  
} )# parallel loop over sims

# Stop cluster
snowfall::sfStop()
