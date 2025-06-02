rm(list = ls())

# Define type of school to be analyzed:
this_type = 'FOB' # or FSC
source('code/1_single/load_libs.R')
mesh_cutoff = 1.5

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path(data_folder, 'weight_data.rds'))
effPoints = readRDS(file = file.path(data_folder, 'effPoints.rds')) # for prediction
load(file.path(data_folder, 'MyGrid.RData'))

# -------------------------------------------------------------------------
# Filter species:
mod_data = weight_data
selsp_data = readRDS(file.path(data_folder, 'model_cat_sp.rds'))
# Tmp data, to make mesh. the selection of sp does not matter since all sp have the same obs locations
tmp_data = mod_data %>% dplyr::filter(sp_name %in% selsp_data$sp_name[1])
n_mod_sp = nrow(selsp_data)

# -------------------------------------------------------------------------
# Make mesh (will be the same for all sp):
prev_mesh = fm_mesh_2d( tmp_data[,c('lon', 'lat')], cutoff = mesh_cutoff )
p1 = ggplot() + geom_fm(data = prev_mesh)
ggsave(paste0('map_mesh', img_type), path = file.path(plot_folder), 
       plot = p1, width = img_width*0.5, height = 80, units = 'mm', dpi = img_res)

# Plot mesh with observations:
mesh_df = data.frame(lon = prev_mesh$loc[,1], lat = prev_mesh$loc[,2]) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
plot_dat = tmp_data %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
p1 = ggplot() + 
  geom_sf(data = plot_dat, color = 'gray80', fill = 'gray80', shape = 21, size = 0.5) +
  geom_sf(data = mesh_df, color = 'black', fill = 'black', shape = 21, size = 0.5) 
p1 = add_sf_map(p1)
ggsave(paste0('map_mesh_obs', img_type), path = plot_folder, 
       plot = p1, width = img_width*0.75, height = 120, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Specify full formula:
full_formula = list(as.formula(bycatch ~ 0 + fyear + quarter + trop_catch + sst),
                    as.formula(bycatch ~ 0 + fyear + quarter + trop_catch + sst))
