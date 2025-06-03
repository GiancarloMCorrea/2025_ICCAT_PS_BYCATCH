rm(list = ls())

# Define type of school to be analyzed:
this_type = 'FOB'
source('code/2_multi/load_libs.R')
mesh_cutoff = 2.5
n_fac_vec = 2:5

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path(data_folder, 'weight_data.rds'))
selsp_data = readRDS(file.path(data_folder, 'model_cat_sp.rds'))
n_sp = 15 # nrow(selsp_data)

# -------------------------------------------------------------------------
# Filter species:
mod_data = weight_data
# Filter to try model:
mod_data = mod_data %>% filter(sp_name %in% selsp_data$sp_name[1:n_sp])
# Define sp ID:
mod_data = mod_data %>% mutate(sp_number = factor(sp_name, levels = selsp_data$sp_name[1:n_sp], labels = 1:n_sp)) %>% 
              mutate(sp_number = as.integer(as.character(sp_number)))
# Factor some columns:
mod_data = mod_data %>% mutate(fsp = factor(sp_number), fyear = factor(year))
mod_data = as.data.frame(mod_data)

# Make mesh
my_mesh = fm_mesh_2d( mod_data[,c('lon', 'lat')], cutoff = mesh_cutoff )
p1 = ggplot() + geom_fm(data = my_mesh)
ggsave(paste0('map_mesh', img_type), path = plot_folder, 
       plot = p1, width = img_width*0.5, height = 80, units = 'mm', dpi = img_res)

# Plot mesh with observations:
mesh_df = data.frame(lon = my_mesh$loc[,1], lat = my_mesh$loc[,2]) %>% 
            st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
plot_dat = mod_data %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
p1 = ggplot() + 
  geom_sf(data = plot_dat, color = 'gray80', fill = 'gray80', shape = 21, size = 0.5) +
  geom_sf(data = mesh_df, color = 'black', fill = 'black', shape = 21, size = 0.5) 
p1 = add_sf_map(p1)
ggsave(paste0('map_mesh_obs', img_type), path = plot_folder, 
       plot = p1, width = img_width*0.75, height = 120, units = 'mm', dpi = img_res)
