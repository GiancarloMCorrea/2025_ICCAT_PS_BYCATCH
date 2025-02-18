rm(list = ls())
require(ggplot2)
require(sf)
require(dplyr)
require(tidyr)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Define type of school to be analyzed:
this_type = 'FOB'

# Define data folder to read data in:
data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_PS_Bycatch'
# Read data in:
load(file.path(data_folder, 'data.weka.ATL.2003_2023_MarineBeacon_enviado_GC_110225.RData'))
# Read all sets data:
allsets_df = read.csv(file.path(data_folder, 'ATL_13a23_setsfleet.csv'))

# -------------------------------------------------------------------------
# Define school type:
save_data_folder = file.path('data', this_type)
dir.create(save_data_folder, recursive = TRUE, showWarnings = FALSE)
plot_dir = file.path('figures', this_type)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Read data ---------------------------------------------------------------
# Explore data:
glimpse(data_weka)

# Define columns to be used:
sp_id_col = grep(pattern = 'tons_|number_|weight_', x = colnames(data_weka))
var_columns = c('ocean_code', 'vessel_code', 'flag_country', 'trip_start_date', 'trip_end_date',
               'observation_date', 'observation_time', 'latitude', 'longitude', 'sst', 
               'school_type', 'sunrise_diference')
var_id_col = which(colnames(data_weka) %in% var_columns)
my_data = data_weka %>% select(all_of(c(var_id_col, sp_id_col)))

# Do some IMPORTANT filtering:
my_data = my_data %>% dplyr::filter(ocean_code == '1', # 1 = ATL
                                    school_type == this_type) # defined above

# Continue..
my_data = my_data %>% mutate(id_set = 1:n(),
                             year = as.integer(format(observation_date, '%Y')),
                             month = format(observation_date, '%m'), .before = 'ocean_code')
glimpse(my_data)
summary(my_data[,var_columns])

# Make plot year and month:
p1 = ggplot(my_data, aes(x = factor(year), fill = month)) + 
  geom_bar() + 
  scale_fill_viridis_d() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = 'bottom') +
  xlab(NULL) + ylab('N obs') +
  guides(fill = guide_legend(title = 'Month', nrow=2))
ggsave(filename = 'obs_num_y-m.png', path = plot_dir, plot = p1, 
       width = img_width, height = 130, units = 'mm', dpi = 300)

# Explore maps:
obsPoints = my_data %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Make plot:
p1 = ggplot(obsPoints) + geom_sf(size = 1, alpha = 0.5)
p1 = add_sf_map(p1)
p1 = p1 + facet_wrap(~ year)
ggsave(paste0('obs_map_sets', img_type), plot = p1, path = plot_dir,
       width = img_width, height = 190, units = 'mm', dpi = img_res) 

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Explore effort data:
glimpse(allsets_df)
eff_data = allsets_df
alt_set_type = this_type
if(this_type == 'FOB') alt_set_type = 'FAD'

# Do some IMPORTANT filtering:
eff_data = eff_data %>% dplyr::filter(oceano == 1, # 1 = ATL
                                      set_type == alt_set_type) # defined above
# Rename some variables:
eff_data = eff_data %>% dplyr::rename(longitude = lon, latitude = lat, year = año,
                                      tuna_catch = suma_de_YFT.BET.SKJ,
                                      no_sets = numero_de_lances)
# Create some columns:
eff_data = eff_data %>% mutate(month = as.integer(format(as.Date(fecha, format = '%Y-%m-%d'), '%m')),
                               .after = 'fecha')
# Remove rows with number of sets = 0
eff_data = eff_data %>% dplyr::filter(no_sets > 0)
# Repeat rows based on the number of sets:
eff_data = eff_data %>% tidyr::uncount(no_sets)

# Make plot year and month:
p1 = ggplot(eff_data, aes(x = factor(year), fill = factor(month))) + 
  geom_bar() + 
  scale_fill_viridis_d() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = 'bottom') +
  xlab(NULL) + ylab('N obs') +
  guides(fill = guide_legend(title = 'Month', nrow=2))
ggsave(filename = 'eff_num_y-m.png', path = plot_dir, plot = p1, 
       width = img_width, height = 130, units = 'mm', dpi = 300)

# Explore maps:
spt_eff_data = eff_data %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Make plot:
p1 = ggplot(spt_eff_data) + geom_sf(size = 1, alpha = 0.5)
p1 = add_sf_map(p1)
p1 = p1 + facet_wrap(~ year)
ggsave(paste0('eff_map_sets', img_type), plot = p1, path = plot_dir,
       width = img_width, height = 150, units = 'mm', dpi = img_res) 

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# Create prediction grid:
# Use effort data:
grid_size = 1 # in degrees

# Point and grid (aggregated):
MyPoints = eff_data %>% mutate(id_set = 1:n()) %>%
            st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
min_lon = floor(min(MyPoints$longitude))
min_lat = floor(min(MyPoints$latitude))
MyGrid = st_make_grid(MyPoints, cellsize = c(grid_size, grid_size), offset = c(min_lon, min_lat)) %>%  # define grid size here
  st_set_crs(4326) %>% st_sf() %>% dplyr::mutate(ID = 1:n())

# Join both Grid and Points:
# Do it in this way to avoid repeated rows when rounded lon/lat fall on 
# grid borders
MyPoints$ID = MyGrid$ID[st_nearest_feature(MyPoints, MyGrid)]

# Also attach ID info to observations:
obsPoints$ID = MyGrid$ID[st_nearest_feature(obsPoints, MyGrid)]
save(obsPoints, file = file.path(save_data_folder, 'obsPoints.RData'))

# Filter based on some criteria:
# Identify grid (and points inside) with some criteria:
my_tab = xtabs(~ ID + year, data = MyPoints) # find frequency of sets per grid and year
my_freq = apply(my_tab, 1, function(x) sum(x > 0)) # find recurrent grids over the years
include_grid = names(my_freq)[which(as.vector(my_freq) >= 0)] # IMPORTANT: include all!
include_grid = as.numeric(include_grid)
# Remove grids:
effPoints = MyPoints %>% dplyr::filter(ID %in% include_grid)
# Save:
save(effPoints, file = file.path(save_data_folder, 'effPoints.RData'))

# Apply filter to grid df and calculate Area km2:
MyGrid = MyGrid %>% dplyr::filter(ID %in% include_grid)
plot(MyGrid)
MyGrid$Area_km2 = as.numeric(st_area(MyGrid))*1e-06 # in km2
save(MyGrid, file = file.path(save_data_folder, 'MyGrid.RData'))

# Save extrapolation grid in VAST format:
user_region = st_centroid(MyGrid) %>% dplyr::mutate(Lon = sf::st_coordinates(.)[,1], Lat = sf::st_coordinates(.)[,2])
st_geometry(user_region) = NULL
saveRDS(user_region, file = file.path(save_data_folder, 'user_region.rds'))

# Calculate number of sets per grid and year:
nsets_df = effPoints %>% st_drop_geometry() %>% group_by(year, ID) %>% summarise(n_sets = n())
all_years = sort(unique(eff_data$year))
n_years = length(all_years)
MyGridSets = purrr::map_dfr(seq_len(n_years), ~user_region)
MyGridSets = MyGridSets %>% dplyr::mutate(year = rep(all_years, each = nrow(user_region)))
MyGridSets = left_join(MyGridSets, nsets_df)
MyGridSets$n_sets[which(is.na(MyGridSets$n_sets))] = 0 # no sets in NAs
saveRDS(MyGridSets, file = file.path(save_data_folder, 'MyGridSets.rds'))
