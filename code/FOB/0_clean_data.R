rm(list = ls())
require(ggplot2)
require(sf)
require(dplyr)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Define type of school to be analyzed:
this_type = 'FOB'

# Define data folder to read data in:
data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_PS_Bycatch'
# Read data in:
load(file.path(data_folder, 'data.weka.ATL_IND.2003_2023_v1.RData'))

# -------------------------------------------------------------------------
# Define school type:
save_data_folder = file.path('data', this_type)
dir.create(save_data_folder, recursive = TRUE, showWarnings = FALSE)
plot_dir = file.path('figures', this_type)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Read data ---------------------------------------------------------------
# Define columns to be used:
sp_id_col = grep(pattern = 'tons_|number_|weight_', x = colnames(data_weka))
var_columns = c('ocean_code', 'vessel_code', 'flag_country', 'trip_start_date', 'trip_end_date',
               'observation_date', 'observation_time', 'latitude', 'longitude', 'sst', 
               'school_type', 'sunrise.sunrise', 'sunrise.sunset', 'sunrise_diference')
var_id_col = which(colnames(data_weka) %in% var_columns)
my_data = data_weka %>% select(all_of(c(var_id_col, sp_id_col)))
my_data = my_data %>% dplyr::filter(ocean_code == '1') # 1 = ATL
my_data = my_data %>% mutate(id_set = 1:n(),
                             year = as.integer(format(observation_date, '%Y')),
                             month = format(observation_date, '%m'), .before = 'ocean_code')
my_data = my_data %>% dplyr::filter(school_type == this_type)
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
ggsave(filename = 'raw_nobs_y-m.png', path = plot_dir, plot = p1, 
       width = img_width, height = 130, units = 'mm', dpi = 300)

# Explore maps:
spt_data = my_data %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Make plot:
p1 = ggplot(spt_data) + geom_sf(size = 1, alpha = 0.5)
p1 = add_sf_map(p1)
p1 = p1 + facet_wrap(~ year)
ggsave(paste0('raw_obs_map', img_type), plot = p1, path = plot_dir,
       width = img_width, height = 190, units = 'mm', dpi = img_res) 

# -------------------------------------------------------------------------
# Create prediction grid:
grid_size = 1 # in degrees

# Point and grid (aggregated):
MyPoints = my_data %>% dplyr::filter(school_type == this_type) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
min_lon = floor(min(MyPoints$longitude))
min_lat = floor(min(MyPoints$latitude))
MyGrid = st_make_grid(MyPoints, cellsize = c(grid_size, grid_size), offset = c(min_lon, min_lat)) %>%  # define grid size here
  st_set_crs(4326) %>% st_sf() %>% dplyr::mutate(ID = 1:n())

# Join both Grid and Points:
joinDF = st_join(MyGrid, left = TRUE, MyPoints) %>% na.omit

# Filter based on some criteria:
# Identify grid (and points inside) with some criteria:
my_tab = xtabs(~ ID + year, data = joinDF) # find frequency of sets per grid and year
my_freq = apply(my_tab, 1, function(x) sum(x > 0)) # find recurrent grids over the years
include_grid = names(my_freq)[which(as.vector(my_freq) >= 3)] # grids to be excluded because infrequent sample (less than X years)
include_grid = as.numeric(include_grid)
# Remove grids:
joinDF = joinDF %>% dplyr::filter(ID %in% include_grid)
# Save:
save(joinDF, file = file.path(save_data_folder, 'joinDF.RData'))

# Apply filter to grid df:
MyGrid = MyGrid %>% dplyr::filter(ID %in% include_grid)
plot(MyGrid)
MyGrid$Area_km2 = as.numeric(st_area(MyGrid))*1e-06 # in km2
save(MyGrid, file = file.path(save_data_folder, 'MyGrid.RData'))

# Save extrapolation grid in VAST format:
user_region = st_centroid(MyGrid) %>% dplyr::mutate(Lon = sf::st_coordinates(.)[,1], Lat = sf::st_coordinates(.)[,2])
st_geometry(user_region) = NULL
saveRDS(user_region, file = file.path(save_data_folder, 'user_region.rds'))
