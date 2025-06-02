rm(list = ls())

# Define type of school to be analyzed:
this_type = 'FOB' # or FSC
source('code/2_multi/load_libs.R')

# -------------------------------------------------------------------------
# Define folder where data is located:
raw_data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_PS_Bycatch'
# Read data in (observations):
load(file.path(raw_data_folder, 'data.weka.ATL.2003_2023_MarineBeacon_enviado_GC_110225.RData'))

# Read data ---------------------------------------------------------------
# Explore data:
glimpse(data_weka)

# Define columns to be used:
sp_id_col = grep(pattern = 'tons_|number_|weight_', x = colnames(data_weka))
var_columns = c('ocean_code', 'vessel_code', 'flag_country', 'trip_start_date', 'trip_end_date',
               'observation_date', 'observation_time', 'latitude', 'longitude',  
               'school_type', 'sunrise_diference')
var_id_col = which(colnames(data_weka) %in% var_columns)
my_data = data_weka %>% select(all_of(c(var_id_col, sp_id_col)))

# Do some IMPORTANT filtering:
my_data = my_data %>% dplyr::filter(ocean_code == '1', # 1 = ATL
                                    school_type == this_type # defined above
                                    ) 

# Continue..
my_data = my_data %>% mutate(id_set = 1:n(),
                             year = as.integer(format(observation_date, '%Y')),
                             quarter = as.factor(ceiling(as.numeric(format(observation_date, '%m'))/3)), 
                             .before = 'ocean_code')
my_data = my_data %>% rename(date = observation_date,
                             trop_catch = tons_target_tuna,
                             lon = longitude,
                             lat = latitude)
# Filter data from 2013: (same as effort data)
my_data = my_data %>% dplyr::filter(year >= 2015) # from 2015

# Explore:
glimpse(my_data)
summary(my_data)

# Delete NA rows:
dim(my_data)
my_data = my_data[complete.cases(my_data), ]
dim(my_data)

# -------------------------------------------------------------------------
# Observations as sf:
obsDF = my_data %>% st_as_sf(coords = c("lon", "lat"), 
                             crs = 4326, remove = FALSE)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Create grid to explore predictions:
grid_size = 1 # in degrees

# Point and grid (aggregated):
min_lon = floor(min(obsDF$lon))
min_lat = floor(min(obsDF$lat))
MyGrid = st_make_grid(obsDF, cellsize = c(grid_size, grid_size), offset = c(min_lon, min_lat)) %>%  # define grid size here
  st_set_crs(4326) %>% st_sf() %>% dplyr::mutate(ID = 1:n())

# Join both Grid and Points (effort data):
obsPoints = st_join(MyGrid, left = TRUE, obsDF) %>% na.omit
# Do it in this way to avoid repeated rows when rounded lon/lat fall on grid borders:
index_dup = obsPoints %>% st_drop_geometry() %>% select(-ID) %>% duplicated # check for duplicates
obsPoints = obsPoints[!index_dup, ]
obsPoints = obsPoints %>% st_drop_geometry()
identical(nrow(obsPoints), nrow(obsDF))

# Filter based on some criteria:
# Identify grid (and points inside) with some criteria:
my_tab = xtabs(~ ID + year, data = obsPoints) # find frequency of sets per grid and year
my_freq = apply(my_tab, 1, function(x) sum(x > 0)) # find recurrent grids over the years
include_grid = names(my_freq)[which(as.vector(my_freq) >= 4)] # IMPORTANT: include > 3 years
include_grid = as.numeric(include_grid)

# Remove grids:
n_obs = nrow(obsPoints)
obsPoints = obsPoints %>% filter(ID %in% include_grid)
cat(round((1-(nrow(obsPoints)/n_obs))*100, 2), '% observations removed\n')
saveRDS(obsPoints, file = file.path(data_folder, 'obsPoints.rds'))

# Apply filter to grid df and calculate Area km2:
MyGrid = MyGrid %>% dplyr::filter(ID %in% include_grid)
plot(MyGrid)
MyGrid$Area_km2 = as.numeric(st_area(MyGrid))*1e-06 # in km2
save(MyGrid, file = file.path(data_folder, 'MyGrid.RData'))
