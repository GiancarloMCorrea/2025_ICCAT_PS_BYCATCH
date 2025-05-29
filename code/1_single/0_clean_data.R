rm(list = ls())

# Define type of school to be analyzed:
this_type = 'FOB' # or FSC
source('code/1_single/load_libs.R')

# -------------------------------------------------------------------------
# Define folder where data is located:
raw_data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_PS_Bycatch'
# Read data in (observations):
load(file.path(raw_data_folder, 'data.weka.ATL.2003_2023_MarineBeacon_enviado_GC_110225.RData'))
# Read all sets (effort) data:
allsets_df = read.csv(file.path(raw_data_folder, 'ATL_13a23_setsfleet.csv'))

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
my_data = my_data %>% dplyr::filter(year >= 2013)

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
# Explore effort data:
glimpse(allsets_df)
eff_data = allsets_df
alt_set_type = this_type
if(this_type == 'FOB') alt_set_type = 'FAD'

# Do some IMPORTANT filtering:
eff_data = eff_data %>% dplyr::filter(oceano == 1, # 1 = ATL
                                      set_type == alt_set_type) # defined above
# Rename some variables:
eff_data = eff_data %>% dplyr::rename(year = año,
                                      trop_catch = suma_de_YFT.BET.SKJ,
                                      no_sets = numero_de_lances)
# Create some columns:
eff_data = eff_data %>% mutate(date = as.Date(fecha), .after = 'fecha')
eff_data = eff_data %>% mutate(quarter = as.factor(ceiling(as.numeric(format(date, '%m'))/3)),
                               .after = 'date')

# Remove rows with number of sets = 0
eff_data = eff_data %>% dplyr::filter(no_sets > 0)
# Distribute tuna catch among no sets:
eff_data = eff_data %>% mutate(trop_catch = trop_catch/no_sets)
# Repeat rows based on the number of sets:
eff_data = eff_data %>% tidyr::uncount(no_sets)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Create grid to explore predictions:
# Use effort data:
grid_size = 5 # in degrees

# Point and grid (aggregated):
effDF = eff_data %>% mutate(id_set = 1:n()) %>%
            st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
min_lon = floor(min(effDF$lon))
min_lat = floor(min(effDF$lat))
MyGrid = st_make_grid(effDF, cellsize = c(grid_size, grid_size), offset = c(min_lon, min_lat)) %>%  # define grid size here
  st_set_crs(4326) %>% st_sf() %>% dplyr::mutate(ID = 1:n())

# Join both Grid and Points (effort data):
effPoints = st_join(MyGrid, left = TRUE, effDF) %>% na.omit
# Do it in this way to avoid repeated rows when rounded lon/lat fall on grid borders:
index_dup = effPoints %>% st_drop_geometry() %>% select(-ID) %>% duplicated # check for duplicates
effPoints = effPoints[!index_dup, ]
effPoints = effPoints %>% st_drop_geometry()
identical(nrow(effPoints), nrow(effDF))
saveRDS(effPoints, file = file.path(data_folder, 'effPoints.rds'))

# Also attach ID info to observations:
# Some points may be removed since obs (all fleets) and eff is SPA
nrow(obsDF)
obsPoints = st_join(MyGrid, left = TRUE, obsDF) %>% na.omit
index_dup = obsPoints %>% st_drop_geometry() %>% select(-ID) %>% duplicated # check for duplicates
obsPoints = obsPoints[!index_dup, ]
nrow(obsPoints)
obsPoints = obsPoints %>% st_drop_geometry()
saveRDS(obsPoints, file = file.path(data_folder, 'obsPoints.rds'))

# Filter based on some criteria:
# Identify grid (and points inside) with some criteria:
my_tab = xtabs(~ ID + year, data = effPoints) # find frequency of sets per grid and year
my_freq = apply(my_tab, 1, function(x) sum(x > 0)) # find recurrent grids over the years
include_grid = names(my_freq)[which(as.vector(my_freq) >= 0)] # IMPORTANT: include all!
include_grid = as.numeric(include_grid)
# Apply filter to grid df and calculate Area km2:
MyGrid = MyGrid %>% dplyr::filter(ID %in% include_grid)
plot(MyGrid)
MyGrid$Area_km2 = as.numeric(st_area(MyGrid))*1e-06 # in km2
save(MyGrid, file = file.path(data_folder, 'MyGrid.RData'))
