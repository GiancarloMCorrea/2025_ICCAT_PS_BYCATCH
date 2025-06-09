rm(list = ls())

# Define type of school to be analyzed:
this_type = 'FOB' 
source('code/2_multi/load_libs.R')

# -------------------------------------------------------------------------
# Load data
my_data = readRDS(file.path(data_folder, 'obsPoints.rds'))

# -------------------------------------------------------------------------
# 'weight' data frame:
del_cols = c('number_')
sel_col = 'weight_'
tmp_data = my_data
tmp_data = tmp_data %>% select(-grep(pattern = paste(del_cols, collapse = "|"), x = colnames(my_data)))
weight_data = tidyr::pivot_longer(data = tmp_data, cols = grep(pattern = sel_col, x = colnames(tmp_data)),
                                  names_to = 'sp_name', values_to = 'value') %>% 
                  mutate(sp_name = gsub(pattern = sel_col, replacement = '', x = sp_name)) 
# Final edits:
weight_data = weight_data %>% dplyr::rename(bycatch = value)
# Shorten species name:
weight_data = weight_data %>% mutate(sp_long_name = sp_name,
                                     sp_name = short_sp_name(sp_long_name))
# Now 'clean' sp short name by hand:
unique(weight_data$sp_name)
weight_data$sp_name[weight_data$sp_name == 'Alopidae'] <- 'Alopiidae'
weight_data$sp_name[weight_data$sp_name == 'Aluterus'] <- 'Aluterus spp.'
weight_data$sp_name[weight_data$sp_name == 'A. spp'] <- 'Auxis spp.'
weight_data$sp_name[weight_data$sp_name == 'bony fish other family'] <- 'Osteichthyes'
weight_data$sp_name[weight_data$sp_name == 'Kyphosus'] <- 'Kyphosus spp.'
weight_data$sp_name[weight_data$sp_name == 'O. balistidae'] <- 'Balistidae'
weight_data$sp_name[weight_data$sp_name == 'O. billfishes'] <- 'Xiphioidea'
weight_data$sp_name[weight_data$sp_name == 'O. Carangidae'] <- 'Carangidae'
weight_data$sp_name[weight_data$sp_name == 'o. Carcharhinidae'] <- 'Carcharhinidae'
weight_data$sp_name[weight_data$sp_name == 'o. Scombridae'] <- 'Scombridae'
weight_data$sp_name[weight_data$sp_name == 'o. shark'] <- 'Elasmobranchii'
weight_data$sp_name[weight_data$sp_name == 's. turtles'] <- 'Chelonioidea'
weight_data$sp_name[weight_data$sp_name == 'Uraspis'] <- 'Uraspis spp.' 

# Combine per id set:
weight_data_clean = weight_data %>% 
        group_by(id_set,ID,year,quarter,vessel_code,sp_name) %>%
        summarise(bycatch=sum(bycatch), lon = mean(lon), lat = mean(lat))

# Delete if any species has zero bycatch:
del_sp = weight_data_clean %>% 
  group_by(sp_name) %>% summarise(tot_bycatch = sum(bycatch, na.rm = TRUE)) %>%
  filter(tot_bycatch == 0) %>% pull(sp_name)
weight_data_clean = weight_data_clean %>% filter(!sp_name %in% del_sp)

# -------------------------------------------------------------------------
# define sp to be deleted from dataset (usually target species)
sort(unique(weight_data_clean$sp_name))
target_sp = c('K. pelamis', 'T. albacares', 'T. obesus',
              'Auxis spp.', 'T. alalunga')
weight_data_clean = weight_data_clean %>% filter(!sp_name %in% target_sp)

# -------------------------------------------------------------------------
# Define modelling categories
n_years = length(unique(weight_data_clean$year))
# Find number of years with presence:
sel_sp_data = weight_data_clean %>% 
  group_by(sp_name) %>% 
  summarise(n_years = n_distinct(year), 
            tot_bycatch = sum(bycatch, na.rm = TRUE)) %>%
  arrange(desc(tot_bycatch))
# Calculate cumulative percentage:
sel_sp_data = sel_sp_data %>% 
  mutate(cum_perc = cumsum(tot_bycatch)/sum(tot_bycatch) * 100)
saveRDS(sel_sp_data, file = file.path(data_folder, 'model_cat_sp.rds'))

# -------------------------------------------------------------------------
# Save created data:
saveRDS(object = weight_data_clean, file = file.path(data_folder, 'weight_data.rds'))
