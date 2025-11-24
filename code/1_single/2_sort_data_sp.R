rm(list = ls())

# Define type of school to be analyzed:
source('code/1_single/load_libs.R')

# -------------------------------------------------------------------------
# Load data
my_data = readRDS(file.path(data_folder, 'obsPoints.rds'))

# -------------------------------------------------------------------------
# # define sp to be deleted from dataset (usually target species)
# del_sp = c('Katsuwonus pelamis', 'Thunnus albacares', 'Thunnus obesus',
#            'Auxis rochei', 'Auxis thazard', 'Auxis thazard, A. rochei',
#            'Euthynnus affinis', 'Euthynnus alletteratus',
#            'Thunnus alalunga')

# -------------------------------------------------------------------------
# 'number' data frame:
del_cols = c('weight_')
sel_col = 'number_'
tmp_data = my_data
tmp_data = tmp_data %>% select(-grep(pattern = paste(del_cols, collapse = "|"), x = colnames(my_data)))
numbers_data = tidyr::pivot_longer(data = tmp_data, cols = grep(pattern = sel_col, x = colnames(tmp_data)),
                                   names_to = 'sp_name', values_to = 'value') %>% 
                  mutate(sp_name = gsub(pattern = sel_col, replacement = '', x = sp_name)) 
                  # dplyr::filter(!(sp_name %in% del_sp))
# Final edits:
numbers_data = numbers_data %>% dplyr::rename(bycatch = value)

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

# Remove species with zero bycatch
zero_sp = weight_data %>% group_by(sp_name) %>% 
  summarise(bycatch = sum(bycatch)) %>% 
  dplyr::filter(bycatch == 0) %>% pull(sp_name)
if(length(zero_sp) > 0) weight_data = weight_data %>% dplyr::filter(!(sp_name %in% zero_sp))

# Now 'clean' sp short name by hand:
unique(weight_data$sp_name)
weight_data$sp_name[weight_data$sp_name == 'Alopidae'] <- 'Alopiidae'
weight_data$sp_name[weight_data$sp_name == 'Aluterus'] <- 'Aluterus spp.'
weight_data$sp_name[weight_data$sp_name == 'A. spp'] <- 'Auxis spp.'
weight_data$sp_name[weight_data$sp_name == 'bony fish other family'] <- 'Osteichthyes'
weight_data$sp_name[weight_data$sp_name == 'Kyphosus'] <- 'Kyphosus spp.'
weight_data$sp_name[weight_data$sp_name == 'O. balistidae'] <- 'Balistidae'
weight_data$sp_name[weight_data$sp_name == 'O. billfishes'] <- 'Istiophoridae' 
weight_data$sp_name[weight_data$sp_name == 'O. Carangidae'] <- 'Carangidae'
weight_data$sp_name[weight_data$sp_name == 'o. Carcharhinidae'] <- 'Carcharhinidae'
weight_data$sp_name[weight_data$sp_name == 'o. Scombridae'] <- 'Scombridae'
weight_data$sp_name[weight_data$sp_name == 'o. shark'] <- 'Elasmobranchii'
weight_data$sp_name[weight_data$sp_name == 's. turtles'] <- 'Chelonioidea'
weight_data$sp_name[weight_data$sp_name == 'Uraspis'] <- 'Uraspis spp.' 
unique(weight_data$sp_name)

# IMPORANT:
# Merge and exclude some sp after discussion with Jon:
weight_data$sp_name[weight_data$sp_name %in% c('Carcharhinidae', 'C. falciformis')] <- 'Carcharhinidae' 
weight_data$sp_name[weight_data$sp_name %in% c('Balistidae', 'C. maculata')] <- 'Balistidae' 
weight_data$sp_name[weight_data$sp_name %in% c('Sphyraenidae', 'S. barracuda')] <- 'Sphyraenidae' 
weight_data$sp_name[weight_data$sp_name %in% c('Carangidae', 'C. crysos')] <- 'Carangidae' 
weight_data = weight_data %>% filter(!(sp_name %in% c('Elasmobranchii', 'Osteichthyes', 'Scombridae', 'T. alalunga')))

# Save data for plotting later:
saveRDS(object = weight_data, file = file.path(data_folder, 'weight_data_all.rds'))

# IMPORTANT:
# Select species for analyses. This will depend on the goal: 'estimation' or 'sample-coverage'
# Do not filter when goal is 'estimation'
if(this_goal == 'sample-coverage') {
  if(this_type == 'FOB') {
    weight_data = weight_data %>% filter(sp_name %in% c('E. bipinnulata', 'Balistidae', 'Coryphaenidae', 'A. solandri', 'Carangidae',
                                                        'Carcharhinidae', 'M. nigricans', 'Sphyrnidae', 'Chelonioidea', 'Mobulidae',
                                                        'Alopiidae', 'Lamnidae', 'P. glauca'))
  }
  if(this_type == 'FSC') {
    weight_data = weight_data %>% filter(sp_name %in% c('Carcharhinidae', 'Mobulidae', 'I. albicans',
                                                        'M. nigricans', 'Sphyrnidae', 'Chelonioidea', 'Molidae',
                                                        'Lamnidae', 'P. glauca'))
  }
}
  
# Combine per id set:
weight_data_clean = weight_data %>% 
        group_by(id_set,ID,year,quarter,vessel_code,sp_name,flag_country,trip_start_date) %>%
        summarise(lon=mean(lon),lat=mean(lat),
                  sst=mean(sst),trop_catch=mean(trop_catch),
                  sst_nonstd = mean(sst_nonstd), 
                  trop_catch_nonstd = mean(trop_catch_nonstd),
                  bycatch=sum(bycatch))

# Delete if any species has zero bycatch:
del_sp = weight_data_clean %>% 
  group_by(sp_name) %>% summarise(tot_bycatch = sum(bycatch, na.rm = TRUE)) %>%
  filter(tot_bycatch == 0) %>% pull(sp_name)
weight_data_clean = weight_data_clean %>% filter(!sp_name %in% del_sp)

# -------------------------------------------------------------------------
# Save created data:
saveRDS(object = numbers_data, file = file.path(data_folder, 'numbers_data.rds'))
saveRDS(object = weight_data_clean, file = file.path(data_folder, 'weight_data.rds'))

# -------------------------------------------------------------------------
# Now select species to be included in model

# Find years with zero presence and Moran I's p-value:
moran_data = weight_data_clean %>% group_by(year, sp_name) %>% 
  summarise(moran_pval = get_moran(bycatch, lon, lat),
            n_zero = length(which(bycatch == 0)),
            porc_zero = n_zero/n() )
moran_data = moran_data %>% mutate(moran_sig = ifelse(moran_pval <= 0.05, 'sig', 'not sig'))
saveRDS(moran_data, file = file.path(data_folder, 'moran_data_weight.rds'))
p1 = ggplot(data = moran_data, aes(x = year, y = moran_pval)) +
  geom_point(aes(color = moran_sig)) +
  scale_x_continuous(breaks = seq(from = 2014, to = 2022, by = 4)) +
  xlab(NULL) + ylab('Moran I p-value') +
  theme(legend.position = 'none') +
  facet_wrap(~ sp_name) 
ggsave(paste0('moran_weight', img_type), path = plot_folder, plot = p1,
       width = img_width*1.5, height = 180, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Define modelling categories
n_years = length(unique(moran_data$year))
sel_sp_data = moran_data %>% filter(!is.na(moran_pval)) %>% 
  group_by(sp_name) %>%
  summarise(n_years = n(),
            prop_zero_per_year = mean(porc_zero),
            n_years_spat_cor = sum(moran_pval <= 0.05)) %>% 
  arrange(desc(n_years_spat_cor) )
sp_quarter = weight_data_clean %>% group_by(quarter, sp_name) %>%
                summarise(bycatch = sum(bycatch)) %>% group_by(sp_name) %>%
                summarise(n_quarters = sum(bycatch > 0),
                          tot_bycatch = sum(bycatch))
sel_sp_data = left_join(sel_sp_data, sp_quarter, by = 'sp_name')
summary(sel_sp_data) # check no NAs
View(sel_sp_data)
sel_sp_data = sel_sp_data %>% arrange(desc(tot_bycatch)) %>% mutate(cum_bycatch = cumsum(tot_bycatch)*100/sum(tot_bycatch))
saveRDS(sel_sp_data, file = file.path(data_folder, 'model_cat_sp.rds'))
