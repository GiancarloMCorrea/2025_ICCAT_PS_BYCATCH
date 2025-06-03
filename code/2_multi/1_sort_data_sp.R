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
weight_data$sp_name[weight_data$sp_name == 'Sphyrnidae'] <- 'Sphyraenidae' # combine in one group
weight_data$sp_name[weight_data$sp_name == 'Uraspis'] <- 'Uraspis spp.' 

# Combine per id set:
weight_data_clean = weight_data %>% 
        group_by(id_set,ID,year,quarter,vessel_code,sp_name,
                 lon,lat) %>%
        summarise(bycatch=sum(bycatch))


# -------------------------------------------------------------------------
# define sp to be deleted from dataset (usually target species)
sort(unique(weight_data_clean$sp_name))
target_sp = c('K. pelamis', 'T. albacares', 'T. obesus',
              'Auxis spp.', 'T. alalunga')
weight_data_clean = weight_data_clean %>% filter(!sp_name %in% target_sp)

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
sel_sp_data = sel_sp_data %>% filter(n_years == n_years & n_years_spat_cor >= 7)
View(sel_sp_data)
saveRDS(sel_sp_data, file = file.path(data_folder, 'model_cat_sp.rds'))


# -------------------------------------------------------------------------

# Include final species with spatial autocorrelation
weight_data_clean = weight_data_clean %>% filter(sp_name %in% sel_sp_data$sp_name)

# -------------------------------------------------------------------------
# Save created data:
saveRDS(object = weight_data_clean, file = file.path(data_folder, 'weight_data.rds'))
