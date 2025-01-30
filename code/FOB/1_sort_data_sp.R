rm(list = ls())
require(ggplot2)
require(sf)
require(dplyr)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Define type of school to be analyzed:
this_type = 'FOB'

# -------------------------------------------------------------------------
# Data and plot folders:
save_data_folder = file.path('data', this_type)
dir.create(save_data_folder, recursive = TRUE, showWarnings = FALSE)
plot_dir = file.path('figures', this_type)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load data
load(file.path(save_data_folder, 'joinDF.RData'))
st_geometry(joinDF) = NULL
my_data = joinDF

# -------------------------------------------------------------------------
# Filter data based on some criteria:
summary(my_data)
my_data = my_data %>% dplyr::filter(year >= 2013)

# -------------------------------------------------------------------------
# define sp to be deleted from dataset (usually target species)
del_sp = c('Katsuwonus pelamis', 'Thunnus albacares', 'Thunnus obesus',
           'Auxis rochei', 'Auxis thazard', 'Auxis thazard, A. rochei',
           'Euthynnus affinis', 'Euthynnus alletteratus',
           'Thunnus alalunga') 

# tons data frame:
del_cols = c('number_', 'weight_')
sel_col = 'tons_'
tmp_data = my_data %>% select(-grep(pattern = paste(del_cols, collapse = "|"), x = colnames(my_data)))
tons_data = tidyr::pivot_longer(data = tmp_data, cols = grep(pattern = sel_col, x = colnames(tmp_data)),
                                names_to = 'sp_code', values_to = 'value') %>% 
  mutate(sp_code = gsub(pattern = sel_col, replacement = '', x = sp_code)) %>%
  dplyr::filter(!(sp_code %in% del_sp))
# number data frame:
del_cols = c('tons_', 'weight_')
sel_col = 'number_'
tmp_data = my_data %>% select(-grep(pattern = paste(del_cols, collapse = "|"), x = colnames(my_data)))
numbers_data = tidyr::pivot_longer(data = tmp_data, cols = grep(pattern = sel_col, x = colnames(tmp_data)),
                                   names_to = 'sp_code', values_to = 'value') %>% 
  mutate(sp_code = gsub(pattern = sel_col, replacement = '', x = sp_code)) %>%
  dplyr::filter(!(sp_code %in% del_sp))
# weight data frame:
del_cols = c('tons_', 'number_')
sel_col = 'weight_'
tmp_data = my_data %>% select(-grep(pattern = paste(del_cols, collapse = "|"), x = colnames(my_data)))
weight_data = tidyr::pivot_longer(data = tmp_data, cols = grep(pattern = sel_col, x = colnames(tmp_data)),
                                  names_to = 'sp_code', values_to = 'value') %>% 
  mutate(sp_code = gsub(pattern = sel_col, replacement = '', x = sp_code)) %>%
  dplyr::filter(!(sp_code %in% del_sp))

# Save created data:
saveRDS(object = tons_data, file = file.path(save_data_folder, 'tons_data.rds'))
saveRDS(object = numbers_data, file = file.path(save_data_folder, 'numbers_data.rds'))
saveRDS(object = weight_data, file = file.path(save_data_folder, 'weight_data.rds'))

# -------------------------------------------------------------------------
# Make figures dominance sp groups:
n_sp = 25 # first N species

# Weight:
plot_data = weight_data %>% group_by(year, sp_code) %>% summarise(value = sum(value))
cumsp_data = plot_data %>% group_by(sp_code) %>% summarise(value = sum(value))
cumsp_data = arrange(cumsp_data, desc(value))
write.csv(cumsp_data, file = file.path(save_data_folder, 'cumsp_data_weight.csv'), row.names = FALSE)
cumsp_data = cumsp_data[1:n_sp, ]
plot_data = plot_data %>% dplyr::filter(sp_code %in% cumsp_data$sp_code) %>% mutate(sp_code = factor(sp_code, levels = cumsp_data$sp_code))

p1 = ggplot(data = plot_data, aes(x = sp_code, y = value)) +
  geom_col() + 
  xlab(NULL) + ylab('Weight (kg)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
  facet_wrap(~ year)
ggsave(paste0('sp_dom_weight', img_type), path = plot_dir, plot = p1,
       width = img_width, height = 140, units = 'mm', dpi = img_res)

# Numbers:
plot_data = numbers_data %>% group_by(year, sp_code) %>% summarise(value = sum(value))
cumsp_data = plot_data %>% group_by(sp_code) %>% summarise(value = sum(value))
cumsp_data = arrange(cumsp_data, desc(value))
write.csv(cumsp_data, file = file.path(save_data_folder, 'cumsp_data_numbers.csv'), row.names = FALSE)
cumsp_data = cumsp_data[1:n_sp, ]
plot_data = plot_data %>% dplyr::filter(sp_code %in% cumsp_data$sp_code) %>% mutate(sp_code = factor(sp_code, levels = cumsp_data$sp_code))

p1 = ggplot(data = plot_data, aes(x = sp_code, y = log(value+1))) +
  geom_col() + 
  xlab(NULL) + ylab('log(Numbers)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
  facet_wrap(~ year)
ggsave(paste0('sp_dom_numbers', img_type), path = plot_dir, plot = p1,
       width = img_width, height = 140, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Make maps by sp group (weight):
save_sp_maps = 'presence_maps'
dir.create(file.path(plot_dir, save_sp_maps), recursive = TRUE, showWarnings = FALSE)
all_sp = unique(weight_data$sp_code)

for(i in seq_along(all_sp)) {
  
  sel_sp = all_sp[i]
  tmp_data = weight_data %>% dplyr::filter(sp_code == sel_sp)
  tmp_data = tmp_data %>% mutate(presence = as.factor(if_else(value > 0, 1, 0)))
  plot_data = tmp_data %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  
  p1 = ggplot(plot_data) + geom_sf(aes(color = presence), size = 1, alpha = 0.5)
  p1 = add_sf_map(p1)
  p1 = p1 + ggtitle(label = sel_sp) + theme(legend.position = 'none') + 
    facet_wrap(~ year)
  ggsave(paste0(gsub(' ', '-', sel_sp), img_type), path = file.path(plot_dir, save_sp_maps), 
         plot = p1, width = img_width, height = 180, units = 'mm', dpi = img_res)
  print(i)
}
