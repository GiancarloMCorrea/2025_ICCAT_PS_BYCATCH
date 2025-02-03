rm(list = ls())
# -------------------------------------------------------------------------
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Select species and school type:
this_type = 'FOB'

# Define folder to save results
plot_folder = file.path('figures', this_type)
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Data folder:
data_folder = file.path('data', this_type)

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path(data_folder, 'weight_data.rds'))
tons_data = readRDS(file = file.path(data_folder, 'tons_data.rds'))
numbers_data = readRDS(file = file.path(data_folder, 'numbers_data.rds'))


# -------------------------------------------------------------------------
# Make figures dominance sp groups:
n_sp = 25 # first N species

# Weight:
plot_data = weight_data %>% group_by(year, sp_code) %>% summarise(value = sum(value))
cumsp_data = plot_data %>% group_by(sp_code) %>% summarise(value = sum(value))
cumsp_data = arrange(cumsp_data, desc(value))
write.csv(cumsp_data %>% select(sp_code), file = file.path(data_folder, 'cumsp_data_weight.csv'), row.names = FALSE)
cumsp_data = cumsp_data[1:n_sp, ]
plot_data = plot_data %>% dplyr::filter(sp_code %in% cumsp_data$sp_code) %>% mutate(sp_code = factor(sp_code, levels = cumsp_data$sp_code))

p1 = ggplot(data = plot_data, aes(x = sp_code, y = value)) +
  geom_col() + 
  xlab(NULL) + ylab('Weight (kg)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
  facet_wrap(~ year)
ggsave(paste0('sp_dom_weight', img_type), path = plot_folder, plot = p1,
       width = img_width, height = 140, units = 'mm', dpi = img_res)

# Numbers:
plot_data = numbers_data %>% group_by(year, sp_code) %>% summarise(value = sum(value))
cumsp_data = plot_data %>% group_by(sp_code) %>% summarise(value = sum(value))
cumsp_data = arrange(cumsp_data, desc(value))
write.csv(cumsp_data %>% select(sp_code), file = file.path(data_folder, 'cumsp_data_numbers.csv'), row.names = FALSE)
cumsp_data = cumsp_data[1:n_sp, ]
plot_data = plot_data %>% dplyr::filter(sp_code %in% cumsp_data$sp_code) %>% mutate(sp_code = factor(sp_code, levels = cumsp_data$sp_code))

p1 = ggplot(data = plot_data, aes(x = sp_code, y = log(value+1))) +
  geom_col() + 
  xlab(NULL) + ylab('log(Numbers)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
  facet_wrap(~ year)
ggsave(paste0('sp_dom_numbers', img_type), path = plot_folder, plot = p1,
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
  ggsave(paste0(gsub(' ', '-', sel_sp), img_type), path = file.path(plot_folder, save_sp_maps), 
         plot = p1, width = img_width, height = 180, units = 'mm', dpi = img_res)
  print(i)
}


# -------------------------------------------------------------------------
# Fishing locations by season/month:
plot_data = weight_data %>% group_by(year, month, id_set) %>% 
              summarise(lon = mean(longitude), lat = mean(latitude))
plot_data = plot_data %>% mutate(month = as.numeric(month),
                                 season = ifelse(month %in% 1:3, '1',
                                                ifelse(month %in% 4:6, '2',
                                                       ifelse(month %in% 7:9, '3', '4'))))
plot_data = plot_data %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Make plot:
p1 = ggplot(plot_data) + geom_sf(aes(color = season), size = 0.5, alpha = 0.5)
p1 = add_sf_map(p1)
p1 = p1 + scale_color_brewer(palette = 'Set1') + 
  theme(legend.position = 'bottom') + facet_wrap(~ year)
ggsave(paste0('fishing_locations_season', img_type), path = plot_folder, plot = p1,
       width = img_width, height = 150, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Number bycatch sp per set
plot_data = weight_data %>% group_by(year, month, id_set) %>% 
  summarise(bycatch_catch = sum(value), n_sp_catch = length(which(value > 0)))
plot_data = plot_data %>% mutate(month = as.numeric(month),
                                 season = ifelse(month %in% 1:3, '1',
                                                 ifelse(month %in% 4:6, '2',
                                                        ifelse(month %in% 7:9, '3', '4'))))

p1 = ggplot(plot_data, aes(x = factor(year), y = n_sp_catch)) + 
  geom_boxplot(aes(fill = season), width = 0.8) +
  scale_fill_brewer(palette = 'Set1') +
  xlab(NULL) + ylab('N sp bycatch') +
  theme(legend.position = 'bottom')
ggsave(paste0('n_sp_per_set', img_type), path = plot_folder, plot = p1,
       width = img_width*0.75, height = 100, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Proportion of sets with positive bycatch:
plot_data = weight_data %>% group_by(year, month, id_set) %>% 
  summarise(n_sp_catch = length(which(value > 0)))
plot_data = plot_data %>% mutate(month = as.numeric(month),
                                 season = ifelse(month %in% 1:3, '1',
                                                 ifelse(month %in% 4:6, '2',
                                                        ifelse(month %in% 7:9, '3', '4'))))
plot_data = plot_data %>% group_by(year, season) %>% 
  summarise(prop_bycatch = (length(which(n_sp_catch > 0))/n())*100)

p1 = ggplot(data = plot_data, aes(x = factor(year), y = prop_bycatch)) +
  geom_col(aes(fill = season), position = 'dodge') +
  scale_fill_brewer(palette = 'Set1') +
  xlab(NULL) + ylab('Sets with bycatch (%)') +
  theme(legend.position = 'bottom')
ggsave(paste0('prop_set_bycatch', img_type), path = plot_folder, plot = p1,
       width = img_width*0.75, height = 100, units = 'mm', dpi = img_res)
