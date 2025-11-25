require(gt)
source('code/3_sampsize/load_libs.R')

# -------------------------------------------------------------------------
# Make diagram sampling experiment:

# Alternative folder to save diagrams:
# This is important in order to make them public and use the Mermaid editor later:
alt_plot_folder = "figures/3_sampsize"

# Create fake data:

# Effort data:
mydat = data.frame(Trip_id = rep(1:6, each = 3), 
                   Set_id = rep(1:3, times = 6),
                   Year = rep(2015:2016, each = 9), 
                   Quarter = rep(c(2,1,4,3,2,1), each = 3), 
                   T_catch = round(exp(rlnorm(n = 18, meanlog = 1, sdlog = 0.2)), digits = 2))
# Simulated data:
mydat2 = mydat %>% mutate(Bycatch = round(exp(rlnorm(n = 18, meanlog = 0.25, sdlog = 0.2)), digits = 2))
# Observers data:
mydat4 = mydat2 %>% filter(Trip_id %in% c(1,3,4,6)) %>% mutate(Bycatch = round(Bycatch*rnorm(n = 12, mean = 1, sd = 0.1), digits = 2))
# Sampled data:
mydat3 = mydat4 %>% filter(Trip_id %in% c(1,4))

# Make figure 
colpal = RColorBrewer::brewer.pal(n = 6, name = 'Set3')

# Make effort data:
mydat %>% gt %>% data_color(
  columns = Trip_id,
  target_columns = everything(),
  palette = colpal
) %>% tab_style(
  style = cell_text(weight = "bold"),
  locations = cells_column_labels(columns = everything())
) %>% gtsave(filename = file.path(alt_plot_folder, 'tab_eff.png'))

# Make simulated data:
mydat2 %>% gt %>% data_color(
  columns = Trip_id,
  target_columns = everything(),
  palette = colpal
) %>% tab_style(
  style = cell_text(weight = "bold"),
  locations = cells_column_labels(columns = everything())
) %>% gtsave(filename = file.path(alt_plot_folder, 'tab_sim.png'))

# Make sampled data:
mydat3 %>% gt %>% data_color(
  columns = Trip_id,
  target_columns = everything(),
  palette = colpal[c(1,3)]
) %>% tab_style(
  style = cell_text(weight = "bold"),
  locations = cells_column_labels(columns = everything())
) %>% gtsave(filename = file.path(alt_plot_folder, 'tab_samp.png'))

# Make simulated data:
mydat4 %>% gt %>% data_color(
  columns = Trip_id,
  target_columns = everything(),
  palette = colpal
) %>% tab_style(
  style = cell_text(weight = "bold"),
  locations = cells_column_labels(columns = everything())
) %>% gtsave(filename = file.path(alt_plot_folder, 'tab_obs.png'))


# -------------------------------------------------------------------------
# Make figure to compare % presence in sets by selected species:

# FOB plot:
fob_df = readRDS(file.path("data/1_single/FOB", "weight_data.rds")) 
plot_data = fob_df %>% group_by(sp_name) %>% summarise(freq = (sum(bycatch > 0)/n())*100)
plot_data = plot_data %>% left_join(fob_sp_df %>% rename(sp_name = sp_levels), by = 'sp_name')
plot_data = plot_data %>% mutate(sp_name = factor(sp_name, levels = fob_sp_df$sp_levels),
                                 sp_type = factor(sp_type, levels = c('Common', 'Special interest', 'Rare')))

p1 = ggplot(data = plot_data, aes(x = sp_name, y = freq)) +
  geom_col(aes(fill = sp_type)) +
  scale_fill_brewer(palette = 'Set1') +
  xlab(NULL) + ylab("Frequency in FOB sets (%)") +
  guides(fill = guide_legend(title = "")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9),
        legend.background = element_rect(fill = NA),
        legend.position = c(0.77, 0.85)) 

# FSC plot:
fsc_df = readRDS(file.path("data/1_single/FSC", "weight_data.rds")) 
plot_data = fsc_df %>% group_by(sp_name) %>% summarise(freq = (sum(bycatch > 0)/n())*100)
plot_data = plot_data %>% left_join(fsc_sp_df %>% rename(sp_name = sp_levels), by = 'sp_name')
plot_data = plot_data %>% mutate(sp_name = factor(sp_name, levels = fsc_sp_df$sp_levels),
                                 sp_type = factor(sp_type, levels = c('Common', 'Special interest', 'Rare')))

p2 = ggplot(data = plot_data, aes(x = sp_name, y = freq)) +
  geom_col(aes(fill = sp_type)) +
  scale_fill_brewer(palette = 'Set1') +
  xlab(NULL) + ylab("Frequency in FSC sets (%)") +
  guides(fill = guide_legend(title = "")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9),
        legend.position = 'none') 

# Merge both plots:
p3 = grid.arrange(p1, p2, ncol = 2)
ggsave(paste0('freq_sp', img_type), plot = p3, path = "figures/3_sampsize", 
       width = img_width , height = 90, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Plot grid by set type:

load('data/1_single/FOB/MyGrid.RData')
p1 = ggplot(MyGrid) + geom_sf(fill = 'white')
p1 = add_sf_map(p1) 
p1 = p1 + ggtitle("FOB")

load('data/1_single/FSC/MyGrid.RData')
p2 = ggplot(MyGrid) + geom_sf(fill = 'white')
p2 = add_sf_map(p2) 
p2 = p2 + ggtitle("FSC")

# Merge both plots:
p3 = grid.arrange(p1, p2, ncol = 2)
ggsave(paste0('grid', img_type), plot = p3, path = "figures/3_sampsize", 
       width = img_width , height = 90, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Map points:

# FOB:
obsPoints = readRDS(file.path("data/1_single/FOB", 'obsPoints.rds'))
obsSF = obsPoints %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
p1 = ggplot(obsSF) + geom_sf(size = 1, alpha = 0.5) + ggtitle('FOB')
p1 = add_sf_map(p1)

# FSC:
obsPoints = readRDS(file.path("data/1_single/FSC", 'obsPoints.rds'))
obsSF = obsPoints %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
p2 = ggplot(obsSF) + geom_sf(size = 1, alpha = 0.5) + ggtitle('FSC')
p2 = add_sf_map(p2)

# Merge both plots:
p3 = grid.arrange(p1, p2, ncol = 2)
ggsave(paste0('sets', img_type), plot = p3, path = "figures/3_sampsize", 
       width = img_width , height = 90, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Number of observations:

# FOB:
obsPoints = readRDS(file.path("data/1_single/FOB", 'obsPoints.rds'))
p1 = ggplot(obsPoints, aes(x = factor(year))) + 
  geom_bar() + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab(NULL) + ylab('Number of observed sets') +
  ggtitle('FOB')

# FSC:
obsPoints = readRDS(file.path("data/1_single/FSC", 'obsPoints.rds'))
p2 = ggplot(obsPoints, aes(x = factor(year))) + 
  geom_bar() + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab(NULL) + ylab('Number of observed sets') +
  ggtitle('FSC')

# Merge both plots:
p3 = grid.arrange(p1, p2, ncol = 2)
ggsave(paste0('nobs', img_type), plot = p3, path = "figures/3_sampsize", 
       width = img_width , height = 90, units = 'mm', dpi = img_res)
