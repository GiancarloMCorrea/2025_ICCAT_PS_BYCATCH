require(gt)
source('code/3_sampsize/load_libs.R')

# -------------------------------------------------------------------------
# Make diagram sampling experiment:

# Alternative folder to save diagrams:
# This is important in order to make them public and use the Mermaid editor later:
alt_plot_folder = 'C:/Use/OneDrive - AZTI/Manuscripts/SampCoverage_PSBycatch'

# Create fake data:

# Effort data:
mydat = data.frame(Trip_id = rep(1:4, each = 3), 
                   Set_id = rep(1:3, times = 4),
                   Year = rep(2015:2016, each = 6), 
                   Quarter = rep(c(2,1,4,5), each = 3), 
                   T_catch = round(exp(rlnorm(n = 12, meanlog = 1, sdlog = 0.2)), digits = 2))
# Simulated data:
mydat2 = mydat %>% mutate(Bycatch = round(exp(rlnorm(n = 12, meanlog = 0.25, sdlog = 0.2)), digits = 2))
# Sampled data:
mydat3 = mydat2 %>% filter(Trip_id %in% c(1,3))

# Make figure 
colpal = RColorBrewer::brewer.pal(n = 4, name = 'Set3')

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
