rm(list = ls())

# Define type of school to be analyzed:
source('code/1_single/load_libs.R')

# -------------------------------------------------------------------------
# Load data
wtData = readRDS(file.path(data_folder, 'weight_data.rds'))
effPoints = readRDS(file.path(data_folder, 'effPoints.rds'))

# -------------------------------------------------------------------------
# Method 1: raised by production (see Amande et al 2010): 10.1051/alr/2011003
# Strata are grid/year/quarter. 

grData = wtData %>% group_by(ID, year, quarter, sp_name) %>% 
                summarise(bycatch = sum(bycatch),
                          production = sum(trop_catch_nonstd))
raiseData = grData %>% mutate(ratio = bycatch/production)
# Replace Inf and NaN with zeros due to production == 0:
# This is a caveat of this method
raiseData = raiseData %>% mutate(ratio = ifelse(is.infinite(ratio), 0, ratio))
raiseData = raiseData %>% mutate(ratio = ifelse(is.nan(ratio), 0, ratio))

# Calculate ratio per year/quarter/species for whole area:
raiseData_tot = grData %>% group_by(year, quarter, sp_name) %>% 
                summarise(ratio_tot = sum(bycatch)/sum(production))

# Summarise effort data:
logbookData = effPoints %>% group_by(ID, year, quarter) %>% 
                summarise(production = sum(trop_catch_nonstd))

# Merge dfs:
mergedData = left_join(logbookData, raiseData[,c('ID', 'year', 'quarter', 'sp_name', 'ratio')])
mergedData = left_join(mergedData, raiseData_tot, by = c('year', 'quarter', 'sp_name'))
# Check if NA comes from grid with no bycatch:
identical( which(is.na(mergedData$sp_name)), which(is.na(mergedData$ratio)) )
# Remove those NA safely:
mergedData = mergedData %>% filter(!is.na(ratio))
# Check again if any NA and replace with tot ratio:
if(any(is.na(mergedData$ratio))) {
  mergedData = mergedData %>% mutate(ratio = ifelse(is.na(ratio), ratio_tot, ratio))
}

# Estimates by year and quarter:
mergedData = mergedData %>% mutate(est_prod = ratio*production)
estimateProd = mergedData %>% group_by(year, quarter, sp_name) %>% 
                summarise(est_prod = sum(est_prod, na.rm = TRUE)) 


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Method 2: raised by effort (number of sets)
# Strata are grid/year/quarter. 

grData = wtData %>% group_by(ID, year, quarter, sp_name) %>% 
  summarise(bycatch = sum(bycatch),
            n_sets = n()) 
raiseData = grData %>% mutate(ratio = bycatch/n_sets)

# Calculate ratio per year/quarter/species for whole area:
raiseData_tot = grData %>% group_by(year, quarter, sp_name) %>% 
  summarise(ratio_tot = sum(bycatch)/n())

# Summarise effort data:
logbookData = effPoints %>% group_by(ID, year, quarter) %>% 
  summarise(n_sets = n())

# Merge dfs:
mergedData = left_join(logbookData, raiseData[,c('ID', 'year', 'quarter', 'sp_name', 'ratio')])
mergedData = left_join(mergedData, raiseData_tot, by = c('year', 'quarter', 'sp_name'))
# Check if NA comes from grid with no bycatch:
identical( which(is.na(mergedData$sp_name)), which(is.na(mergedData$ratio)) )
# Remove those NA safely:
mergedData = mergedData %>% filter(!is.na(ratio))
# Check again if any NA and replace with tot ratio:
if(any(is.na(mergedData$ratio))) {
  mergedData = mergedData %>% mutate(ratio = ifelse(is.na(ratio), ratio_tot, ratio))
}

# Estimates by year and quarter:
mergedData = mergedData %>% mutate(est_sets = ratio*n_sets)
estimateSets = mergedData %>% group_by(year, quarter, sp_name) %>% 
  summarise(est_sets = sum(est_sets, na.rm = TRUE)) 


# -------------------------------------------------------------------------
# Save ratio estimates:
outData = left_join(estimateProd, estimateSets, by = c('year', 'quarter', 'sp_name'))
saveRDS(outData, file = file.path(data_folder, 'ratio_estimates.rds'))

# -------------------------------------------------------------------------
# Compare plot by year:
plot_data = outData %>% 
  pivot_longer(cols = c(est_prod, est_sets), names_to = 'est_type', values_to = 'est_value') %>%
  mutate(est_type = factor(est_type, levels = c('est_prod', 'est_sets'), 
                           labels = c('Production', 'Effort (sets)')))
plot_data = plot_data %>% group_by(year, sp_name, est_type) %>%
  summarise(est_value = sum(est_value, na.rm = TRUE))

ggplot(plot_data, aes(x = year, y = est_value, color = est_type)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(from = 2014, to = 2022, by = 4)) +
  facet_wrap(~ sp_name, scales = 'free_y') +
  labs(x = 'Year', y = 'Estimated bycatch (tons)', color = 'Estimation type') +
  theme_classic() +
  theme(legend.position = c(0.8, 0.03))
ggsave(file.path(plot_folder, 'ratio_estimates_by_year.png'), width = img_width*1.5, 
       height = 180, units = 'mm', dpi = img_res)
