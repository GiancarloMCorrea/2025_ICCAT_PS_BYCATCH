rm(list = ls())

# Define type of school to be analyzed:
source('code/3_sampsize/load_libs.R')

# Read effort data:
effPoints = readRDS(file = file.path('data/1_single/FOB', 'effPoints.rds')) # for prediction

# -------------------------------------------------------------------------

# Find species with fitted model:
all_sp = list.files(file.path("model/1_single/FOB"))
these_sp = numeric(0)
for(k in seq_along(all_sp)) {
  this_file = file.path("model/1_single/FOB", all_sp[k], 'mod.RData')
  if(file.exists(this_file)) these_sp = c(these_sp, all_sp[k])
}

# Make seeds:
set.seed(8675309)
seeds = sample(x = (-1e9):(1e9), size = 1000, replace = FALSE)
saveRDS(seeds, file.path(data_folder, "seeds.RDS"))

# Sample effort data:
frac_vector = seq(from = 1, to = 0.1, by = -0.05)
i = 1 # loop over sample fraction vector
sampled_data = effPoints %>% group_by(year, quarter) %>% sample_frac(frac_vector[i]) # fraction

# Make simulations
nSims = 100
isp = 1 # loop over these_sp

rm(mod_upd)
load(file.path("model/1_single/FOB", these_sp[isp], 'mod.RData'))
mod_dat = mod_upd$data
yr_keep = unique(mod_dat$year)
qt_keep = unique(mod_dat$quarter)
pred_data = sampled_data %>% filter(year %in% yr_keep, quarter %in% qt_keep)
pred_data = pred_data %>% mutate(fyear = factor(year, levels = sort(unique(pred_data$year))))
sim_data = simulate(mod_upd, nsim = nSims, seed = seeds[1:nSims], newdata = pred_data)



dat = mod_upd$data
s1 = simulate(mod_upd, nsim = 50)
dim(s1)

sum(s1 == 0)/length(s1)
sum(dat$bycatch == 0) / length(dat$bycatch)

pos = which(dat$year == 2015 & dat$quarter == 3)
plot(dat$lon[pos], dat$lat[pos], pch = 1, cex = dat$bycatch[pos]*5)
points(dat$lon[pos], dat$lat[pos], pch = 1, col = "red", cex = s1[pos, 1]*5)
points(dat$lon[pos], dat$lat[pos], pch = 1, col = "green", cex = s1[pos, 2]*5)

pred_dat = data.frame(year = 2016, quarter = 1, 
                      lon = dat$lon[pos[1:5]], lat = dat$lat[pos[1:5]],
                      trop_catch = 0, sst = 0)
pred_dat = pred_dat %>% mutate(fyear = as.factor(year), quarter = as.factor(quarter))
s2
points(dat$lon[pos[1:5]], dat$lat[pos[1:5]], pch = 1, col = "blue", cex = s2[, 1]*2)
