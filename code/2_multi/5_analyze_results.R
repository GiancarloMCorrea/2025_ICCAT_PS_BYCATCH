rm(list = ls())

# Define type of school to be analyzed:
this_type = 'FOB' 
source('code/2_multi/load_libs.R')
sel_nfac = '2'

# -------------------------------------------------------------------------
# Make main plot for each component:

sel_comp = '1'
# For component 1 Load rotated loading matrices:
rot_load = readRDS(file.path(model_folder, sel_nfac, paste0('rotated_loading', sel_comp,'.rds')))
rot_load = rot_load[,c(1,2,ncol(rot_load))]
axis_names = colnames(rot_load)[1:2]
colnames(rot_load) = c('fac1', 'fac2', 'species')
p1 = ggplot(data = rot_load, aes(x = fac1, y = fac2, label = species)) +
  geom_point() +
  geom_text(hjust=0, vjust=0, size = 2.5) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  coord_cartesian(xlim = c(min(rot_load$fac1), max(rot_load$fac1) + 0.22), 
                  ylim = range(rot_load$fac2)) +
  labs(x = axis_names[1], y = axis_names[2]) 

# Plot cluster:
load(file.path(model_folder, sel_nfac, paste0('loading', sel_comp,'.RData')))
rownames(Lhat_cf) = rot_load$species
Dist = dist(Lhat_cf, diag=TRUE, upper=TRUE)
hc = hclust(Dist)
pamk_best = pamk(Dist)
h1 = fviz_dend(hc, k = pamk_best$nc, 
          cex = 0.6, 
          labels_track_height = 0.25,
          main = NULL,
          rect = TRUE,
          k_colors = c("#2E9FDF","#00AFBB", "#E7B800" ,"#FC4E07"),
          color_labels_by_k = TRUE, 
          ggtheme = theme_classic() 
)

# Plot omegas:
omega = readRDS(file.path(model_folder, sel_nfac, paste0('omega', sel_comp,'.rds')))
plot_dat = omega %>% filter(type_fac == 'factor_1')
plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
m1 = ggplot(plot_dat) + geom_sf(aes(color = omega), size = 2) + 
  scale_colour_gradient2(low = muted("blue"), high = muted("red")) + 
  labs(title = expression(omega*' (Factor 1)'), color = NULL) +
  theme(legend.position = c(0.85, 0.78),
        legend.background = element_rect(fill = 'transparent', color = 'transparent')) 
m1 = add_sf_map(m1)
plot_dat = omega %>% filter(type_fac == 'factor_2')
plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
m2 = ggplot(plot_dat) + geom_sf(aes(color = omega), size = 2) + 
  scale_colour_gradient2(low = muted("blue"), high = muted("red")) + 
  labs(title = expression(omega*' (Factor 2)'), color = NULL) +
  theme(legend.position = c(0.85, 0.78),
        legend.background = element_rect(fill = 'transparent', color = 'transparent')) 
m2 = add_sf_map(m2)

# Merge plots:
p_comb = grid.arrange(p1, h1, m1, m2,  ncol = 2)
ggsave(filename = paste0('Figure_summ', sel_comp, img_type), path = plot_folder, plot = p_comb, 
       width = img_width, height = 170, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Make predictions by species:

# Load Grid:
load(file.path(data_folder, 'MyGrid.RData'))

# Read joint predictions:
joint_data = readRDS(file.path(model_folder, sel_nfac, 'predictions.rds'))
joint_data$type_model = 'Joint'

# Read single predictions:
list_sp = list.files(file.path(model_folder, 'single'))
save_pred = list()
for(i in 1:length(list_sp)) {
  save_pred[[i]] = readRDS(file.path(model_folder, 'single', list_sp[i], 'predictions.rds'))
}
save_pred = bind_rows(save_pred)
save_pred$type_model = 'Single'

# Make plot:
plot_data = bind_rows(joint_data, save_pred)
all_sp = unique(plot_data$sp_name)
save_plot = list()
for(k in 1:length(all_sp)) {
  tmp_dat = plot_data %>% filter(sp_name == all_sp[k])
  grid_dat = left_join(MyGrid, tmp_dat, by = 'ID')
  p1 = ggplot(grid_dat) + geom_sf(aes(color = pred, fill = pred)) + 
    scale_colour_viridis() + scale_fill_viridis() + 
    guides(color = 'none') + 
    theme_void() +
    labs(fill = NULL) +
    theme(legend.background = element_rect(fill = "transparent",
                                           color="transparent"),
          legend.position = c(0.15, 0.55),
          legend.key.height = unit(0.25, 'cm'),
          legend.key.width = unit(0.25, 'cm'),
          legend.text = element_text(size=5),
          strip.text.y = element_text(angle = 90)) 
  p1 = add_sf_map(p1)
  p1 = p1 + facet_grid(rows = vars(type_model), cols = vars(sp_name), switch = "y")
  save_plot[[k]] = p1
}

merged_1 = do.call("grid.arrange", c(save_plot[1:6], nrow = 1))
merged_2 = do.call("grid.arrange", c(save_plot[7:12], nrow = 1))
merged_3 = do.call("grid.arrange", c(save_plot[13:18], nrow = 1))
merged_plot = grid.arrange(merged_1, merged_2, merged_3, ncol = 1)
ggsave(filename = paste0('predicted', img_type), path = plot_folder, plot = merged_plot, 
       width = img_width, height = 180, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Combine loading matrices:
# 
# cp = 1
# n_sp = 17
# n_fac = 3
# theta_slot_vec = c('theta_z', 'theta2_z')
# theta_slot = theta_slot_vec[cp]
# 
# L1 = matrix(0, nrow = n_sp, ncol = n_fac)
# theta_vec = as.list(jtVModel$sdrep, what="Estimate")[['theta_z']]
# L1[lower.tri(L1, diag=TRUE)] = theta_vec
# plot(L1[,1], type = 'h', lwd = 4)
# abline(h = 0, lty = 2, col = 'grey')
# 
# L2 = matrix(0, nrow = n_sp, ncol = n_fac)
# theta_vec = as.list(jtVModel$sdrep, what="Estimate")[['theta2_z']]
# L2[lower.tri(L2, diag=TRUE)] = theta_vec
# lines(L2[,1], type = 'h', col = 2, lwd = 3)
# 
# V_tot = L1 %*% t(L1) + L2 %*% t(L2)
# out = svd(V_tot)
# Lhat_cf_rot = rotate_pca( L_tf = out$v[,1:2], order = 'decreasing' )$L_tf
# lines(Lhat_cf_rot[,1], type = 'h', col = 3, lwd = 2)


# -------------------------------------------------------------------------
# Make pretty map with species silhouettes:
sel_comp = '1'
# Get fish images:
s_barracuda = get_uuid(name = "Sphyraena", n = 1)
s_barracuda = get_phylopic(uuid = s_barracuda)
a_solandri = get_uuid(name = "Acanthocybium solandri", n = 1)
a_solandri = get_phylopic(uuid = a_solandri)
sphyrnidae = get_uuid(name = "Sphyrnidae", n = 1)
sphyrnidae = get_phylopic(uuid = sphyrnidae)
p_glauca = get_uuid(name = "Prionace glauca", n = 1)
p_glauca = get_phylopic(uuid = p_glauca)
chelonioidea = get_uuid(name = "Chelonioidea", n = 1)
chelonioidea = get_phylopic(uuid = chelonioidea)
molidae = get_uuid(name = "Molidae", n = 1)
molidae = get_phylopic(uuid = molidae)
lamnidae = get_uuid(name = "Lamnidae", n = 1)
lamnidae = get_phylopic(uuid = lamnidae)

# Read omega data:
omega = readRDS(file.path(model_folder, sel_nfac, paste0('omega', sel_comp,'.rds')))

# Make plot:
text_size = 1.5
alpha_img = 0.7
plot_dat = omega %>% filter(type_fac == 'factor_1') %>% 
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
m1 = ggplot(plot_dat) + geom_sf(aes(color = omega), size = 3) + 
  scale_colour_gradient2(low = alpha("blue", alpha=0.1), high = alpha("red", alpha=0.1)) + 
  theme_classic() +
  labs(color = NULL) + theme(legend.position = 'none') +
  geom_sf(data = worldmap, fill = "gray60", color = "gray60") +
  coord_sf(expand = FALSE, xlim = c(-40, 20), ylim = c(-25, 30)) +
  xlab(NULL) + ylab(NULL) +
  scale_y_continuous(breaks = yBreaks) + scale_x_continuous(breaks = xBreaks)
m1 = m1 + 
  add_phylopic(x = -10, y = 0, img = s_barracuda, width = 14, alpha = alpha_img) +
  add_phylopic(x = -20, y = 3, img = a_solandri, width = 14, alpha = alpha_img) +
  add_phylopic(x = 2, y = 23, img = sphyrnidae, width = 14, alpha = alpha_img) +
  add_phylopic(x = 6, y = 26, img = chelonioidea, width = 7, angle = 90, alpha = alpha_img) +
  add_phylopic(x = 6, y = 18, img = molidae, width = 7, alpha = alpha_img) +
  add_phylopic(x = 13, y = 23, img = p_glauca, width = 11, alpha = alpha_img) +
  add_phylopic(x = 13, y = 20, img = lamnidae, width = 11, alpha = alpha_img)
m1 = m1 + annotate("text", x = -10, y = 0, label = "italic(S.~barracuda)", 
                   size = text_size, color = "white", parse = TRUE) +
  annotate("text", x = -20, y = 3, label = "italic(A.~solandri)", 
           size = text_size, color = "white", parse = TRUE) +
  annotate("text", x = 2, y = 23, label = "Sphyrnidae", 
                   size = text_size, color = "white") +
  annotate("text", x = 6, y = 26, label = "Chelonioidea", 
                   size = text_size, color = "white") +
  annotate("text", x = 6, y = 18, label = "Molidae", 
                   size = text_size, color = "white") +
  annotate("text", x = 13, y = 23, label = "italic(P.~glauca)", 
                   size = text_size, color = "white", parse = TRUE) +
  annotate("text", x = 13, y = 20, label = "Lamnidae", 
                   size = text_size, color = "white")
m1 = m1 + geom_ellipse(aes(x0 = -15, y0 = 1, a = 14, b = 5, angle = 0)) +
          geom_ellipse(aes(x0 = 7, y0 = 21, a = 12, b = 10, angle = 0))
m1 = m1 + geom_segment(aes(x = -5, y = 20, xend = -20, yend = 20)) +
          geom_segment(aes(x = 8, y = 11, xend = 10, yend = -10)) +
          geom_point(aes(x = -20, y = 20), size = 2, color = 'black') + 
          geom_point(aes(x = 10, y = -10), size = 2, color = 'black') 
ggsave(filename = paste0('Figure_map_sp', sel_comp, img_type), path = plot_folder, plot = m1, 
       width = img_width*0.6, height = 90, units = 'mm', dpi = img_res)
