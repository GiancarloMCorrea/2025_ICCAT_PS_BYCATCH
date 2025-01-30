# These parameters will be used in all R scripts for plotting

# For maps, define limits:
xLim = c(-55, 15)
yLim = c(-25, 30)

# Read land map information:
worldmap = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Define breaks in maps:
xBreaks = seq(from = -60, to = 20, by = 20)
yBreaks = seq(from = -40, to = 40, by = 20)

# Define image parameters:
img_type = '.png'
img_res = 300 # dpi
img_width = 170 # max width image

# Theme for ggplot:
theme_set(theme_bw())

# Folder in sharepoint where to save plots:
# plot_dir = 'figures'