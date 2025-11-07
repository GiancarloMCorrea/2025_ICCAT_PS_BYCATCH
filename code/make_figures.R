require(gt)

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

# -------------------------------------------------------------------------
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
) %>% gtsave(filename = file.path('figures/3_sampsize/tab_eff.png'))

# Make simulated data:
mydat2 %>% gt %>% data_color(
  columns = Trip_id,
  target_columns = everything(),
  palette = colpal
) %>% tab_style(
  style = cell_text(weight = "bold"),
  locations = cells_column_labels(columns = everything())
) %>% gtsave(filename = file.path('figures/3_sampsize/tab_sim.png'))

# Make sampled data:
mydat3 %>% gt %>% data_color(
  columns = Trip_id,
  target_columns = everything(),
  palette = colpal[c(1,3)]
) %>% tab_style(
  style = cell_text(weight = "bold"),
  locations = cells_column_labels(columns = everything())
) %>% gtsave(filename = file.path('figures/3_sampsize/tab_samp.png'))


# -------------------------------------------------------------------------
# Make flowchart:
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
mg <- mermaid("graph LR;
    A--> B;
    B--> A;
    C--> A;
    D--> A;
    E--> B;")
mg <- mermaid("graph TD
  DIR('<img src='https://iconscout.com/ms-icon-310x310.png'; width='30' />')
              ")
# Save:
DPI = 500
WidthCM = 17
HeightCM = 10
mg %>% export_svg %>% charToRaw %>% 
  rsvg(width = WidthCM *(DPI/2.54), height = HeightCM *(DPI/2.54)) %>% 
  png::writePNG('figures/3_sampsize/flowchart.png')
