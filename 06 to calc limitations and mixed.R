

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, colourpicker, pROC, pvclust, RColorBrewer, colourpicker, dendextend, randomForest, ggspatial, ggdendro, usdm, sf, glue, rgeos, gtools, tidyverse, rnaturalearthdata, rnaturalearth, dismo, raster)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

# Spatial data
zone <- st_read('gpkg/westafrica.gpkg') %>% filter(sov_a3 %in% c('GHA', 'CIV'))
wrld <- ne_countries(returnclass = 'sf', scale = 50)

# Random forest results
fles <- dir_ls('rf/output/run_1/results/raw', regexp = '.asc$')
clst <- grep('Clust', fles, value = T) %>% terra::rast()
prob <- grep('Prob', fles, value = T) %>% terra::rast()
uncr <- grep('Unc', fles, value = T) %>% terra::rast()

# Load presences
load(file = 'rData/run_1/clustereddata.rData')
occr <- clusteredpresdata
rm(clusteredpresdata, labelRF, prsc)

# Extract the values for the prob and uncertainty -------------------------
occr <- as_tibble(cbind(occr, prob = terra::extract(prob, occr[,c('x', 'y')])[,2], uncr = terra::extract(uncr, occr[,c('x', 'y')])[,2]))
qntl.prob <- setNames(rownames_to_column(as.data.frame(quantile(pull(occr, prob), seq(0, 1, 0.05)))), c('percentile', 'value'))
qntl.uncr <- setNames(rownames_to_column(as.data.frame(quantile(pull(occr, uncr), seq(0, 1, 0.05)))), c('percentile', 'value'))
thrs.prob <- filter(qntl.prob, percentile == '5%') %>% pull(2)
thrs.uncr <- filter(qntl.uncr, percentile == '5%') %>% pull(2)

# To reclassify the raster ------------------------------------------------

# Limitations -------------------------------------------------------------
prob.rclf <- terra::ifel(prob > thrs.prob, 2, 0)
clst.rclf <- terra::ifel(clst < 2.5, 0, 1)

# To make the difference 
dffr <- prob.rclf - clst.rclf
rslt <- clst
rslt[dffr == -1] <- 8
rslt[dffr == 2] <- 8

terra::writeRaster(x = rslt, filename = 'rf/output/run_1/results/process/RF5_Clust_lim_bsl.tif', overwrite = TRUE)

# Mixed -------------------------------------------------------------------
fnal <- rslt
fnal[uncr < thrs.uncr & prob > thrs.prob] <- 9
terra::writeRaster(x = fnal, filename = 'rf/output/run_1/results/process/RF5_clust_unc_bsl.tif', overwrite = T)

# To make the map  --------------------------------------------------------
tble <- terra::as.data.frame(fnal, xy = T)
tble <- as_tibble(tble)
tble <- setNames(tble, c('x', 'y', 'value'))

lbls <- tibble(value = 1:9, class = c(rep('Unsuit', 2), 'Moderate - Very dry', 'Hot - Dry', 'Hot - Wet', 'Cold - Wet', 'Very hot - Very wet', 'Limitations', 'Mixed'))
tble <- inner_join(tble, lbls, by = 'value')
tble <- mutate(tble, class = factor(class, levels = unique(lbls$class)))

colourWidget()

g1 <- ggplot() + 
  geom_tile(data = tble, aes(x = x , y = y, fill = class)) + 
  scale_fill_manual(values = c('#F5F0F0', '#7A772A', '#DE8E16', '#3C8548', '#14AFC7', '#D929D0', '#878787', '#D1C98A')) + 
  geom_sf(data = zone, fill = NA, col = 'grey30') + 
  geom_sf(data = wrld, fill = NA, col = 'grey70') +
  coord_sf(xlim = ext(zone)[1:2], ylim = ext(zone)[3:4]) + 
  ggtitle(label = "Agroclimatic zones for cocoa in Ghana and Côte dIvoire") +
  labs(x = 'Lon', y = 'Lat', fill = 'AEZ') + 
  theme_minimal() + 
  theme(legend.position = 'bottom', 
        plot.title = element_text(face = 'bold', hjust = 0.5),
        text = element_text(family = 'Barlow')) + 
  annotation_scale(location =  "bl", width_hint = 0.5, text_family = 'Barlow', text_col = 'grey60', bar_cols = c('grey60', 'grey99'), line_width = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_family = 'Barlow', text_col = 'grey40', line_col = 'grey60', fill = c('grey60', 'grey99'))) 

ggsave(plot = g1, filename = 'png/maps/cocoa_rf_bsl_run1.png', units = 'in', width = 9, height = 8, dpi = 300)
  

