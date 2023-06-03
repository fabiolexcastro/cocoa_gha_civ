
# Data source: https://nk.users.earthengine.app/view/cocoa-map

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, colourpicker, sf, tidyverse, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
wrld <- ne_countries(returnclass = 'sf', scale = 50)
rstr <- rast('tif/cocoa_map_th_v1.tif')
prob <- rast('tif/cocoa_map.tif')
zone <- filter(wrld, sov_a3 %in% c('GHA', 'CIV'))
zone <- vect(zone)

# To extract by mask  -----------------------------------------------------
rstr <- terra::crop(rstr, zone)
rstr <- terra::mask(rstr, zone)
zone <- st_as_sf(zone)

# Raster to table ---------------------------------------------------------
tble <- as_tibble(terra::as.data.frame(rstr, xy = T))
tble <- mutate(tble, class = 'Cocoa')

# To make the map  --------------------------------------------------------
g1 <- ggplot() + 
  geom_tile(data = tble, aes(x = x, y = y, fill = class)) + 
  scale_fill_manual(values = '#753F13') +
  geom_sf(data = wrld, fill = NA, col = 'grey50') + 
  geom_sf(data = zone, fill = NA, col = 'grey20') + 
  coord_sf(xlim = ext(zone)[1:2], ylim = ext(zone)[3:4]) + 
  labs(fill = '') +
  theme_void() + 
  theme(legend.position = 'bottom')

ggsave(plot = g1, filename = 'png/maps/cocoa_th_5km.png', units = 'in', width = 9, height = 7, dpi = 300)

write.csv(tble, 'tbl/points/points_v1.csv', row.names = FALSE)

