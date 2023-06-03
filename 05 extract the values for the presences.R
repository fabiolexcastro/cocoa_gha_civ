

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, colourpicker, pROC, pvclust, RColorBrewer, dendextend, randomForest, ggdendro, usdm, sf, glue, rgeos, gtools, tidyverse, rnaturalearthdata, rnaturalearth, dismo, raster)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

# Random forest results
fles <- dir_ls('rf/output/run_1/results/raw', regexp = '.asc$')
prob <- terra::rast(grep('Prob', fles, value = T))
clst <- terra::rast(grep('Clust', fles, value = T))
uncr <- terra::rast(grep('Unc', fles, value = T))

# Presence data
load(file = 'rData/run_1/clustereddata.rData')
prsc <- clusteredpresdata

# Climate data 
clma <- terra::rast('tif/climate/westafrica_baseline/bioc_allv.tif')
vars <- readRDS(file = 'rds/run_1/vars_vif.rds')
clma <- clma[[grep(paste0(paste0(vars, '$'), collapse = '|'), names(clma), value = T)]]

# To make the graph  ------------------------------------------------------
prsc <- mutate(prsc, cluster = factor(cluster, levels = as.character(1:5)))
dfrm <- prsc %>% 
  gather(var, value, -x, -y, -cluster, -pb) %>% 
  mutate(var = factor(var, levels = vars))

g1 <- ggplot(data = dfrm, aes(y = value, x = cluster)) + 
  facet_wrap(.~var, scales = 'free') +
  geom_boxplot() + 
  labs(x = '', y = 'Value') +
  ggtitle(label = 'Bioclimatic values for each cluster') + 
  theme_minimal() + 
  theme(strip.text = element_text(face = 'bold'), 
        axis.text.x = element_blank(), 
        plot.title = element_text(face = 'bold', hjust = 0.5),
        text = element_text(family = 'Barlow'))

ggsave(plot = g1, filename = 'png/graphs/boxplot_run1.png', units = 'in', width = 9, height = 9, dpi = 300)


