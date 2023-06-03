
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, colourpicker, usdm, sf, glue, rgeos, gtools, tidyverse, rnaturalearthdata, rnaturalearth, dismo, raster)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
bioc <- terra::rast('tif/climate/westafrica_baseline/bioc_allv.tif')
pnts <- read_csv('tbl/points/points_swd_sub.csv')

# Extract the values from all the points ----------------------------------
vles <- terra::extract(bioc, pnts[,c(1, 2)])
pnts <- cbind(pnts, vles[,2:ncol(vles)])
pnts <- as_tibble(pnts)
pnts <- mutate(pnts, pb = 1)

# To make the VIF  --------------------------------------------------------
mtrx <- pnts[,5:ncol(pnts)]
mtrx <- as.data.frame(mtrx)

vif.res <- vifstep(x = mtrx, th = 5)
vars <- as.character(vif.res@results$Variables)
dir.create('rds/run_1', recursive = T)
saveRDS(object = vars, file = 'rds/run_1/vars_vif.rds')

# Pseudo-absences ---------------------------------------------------------
mask <- bioc[[1]] * 0 + 1
cell <- terra::extract(mask, pnts[,1:2], cell = T)$cell
mask[cell] <- NA
back <- terra::as.data.frame(mask, xy = T)
back <- sample_n(tbl = back, size = nrow(pnts), replace = FALSE)
back <- back[,1:2]

# To extract the value for the presences 
back <- as_tibble(cbind(back, terra::extract(bioc, back[,1:2])))
back <- dplyr::select(back, -ID)
back <- mutate(back, pb = 0)
back <- dplyr::select(back, x, y, pb, vars)

# To join all the dataframes into only one  -------------------------------
pnts <- dplyr::select(pnts, x, y, pb, vars)
allp <- rbind(pnts, back)
write.csv(allp, './tbl/points/points_back.csv', row.names = FALSE)







