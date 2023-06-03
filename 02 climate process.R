
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, colourpicker, sf, glue, rgeos, gtools, tidyverse, rnaturalearthdata, rnaturalearth, dismo, raster)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions ---------------------------------------------------------------
cumTemp <- function(x) {
  
  p <- matrix(nrow = 1, ncol = 4)
  colnames(p) <- paste('bio', 21:24, sep = '')
  
  w <- x[25:36] ### tmax
  y <- x[13:24] ### tmean
  x <- x[1:12]  ### Prec-PET
  z <- x
  
  ### if the values are NA the bios are NA
  if(all(is.na(x))) {
    p[,'bio21'] <- NA
    p[,'bio22'] <- NA
    p[,'bio23'] <- NA
    p[,'bio24'] <- NA
  } else {
    
    ## cumulative deficit to determine dry season (=Bio22)
    
    # print('Bio 22...')
    
    x <- z
    lng <- length(x)
    x <- c(x, x[1:12])
    x[x>0] <- NA
    cumdef <- matrix(ncol = 12, nrow = lng)
    for (i in 1:12) {
      cumdef[, i] <- x[i:(lng + i - 1)]
    }
    p[,'bio22'] <- min(c(0,apply(cumdef, MARGIN = 1, FUN = cumsum)),na.rm=T)
    
    ## cumulative surplus to determine growing season
    x <- z
    lng <- length(x)
    x <- c(z, z[1:12])
    x[x<0] <- NA
    cumplus <- matrix(ncol = 12, nrow = lng)
    
    for (i in 1:12) {
      
      cumplus[, i] <- x[i:(lng + i - 1)]
      
    }
    
    ### If there is no dry season
    ### the length becomes 0
    ### the growing season temp is the mean of monthly mean temp
    ### the dry season max temp is the max temp of the driest month 
    
    if(p[,'bio22']==0){
      
      p[,'bio21'] <- 0
      p[,'bio23'] <- mean(y)
      p[,'bio24'] <- w[which.min(z)]
      
    } else {
      
      ### the mean temperatures for all possible seasons
      y <- c(y, y[1:12])
      n <- matrix(ncol = 12, nrow = lng)
      for (i in 1:12) {
        
        n[, i] <- y[i:(lng + i - 1)]
        
      }
      
      meantemp <- apply(n, MARGIN = 1, FUN = cumsum)
      
      ### the max temperatures for all possible seasons
      w <- c(w, w[1:12])
      n <- matrix(ncol = 12, nrow = lng)
      
      for (i in 1:12) {
        
        n[, i] <- w[i:(lng + i - 1)]
        
      }
      maxtemp <- apply(n, MARGIN = 1, FUN = cumsum)
      
      ### Consecutive months with Prec<PET (=bio21)
      x <- z
      x <- c(x, x[1:12])
      x[x>0] <- NA
      x[x<0] <- 1
      o <- matrix(ncol = 12, nrow = lng)
      
      for (i in 1:12) {
        
        o[, i] <- x[i:(lng + i - 1)]
        
      }
      
      con_months <- max(apply(o,1,cumsum),na.rm=T)
      p[,'bio21'] <- con_months
      
      ### if the dry season is 12 months the growing season mean is the mean of the wettest month
      
      if(con_months==12){
        
        p[,'bio23'] <- y[which.max(z)]
        
      } else { 
        
        ### The meantemp of the wettest season
        p[,'bio23'] <- meantemp[which.max(apply(cumplus, MARGIN = 1, FUN = cumsum))]/(12-con_months)
        
      }
      ### The mean maxtemp of the driest season
      
      p[,'bio24'] <- maxtemp[which.min(apply(cumdef, MARGIN = 1, FUN = cumsum))]/con_months    
      
    }
    
  }
  
  return(p)
  
}
etpvars <- function(x){
  
  p <- matrix(nrow = 1, ncol = 9)
  colnames(p) = paste('bio', 25:33, sep = '')
  
  tavg <- x[25:36] ### Temp
  prec <- x[13:24] ### PREC
  pet <- x[1:12]  ### PET
  
  ### if the values are NA the bios are NA
  if(all(is.na(x))) { 
    return(p)
  } else {
    
    window <- function(x)  { 
      lng <- length(x)
      x <- c(x,  x[1:3])
      m <- matrix(ncol = 3, nrow = lng)
      for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
      apply(m, MARGIN = 1, FUN = sum)
    }
    
    ### BIO_23: Annual PET
    p[,1] <- sum(pet)
    ### BIO_24: PET seasonality (Coefficient of Variation)
    p[,2] <- cv(pet)
    ### BIO_25: MAX PET
    p[,3] <- max(pet)
    ### BIO_26: Min PET
    p[,4] <- min(pet)
    ### BIO_27: Range of PET (PETmax-PETmin)
    p[,5] <- p[,3]-p[,4]
    
    wet <- window(prec)
    hot <- window(tavg)/3
    pet2 <- c(pet,pet[1:2])
    
    ### BIO_28: PET of wettest quarter
    p[,6] <- sum(pet2[c(which.max(wet):(which.max(wet)+2))])
    ### BIO_29:	PET of driest quarter
    p[,7] <- sum(pet2[c(which.min(wet):(which.min(wet)+2))])
    ### BIO_30:	PET of warmest quarter
    p[,8] <- sum(pet2[c(which.max(hot):(which.max(hot)+2))])
    ### BIO_31:	PET of coldest quarter
    p[,9] <- sum(pet2[c(which.min(hot):(which.min(hot)+2))])
    
  }
  
  round(p, digits = 2)
  return(p)
  
} 

# Load data ---------------------------------------------------------------
wrld <- ne_countries(returnclass = 'sf', scale = 50)
pnts <- suppressMessages(read_csv('tbl/points/points_v1.csv'))

# Baseline ----------------------------------------------------------------
bsln <- as.character(dir_ls('tif/climate/westafrica_baseline', regexp = '.tif$'))

isos <- c('GHA', 'CIV')
zone <- filter(wrld, sov_a3 %in% isos)
st_write(zone, 'gpkg/westafrica.gpkg')

# Read as rasters
prec <- rast(grep('prec', bsln, value = T))
tmax <- rast(grep('tmax', bsln, value = T))
tmin <- rast(grep('tmin', bsln, value = T))
tavg <- rast(grep('tavg', bsln, value = T))
etps <- rast(grep('etps', bsln, value = T))

# To extract by nmask  ----------------------------------------------------
prec <- terra::mask(terra::crop(prec, vect(zone)), vect(zone))
tmax <- terra::mask(terra::crop(tmax, vect(zone)), vect(zone))
tmin <- terra::mask(terra::crop(tmin, vect(zone)), vect(zone))
tavg <- terra::mask(terra::crop(tavg, vect(zone)), vect(zone))
etps <- terra::mask(terra::crop(etps, vect(zone)), vect(zone))

# Clean NAs ---------------------------------------------------------------
stck <- c(prec, tmax, tmin, tavg, etps) 
rstr <- terra::as.data.frame(stck, xy = T) %>% as_tibble() %>% drop_na()
rstr <- terra::rast(rstr, type = 'xyz')

prec <- rstr[[grep('prec', names(rstr), value = FALSE)]]
tmax <- rstr[[grep('tmax', names(rstr), value = FALSE)]]
tmin <- rstr[[grep('tmin', names(rstr), value = FALSE)]]
tavg <- rstr[[grep('tavg', names(rstr), value = FALSE)]]
etps <- rstr[[grep('etps', names(rstr), value = FALSE)]]

# To create bioclimatic variables -----------------------------------------
bioc <- rast(dismo::biovars(prec = stack(prec), tmax = stack(tmax), tmin = stack(tmin)))
terra::writeRaster(x = bioc, filename = 'tif/climate/westafrica_baseline/bioc_bsln.tif', overwrite = TRUE)

# Filtering the study zone ------------------------------------------------

# To extract the value  ---------------------------------------------------
plot(st_geometry(zone))
points(pnts$x, pnts$y, pch = 16, col = 'red')

vles <- terra::extract(bioc, pnts[,c('x', 'y')])
pnts <- cbind(pnts[,c(1, 2, 3, 4)], vles[,2:ncol(vles)])
pnts <- as_tibble(pnts)
write.csv(pnts, 'tbl/points/points_swd_all.csv', row.names = FALSE)

# To make a sample --------------------------------------------------------
subs <- sample_n(tbl = pnts, size = nrow(pnts) * 0.1, replace = FALSE)

write.csv(subs[,1:4], 'tbl/points/points_swd_sub.csv', row.names = FALSE)

# To make a map  ----------------------------------------------------------
dout <- 'png/maps'; dir_create(dout)
g1 <- ggplot() +
  geom_point(data = pnts, aes(x = x, y = y, col = 'red'), size = 0.1) + 
  geom_sf(data = wrld, fill = NA, col = 'grey80') +
  geom_sf(data = zone, fill = NA, col = 'grey30') + 
  coord_sf(xlim = ext(zone)[1:2], ylim = ext(zone)[3:4]) + 
  labs(x = 'Lon', y = 'Lat') +
  theme_minimal() +
  theme(axis.text.x = element_text(face = 'bold', hjust = 0.5, size = 7), 
        axis.text.y = element_text(face = 'bold', hjust = 0.5, size = 7, angle = 90), 
        axis.title.x = element_text(face = 'bold'), 
        axis.title.y = element_text(face = 'bold'), 
        legend.position = 'none',
        text = element_text(family = 'Barlow'))

g2 <- ggplot() +
  geom_point(data = subs, aes(x = x, y = y, col = 'red'), size = 0.1) + 
  geom_sf(data = wrld, fill = NA, col = 'grey80') +
  geom_sf(data = zone, fill = NA, col = 'grey30') + 
  coord_sf(xlim = ext(zone)[1:2], ylim = ext(zone)[3:4]) + 
  labs(x = 'Lon', y = 'Lat') +
  theme_minimal() +
  theme(axis.text.x = element_text(face = 'bold', hjust = 0.5, size = 7), 
        axis.text.y = element_text(face = 'bold', hjust = 0.5, size = 7, angle = 90), 
        axis.title.x = element_text(face = 'bold'), 
        axis.title.y = element_text(face = 'bold'), 
        legend.position = 'none',
        text = element_text(family = 'Barlow'))


ggsave(plot = g1, filename = glue('{dout}/presences_v1.png'), units = 'in', width = 9, height = 7, dpi = 300)
ggsave(plot = g2, filename = glue('{dout}/presences_v2.png'), units = 'in', width = 9, height = 7, dpi = 300)

# Create bioclimatic variables --------------------------------------------

# Bio 20  -----------------------------------------------------------------
prec.rclf <- terra::rast(raster::reclassify(stack(prec), c(-Inf, 100, 1, 100, Inf, NA)))
prec.two  <- c(prec.rclf, prec.rclf)
prec.two  <- stack(prec.two)
allperiods <- stack()
for(i in 1:12){
  oneyear <- prec.two[[i:(i+11)]]
  drymonths <- cumsum(oneyear)
  maxnumber <- max(drymonths, na.rm = T)
  allperiods <- addLayer(allperiods, maxnumber) 
}
bio20 <- max(allperiods, na.rm = T)
bio20[is.na(bio20)] <- 0
bio20 <- rast(bio20) %>% mask(., vect(zone))
msk <- bio20 * 0 
terra::writeRaster(x = bio20, filename = glue('tif/climate/westafrica_baseline/bioc20_bsln.tif'), overwrite = TRUE)

# Balance variables -------------------------------------------------------
defc <- prec - etps
DefAndTemp <- cbind(as.matrix(defc), as.matrix(tavg), as.matrix(tmax))
biovalues  <- t(apply(DefAndTemp, 1, cumTemp))
nms <- paste0('bioc', 30:33)
map(.x = 1:length(nms), .f = function(i){
  mask <- bio20 * 0 
  terra::values(mask) <- biovalues[,i]
  terra::writeRaster(x = mask, filename = glue('tif/climate/westafrica_baseline/{nms[i]}_bsln.tif'), overwrite = TRUE)
})

# ETP Bioclimatic ---------------------------------------------------------
names(etps) <- paste0('etp_', 1:12)
names(prec) <- paste0('prec_', 1:12)
names(tavg) <- paste0('tmean_', 1:12)

ETPAndPrec <- cbind(as.matrix(etps),as.matrix(prec),as.matrix(tavg))
etpbios    <- t(apply(ETPAndPrec, 1, etpvars))
nms <- paste0('bioc', 21:29)
map(.x = 1:ncol(etpbios), .f = function(i){
  print(i)
  rsl <- msk
  terra::values(rsl) <- etpbios[,i]
  terra::writeRaster(x = rsl, filename = glue('tif/climate/westafrica_baseline/{nms[i]}_bsln.tif'), overwrite = TRUE)
})

# Stacking all variables into only one raster -----------------------------
fles <- dir_ls('tif/climate/westafrica_baseline', regexp = '.tif$') %>% as.character() %>% mixedsort() %>% grep('bioc', ., value = T)
fles <- fles[c(15, 1:14)]
rstr <- terra::rast(fles)
names(rstr) <- glue('bioc_{1:33}')
terra::writeRaster(x = rstr, filename = 'tif/climate/westafrica_baseline/bioc_allv.tif', overwrite = T)



