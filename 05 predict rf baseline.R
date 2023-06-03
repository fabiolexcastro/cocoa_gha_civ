

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, colourpicker, pROC, pvclust, RColorBrewer, dendextend, randomForest, ggdendro, usdm, sf, glue, rgeos, gtools, tidyverse, rnaturalearthdata, rnaturalearth, dismo, raster)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions to use --------------------------------------------------------
source('functionsRFclustering.R')
rf.clust <- function(occ, nforest, ntrees, nVars, nclasses){
  # occ = back_swd; nforest = 50; ntrees = 500; nVars = 8; nclasses = 2
  datRF_presences <- occ[,3:ncol(occ)] %>% as.data.frame()
  print(nrow(datRF_presences))
  attach(datRF_presences)
  no.forests <- nforest
  no.trees <- ntrees
  distRF_presences <- RFdist(datRF_presences, mtry1 = nVars, no.trees, no.forests, addcl1 = T, addcl2 = F, imp = T, oob.prox1 = T)
  no.presencesclasses <- nclasses
  labelRF <- pamNew(distRF_presences$cl1, no.presencesclasses)
  print(table(labelRF))
  clusterdata <- hclust(as.dist(distRF_presences$cl1), method = 'ward.D2')
  return(list(labelRF, clusterdata))
}

# Load data ---------------------------------------------------------------
clma <- rast('tif/climate/westafrica_baseline/bioc_allv.tif')
vars <- readRDS(file = 'rds/run_1/vars_vif.rds')
zone <- st_read('gpkg/westafrica.gpkg') %>% filter(sov_a3 %in% c('GHA', 'CIV'))

# Load clustered data
load(file = 'rData/run_1/clustereddata.rData')

# Bias raster - process ---------------------------------------------------
back <- clma[[1]] * 0 + 1
names(back) <- 'mask'
cell <- terra::extract(back, clusteredpresdata[,c('x', 'y')], cell = TRUE)$cell
back[cell] <- NA
samplesize <- round(min(summary(as.factor(clusteredpresdata$cluster))) / 2, 0) 
NumberOfClusters <- 5
numberofpresences <- nrow(clusteredpresdata) 
back <- randomPoints(raster(back), 1*numberofpresences) %>% as_tibble()
back_swd  <- raster::extract(clma, back[,c(1, 2)]) %>% cbind(back, .)
nrow(back_swd) == nrow(back_swd[complete.cases(back_swd),])
back_swd <- dplyr::select(back_swd, x, y, vars)
write.csv(back_swd, 'tbl/points/run_1/back_swd.csv', row.names = FALSE)
write.csv(prsc, 'tbl/points/run_1/occ_swd.csv', row.names = FALSE)

# Cluster analysis to pseudoabsences
bckclust <- rf.clust(occ = back_swd, nforest = 50, ntrees = 500, nVars = 8, nclasses = 2)
datRF <- as.data.frame(back_swd[,3:ncol(back_swd)])
attach(datRF)
no.forests <- 50#raw = 25
no.trees <- 500
distRF <- RFdist(datRF, mtry1 = 8, no.trees, no.forests, addcl1 = T, addcl2 = F, imp =T, oob.prox1 = T)# mtry1 = 4 raw  # es la cantidad de variables a utilizar en cada no
no.absenceclasses <- 2
labelRF <- pamNew(distRF$cl1, no.absenceclasses)
detach(datRF)
classdata <- cbind(pb = as.factor(labelRF), back_swd[,3:ncol(back_swd)])

no.absenceclasses <- 2

presvalue_swd  <- clusteredpresdata[,3:ncol(clusteredpresdata)] %>%
  mutate(pb = cluster + no.absenceclasses) %>% 
  na.omit() %>%
  as.data.frame() %>%
  mutate(cluster = cluster + no.absenceclasses)
presvalue_swd <- dplyr::select(presvalue_swd, pb, vars)
presvalue_swd <- mutate(presvalue_swd, pb = as.factor(pb))
classdata_2 <- cbind(pb = as.data.frame(classdata)$pb, classdata[,2:ncol(classdata)]) # Background

dim(classdata_2); dim(presvalue_swd)
presvalue_swd <- presvalue_swd %>% dplyr::select(-cluster)

allclasses_swd <- rbind(classdata_2, presvalue_swd[,1:ncol(classdata_2)])
unique(allclasses_swd$pb)
write.csv(allclasses_swd, 'tbl/points/run_1/all_classes_swd.csv', row.names = FALSE)

# To make the random forest analysis --------------------------------------
vrs <- vars
vrs <- gsub('.asc', '', vrs) 
vrs <- gsub('\\$', '', vrs)
model1 <- as.formula(paste('factor(pb) ~', paste(paste(vrs), collapse = '+', sep =' ')))
rflist <- vector('list', 50) 
auc <- vector('list', 50)

for(repe in 1:50){ # 50 bosques
  
  print(repe)
  pressample <- list()
  
  for (i in 1:(NumberOfClusters+no.absenceclasses)){
    
    if(any(i==c(1:no.absenceclasses))) { 
      
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), 
                     size = samplesize*NumberOfClusters/2/no.absenceclasses)
    } else {
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), size=samplesize)
    }
    pressample[[i]] <- allclasses_swd[rows,] 
  }
  
  species <- na.omit(do.call(rbind, pressample)) 
  head(species)
  Samplesplit <- sample(rownames(species)) 
  
  envtrain <- species[Samplesplit[1:(0.8*nrow(species))],] 
  envtest <- species[Samplesplit[(0.8*nrow(species)):nrow(species)],] 
  
  rfmodel <- randomForest(model1, data = envtrain, ntree = 500, na.action = na.omit, nodesize = 2) 
  
  save(rfmodel, file = paste('rf/output/run_1/models/rf', NumberOfClusters, 'Prob_' , 'rep_' ,repe, '.rdata' ,sep=''))
  rflist[[repe]] <- rfmodel
  
  # AUC 
  predicted <- as.numeric(predict(rfmodel, envtest))
  observed <- as.vector(envtest[,'pb'])
  auc[[repe]] <- auc(observed, predicted) 
  rm(rfmodel)
  
  cat(auc[[repe]] ,'\n')
  
}

auc <- unlist(auc)
rff <- do.call(randomForest::combine, rflist)
importance <- as.data.frame(rff$importance)

run <- 'run_1'

save(rflist, file = paste('./rData/', run, '/rflist_', NumberOfClusters, '.rdata', sep = ''))
save(importance, file = paste0('./rData/', run, '/importanceRF.rData'))
save(auc, file = paste0('./rData/', run, '/aucRF_dist.rData'))
save(rff, file = paste0('./rData/', run, '/rff_dist.rData'))

# Predict modell
clma <- clma[[grep(paste0(paste0(vars, '$'), collapse = '|'), names(clma), value = F)]]
climatevalues  <- data.frame(values(clma))
NumberOfClusters <- 5

rasterProbs <- predict(rff, climatevalues, type = 'prob') # proximity = T
rasterProbs_na <- na.omit(rasterProbs)
sum_rasterProbs_na <- apply(rasterProbs_na, 1, sum)

rasterRF <- rowSums(rasterProbs[,c(3:(NumberOfClusters+2))])
uncertainty <- apply(rasterProbs, 1, max)  

rasterRFprob <- clma[[1]]
values(rasterRFprob) <- rasterRF 

rasterRFuncertainty <- clma[[1]]
values(rasterRFuncertainty) <- uncertainty 

rasterRF <- max.col(rasterProbs, 'first')
rasterRFclass <- clma[[1]]
values(rasterRFclass) <- rasterRF

terra::writeRaster(rasterRFclass, paste0('rf/output/run_1/results/raw/RF_5Clust_current.asc'), overwrite = T)
terra::writeRaster(rasterRFprob, paste0('rf/output/run_1/results/raw/RF_5Prob_current.asc'), overwrite = T)
terra::writeRaster(rasterRFuncertainty, paste0('rf/output/run_1/results/raw/RF_5Unc_current.asc'), overwrite = T)










