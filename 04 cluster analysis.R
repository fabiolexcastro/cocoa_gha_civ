
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, colourpicker, pvclust, RColorBrewer, dendextend, randomForest, ggdendro, usdm, sf, glue, rgeos, gtools, tidyverse, rnaturalearthdata, rnaturalearth, dismo, raster)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions ---------------------------------------------------------------
source('functionsRFclustering.R')
source('dendogram.R')
rf.clust <- function(occ, nforest, ntrees, nVars, nclasses){
  
  datRF_presences <- occ
  print(nrow(datRF))
  
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
pnts <- suppressMessages(read_csv('tbl/points/points_back.csv'))
bioc <- terra::rast('tif/climate/westafrica_baseline/bioc_allv.tif')
vars <- dplyr::select(pnts, contains('bioc')) %>% colnames()
pnts <- drop_na(pnts)

# Clustering using random forest ------------------------------------------
prsc <- filter(pnts, pb == 1)
datRF <- as.data.frame(prsc[,4:ncol(prsc)])
d <- dist(datRF, method = 'euclidean')
no.clusters <- 5
rfClust <- rf.clust(occ = prsc[,4:ncol(prsc)], nforest = 25, ntrees = 100, nVars = 5, nclasses = 5)
labelRF <- rfClust[[1]]
clusterdata <- rfClust[[2]]

classdata <- cbind(pb = as.factor(labelRF), prsc[,3:ncol(prsc)])
clusteredpresdata <- cbind(prsc, cluster = labelRF) %>% na.omit() %>% tbl_df()

# To save the objects 
dout <- 'rData/run_1'
save(datRF, file = paste0(dout, '/datRF.rData'))
save(clusterdata, file = paste0(dout, '/clusterdata.rData'))
save(prsc, clusteredpresdata, no.clusters, labelRF, file = paste0(dout, '/clustereddata.rData'))

# Dendogram way 1 ---------------------------------------------------------
clusterdata
dend <- as.dendrogram(clusterdata) %>%
  set("branches_k_color", k = 5) %>% 
  set("branches_lwd", 0.8) %>%
  set("labels_colors", brewer.pal(n = 5, name = 'Set2')) %>% 
  set("labels_cex", c(.9,1.2)) %>% 
  set("leaves_pch", 19) %>% 
  set("leaves_col", c("blue", "red")) 
ggdn <- ggplot(as.ggdend(dend)) + 
  # scale_color_manual(values = brewer.pal(n = 5, name = 'Set2')) +
  theme_minimal() + 
  ggtitle(label = 'Dendogram') + 
  theme(axis.text.x = element_blank(), 
        axis.title = element_blank(), 
        axis.text.y = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5))
ggsave(plot = ggdn, filename = 'png/graphs/dend_run1.png', units = 'in', width = 9, height = 7, dpi = 300)

# Dendogram way 2 ---------------------------------------------------------
clusterdata
ED.clusters <- 5

## Make a nice Dendrogram
ppar          <- par(no.readonly=TRUE)
clusterlabels <- paste("Cluster",1:ED.clusters)

png(filename = 'png/graphs/dend_run1_v2.png', width = 9, height = 6, units = "in", res = 300)
opar <- par(no.readonly=TRUE)
A2Rplot.hclust(clusterdata, k = ED.clusters, boxes = FALSE, col.up = "gray50", #clusterdata
               col.down = c("#4bd14b", "#ce53ed", "#dbc24f","#a87000","#57debe"),
               #        main="Dendrogram of agglomerative clustering",
               main=NULL,
               ylab="Height",
               mtext(seq(0, 1000000, 10000), side = 2, at = seq(0, 1000000, 10000), # seq(0,10000,5000), at = seq(0, 5000, 1000),
                     line = 1, las = 1),
               hang=-1,axes=T,show.labels=F,lty.up = 2,lty.down = 1)
axis(side = 2, at = seq(0, 5000, 1000), labels = F, lwd = 1, line = -4.5, tick = T) # at = seq(0, 2, 0.36), at = seq(0, 5000, 1000)
title(ylab = "RF Clustering", line = -2.5)# 'Euclidenn distance'
par(opar)

# legend coordinates depend on Euclidean distance (y) and number of cases (x)
legend("topright", legend = c('Type 1', 'Type 2', 'Type 3', 'Type 4', 'Type 5'), col = c("#4bd14b", "#ce53ed", "#dbc24f", '#a87000', '#57debe'),
       lty = 1, lwd = 2) # legend = clusterlabels, y.intersp = 1.1
## Dendrogram Ende
dev.off()



