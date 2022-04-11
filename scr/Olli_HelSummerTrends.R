# Half-century trends in alpha and beta diversity of phytoplankton summer communities in the Helsinki Archipelago, the Baltic Sea
# Kalle Olli, Emil Nyman, Timo Tamminen
# 
# Olli_HelSummerTrends.R

#  Script summary
#  1. Load the main data files
#  2. Split the sampling stations to inner and outer arhipleago based on NMDS ordination
#  3. Map of the study area (Olli_HelSummerTrends_Figure1)
#  4. NMDS ordination of inner and outer archipelago, superimposed by envfit, kmeans clusters 
#  5. Diveristy profiles of communities
#  6. vegan::envfit profiles
#  7.Save NMDS ordination with envfit vectors: save Olli_HelSummerTrends_Figure2
#  8. NMDS1 and species richness trends, save Olli_HelSummerTrends_Figure3
#  9. Community evenness trends, save Olli_HelSummerTrends_FigureS6
# 10. Species cumulative likelihood of occurrence, saves file Olli_HelSummerTrends_Figure4
# 11. Most frequent species, Supplementary tables S2 and S3
# 12. Beta diversity calculations
# 13. Trends in temporal beta diversity; save Olli_HelSummerTrends_Figure5
# 14. Trends in spatial beta diversity; save Olli_HelSummerTrends_Figure6
# 15. Beta diverstity decomposition; save Olli_HelSummerTrends_Figure7
# 16. Distance decay; save Olli_HelSummerTrends_Figure8
# 17. GAM effect size, save Olli_HelSummerTrends_FigureS2.pdf
# 18. The effect of microscopists
# 19. Mantel and adonis block for Table 1 and Fig S3
# 20. Table 1
# 21. Hierarchical GAM, Figs S4, S5, Table S1
# 22. Nutrient block, saves Olli_HelSummerTrends_FigureS7.pdf

# load libraries
if(T){
  library(tidyverse)
  library(mgcv)
  library(nlme)
  library(vegan)
  library(adespatial) # beta.div (Var(Y) based beta: var.div.comp - decomposition of beta)
  library(proxy) # for row.dist and col.dist to get row and column indexes of dist object
  library(sp)
  library(sf)
  library(parallelDist)
  library(parallel)
  library(entropart) # entropart::Dqz
  library(lubridate) # lubridate::decimal_date
  library(cowplot)
  library(stargazer) # for visualizing tables in R
  
  # required external files:
  # ./Olli_HelSummerTrends_meta.txt
  # ./Olli_HelSummerTrends_species.txt'
  # ./Olli_HelSummerTrends_renyi2_renyi2.R'  # modified vegan::renyi function to allow negative orders of diversity
  # ./Olli_HelSummerTrends_bscl.txt # taxon classification according to WORMS
  # ./Olli_HelSummerTrends_bstree.txt # crude phylogeny from WORMS classification
  
}

## LOAD MAIN DATA TABLES ####

if(T){ # provides data tables dat, meta. bscl
  # load sample metadata
  meta <- read.table('Olli_HelSummerTrends_meta.txt', header = TRUE, sep = '\t') # 4630 samples
  
  # load long format species data
  dat <- read.table('Olli_HelSummerTrends_species.txt', header = TRUE, sep = '\t') # 155301 × 4
  
  # load worms classification table
  bscl <- read.table('./Olli_HelSummerTrends_bscl.txt', header = TRUE, sep = '\t')
  
  # add taxon name to dat
  dat <- left_join(dat, select(bscl, id, valid_name), by = c('valid_AphiaID' = 'id'))
  
  # make global community matrix - species in columns, samples in rows
  # this matrix is used for NMDS ordination, but only to split the sampling stations into inner bays and outer archipelago
  tmp <- group_by(dat, sampleID, valid_name) %>% summarise(ww = mean(ww, na.rm = T)) %>% pivot_wider(names_from = valid_name, values_from = ww)
  datcm  <- data.matrix(tmp[ , -1])# community matrix on dat 4630 samples × 620 taxa
  rownames(datcm) <- tmp[[1]] # add sampleID as rownames
  datcm[is.na(datcm)] <- 0 # change all NA's to zeros
  
  # supplement meta table with additional time and diversity variables
  meta <- mutate(meta,
                 obstime = as.Date(obstime), # convert data to Date format
                 year = as.numeric(format(obstime, '%Y')), # extract sampling year
                 mon = as.numeric(format(obstime, '%m')), # extract sampling month
                 jul = as.numeric(format(obstime, '%j')), # extract day of the year
                 N0 = rowSums(datcm[as.character(meta$sampleID), ] > 0),  # Species richness
                 H = vegan::diversity(datcm[as.character(meta$sampleID), ]), # Shannon entropy (base e)
                 N1 = exp(H), # Shannon diversity (base e)
                 N2 = vegan::diversity(datcm[as.character(meta$sampleID), ], "inv"), # Simpson diversity
                 J = H / log(N0), # Pielou evenness
                 E10 = N1 / N0,   # Shannon evenness (Hill's ratio)
                 fstn = factor(stn), # station as a factor
                 time = lubridate::decimal_date(obstime) # continuous time in year units
  ) 
  # use sampleID as rownames
  rownames(meta) <- meta$sampleID
} # provides data tables dat, meta, bscl

# run Nonmetric Multidimensional Scaling (vegan::metaMDS) on the global community matrix.
# we use NMDS ordination to split the sampling stations into two groups: (i) inner bays and (ii) outer archipelago

## SPLIT THE STATIONS TO INNER AND OUTER ARCHIPELAGO ####

if(T){ #  provides datmds, stni, stno, updates meta with inner/outer split
  if(F){ # NB! long calculation
    datmds <- metaMDS(datcm, k = 2) # may take some time, ca 30 min on my laptop
    save(datmds, file = 'Olli_HelSummerTrends_datmds.rda')
  } # NB! long calculation
  load(file = 'Olli_HelSummerTrends_datmds.rda')
  
  # the fist global nmds ordination axis correlates with time; the second axis correlates with the coastal-pelagic gradient. 
  # we rotate the ordination to make station depth (proxy of the coastal-pelagic gradient) parallel to the 1st ordination axis
  datmdsr <- MDSrotate(datmds, meta$stnDepth)
  
  # Calculate y-centers of gravity for stations
  
  # extract the scores of the ordination and combine them with the station variables
  mdsGlobalr   <- scores(datmdsr)[rownames(meta), ] %>% data.frame() %>% cbind(select(meta, stn, dist, lat, lon,stnDepth))
  #  calculate y-centers of gravity 
  mdsGlobalrAg <- group_by(mdsGlobalr, stn, lat, lon, dist, stnDepth) %>% summarize(n = n(), mds = mean(NMDS1))
  
  plot(sort(mdsGlobalrAg$mds), ylim = c(-0.5, 0.53), ylab = 'Community structure index')
  
  # there seem to be a few gaps in the community structure along the coastal-pelagic gradient, which we can use to separate the coastal bays from the outer archipelago
  # a reasonable split is at 0.07, which we set to categorize the stations
  abline(h = 0.07)
  
  mdsGlobalrAg$pel <- ifelse(mdsGlobalrAg$mds > 0.07, 0, 1) # 0 = 'inner bays'; 1 - outer archipelago
  
  # insert station categories to meta table
  meta$pel <- mdsGlobalrAg[match(meta$stn, mdsGlobalrAg$stn), ]$pel
  
  # inner and outer station names
  stno <- filter(meta, pel == 1)$stn %>% unique()
  stni <- filter(meta, pel == 0)$stn %>% unique()
  
  # stno <- c(68, 71, 75, 44, 36, 18, 111, 125, 107, 114, 113, 124, 122, 117, 118, 120, 142, 147, 62, 154, 166, 149, 168, 148, 174, 180)
  # stni <- c(87, 81, 25, 92, 4, 140, 23, 146, 97, 94)
}  #  provides datmds, stni, stno, updates meta with inner/outer split

## Fig 1S, global NMDS ##

if(T){ # saves file Olli_HelSummerTrends_Figure1S
  # Fig 1S global nmds ordination of samples, showing the split between inner bay and pelagic stations, and the station centres of gravity along the coasta-pelagic gradient
  
  mdsGlobal   <- scores(datmds)[rownames(meta), ] %>% data.frame() %>% cbind(select(meta, stn, dist, lat, lon, stnDepth, pel))
  mdsGlobalAg <- group_by(mdsGlobal, stn, lat, lon) %>% summarize(n = n(), mds1 = mean(NMDS1), mds2 = mean(NMDS2), pel = mean(pel))
  
  Fig1S <- ggplot(mdsGlobal, aes(x = NMDS1, y = NMDS2, col = factor(pel))) + 
    geom_point(alpha = 0.5) +
    scale_colour_manual(values =(mdsGlobal$pel+3), name = 'Stations', labels = c('inner bays','outer archipelago')) +
    geom_point(data = mdsGlobalAg, aes(x=mds1, y=mds2), fill = (mdsGlobalAg$pel + 3), col = 1, shape = 21, size = 4, show.legend  = FALSE) + 
    theme(legend.position = c(0.75, 0.85), legend.background = element_rect(fill = 0)) +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),legend.text = element_text(size = 14), legend.title = element_text(size = 14)) 
  
   ggsave2(Fig1S, file = './Olli_HelSummerTrends_FigureS1.pdf', width = 6, height = 6)
   ggsave2(Fig1S, file = './Olli_HelSummerTrends_FigureS1.svg', width = 6, height = 6)
} # saves file Olli_HelSummerTrends_Figure1S

# Fig 1 STATION MAP ####

if(T){ # saves file Olli_HelSummerTrends_Figure1.pdf
  library(maptools); gpclibPermit();  gpclibPermitStatus(); rgeosStatus() 
  
  # for the coastlines we need gshhs_f.b and gshhs_l.b binary files, freely available from http://www.soest.hawaii.edu/pwessel/gshhg/
  # put the files in the maptools/share folder, or point the two commands below wherever you have the coastline files
  
  BZ <- getRgshhsMap(system.file("share/gshhs_f.b", package="maptools"), xlim=c(24.55, 25.3), ylim=c(59.9, 60.25),level=1) %>% st_as_sf() # Helsinki archipelago
  BS <- getRgshhsMap(system.file("share/gshhs_l.b", package="maptools"), xlim=c(10, 31), ylim=c(53, 67),level=1) %>% st_as_sf() # Baltic Sea insert coastline 
  
  # add projection
  df <- st_as_sf(x = mdsGlobalAg, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  # Baltic Sea insert map
  ggm1 <- ggplot(data = BS) + geom_sf(fill = 'whitesmoke') + theme_void() + geom_sf(data = st_as_sfc(st_bbox(BZ_sf)), fill = NA, color = 'red', size = 0.8 )
  ggm1 <- ggplot(data = BS) + geom_sf(fill = 'whitesmoke') + theme_void() + geom_sf(data = st_as_sfc(st_bbox(BZ)), fill = NA, color = 'red', size = 0.8 )
  
  library(ggsflabel)
  library(ggspatial) # scales and N arrows
  library(gghighlight)
  
  ggm2 <- ggplot(data = BZ) + # Helsinki archipelago map
    geom_sf(fill= "whitesmoke", color = gray(.6)) + 
    coord_sf(xlim = c(24.55, 25.3), ylim = c(59.935, 60.25), expand = FALSE) +
    geom_point(data = df, aes(geometry = geometry, stat(X), stat(Y), fill = I(pel+3), size = sqrt(n)), stat = StatSfCoordinates, fun.geometry = identity, show.legend = FALSE, shape = 21) +
    geom_sf_label_repel(data = df, aes(label = stn), label.size = NA, label.padding = unit(0.1, 'lines'), alpha = 0.7, fontface = 'bold') +
    xlab('Longitude') + ylab('Latitude') +
    annotation_scale(location = "bl", width_hint = 0.3, pad_x = unit(0.5, "in"), pad_y = unit(0.3, "in")) +
    annotation_north_arrow(location = "tl", which_north = "true", pad_x = unit(0.5, "in"), pad_y = unit(0.3, "in"), style = north_arrow_fancy_orienteering) +
    theme(panel.grid.major = element_line(color = gray(.5), linetype = 'dashed', size = 0.2), panel.background = element_rect(fill = 'aliceblue'), panel.border = element_rect(colour = "black", fill=NA, size=1)) 
  
  # finalize and save Figure 1
  
  Figure1 <- ggdraw() +
    draw_plot(ggm2) +
    draw_plot(ggm1, scale = 0.4, halign = 1.11, valign = 0.145)
  
  ggsave2(Figure1, file = './Olli_HelSummerTrends_Figure1.pdf', width = 7.5, height = 6.5)
  ggsave2(Figure1, file = './Olli_HelSummerTrends_Figure1.svg', width = 7.5, height = 6.5, device = 'svg')
  
} # saves file Olli_HelSummerTrends_Figure1

## MAIN CALCULATION ####

# 1. phytoplankton community NMDS ordination in inner bays and outer archipelago
# 2. k-means clustering of communities, based on NMDS scores 

# 4. external variable correlation with the NMDS community ordination, by using vegan::envfit
  
  if(T){ # community ordination
    # We are ready to split the phytoplankton community matrix into inner bays and outer archipelago subsets, based on the newly created pel category
    MXl <- list()
    MXl[['out']] <- datcm[as.character(filter(meta, pel == 1)$sampleID), ] # 1115  619
    MXl[['in']]  <- datcm[as.character(filter(meta, pel == 0)$sampleID), ] # 3515  619
    
    # remove columns with no species occurrences
    MXl <- lapply(MXl, function(x){x <- x[, colSums(x) > 0]})
    
    sapply(MXl, dim) # remaining species richness: 546 and  505 taxa in outer archipelago and inner bays, respectively.
    
    # NMDS COMMUNITY ORDINATION ####
    
    # 1. perform ndms on outer archipelago and inner bay phytoplankton communities separately
    # NB long calculation
    MDSl <- list()
    MDSl[['in']]  <- metaMDS(MXl[['in']],  k = 2, trymax = 50)
    MDSl[['out']] <- metaMDS(MXl[['out']], k = 2, trymax = 50)
    
    # 2. categorize the samples based on nmds ordination
    
    KMEANSl <- mclapply(MDSl,  function(x){kmeans(scores(x), centers = 3, nstart = 100)})
    
    clus = c(KMEANSl[[1]]$cluster, KMEANSl[[2]]$cluster)
    meta$clus <- clus[rownames(meta)]
    
    mdscores <- rbind(scores(MDSl[[1]]), scores(MDSl[[2]]))
    meta <- cbind(meta, mdscores[rownames(meta),])
    
  } # community ordination; provides 2-element lists: MXl, MDSl, KMEANSl

 # DIVERSITY PROFILES ####

# Diversity profiles of communities
#  we use diversity profiles to reveal which order of diversity associates best with the community ordination
#  we test species diversity profiles (renyi diversity) and phylogenetic diversity, which accounts for the inherent non-independence of species due to their common evolutionary history. For the latter we use similarity based community diversity function entropart::Dqz, using the (cophenetic) phylogenetic similarity matrix
 
 if(T){ # diversity profiles; provides DqzPDiv and RenyiD, and updates meta with diversity profiles
   # we load a modified vegan::renyi function, which allows for negative orders of diversity (overweighting rare species)
   source('./Olli_HelSummerTrends_renyi2.R')
   
   ## Helper functions
   nrm <- function(x){x/sum(x)} # normalization, sum to 1
   d2s <- function(x){return(1-x/max(x))} # distance 2 similarity
   
   # to speed up parallel calculations we make probability vector list of communities
   datcm.lst <- mclapply(1:nrow(datcm), function(x){x <- nrm(datcm[x, which(datcm[x, ] > 0), drop = FALSE]); cnames <- colnames(x); x <- c(x); names(x) <- cnames; return(x)}, mc.cores = 12)
   
   # load species phylogeny
   bstree <- ape::read.tree(file = "./Olli_HelSummerTrends_bstree.txt")
   bstree$tip.label <- gsub('_',' ', bstree$tip.label)
   
   # calculate phylogenetic distance matrix
   PDist <- cophenetic(bstree)[colnames(datcm), colnames(datcm)]
   
   # transform phylogenetic similarity to distance matrix
   PSimil <- d2s(PDist)
   
   # define a vector of the orders of diversity; length 26
   q0 <- seq(-0.5, 2, by = 0.1) # diversity profile; length 26
   
   # calculate phylogenetic diversity profiles
   DqzPDiv <- list()
   for(i in 1:length(q0)){
     DqzPDiv[[i]] <- mclapply(datcm.lst, function(x){Dqz(x, q0[i], as.matrix(PSimil[names(x), names(x)]))}, mc.cores = detectCores()) %>% unlist()
   }
   # phylogenetic diversity profiles list to matrix
   DqzPDiv <- sapply(DqzPDiv, rbind)
   
   # calculate renyi diversity profiles (hill numbers)
   RenyiD <- renyi2(datcm,  hill = TRUE, scales = q0)
   
   # add column and rownames when missing
   colnames(DqzPDiv) <- sub('-','n', paste('DqzPD', q0, sep = ''))
   colnames(RenyiD) <- sub('-','n', paste('RenyiD', q0, sep = ''))
   
   rownames(DqzPDiv) <- names(datcm.lst)
   
   # add the diversity profiles to meta table
   meta <- cbind(meta, cbind(DqzPDiv, RenyiD)[as.character(meta$sampleID),])
 } # diversity profiles; provides DqzPDiv and RenyiD, and updates meta with diversity profiles

# VEGAN::ENVFIT PROFILES ####

# How external variables correlate with the community ordination

if(T){# calculate vegan::envfit 
  # variables with NA. These have to be treated individually with envfit, excluding NAs. Otherwise envfit will remove all records, where any variable is NA
  ext_var <- names(meta)[sapply(meta, anyNA)]
  
  # vegan::envfit - fit external variables onto NMDS ordination
  
  ENVFITl.lst <- mclapply(MDSl, function(x){
    # when NAs are present, run one variable at a time
    tmp.lst <- list()
    for(i in 1:length(ext_var)){
      tmp.lst[[i]] <- envfit(x, filter(meta, sampleID %in% rownames(scores(x))) %>%  select(ext_var[i]), na.rm=TRUE)
    }
    # run diversity profiles in batch against NMDS ordinations
    
    res <- c(tmp.lst, list(envfit(x, filter(meta, sampleID %in% rownames(scores(x))) %>% select(!all_of(c(ext_var, 'sampleID','stn','obstime')))) ) )
    return(res)
  }, mc.cores = 2)
  
  
  # extract the correlation coefficients between the variable and ordination form th envfit objects
  envfitl.r <- lapply(ENVFITl.lst, function(x){sapply(x, '[', 'vectors') %>% sapply('[', 'r') %>% unlist()}) # length 159
} # calculate vegan::envfit; provides ENVFITl.lst list and envfitl.r summaries

#### Fig 2 NMDS ordination ####

# community NMDS ordination plots. Samples will be colored according to sampling time
# sample symbol will correspond to clusterswe add ellipses to highlight the clusters
# we also add envfit vectors for selected environmental variables

if(T){  #  NMDS plot as Olli_HelSummerTrends_Figure2.pdf
  
  # prepare arrows
  tmp <- sapply(ENVFITl.lst[[1]], '[', 'vectors') %>% sapply('[', 'arrows')
  ar <- do.call(rbind, lapply(tmp, rbind)) %>% data.frame
  ar$r <- sapply(ENVFITl.lst[[1]], '[', 'vectors') %>% sapply('[', 'r') %>% unlist()
  ar <- mutate(ar, NMDS1 = NMDS1 * 1.5 * r, NMDS2 = NMDS2 * 1.5 * r)[c('year','N0','E20','jul'),]
  ar$var <- c('Year', 'NSP', 'E', 'J')
  
  
  Fig2a <- ggplot(filter(meta, pel == 1) , aes(x=NMDS1, y=NMDS2, fill = time, shape = factor(clus))) +
    geom_point(size = 3) +
    scale_shape_manual(values=c(22, 21, 23), guide = 'none') +
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'Year') +
    scale_x_continuous(expand = c(0,0), limits = c(-1.5, 1.7)) +
    scale_y_continuous(expand = c(0,0), limits = c(-1.5, 1.4)) + 
    theme(legend.position = c(0.94, 0.15), legend.background = element_rect(fill = 0)) +
    stat_ellipse(size = 1) +  
    annotate('segment',x = 0, y = 0, xend = ar[, 'NMDS1'], yend = ar[,'NMDS2'], arrow = arrow(length = unit(0.5, "cm")), size = 1) +
    annotate('label', x = ar[, 'NMDS1'], y = ar[,'NMDS2'], label = ar[,'var'], size = 6, fontface = 'bold',  vjust = "outward", hjust = "outward", fill = 'white', alpha =  0.7, label.size = NA, label.padding = unit(0.1, "lines")) + 
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), legend.text = element_text(size = 14), legend.title = element_text(size = 14))
  
  tmp <- sapply(ENVFITl.lst[[2]], '[', 'vectors') %>% sapply('[', 'arrows')
  ar <- do.call(rbind, lapply(tmp, rbind)) %>% data.frame
  ar$r <- sapply(ENVFITl.lst[[2]], '[', 'vectors') %>% sapply('[', 'r') %>% unlist()
  ar <- mutate(ar, NMDS1 = NMDS1 * 1.5 * r, NMDS2 = NMDS2 * 1.5 * r)[c('year','N0','pH','chla','PO4','Ptot'),]
  ar$var <- c('Year', 'NSP','pH','chla', 'PO4','Ptot')   
  
  Fig2b <- ggplot(filter(meta, pel == 0), aes(x=NMDS1, y=NMDS2, fill = time, shape = factor(clus))) +
    geom_point(size = 3) +
    scale_shape_manual(values=c(22, 21, 23), guide = 'none') +
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'Year') +
    scale_x_continuous(expand = c(0,0), limits = c(-1.6, 2.2)) +
    scale_y_continuous(expand = c(0,0), limits = c(-1.5, 1.4)) + 
    theme(legend.position = c(0.50, 0.15), legend.background = element_rect(fill = 0)) +
    stat_ellipse(size = 1) +
    annotate('segment',x = 0, y = 0, xend = ar[, 'NMDS1'], yend = ar[,'NMDS2'], arrow = arrow(length = unit(0.5, "cm")), size = 1) +
    annotate('label', x = ar[, 'NMDS1'], y = ar[,'NMDS2'], label = ar[,'var'], size = 6, fontface = 'bold',  vjust = "outward", hjust = "outward", fill = 'white', alpha =  0.7, label.size = NA, label.padding = unit(0.1, "lines")) + 
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), legend.text = element_text(size = 14), legend.title = element_text(size = 14))  
  
  
  Figure2 <- plot_grid(Fig2a, Fig2b, labels = "AUTO", label_size = 24) # glue panels together
  cowplot::save_plot("Olli_HelSummerTrends_Figure2.pdf", Figure2, ncol = 2, base_height = 6, base_width = 6)
  cowplot::save_plot("Olli_HelSummerTrends_Figure2.svg", Figure2, ncol = 2, base_height = 6, base_width = 6)

} # saves Fig 2 NMDS plot as Olli_HelSummerTrends_Figure2

#### Fig 3 NMDS1 and species richness trends ####

# 1st NMDS axis is strongly associated with time trend
# we visualize this by plotting the NMDS scores against time. 
# symbols correspond to sample clusters, color corresponds to species richness
# the GAM smooth trend line is split into sections of significant increase (red) or decrease (black) - i.e. the derivative of the slope is significantly different from zero

# Species richness correlates well with the 1st NMDS axis, and has a conspicuous increasing long-term trend
# we plot the species richness trend
# symbols correspond to sample clusters, color corresponds to NMDS1 score values
# the GAM smooth trend line is split into sections of significant increase (red) or decrease (black) - i.e. the derivative of the slope is significantly
# the variance explained is notably lower than with NMDS1 scores

if(T){ # saves Fig 3 NMDS and richness trends as Olli_HelSummerTrends_Figure3.pdf
  
  ctrl <- list(niterEM = 0, msVerbose = FALSE, optimMethod="L-BFGS-B")
  # m0 model
  m0a  <- gamm(NMDS1 ~ s(time, bs = 'tp') + s(jul),  data = filter(meta, pel == 1), correlation = corCAR1(form = ~ 1|year), control = ctrl, method="REML") # C
  preda <- tidymv::predict_gam(m0a$gam, exclude_terms = c('s(jul)')) 
  der <- gratia::derivatives(m0a$gam, term = 's(time)', newdata = preda)
  # add increasing and decreasing
  preda$incr <- ifelse(der$derivative > 0 & der$lower > 0, preda$fit, NA)
  preda$decr <- ifelse(der$derivative < 0 & der$upper < 0, preda$fit, NA)
  
  Fig3A <- ggplot(filter(meta, pel == 1), aes(x=time, y=NMDS1)) + 
    geom_point(aes(shape = factor(clus), fill = N0), size = 3) +
    scale_shape_manual(values=c(22, 21, 23), guide = FALSE) +
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'NSP') +
    scale_y_continuous(expand = c(0,0), limits = c(-1.5, 1.9)) + 
    scale_x_continuous(breaks = seq(1970, 2020, by = 10), limits = c(1965, 2020)) + 
    theme(legend.position = c(0.85, 0.15), legend.background = element_rect(fill = 0)) +
    geom_line(data = preda, aes(x = time, y = fit), col = 'white', size = 1) +
    geom_line(data = preda, aes(x = time, y = I(fit+2*se.fit)), col = 'gray',size = 1) +
    geom_line(data = preda, aes(x = time, y = I(fit-2*se.fit)), col = 'gray',size = 1) +
    geom_line(data = preda, aes(x = time, y=incr), col = 'red', size = 1.5) + # slight difference in smooth compared to mg model
    geom_line(data = preda, aes(x = time, y=decr), col = 'black', size = 1.5) +
    xlab('Year')  + 
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), legend.text = element_text(size = 14), legend.title = element_text(size = 14))  
  

  m0b  <- gamm(NMDS1 ~ s(time, bs = 'tp') + s(jul),  data = filter(meta, pel == 0), correlation = corCAR1(form = ~ 1|year), control = ctrl, method="REML") # C
  predb <- tidymv::predict_gam(m0b$gam, exclude_terms = c('s(jul)')) # works not, too many terms to exclude
  der <- gratia::derivatives(m0b$gam, term = 's(time)', newdata = predb)
  # add increasing and decreasing
  predb$incr <- predb$decr <- 0
  predb$incr <- ifelse(der$derivative > 0 & der$lower > 0, predb$fit, as.numeric(NA))
  predb$decr <- ifelse(der$derivative < 0 & der$upper < 0, predb$fit, as.numeric(NA))
  
  Fig3B <- ggplot(filter(meta, pel == 0), aes(x = time, y = NMDS1)) + 
    geom_point(aes(shape = factor(clus), fill = N0), size = 3) +
    scale_shape_manual(values=c(22, 21, 23), guide = FALSE) +
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'NSP') +
    scale_y_continuous(expand = c(0,0), limits = c(-1.5, 2.2)) + 
    scale_x_continuous(breaks = seq(1970, 2020, by = 10), limits = c(1965, 2020)) + 
    theme(legend.position = c(0.85, 0.15), legend.background = element_rect(fill = 0)) +
    geom_line(data = predb, aes(x = time, y = fit), col = 'white', size = 1) +
    geom_line(data = predb, aes(x = time, y = I(fit+2*se.fit)), col = 'gray',size = 1) +
    geom_line(data = predb, aes(x = time, y = I(fit-2*se.fit)), col = 'gray',size = 1) +
    geom_line(data = predb, aes(x = time, y=incr), col = 'red', size = 1.5) + # slight difference in smooth compared to mg model
    geom_line(data = predb, aes(x = time, y=decr), col = 'black', size = 1.5) +
    xlab('Year') + 
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), legend.text = element_text(size = 14), legend.title = element_text(size = 14))  
    
  
  m0c  <- gamm(N0 ~ s(time, bs = 'tp') + s(jul),  data = filter(meta, pel == 1), correlation = corCAR1(form = ~ 1|year), control = ctrl, method="REML")
  predc <- tidymv::predict_gam(m0c$gam, exclude_terms = c('s(jul)'))
  der <- gratia::derivatives(m0c$gam, term = 's(time)', newdata = predc)
  # add increasing and decreasing segments
  predc$incr <- ifelse(der$derivative > 0 & der$lower > 0, predc$fit, NA)
  predc$decr <- ifelse(der$derivative < 0 & der$upper < 0, predc$fit, NA)
  
  Fig3C <- ggplot(filter(meta, pel == 1), aes(x = time, y=N0)) + 
    geom_point(aes(shape = factor(clus), fill = NMDS1), size = 3) +
    scale_shape_manual(values=c(22, 21, 23), guide = FALSE) +
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'NMDS1') +
    scale_y_continuous(expand = c(0,0), limits = c(0, 70)) + 
    scale_x_continuous(breaks = seq(1970, 2020, by = 10), limits = c(1965, 2020)) + 
    theme(legend.position = c(0.9, 0.15), legend.background = element_rect(fill = 0)) +
    geom_line(data = predc, aes(x = time, y = fit), col = 'white', size = 1) +
    geom_line(data = predc, aes(x = time, y = I(fit+2*se.fit)), col = 'gray',size = 1) +
    geom_line(data = predc, aes(x = time, y = I(fit-2*se.fit)), col = 'gray',size = 1) +
    geom_line(data = predc, aes(x = time, y=incr), col = 'red', size = 1.5) + 
    geom_line(data = predc, aes(x = time, y=decr), col = 'black', size = 1.5) +
    xlab('Year') + ylab('Species richness') +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) 
  ###
  
  m0d  <- gamm(N0 ~ s(time, bs = 'tp') + s(jul),  data = filter(meta, pel == 0), correlation = corCAR1(form = ~ 1|year), control = ctrl, method="REML")
  predd <- tidymv::predict_gam(m0d$gam, exclude_terms = c('s(jul)')) 
  der <- gratia::derivatives(m0d$gam, term = 's(time)', newdata = predd)
  # add increasing and decreasing segments
  predd$incr <- predd$decr <- 0
  predd$incr <- ifelse(der$derivative > 0 & der$lower > 0, predd$fit, as.numeric(NA))
  predd$decr <- ifelse(der$derivative < 0 & der$upper < 0, predd$fit, as.numeric(NA))
  
  Fig3D <- ggplot(filter(meta, pel == 0), aes(x = time, y=N0)) + 
    geom_point(aes(shape = factor(clus), fill = NMDS1), size = 3) +
    scale_shape_manual(values=c(22, 21, 23), guide = FALSE) +
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'NMDS1') +
    scale_y_continuous(expand = c(0,0), limits = c(0,70)) + 
    scale_x_continuous(breaks = seq(1970, 2020, by = 10), limits = c(1965, 2020)) + 
    theme(legend.position = c(0.9, 0.15), legend.background = element_rect(fill = 0)) +
    geom_line(data = predd, aes(x = time, y = fit), col = 'white', size = 1) +
    geom_line(data = predd, aes(x = time, y = I(fit+2*se.fit)), col = 'gray',size = 1) +
    geom_line(data = predd, aes(x = time, y = I(fit-2*se.fit)), col = 'gray',size = 1) +
    geom_line(data = predd, aes(x = time, y = incr), col = 'red', size = 1.5) + # slight difference in smooth compared to mg model
    geom_line(data = predd, aes(x = time, y = decr), col = 'black', size = 1.5) +
    xlab('Year') + ylab('Species richness') +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) 

  Figure3 <- plot_grid(Fig3A, Fig3B, Fig3C, Fig3D, labels = "AUTO", label_size = 18) # glue panels together
  cowplot::save_plot("./Olli_HelSummerTrends_Figure3.pdf", Figure3, ncol=2, base_height = 12, base_width = 8)
  cowplot::save_plot("./Olli_HelSummerTrends_Figure3.svg", Figure3, ncol=2, base_height = 12, base_width = 8)
  
} # saves Fig 3 NMDS1 and richness trends as Olli_HelSummerTrends_Figure3.pdf

#### Fig S6 evenness trends ####

# Species evenness (Shannon evenness or Hill's ratio) is another essential facet of diversity. It does have a long-term trend, but the variance explained is notably lower than with species richness
# we plot the species richness trend
# symbols correspond to sample clusters, color corresponds to NMDS1 score values
# the GAM smooth trend line is split into sections of significant increase (red) or decrease (black) - i.e. the derivative of the slope is significantly

if(T){ # saves Fig S6 richness trends as Olli_HelSummerTrends_FigureS6.pdf
  ctrl <- list(niterEM = 0, msVerbose = FALSE, optimMethod="L-BFGS-B")
  # m0 model
  m0  <- gamm(E10 ~ s(time, bs = 'tp') + s(jul),  data = filter(meta, pel == 1), correlation = corCAR1(form = ~ 1|year), control = ctrl, method="REML")
  pred <- tidymv::predict_gam(m0$gam, exclude_terms = c('s(jul)'))
  der <- gratia::derivatives(m0$gam, term = 's(time)', newdata = pred)
  # add increasing and decreasing segments
  pred$incr <- ifelse(der$derivative > 0 & der$lower > 0, pred$fit, NA)
  pred$decr <- ifelse(der$derivative < 0 & der$upper < 0, pred$fit, NA)
  
  FigS6A <- ggplot(filter(meta, pel == 1), aes(x = time, y=E10)) + 
    geom_point(aes(shape = factor(clus), fill = NMDS1), size = 3) +
    scale_shape_manual(values=c(22, 21, 23), guide = FALSE) +
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'NMDS1') +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.05)) + 
    scale_x_continuous(breaks = seq(1970, 2020, by = 10), limits = c(1965, 2020)) + 
    theme(legend.position = c(0.85, 0.85), legend.background = element_rect(fill = 0)) +
    geom_line(data = pred, aes(x = time, y = fit), col = 'white', size = 1) +
    geom_line(data = pred, aes(x = time, y = I(fit+2*se.fit)), col = 'gray',size = 1) +
    geom_line(data = pred, aes(x = time, y = I(fit-2*se.fit)), col = 'gray',size = 1) +
    geom_line(data = pred, aes(x = time, y=incr), col = 'red', size = 1.5) +
    geom_line(data = pred, aes(x = time, y=decr), col = 'black', size = 1.5) +
    xlab('Year') + ylab('Species evenness') +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),legend.text = element_text(size = 14), legend.title = element_text(size = 14)) 
  
  m0  <- gamm(E10 ~ s(time, bs = 'tp') + s(jul),  data = filter(meta, pel == 0), correlation = corCAR1(form = ~ 1|year), control = ctrl, method="REML") # C
  pred <- tidymv::predict_gam(m0$gam, exclude_terms = c('s(jul)'))
  der <- gratia::derivatives(m0$gam, term = 's(time)', newdata = pred)
  # add increasing and decreasing segments
  pred$incr <- pred$decr <- 0
  pred$incr <- ifelse(der$derivative > 0 & der$lower > 0, pred$fit, as.numeric(NA))
  pred$decr <- ifelse(der$derivative < 0 & der$upper < 0, pred$fit, as.numeric(NA))
  
  FigS6B <- ggplot(filter(meta, pel == 0), aes(x = time, y=E10)) + 
    geom_point(aes(shape = factor(clus), fill = NMDS1), size = 3) +
    scale_shape_manual(values=c(22, 21, 23), guide = FALSE) +
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'NSP') +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.05)) + 
    scale_x_continuous(breaks = seq(1970, 2020, by = 10), limits = c(1965, 2020)) + 
    theme(legend.position = c(0.85, 0.85), legend.background = element_rect(fill = 0)) +
    geom_line(data = pred, aes(x = time, y = fit), col = 'white', size = 1) +
    geom_line(data = pred, aes(x = time, y = I(fit+2*se.fit)), col = 'gray',size = 1) +
    geom_line(data = pred, aes(x = time, y = I(fit-2*se.fit)), col = 'gray',size = 1) +
    geom_line(data = pred, aes(x = time, y=incr), col = 'red', size = 1.5) + 
    geom_line(data = pred, aes(x = time, y=decr), col = 'black', size = 1.5) +
    xlab('Year') + ylab('Species evenness') +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),legend.text = element_text(size = 14), legend.title = element_text(size = 14)) 
  
  FigureS6 <- plot_grid(FigS6A, FigS6B, labels = "AUTO", label_size = 24) # glue panels together
  cowplot::save_plot("Olli_HelSummerTrends_FigureS6.pdf", FigureS6, ncol = 2, base_height = 6, base_width = 6)
  cowplot::save_plot("Olli_HelSummerTrends_FigureS6.svg", FigureS6, ncol = 2, base_height = 6, base_width = 6)
} # saves Fig S6 richness trends as Olli_HelSummerTrends_FigureS6


#### Fig 4 SPECIES CUMULATIVE LIKELIHOOD BLOCK ####

# What has changed in the community composition over the decades? Which groups have increased, which ones have decreased?
# We calculate the likelihood of presence of species over the 52 year span by fitting a binomial GAM model
# We order the species along time by calculating the the center of gravity of occurrence likelihood along the time axis


if(T){  # saves species likelihood figure Olli_HelSummerTrends_Figure4.pdf 
  
  # Species likelihood matrices ##
  # We use binary response (presence/absence) and delete rare species (5 or less occurrences) and at least genus level or lower taxonomic level
  datcm0 <-  decostand(MXl[[2]], 'pa') %>% .[, colSums(.) > 5] %>% .[, colnames(.) %in% filter(bscl, !is.na(Genus))$valid_name] # 1115  297
  datcm1 <-  decostand(MXl[[1]], 'pa') %>% .[, colSums(.) > 5] %>% .[, colnames(.) %in% filter(bscl, !is.na(Genus))$valid_name] # 3515  345
  
  # define time sequence
  dtdf <- data.frame(time = seq(1966, 2019, by =  0.25))
  
  # initialize cumulative likelihood matrices for inner bays (cumlik0) and outer archipelago (cumlik1)
  predcumul0 <- matrix(0, ncol = ncol(datcm0), nrow = nrow(dtdf), dimnames = list(as.character(dtdf$time), colnames(datcm0))) # 213 296
  predcumul1 <- matrix(0, ncol = ncol(datcm1), nrow = nrow(dtdf), dimnames = list(as.character(dtdf$time), colnames(datcm1))) # 213 344
  
  # Fill the matrices with binomial GAM predictions, one species at a time
  # NB long calculation
  for(i in 1:ncol(predcumul1)){
    predcumul1[, i] <- predict(gam(datcm1[, i] ~ s(time), data = meta[rownames(datcm1), ], family = binomial, method = 'REML'), dtdf, type  = 'response')
  }
  # NB long calculation
  for(i in 1:ncol(predcumul0)){
    predcumul0[, i] <- predict(gam(datcm0[, i] ~ s(time), data = meta[rownames(datcm0), ], family = binomial, method = 'REML'),  dtdf, type  = 'response')
  }
  
  # species temporal optimal along the time axis
  dk1.opt <- array(dim=ncol(predcumul1))
  for(i in 1:ncol(predcumul1)){
    pred <- predcumul1[,i]
    dk1.opt[i] <- weighted.mean(dtdf$time, predcumul1[, i])
  }
  names(dk1.opt) <- colnames(predcumul1)
  
  dk0.opt <- array(dim=ncol(predcumul0))
  for(i in 1:ncol(predcumul0)){
    pred <- predcumul0[,i]
    dk0.opt[i] <- weighted.mean(dtdf$time, predcumul0[, i])
  }
  names(dk0.opt) <- colnames(predcumul0)
  
  # re-order taxa according to their occurrence optima
  
  dk0.opt.ord <- sort(dk0.opt)
  dk1.opt.ord <- sort(dk1.opt)
  
  predcumul1.ord <- predcumul1[, names(dk1.opt.ord)] 
  predcumul0.ord <- predcumul0[, names(dk0.opt.ord)]
  
  # split the taxa into thee groups, accordin to their optimal time of occurrence
  tg0 <- data.frame(tax = names(dk0.opt), optloc = dk0.opt,  group = ntile(dk0.opt, 3)) %>% arrange(optloc)
  tg1 <- data.frame(tax = names(dk1.opt), optloc = dk1.opt,  group = ntile(dk1.opt, 3)) %>% arrange(optloc)
  
  
  # data frames for cumulative likelihood plot
  predcumul0.df <- data.frame(value = as.vector(predcumul0.ord), 
                              time = rep(as.numeric(rownames(predcumul0.ord)), ncol(predcumul0.ord)), 
                              tax = (rep(colnames(predcumul0.ord), each = nrow(predcumul0.ord))))
  
  predcumul1.df <- data.frame(value = as.vector(predcumul1.ord), 
                              time = rep(as.numeric(rownames(predcumul1.ord)), ncol(predcumul1.ord)), 
                              tax = (rep(colnames(predcumul1.ord), each = nrow(predcumul1.ord))))
  
  # add order for plotting the area graph
  predcumul0.df <- left_join(predcumul0.df, data.frame(tax = names(dk0.opt.ord), ord = dk0.opt.ord))
  predcumul1.df <- left_join(predcumul1.df, data.frame(tax = names(dk1.opt.ord), ord = dk1.opt.ord))
  
  # data frame for discrete cumulative lines
  
  linescumul0 <- data.frame(time = as.numeric(rownames(predcumul0.ord)), val1 = rowSums(predcumul0.ord), 
                            val2 = rowSums(predcumul0.ord[, rownames(filter(tg0, group == 3))]),  
                            val3 = rowSums(predcumul0.ord[, rownames(filter(tg0, group != 1))]) )
  
  linescumul1 <- data.frame(time = as.numeric(rownames(predcumul1.ord)), val1 = rowSums(predcumul1.ord), 
                            val2 = rowSums(predcumul1.ord[, rownames(filter(tg1, group == 3))]),  
                            val3 = rowSums(predcumul1.ord[, rownames(filter(tg1, group != 1))]) )
  
  
  
  Fig4A <- ggplot(predcumul0.df, aes(x = time, y = value)) +  
    geom_area(show.legend = FALSE, aes(fill = factor(ord))) +
    scale_fill_discrete(type = colorRamps::matlab.like(ncol(predcumul0.ord))) +
    geom_line(data = linescumul0, aes(x = time, y = val1), col = 1) +
    geom_line(data = linescumul0, aes(x = time, y = val2), col = 1) +
    geom_line(data = linescumul0, aes(x = time, y = val3), col = 1) +
    xlab('Year') + 
    ylab('Cumulative likelihood of presence') + 
    scale_y_continuous(expand = c(0,0), limits = c(0, 40), breaks = seq(0, 30, by = 10)) + 
    scale_x_continuous(expand = c(0,0), breaks = seq(1970, 2020, by = 10), limits = c(1966, 2015)) +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18)) 
  
  Fig4B <- ggplot(predcumul1.df, aes(x = time, y = value)) +  
    geom_area(show.legend = FALSE, aes(fill = factor(ord))) +
    scale_fill_discrete(type = colorRamps::matlab.like(ncol(predcumul1.ord))) +
    geom_line(data = linescumul1, aes(x = time, y = val1), col = 1) +
    geom_line(data = linescumul1, aes(x = time, y = val2), col = 1) +
    geom_line(data = linescumul1, aes(x = time, y = val3), col = 1) +
    xlab('Year') + 
    ylab('Cumulative likelihood of presence') + 
    scale_y_continuous(expand = c(0,0), limits = c(0, 40), breaks = seq(0, 30, by = 10)) + 
    scale_x_continuous(expand = c(0,0), breaks = seq(1970, 2020, by = 10), limits = c(1966, 2015)) +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18)) 
  
  Fig4 <- plot_grid(Fig4A, Fig4B, labels = "AUTO", label_size = 24) # glue panels together
  cowplot::save_plot("./Olli_HelSummerTrends_Figure4.pdf", Fig4, ncol=2, base_height = 6, base_width = 6)
  cowplot::save_plot("./Olli_HelSummerTrends_Figure4.svg", Fig4, ncol=2, base_height = 6, base_width = 6)
  
} # saves species likelihood figure Olli_HelSummerTrends_Figure4.pdf 

### Tables S2 and S3  ####

if(T){ # species likelihood tables
  
  # cumulative table-S
  tblS2a <- data.frame(valid_name =  names(dk0.opt), optloc = dk0.opt, freq = colSums(datcm0), ntile = ntile(dk0.opt, 3))
  tblS2a <- left_join(tblS2a, select(bscl, valid_name, Kingdom, Class))
  TableS2a <- table(tblS2a$Class, tblS2a$ntile)[c('Bacillariophyceae', 'Chlorophyceae','Trebouxiophyceae', 'Cyanophyceae', 'Dinophyceae'),] 
  # diatoms and Chlorophyceae decrease, Cyanobacteria and Dinoflagellates increase
  
  # ntile boundaries
  arrange(tblS2a, optloc) %>% group_by( ntile) %>% summarise(beg=first(optloc), end=last(optloc)) %>% data.frame()
  
  tblS2b <- data.frame(valid_name =  names(dk1.opt), optloc = dk1.opt, freq = colSums(datcm1), ntile = ntile(dk1.opt, 3))
  tblS2b <- left_join(tblS2b, select(bscl, valid_name, Kingdom, Class))
  TableS2b <- table(tblS2b$Class, tblS2b$ntile)[c('Bacillariophyceae', 'Chlorophyceae','Trebouxiophyceae', 'Cyanophyceae', 'Dinophyceae'),]
  
  # ntile boundaries
  arrange(tblS2b, optloc) %>% group_by( ntile) %>% summarise(beg=first(optloc), end=last(optloc)) %>% data.frame()
  
  TableS2 <- cbind(TableS2a, TableS2b)
  colnames(TableS2) <- paste(rep(c('Outer', 'Inner'), each = 3), colnames(t1))
  
  stargazer(TableS2 , type = 'text', title = 'Table S2')
  
  
  # what about species?
  tblS3a <- group_by(tblS2a, ntile, Class) %>% arrange(desc(freq), by_group = TRUE) %>%  summarise(Taxon = head(valid_name, 3), freq = head(freq, 3)) %>% filter(Class %in% c('Bacillariophyceae', 'Chlorophyceae','Trebouxiophyceae', 'Cyanophyceae', 'Dinophyceae')) %>% arrange(Class, ntile, desc(freq)) %>% data.frame()
  
  tblS3b <- group_by(tblS2b, ntile, Class) %>% arrange(desc(freq), by_group = TRUE) %>%  summarise(Taxon = head(valid_name, 3), freq = head(freq, 3)) %>% filter(Class %in% c('Bacillariophyceae', 'Chlorophyceae','Trebouxiophyceae', 'Cyanophyceae', 'Dinophyceae')) %>% arrange(Class, ntile, desc(freq)) %>% data.frame()
  
  TableS3 <- cbind(tblS3a, tblS3b[, 3:4])
  
  # most frequent taxa per Class and time period
  stargazer(as.matrix(TableS3) , type = 'text', title = 'Table S3')
  
  
} # most frequent species tables S2 and S3


#### BETA DIVERSITY BLOCK ####

# we use Bray-Curtis dissimilarity between a pair of samples as a beta diversity metrics.
# first we build a table of pairwise dissimilarities. As we have a total of 4630 samples, this gives us 4630*4629/2 = 10,716,135 pairwise dissimilarities
# to visualize the trends in beta diversity, we calculate mean dissimilarities per year, constraining the sample pairs to the same year, month (to avoid seasonal confounding effect), and station (to avoid spatial confounding effect)

if(TRUE){ # make df - the beta table
  # sqrt transform the community raw community matrix
  datcmsqrt <- sqrt(datcm[as.character(meta$sampleID), ])
  
  if(F){ # NB long calculation
    beta.decomposition <- adespatial::beta.div.comp(datcmsqrt, coef = "S", quant = T) # D%diff, Sørensen aka Bray
    save(beta.decomposition, file = 'Olli_HelSummerTrends_beta_decomposition.rda')
  } # NB long calculation
  load(file = 'Olli_HelSummerTrends_beta_decomposition.rda')
  # provides S1 = mean(S1d), S1rich = mean(S1rich), S1repl = mean(S1repl)
  
  # construct pairwise distance table
  Dyear <- stats::dist(meta$year) # euclidean distance between sampling years; units = year; l = 10716135
  idxs <- row.dist(Dyear) # start year index, min 2, max 4630
  idxe <- col.dist(Dyear) # end year index; min 1, max 4629
  Dobs = as.numeric(stats::dist(meta$time))
  
  # construct the betadf data frame
  betadf <- data.frame(Dyear = as.numeric(Dyear), Dobs = Dobs)
  
  # add sampling years of the first (years) and second sample (yeare)
  betadf <- mutate(betadf, years = meta[idxs,]$year, yeare = meta[idxe,]$year)
  
  # add sampling station of the first (stns) and second sample (stne)
  betadf <- mutate(betadf, stns = meta[idxs,]$stn, stne = meta[idxe,]$stn)
  
  # add sampling month of the first (mons) and second sample (mone)
  betadf <- mutate(betadf, mons = meta[idxs,]$mon, mone = meta[idxe,]$mon )
  
  # add beta diversity decomposition outcome
  betadf <- mutate(betadf, S = as.numeric(beta.decomposition$D), Srich = as.numeric(beta.decomposition$rich), Srepl = as.numeric(beta.decomposition$repl))
  
  # object.size(betadf) # 8.6 Gb
  
} # provides betadf - a 8.6 Gb beta diversity data frame

# calculates some beta statistics, which are used in the manuscript text

if(T){ # beta statistics
  # Results::Beta diversity decomposition
  
  # "The mean pairwise beta diversity across the whole investigated period, constrained within the same year, month, and station, was 0.46 and 0.39 in the outer archipelago and inner bays, respectively."
# "The mean species turnover components were 0.30 and 0.25, respectively, and the richness difference components 0.16 and 0.14."
  # "the variation in species richness accounted for no more than 35-36% of the mean between sample beta diversity"
  filter(betadf,  mone == mons, stns == stne, stns %in% stni,  Dyear == 0)  %>% dplyr::summarise(S = mean(S), Srich = mean(Srich), Srepl = mean(Srepl)) %>% round(2)
  filter(betadf,  mone == mons, stns == stne, stns %in% stno,  Dyear == 0)  %>% dplyr::summarise(S = mean(S), Srich = mean(Srich), Srepl = mean(Srepl)) %>% round(2)
  
} # beta statistics

### Fig. 5 TRENDS IF BETA DIVERSITY ####

if(T){ # saves Fig 5 beta diversity trends as Olli_HelSummerTrends_Figure5.pdf
  F5A <- filter(betadf, mone == mons, stns == stne, stns %in% stni,  Dyear < 4) %>% group_by(years, Dyear) %>% dplyr::summarise(S = mean(S)) # 150 x 5
  # outer stations
  F5B <- filter(betadf, mone == mons, stns == stne, stns %in% stno,  Dyear < 4) %>% group_by(years, Dyear) %>%  dplyr::summarise(S = mean(S)) # 206 x 5
  
  Fig5A <- ggplot(F5A, aes(years, S)) + # S1 is Sørensen aka Bray aka D%Diff, J1 is Jaccard aka Ružička 
    geom_line(aes(size = factor(Dyear))) + 
    geom_smooth(aes(linetype = factor(Dyear)), col = 1) +
    scale_colour_manual(values = c("black", "black", "black", "black"),  guide = 'none') +
    scale_size_manual(values = c(.75, 0.5,0.5,0.5), guide = 'none') + 
    scale_x_continuous(breaks = seq(1970, 2020, by = 10), limits = c(1965, 2020)) + 
    scale_linetype_manual(values = c(1, 5,2,3), name = 'Year distance', labels = c('same year','1 year apart', '2 years apart', '3 years apart')) + 
    theme(legend.position = c(0.2, 0.88), legend.background = element_rect(fill = 0)) +
    xlab('Years') + 
    ylab('Bray-Curtis dissimilarity') +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), legend.text = element_text(size = 14), legend.title = element_text(size = 14), legend.key.width= unit(1.5, 'cm')) 
  
  Fig5B <- ggplot(F5B, aes(years, S)) + # S1 is Sørensen aka Bray aka D%Diff, J1 is Jaccard aka Ružička 
    geom_line(aes(size = factor(Dyear))) + 
    geom_smooth(aes(linetype = factor(Dyear)), col = 1) +
    scale_colour_manual(values = c("black", "black", "black", "black"),  guide = 'none') +
    scale_size_manual(values = c(.75, 0.5,0.5,0.5), guide = 'none') + 
    scale_x_continuous(breaks = seq(1970, 2020, by = 10), limits = c(1965, 2020)) + 
    scale_linetype_manual(values = c(1, 5,2,3), name = 'Year distance', labels = c('same year','1 year apart', '2 years apart', '3 years apart')) + 
    theme(legend.position = c(0.2, 0.12), legend.background = element_rect(fill = 0)) +
    xlab('Years') + 
    ylab('Bray-Curtis dissimilarity') +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),legend.text = element_text(size = 14), legend.title = element_text(size = 14), legend.key.width= unit(1.5, 'cm')) 
  
  Figure5 <- plot_grid(Fig5B, Fig5A, labels = "AUTO", label_size = 24) # glue panels together
  cowplot::save_plot("./Olli_HelSummerTrends_Figure5.pdf", Figure5, ncol = 2, base_height = 6, base_width = 6)
  cowplot::save_plot("./Olli_HelSummerTrends_Figure5.svg", Figure5, ncol = 2, base_height = 6, base_width = 6)
  
} # saves Fig 5 beta diversity trends as Olli_HelSummerTrends_Figure5.pdf

### Fig. 6 TRENDS IN SPATIAL BETA DIVERSITY ####
## aka inter-station beta

if(T){ # saves Fig 6 beta diversity trends as Olli_HelSummerTrends_Figure6.pdf
  # inner stations
  F6A <- filter(betadf, mone == mons, stns != stne, stns %in% stni,  Dyear < 1) %>% group_by(years, Dyear) %>% dplyr::summarise(S = mean(S)) # 254 x  5
  # outer stations
  F6B <- filter(betadf, mone == mons, stns != stne, stns %in% stno,  Dyear < 1) %>% group_by(years, Dyear) %>% dplyr::summarise(S = mean(S))

  F6A$pel <- 'Coastal'
  F6B$pel <- 'Pelagic'
  F6 <- rbind(F6A, F6B)
  
  Figure6 <- ggplot(F6, aes(years, S,  linetype = pel)) + 
    geom_line() + 
    geom_smooth(col = 1) + xlab('Years') + 
    ylab('Bray-Curtis dissimilarity') + 
    theme(legend.position = c(0.2, 0.2), legend.background = element_rect(fill = 0)) + 
    theme(legend.title=element_blank()) +
    scale_x_continuous(breaks = seq(1970, 2020, by = 10), limits = c(1965, 2020)) + 
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),legend.text = element_text(size = 14)) 
  
  ggsave('./Olli_HelSummerTrends_Figure6.pdf', Figure6, width = 6, height = 6)
  ggsave('./Olli_HelSummerTrends_Figure6.svg', Figure6, width = 6, height = 6)
  
} # saves Fig 6 beta diversity trends as Olli_HelSummerTrends_Figure6.pdf

### Fig. 7 BETA DIVERSITY DECOMPOSITION ####

if(T){ # saves Fig 7 beta diversity trends as Olli_HelSummerTrends_Figure7.pdf
  
  # month and station constrained beta decomposition with 0 and 3 year lag
  # inner stations
  F7A <- filter(betadf, mone == mons, stns == stne, stns %in% stni, Dyear %in% c(0, 3)) %>% group_by(years, Dyear) %>% dplyr::summarise(Srich = mean(Srich)) 
  # outer stations
  F7B <- filter(betadf, mone == mons, stns == stne, stns %in% stno, Dyear %in% c(0, 3)) %>% group_by(years, Dyear) %>% dplyr::summarise(Srich = mean(Srich))
  
  Fig7A <- ggplot(F7A, aes(years, Srich, col = factor(Dyear))) +
    geom_line(aes(size = factor(Dyear))) + 
    geom_smooth(aes(linetype = factor(Dyear))) +
    scale_colour_manual(values = c("black", "black"), guide = 'none' ) +
    scale_size_manual(values = c(.75, 0.5), guide = 'none') + 
    scale_linetype_manual(values = c(1, 5), guide = 'none') + 
    scale_x_continuous(breaks = seq(1970, 2020, by = 10), limits = c(1965, 2020)) + 
    xlab('Years') + 
    ylab('Relativized abundance difference')+
    theme(legend.position = c(0.1, 0.2), legend.background = element_rect(fill = 0)) + 
    theme(legend.title=element_blank(), axis.title = element_text(size = 18), axis.text = element_text(size = 18))
  
  Fig7B <- ggplot(F7B , aes(years, Srich, col = factor(Dyear))) +
    geom_line(aes(size = factor(Dyear))) + 
    geom_smooth(aes(linetype = factor(Dyear))) +
    scale_colour_manual(values = c("black", "black"),  guide = 'none') +
    scale_size_manual(values = c(.75, 0.5), guide = 'none') + 
    scale_linetype_manual(values = c(1, 5),guide = 'none') + 
    scale_x_continuous(breaks = seq(1970, 2020, by = 10), limits = c(1965, 2020)) + 
    xlab('Years') + 
    ylab('Relativized abundance difference') +
    theme(legend.title=element_blank(), axis.title = element_text(size = 18), axis.text = element_text(size = 18))
  
  F7 <- plot_grid(Fig7B, Fig7A, labels = "AUTO", label_size = 24) # glue panels together
  
  cowplot::save_plot("./Olli_HelSummerTrends_Figure7.pdf", F7, ncol = 2, base_height = 6, base_width = 6)
  cowplot::save_plot("./Olli_HelSummerTrends_Figure7.svg", F7, ncol = 2, base_height = 6, base_width = 6)
  
} # # saves Fig 7 beta diversity trends as Olli_HelSummerTrends_Figure7.pdf

### Fig. 8 DISTANCE DECAY ####

if(T){ # saves Fig 9 beta diversity trends as Olli_HelSummerTrends_Figure9.pdf
  
  DisDec <- filter(betadf, stns == stne) %>% select(Dobs, stne, S)
  DisDec <- mutate(DisDec, S = 1-S) # distance to similarity
  
  # distance decay glm models and statistics
  brayi.glm <- glm(S ~ Dobs, family = binomial(link = log), data = filter(DisDec, stne %in% stni))
  brayo.glm <- glm(S ~ Dobs, family = binomial(link = log), data = filter(DisDec, stne %in% stno))
  
  # bray distance decay with geom_hex()
  # predict for plotting the model fit
  newd <- data.frame(Dobs = 0:52) # prediction years
  newd$fito <- exp(predict(brayo.glm, newdata = newd))
  newd$fiti <- exp(predict(brayi.glm, newdata = newd))
  
  # inner station halving distance
   -log(2)/coef(brayi.glm)[2] # 11.3 y
  # outer station halving distance
   -log(2)/coef(brayo.glm)[2] # 23.6 y
  # 
  # initial beta similarity inner bays
  exp(coef(brayi.glm)[1]) # 0.46
  # initial beta similarity outer archipelago
  exp(coef(brayo.glm)[1]) # 0.38
  
  # Figure 9
  #  outer station panel
  F8B <- ggplot(data = filter(DisDec, stne %in% stno), aes(Dobs, y=S)) + 
    geom_hex(bins = 30) + 
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'Count') + 
    geom_line(data = newd, aes(Dobs, fito), size = 1.3, col = 'black') +
    geom_vline(xintercept = -log(2)/coef(brayo.glm)[2], size = 1, col = 'red') +
    theme(legend.position = c(0.85, 0.77), legend.background = element_rect(fill = 0)) +
    ylim(0, 1) + 
    xlab("Years") + 
    ylab("Bray-Curtis similarity") +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),legend.text = element_text(size = 14), legend.title = element_text(size = 14)) 
  
  #  inner station panel
  F8A <- ggplot(data = filter(DisDec, stne %in% stni), aes(Dobs, y=S)) + 
    geom_hex(bins = 30) + 
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'Count') + 
    geom_line(data = newd, aes(Dobs, fiti), size = 1.3, col = 'black') +
    geom_vline(xintercept = -log(2)/coef(brayi.glm)[2], size = 1, col = 'red') +
    theme(legend.position = c(0.85, 0.77), legend.background = element_rect(fill = 0)) +
    ylim(0, 1) + 
    xlab("Years") + 
    ylab("Bray-Curtis similarity") +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),legend.text = element_text(size = 14), legend.title = element_text(size = 14)) 
  
  Fig8 <- plot_grid(F8B, F8A, labels = "AUTO", label_size = 24) # glue panels together
  
  cowplot::save_plot("./Olli_HelSummerTrends_Figure8.pdf", Fig8, ncol = 2, base_height = 6, base_width = 6)
  cowplot::save_plot("./Olli_HelSummerTrends_Figure8.svg", Fig8, ncol = 2, base_height = 6, base_width = 6)
}  # saves Fig 8 beta diversity trends as Olli_HelSummerTrends_Figure8.pdf


###  Fig. S2 GAM EFFECT SIZE ####

if(FALSE){ # saves Fig S2 beta diversity trends as Olli_HelSummerTrends_FigureS2.pdf
  
  ctrl <- list(niterEM = 0, msVerbose = FALSE, optimMethod="L-BFGS-B")
  FigS2_meta <- filter(meta, pel == 1) %>% rename('doy' = jul)
  
  mod <- gamm(N0 ~ s(time, bs = 'tp') + s(doy),  data = FigS2_meta, correlation = corCAR1(form = ~ 1|year), control = ctrl, method="REML")
  
  FigS2 <- gratia::draw(mod$gam)
  cowplot::save_plot("./Olli_HelSummerTrends_FigureS2.pdf", FigS2, base_height = 3, base_width = 6)
  cowplot::save_plot("./Olli_HelSummerTrends_FigureS2.svg", FigS2, base_height = 3, base_width = 6)
  
} # saves Fig S2 beta diversity trends as Olli_HelSummerTrends_FigureS2.pdf


#### The effect of microscopists ####

# to what extent the change of microscopists affects the community structure?
# we fit a gam model to the data and regress the microscopist on the residuals of the gam model. This tells how much of the residual variance can be explained by the subjective judgement of the microscopists. 

if(T){ # switching the microscopist - did it change community structure?
  # amend with microscopist code - variable Det
  meta$Det <- ifelse(meta$year %in% c(1970, 1972), 'A', NA) # 1970, 1972 have been counted 100% by one specific person
  meta$Det <- ifelse(meta$year == 1983, 'B', meta$Det) # 1983 that was counted 100 % by another person
  meta$Det <- ifelse(meta$year %in% 1970:1982 & is.na(meta$Det), 'C', meta$Det) # 1983-1987
  meta$Det <- ifelse(meta$year %in% 1984:1987, 'D', meta$Det) # 1983-1987
  meta$Det <- ifelse(meta$year %in% 1988:1990, 'E', meta$Det) #  1988-1990
  meta$Det <- ifelse(meta$year %in% 1991:1992, 'F', meta$Det) #  1991-1992
  meta$Det <- ifelse(meta$year == 1993, 'G', meta$Det) #
  meta$Det <- ifelse(meta$year >= 1994, 'H', meta$Det) # from 1994 onwards most samples have been counted by one specific person
  meta$Det <- ifelse(meta$year >= 2008 & meta$stn %in% c(4, 25, 87), 'I', meta$Det) # with the most notable exceptions being stations 4, 25 and 87 that have been counted from approximately 2008 onwards by another specific person.
  meta$Det <- ifelse(meta$year %in% 1966:1969 & is.na(meta$Det), 'J', meta$Det) #  samples from 60's were counted by one person
  
  # re-code microscopist (Det) and station to factors
  meta$fDet <- factor(meta$Det)
  meta$fstn <- factor(meta$stn)
  
  gmod1 <- gamm(NMDS1 ~ s(time, bs = 'tp') + s(jul) + s(fstn, bs = 're'),  data = filter(meta, pel == 1),  method = "REML")
  gmodResid1 <- residuals(gmod1$gam)
  lm(gmodResid1~Det,  data = filter(meta, pel == 1)) %>% summary()
  # in the outer archipelago R-squared is < 5%
  
  gmod0 <- gamm(NMDS1 ~ s(time, bs = 'tp') + s(jul) + s(fstn, bs = 're'),  data = filter(meta, pel == 0),  method = "REML")
  gmodResid0 <- residuals(gmod0$gam)
  lm(gmodResid0~Det,  data = filter(meta, pel == 0)) %>% summary()
  # in the inner bays the adjusted R-squared is ca 0.6% 

} # switching the microscopist - did it change community structure?


#### MANTEL AND ADONIS BLOCK ####

# mantel and adonis profiles as predictors of community structure #
# needed for Table 1 and Fig S3

if(T){ # mantel/adonis block for Table 1 and Fig S3
  # mantel/adonis block for Table 1 and Fig S3
  # outer archipelago and inner bay phytoplankton communities are in datcm1 [3515 x 547], datcm0 [1115 x 506]
  # env variables in meta
  # dissimilarities of inner and outer communities
  
  # Table 1 env.var needs: 
  env.var <- c('year','jul','Ptot','Ntot','PO4','NO3','NH4','pH','Salin','Temp','chla','N0','E10', grep('Dqz|RenyiD', names(meta), value = T))
  
  if(F){ # long calculations
    MANTELl.lst <- list()
    for(i in 1:length(MXl)){
      com <- MXl[[i]]
      print(i)
      MANTELl.lst[[i]] <- mclapply(meta[rownames(com), env.var], function(x){ # takes env dataframe columns one at a time as list elements
        idx <- which(!is.na(x)) # check if NAs exist
        env.dis <- x[idx] %>% decostand(method = 'standardize') %>% stats::dist()
        com <- com[idx,] # use only samples, with non-NA environmental record
        com <- com[ ,colSums(com) > 0] # delete species with no occurrence
        com.dis <- vegdist(wisconsin(sqrt(com)),"bray") # use bray distance similar to NMDS
        mantel(com.dis, env.dis, perm = 1, parallel = 1)
      }, mc.cores = detectCores())
    } # 230 sek
    
    
    ADONIS.lst <- list()
    for(i in 1:length(MXl)){
      com <- MXl[[i]]
      print(i)
      ADONIS.lst[[i]] <- mclapply(meta[rownames(com), env.var], function(x){ # takes env dataframe columns one at a time as list elements
        idx <- which(!is.na(x)) # check if NAs exist
        env.df <- data.frame(var = x[idx])
        com <- com[idx, ] # use only samples, with non-NA environmental record
        com <- com[ ,colSums(com) > 0] # delete species with no occurrence
        com.dis <- vegdist(wisconsin(sqrt(com)),"bray") # use bray distance similar to NMDS
        adonis2(com.dis~var, data = env.df, perm = 1, parallel = 1)
      }, mc.cores = detectCores())
    }
    
    save(MANTELl.lst, ADONIS.lst , file = 'Olli_HelSummerTrends_mantel_adonis.rda')
    
  } # long calculations
  load(file = 'Olli_HelSummerTrends_mantel_adonis.rda')

} ## mantel/adonis block for Table 1 and Fig S3

#### Fig. S3 DIVERSITY PROFILES AS CORRELATES WIHT ORDINATION ####

if(T){ # saves Fig S3 beta diversity trends as Olli_HelSummerTrends_FigureS3.pdf
  
  # extract diversity profile statistics
  envfitl.r <- lapply(ENVFITl.lst, function(x){sapply(x, '[', 'vectors') %>% sapply('[', 'r') %>% unlist()}) # length 159
  mantel.r <- lapply(MANTELl.lst, function(x){sapply(x, '[', 'statistic') %>% sapply('[', 1) %>% unlist()})
  adonis.R2 <- lapply(ADONIS.lst, function(x){sapply(x, '[', 'R2') %>% sapply('[', 1) %>% unlist()})
  
  # ENVFIT panels
  ef <- rbind(data.frame(value = c(envfitl.r[[1]][grep('RenyiD', names(envfitl.r[[1]]))],
              envfitl.r[[1]][grep('DqzPD', names(envfitl.r[[1]]))]),
              io = 'Outer Archipelago',
              metric = factor(rep(c(1,2), each = length(q0)), ordered = TRUE, labels = c('Species diversity','Phylogenetic diversity'))),
              data.frame(value = c(envfitl.r[[2]][grep('RenyiD', names(envfitl.r[[2]]))],
              envfitl.r[[2]][grep('DqzPD', names(envfitl.r[[2]]))]),
              io = 'Inner Bays',
              metric = factor(rep(c(1,2), each = length(q0)), ordered = TRUE, labels = c('Species diversity','Phylogenetic diversity'))))
  
  ef$q0 <- rep(q0, nrow(ef)/length(q0))
  
  enpl <- ggplot(data = ef, aes(x = q0, y = value,  group = io)) + 
    geom_vline(xintercept = 0, alpha = 0.5) + 
    geom_line(aes(colour = io), size = 1) + 
    facet_wrap(vars(metric), scales = 'free_y') + 
    xlab("Order of diversity (q)") + 
    ylab("Coef. of determination") +
    theme(legend.position = c(0.8, 0.17), legend.background = element_rect(fill = 0)) + 
    labs(color = "")
  
  # MANTEL panels
 
  ma <- rbind(data.frame(value = c(mantel.r[[1]][grep('RenyiD', names(mantel.r[[1]]))],
                DqzPD = mantel.r[[1]][grep('DqzPD', names(mantel.r[[1]]))]),
                io = 'Outer Archipelago',
                metric = factor(rep(c(1,2), each = length(q0)), ordered = TRUE, labels = c('Species diversity','Phylogenetic diversity'))),
                data.frame(value = c(mantel.r[[2]][grep('RenyiD', names(mantel.r[[2]]))],
                DqzPD = mantel.r[[2]][grep('DqzPD', names(mantel.r[[2]]))]),
                io = 'Inner Bays',
                metric = factor(rep(c(1,2), each = length(q0)), ordered = TRUE, labels = c('Species diversity','Phylogenetic diversity'))))
  ma$q0 <- rep(q0, nrow(ma)/length(q0))
  
  mapl <- ggplot(data = ma, aes(x = q0, y = value,  group = io)) + 
    geom_vline(xintercept = 0, alpha = 0.5) + 
    geom_line(aes(colour = io), size = 1) + 
    facet_wrap(vars(metric), scales = 'free_y') + 
    xlab("Order of diversity (q)") + 
    ylab("Pearson correlation") +
    theme(legend.position = c(0.8, 0.17), legend.background = element_rect(fill = 0)) + 
    labs(color = "")
  
  # ADONIS panesl
 
  ad <- rbind(data.frame(value = c(adonis.R2[[1]][grep('RenyiD', names(adonis.R2[[1]]))],
              DqzPD = adonis.R2[[1]][grep('DqzPD', names(adonis.R2[[1]]))]),
              io = 'Outer Archipelago',
              metric = factor(rep(c(1,2), each = length(q0)), ordered = TRUE, 
              labels = c('Species diversity','Phylogenetic diversity'))),
              data.frame(value = c(adonis.R2[[2]][grep('RenyiD', names(adonis.R2[[2]]))],
              DqzPD = adonis.R2[[2]][grep('DqzPD', names(adonis.R2[[2]]))]),
              io = 'Inner Bays',
              metric = factor(rep(c(1,2), each = length(q0)), ordered = TRUE, labels = c('Species diversity','Phylogenetic diversity'))))
  
  ad$q0 <- rep(q0, nrow(ad)/length(q0))
  
  adpl <- ggplot(data = ad, aes(x = q0, y = value,  group = io)) + geom_vline(xintercept = 0, alpha = 0.5) + 
    geom_line(aes(colour = io), size = 1) + 
    facet_wrap(vars(metric), scales = 'free_y') + 
    xlab("Order of diversity (q)") + 
    ylab("Coef. of determination") +
    theme(legend.position = c(0.8, 0.17), legend.background = element_rect(fill = 0)) + 
    labs(color = "")
  
  FigS3 <- plot_grid(enpl, mapl, adpl, labels = "AUTO", label_size = 18,  ncol = 1)
  
  cowplot::save_plot("./Olli_HelSummerTrends_FigureS3.pdf", FigS3, base_height = 8, base_width = 5)
  cowplot::save_plot("./Olli_HelSummerTrends_FigureS3.svg", FigS3, base_height = 8, base_width = 5)
  
} # saves Fig S3 beta diversity trends as Olli_HelSummerTrends_FigureS3.pdf

### TABLE 1 ####

if(T){ # Table 1
  # extract all statistics from envfit, mantel and adonis lists
  envfitl.r <- lapply(ENVFITl.lst, function(x){sapply(x, '[', 'vectors') %>% sapply('[', 'r') %>% unlist()}) 
  mantel.r <- lapply(MANTELl.lst, function(x){sapply(x, '[', 'statistic') %>% sapply('[', 1) %>% unlist()})
  adonis.R2 <- lapply(ADONIS.lst, function(x){sapply(x, '[', 'R2') %>% sapply('[', 1) %>% unlist()})
  
  # external variable names of interest
  env_var <- c('year', 'jul', 'Ptot', 'Ntot','chla', 'PO4', 'NO3', 'NH4', 'pH', 'Salin', 'Temp', 'N0','E10','DqzPD0')
  
  # extract the variable of interest and round to two significant numbers
  T1e <- envfitl.r %>% sapply ('[') %>% .[paste('vectors.r.',env_var, sep = ''),] %>% signif(2)
  T1m <- mantel.r %>% sapply ('[') %>% .[paste(env_var, '.statistic', sep = ''), ] %>% signif(2)
  T1a <- adonis.R2 %>% sapply ('[') %>% .[paste(env_var, '.R2', sep = ''), ] %>% signif(2)
  
  # combine the results forom the three different methods (envfit, mantel, adonis) into one final table, add row and column names   
  T1 <- cbind(T1e, T1m, T1a) 
  dimnames(T1) <- list(c('Year', 'Season', 'Ptot', 'Ntot','Chl a', 'PO4', 'NO3', 'NH4', 'pH', 'Salinity', 'Temperature', 'Species richness','Species evenness','Phylogenetic diversity'), c('Envfit outer','Envfit inner', 'Mantel outer','Mantel inner' ,'Adonis outer','Adonis inner'))
  
  # The final Table 1
  stargazer(T1, type = 'text', title = 'Table 1')
  
} # Table 1

#### HIERARCHICAL GAMM BLOCK; Figs S4, S5, Table S1 ####

if(T){ # Hierarchical GAM; provides Table S1, saves Olli_HelSummerTrends_FigureS4.pdf, Olli_HelSummerTrends_FigureS5.pdf
  
  ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
  # easier to split the inner bay and outer archipelago part of meta into list elements
  META.lst <- list()
  META.lst[['in']] <- filter(meta, pel == 0) %>% mutate(fstn = factor(stn))
  META.lst[['out']] <- filter(meta, pel == 1) %>% mutate(fstn = factor(stn))
  
  XXgamm <- function(var='NMDS1'){gm <- res <- list()
  for(i in 1:2){
    print(var)
    gm[['0']] <-  eval(parse(text=paste("(gamm(",var," ~ s(time, bs = 'tp') + s(jul),  data = META.lst[[i]], control = ctrl, method='REML'))", sep = ''))) # null model
    gm[['C']]   <- eval(parse(text=paste("(gamm(",var," ~ s(time, bs = 'tp') + s(jul), data = META.lst[[i]], correlation = corCAR1(form = ~ 1|year), control = ctrl, method='REML'))", sep = ''))) # CAR1
    gm[['G']] <- eval(parse(text=paste("(gamm(",var," ~ s(time, bs = 'tp') + s(jul) + s(fstn, bs = 're'), data = META.lst[[i]], correlation = corCAR1(form = ~ 1|year), control = ctrl, method='REML'))", sep = ''))) # G + CAR1
    gm[['GS']] <- eval(parse(text=paste("(gamm(",var," ~ s(time, bs = 'tp', m = 2) + s(jul) + t2(time, fstn, bs = c('cr','re'), m = 2), data = META.lst[[i]], correlation = corCAR1(form = ~ 1|year), control = ctrl, method='REML'))", sep = ''))) # GS + CAR
    gm[['GI']] <- eval(parse(text=paste("(gamm(",var," ~ s(time, bs = 'tp', m = 2) + s(jul) + t2(time, by = fstn, m = 1, bs ='tp') + s(fstn, bs = 're'), data = META.lst[[i]], correlation = corCAR1(form = ~ 1|year), control = list(niterEM = 0, msVerbose = TRUE, optimMethod='L-BFGS-B', msMaxIter = 200), method='REML'))", sep = ''))) # GI + CAR
    
    res[[i]] <- gm
  }
  return(res)
  }
  
  # splits a lot of output
hgam_N0 <- XXgamm('N0')
hgam_NMDS <- XXgamm('NMDS1')
hgam_E10 <- XXgamm('E10')
  
## TABLE S1 ####

# for Table S1 should extract AIC, RMSE, and R2
  hgmf <- function(x){data.frame(
    AIC = sapply(x, '[', 1) %>% sapply(AIC) %>% round(),
    RMSE = sapply(x, function(x){round(sqrt(sum(x$gam$residuals^2)),2)}),
    R2 = sapply(x, '[', 2) %>% lapply(summary) %>% sapply('[', 'r.sq') %>% unlist() %>% signif(3)
  )}

  tableS1 <- rbind(cbind(hgmf(hgam_NMDS[[1]]), hgmf(hgam_NMDS[[2]])), cbind(hgmf(hgam_N0[[1]]), hgmf(hgam_N0[[2]])), cbind(hgmf(hgam_E10[[1]]), hgmf(hgam_E10[[2]])))
  
  # format column and row names
  colnames(tableS1) <- paste(rep(c('Outer', 'Inner'), each = 3), colnames(tableS1))
  rownames(tableS1) <- paste(rep(c('NMDS1','NSP','E10'), each = 5), rep(c('0', 'C', 'G', 'GS', 'GI'), 3), sep = '_')
  
  # Table S1
  stargazer(as.matrix(tableS1) , type = 'text', title = 'Table S1')
  
  # hypothesis testing - compare the hierarchical GAM models with anova
  i <- 1 # inner bays; chagne i <- 2 for outer archipelago 
  anova(hgam_NMDS[[i]][['0']]$lme, hgam_NMDS[[i]][['C']]$lme, hgam_NMDS[[i]][['G']]$lme,hgam_NMDS[[i]][['GS']]$lme, hgam_NMDS[[i]][['GI']]$lme)
  anova(hgam_N0[[i]][['0']]$lme, hgam_N0[[i]][['C']]$lme, hgam_N0[[i]][['G']]$lme,hgam_N0[[i]][['GS']]$lme, hgam_N0[[i]][['GI']]$lme)
  anova(hgam_E10[[i]][['0']]$lme, hgam_E10[[i]][['C']]$lme, hgam_E10[[i]][['G']]$lme,hgam_E10[[i]][['GS']]$lme, hgam_E10[[i]][['GI']]$lme)
  
  ## FIG S4 ####
  
  # Supplementary Fig S4 - model GI faceted gam fit NMDS1 scores
  
  dat1 <- tidymv::predict_gam(hgam_NMDS[[2]][['GS']]$gam, exclude_terms = "s(jul)", length_out = 200, values = list(jul = NULL))
  dat2 <- hgam_NMDS[[2]][['GS']]$gam$model
  
  Fig4S <- ggplot(data = dat1, aes(x = time, y = fit, group = fstn)) +
    geom_point(data = dat2, aes(x=time, y=NMDS1, group = fstn), alpha = 0.3) + 
    geom_line(aes(x = time, y = I(fit+2*se.fit)), col = 'gray') +
    geom_line(aes(x = time, y = I(fit-2*se.fit)), col = 'gray') +
    geom_line( color = 'red')  +
    facet_wrap(~fstn) + xlab('Year') + ylab('NMDS1 scores') +
    scale_x_continuous(expand = c(0.1,0.1)) # , limits = c(1965, 2020)
  
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18)) 
  
  ggsave('./Olli_HelSummerTrends_FigureS4.pdf', Fig4S, width = 6, height = 6)
  ggsave('./Olli_HelSummerTrends_FigureS4.svg', Fig4S, width = 6, height = 6)
  
  ## FIG S5 ####
  # Supplementary Fig S5 - model GI faceted gam fit species richness
  
  dat1 <- tidymv::predict_gam(hgam_N0[[2]][['GS']]$gam, exclude_terms = "s(jul)", length_out = 200, values = list(jul = NULL))
  dat2 <- hgam_N0[[2]][['GS']]$gam$model
  
  Fig5S <- ggplot(data = dat1, aes(x = time, y = fit, group = fstn)) +
    geom_point(data = dat2, aes(x = time, y=N0, group = fstn), alpha = 0.3) + 
    geom_line( color = 'red')  +
    geom_line(aes(x = time, y = I(fit+2*se.fit)), col = 'gray') +
    geom_line(aes(x = time, y = I(fit-2*se.fit)), col = 'gray') +
    facet_wrap(~fstn) + xlab('Year') + ylab('Species richness')
  
  ggsave('./Olli_HelSummerTrends_FigureS5.pdf', Fig5S, width = 6, height = 6)
  ggsave('./Olli_HelSummerTrends_FigureS5.svg', Fig5S, width = 6, height = 6)
  
} # Hierarchical GAM; provides Table S1, saves Olli_HelSummerTrends_FigureS4.pdf, Olli_HelSummerTrends_FigureS5.pdf


#### Fig S7  NUTRIENT BLOCK ####

if(TRUE){ # saves nutreint concentration figure figure Olli_HelSummerTrends_FigureS7.pdf 
  nuts <- read.table(file = 'Olli_HelSummerTrends_Nutreints.txt', head = TRUE, sep = '\t')

  # actual plots
  FigS7_N0 <- ggplot(filter(nuts, var == 'ntot', pel == 0), aes(x=time, y=value)) + 
    geom_point() + 
    geom_smooth() + 
    xlab('Year') + 
    ylab(expression(Total~nitrogen~(mu*gL**-1))) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(1970, 2020, by = 10), limits = c(1971, 2023)) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12)) 
  
  FigS7_N1 <- ggplot(filter(nuts, var == 'ntot', pel == 1), aes(x=time, y=value)) + 
    geom_point() + geom_smooth() + 
    xlab('Year') + 
    ylab(expression(Total~nitrogen~(mu*gL**-1))) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(1970, 2020, by = 10), limits = c(1971, 2023)) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12))
  
  FigS7_P0 <- ggplot(filter(nuts, var == 'ptot', pel == 0), aes(x=time, y=value)) + 
    geom_point() + geom_smooth() + 
    xlab('Year') + 
    ylab(expression(Total~phosphorus~(mu*gL**-1))) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(1970, 2020, by = 10), limits = c(1971, 2023)) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12))
  
  FigS7_P1 <- ggplot(filter(nuts, var == 'ptot', pel == 1), aes(x=time, y=value)) + 
    geom_point() + 
    geom_smooth() + 
    xlab('Year') + 
    ylab(expression(Total~phosphorus~(mu*gL**-1))) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(1970, 2020, by = 10), limits = c(1971, 2023)) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12))
  
  FigS7 <- plot_grid(FigS7_N0, FigS7_P0, FigS7_N1, FigS7_P1, labels = "AUTO", label_size = 18) # glue panels together
  cowplot::save_plot("./Olli_HelSummerTrends_FigureS7.pdf", FigS7, ncol=2, base_height = 6, base_width = 4)
  cowplot::save_plot("./Olli_HelSummerTrends_FigureS7.svg", FigS7, ncol=2, base_height = 6, base_width = 4)
} # saves nutreint concentration figure figure Olli_HelSummerTrends_FigureS7.pdf 
