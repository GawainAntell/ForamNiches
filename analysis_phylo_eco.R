library(phytools)
library(paleotree)
library(ape)
library(geiger)
library(nlme)
library(xtable)

ss <- TRUE
day <- as.Date(date(), format="%a %b %d %H:%M:%S %Y")

# Niche data --------------------------------------------------------------

spAttr <- read.csv('Data/foram-spp-data_2020-07-21.csv')
if (ss){
  # using RT bandwidth gives lambda = 0; SJ, 1
 # df <- read.csv('Data/niche-sumry-metrics_nrd0_SS_2020-07-22.csv')
  df <- read.csv('Data/niche-sumry-metrics_SJ-ste_SS_2020-07-24.csv')
} else {
 #  df <- read.csv('Data/niche-sumry-metrics_nrd0_hab_2020-07-22.csv')
  df <- read.csv('Data/niche-sumry-metrics_SJ-ste_hab_2020-07-24.csv')
}
ordr <- order(df$bin, decreasing = TRUE)
df <- df[ordr,]
bins <- unique(df$bin)
spp <- unique(df$sp)

# Phylogenetic tree -------------------------------------------------------

data(macroperforateForam)

# 2 species are assigned to a different genus in the phylogeny
foramAMb$Species_name <- gsub('Globigerinoides_sacculifer', 'Trilobatus_sacculifer', 
                              foramAMb$Species_name)
foramAMb$Species_name <- gsub('Globigerinoides_trilobus', 'Trilobatus_trilobus', 
                              foramAMb$Species_name)

# from paleotree example: convert Aze et al's suppmat to paleotree-readable format
createTaxaData <- function(table){
  #reorder table by first appearance
  timetable <- table[order(-as.numeric(table[,3])),]
  ID <- 1:nrow(table)
  anc <- sapply(table[,2], function(x){
    if(!is.na(x)){
      which(x == table[,1])
    } else { NA }
  })
  stillAlive <- as.numeric(table[,4] == 0)
  ages <- cbind(as.numeric(table[,3]),
                as.numeric(table[,4]))
  res <- cbind(ID,anc,ages,stillAlive,ID)
  colnames(res) <- c('taxon.id','ancestor.id','orig.time',
                     'ext.time','still.alive','looks.like')
  rownames(res) <- table[,1]
  return(res)
}

# warning: slow!
taxaAMb <- createTaxaData(foramAMb)
treeAMb <- taxa2phylo(taxaAMb)

# drop anagenetic zero-length-terminal-edge ancestors
cleanAMb <- dropZLB(treeAMb)

# rename tips from alphanumeric codes to species names
sppCodes <- cleanAMb$tip.label
rowOrdr <- match(sppCodes, foramAMb$Species_code)
cleanAMb$tip.label <- foramAMb$Species_name[rowOrdr]

# Summarise H for each species --------------------------------------------

spSmry <- function(s, dat){
  sBool <- dat$sp==s
  sDf <- dat[sBool,]
  numCol <- ! colnames(dat) %in% c('bin','sp')
  out <- apply(sDf[,numCol], 2, mean, na.rm=TRUE)
  out <- data.frame(t(out))
  out$nInt <- nrow(sDf)
  out$species <- s
  out
}
spL <- lapply(spp, spSmry, dat=df)
spDf <- do.call(rbind, spL)

# Combine with species attribute data
spDf <- merge(spDf, spAttr, by = 'species')
# set open ocean thermocline species as the reference level
spDf$eco <- factor(spDf$eco, levels=c(3,1,2,4,5))
ecoLbl <- c('Thermocline','Mix layer symbiotic','Mix layer heterotroph',
            'Subthermocline','High latitude')
# From Aze et al. 2011, ecotype codes are:
# 1 = open ocean, mixed layer, trop/subtrop, w symbionts
# 2 = open ocean, mixed layer, trop/subtrop, w/o symbionts
# 3 = open ocean, thermocline
# 4 = open ocean, sub-thermocline
# 5 = high lat
# 6 = upwelling/high productivity

# Supplemental table ------------------------------------------------------

alph <- order(spDf$species)
cols2export <- c('species','eco','ecoRef','DepthHabitat','ref','h','pe')
tbl <- spDf[alph,cols2export]

tbl$ecoRef <- gsub('\\(','',tbl$ecoRef)
tbl$ecoRef <- gsub(')','',tbl$ecoRef)
tbl$ecoRef <- gsub('&','and',tbl$ecoRef)
tbl$ecoRef <- gsub(';',',',tbl$ecoRef)
tbl$ecoRef <- gsub('This study','Fenton et al. in prep',tbl$ecoRef)
tbl$ecoRef <- gsub('2001b','2001',tbl$ecoRef)
tbl$ecoRef <- gsub('1981.','1981',tbl$ecoRef)

tbl$DepthHabitat <- gsub('Surface.subsurface','Surface-subsurf',tbl$DepthHabitat)
code2habitat <- function(x){
  switch(paste(x), '1' = 'Mix layer symbiotic', '2' = 'Mix layer heterotroph', 
         '3' = 'Thermocline', '4' = 'Subthermocline', '5' = 'High latitude', 'NA' = '-')
} 
tbl$eco <- sapply(tbl$eco, code2habitat)

tbl$h <- sprintf('%.2f', round(tbl$h, 2))
tbl$pe <- sprintf('%.1f', round(tbl$pe, 1))
colnames(tbl) <- c('Species','Ecotype','Eco. ref.','Depth','Depth ref.','H','Optimum')
tblx <- xtable(tbl, align=c(rep('l',6),'r','r'))
tblNm <- paste0('Figs/species-attribute-table_SS_',day,'.tex')
if (ss){
  print(tblx, file=tblNm, include.rownames = FALSE, 
        floating = FALSE, tabular.environment = 'longtable')
}

# Phylogenetic comparison -------------------------------------------------
# http://www.phytools.org/Cordoba2017/ex/4/PGLS.html

# Ensure uniform species name set in data and tree
spDf$species <- gsub(' ', '_', spDf$species)
toss <- name.check(cleanAMb, spDf$species, data.names = spDf$species)
noPhy <- toss$data_not_tree
if (length(noPhy) > 0){
  rows2toss <- spDf$species %in% noPhy
  spDf <- spDf[!rows2toss,]
  paste(paste(noPhy, collapse = ' '), 'removed')
} 
trTrim <- drop.tip(cleanAMb, toss$tree_not_data)
#  dropPaleoTip(treeAMb, toss$tree_not_data)
name.check(trTrim, spDf$species, data.names = spDf$species)

# TODO consider exporting as a supplemental figure
plot(trTrim, main ='timetreeAMb', show.tip.label = T)

# use nlme method - caper oddly drops 8 tree tips

pgl <- corPagel(1, form = ~ species, trTrim)
pglMod <- gls(h ~ eco, data = spDf, correlation = pgl)
summary(pglMod)

# ANOVA

h <- spDf$h
eco <- spDf$eco
names(eco) <- names(h) <- spDf$species
aovRes <- aov.phylo(h ~ eco, trTrim, nsim = 10000)

coefs <- summary(aovRes)$coefficients
ecoEsts <- sprintf('%.3f', round(coefs[-1, 'Estimate'],  3))
ecoErr <-  sprintf('%.3f', round(coefs[-1,'Std. Error'], 3))
intrcpt <- coefs['(Intercept)', c('Estimate','Std. Error')]
intrcptRnd <- sprintf('%.3f', round(intrcpt, 3))
interceptRow <- cbind('Intercept', intrcptRnd[1], intrcptRnd[2])
ecoRows <- cbind(Ecotype = ecoLbl, Coef = c(0, ecoEsts), SE = c('-', ecoErr))
aovTbl <- rbind(interceptRow, ecoRows)
aovX <- xtable(aovTbl)
if (ss){
  aovNm <- paste0('Figs/phylo-ANOVA-table_SS_',day,'.tex')
} else {
  aovNm <- paste0('Figs/phylo-ANOVA-table_hab_',day,'.tex')
}
print(aovX, file = aovNm, include.rownames = FALSE)
