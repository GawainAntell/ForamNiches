library(phytools)
library(paleoPhylo)
library(ape)
library(geiger)
library(nlme)
library(xtable)

ss <- TRUE

# Data prep ---------------------------------------------------------------

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y-%m-%d')

# Foram niche data
spAttr <- read.csv('Data/foram_spp_data_20-04-05.csv', stringsAsFactors=FALSE)
if (ss){
  df <- read.csv('Data/foram_niche_sumry_metrics_0m_20-04-05.csv', stringsAsFactors=FALSE)
} else {
  df <- read.csv('Data/foram_niche_sumry_metrics_20-04-05.csv', stringsAsFactors=FALSE)
}
ordr <- order(df$bin, decreasing = TRUE)
df <- df[ordr,]
bins <- unique(df$bin)
spp <- unique(df$sp)
nspp <- length(spp)

# Phylogenetic tree
trRaw <- read.csv('Data/Aze_et_al_2011_bifurcating_tree_data.csv', stringsAsFactors=FALSE)
trRaw$Species.name <- gsub('Globigerinoides sacculifer', 'Trilobatus sacculifer', trRaw$Species.name)
trRaw$Species.name <- gsub('Globigerinoides trilobus', 'Trilobatus trilobus', trRaw$Species.name)
trPaleo <- with(trRaw, 
                 as.paleoPhylo(Species.code, Ancestor.code, Start.date, End.date)
)
tr <- buildApe(trPaleo)
sppCodes <- tr$tip.label
rowOrdr <- match(sppCodes, trRaw$Species.code)
tr$tip.label <- trRaw$Species.name[rowOrdr]

# Summarise niche data for each species as the median among intervals, to plot
spSmry <- function(s, dat){
  sBool <- dat$sp==s
  sDf <- dat[sBool,]
  numCol <- ! colnames(dat) %in% c('bin','sp')
  out <- apply(sDf[,numCol], 2, median, na.rm=TRUE)
  out <- data.frame(t(out))
  out$nInt <- nrow(sDf)
  out$species <- s
  out
}
spL <- lapply(spp, spSmry, dat=df)
spDf <- do.call(rbind, spL)

# Combine with species attribute data
spDf <- merge(spDf, spAttr, by='species')
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
tblNm <- paste0('Figs/species_attribute_table_0m_',day,'.tex')
if (ss){
  print(tblx, file=tblNm, include.rownames=FALSE, 
        floating = FALSE, tabular.environment = 'longtable')
}

# Phylogenetic comparison -------------------------------------------------

# Ensure uniform species name set in data and tree
toss <- name.check(tr, spDf$species, data.names=spp)
noPhy <- toss$data_not_tree
if (length(noPhy) > 0){
  rows2toss <- spDf$species %in% noPhy
  spDf <- spDf[!rows2toss,]
  paste(paste(noPhy, collapse = ' '), 'removed')
  spp <- setdiff(spp, noPhy)
} 
trTrim <- drop.tip(tr, toss$tree_not_data)
name.check(trTrim, spDf$species, data.names=spp)

# use nlme method - caper oddly drops 8 tree tips

row.names(spDf) <- spDf$species
spDf <- spDf[trTrim$tip.label,]

bm <- corBrownian(1, trTrim)
bmMod <- gls(h~eco, data=spDf, correlation=bm)
summary(bmMod)

pgl <- corPagel(1, trTrim)
pglMod <- gls(h~eco, data=spDf, correlation = pgl)
summary(pglMod)

# ANOVA

h <- spDf$h
eco <- spDf$eco
names(eco) <- names(h) <- spDf$species
aovRes <- aov.phylo(h ~ eco, trTrim, nsim=10000)

coefs <- summary(aovRes)$coefficients
ecoEsts <- sprintf('%.3f', round(coefs[-1,'Estimate'], 3))
ecoErr <- sprintf('%.3f', round(coefs[-1,'Std. Error'], 3))
intrcpt <- coefs['(Intercept)',c('Estimate','Std. Error')]
intrcptRnd <- sprintf('%.3f', round(intrcpt, 3))
interceptRow <- cbind('Intercept', intrcptRnd[1], intrcptRnd[2])
ecoRows <- cbind(Ecotype=ecoLbl, Coef=c(0, ecoEsts), SE=c('-', ecoErr))
aovTbl <- rbind(interceptRow, ecoRows)
aovX <- xtable(aovTbl)
if (ss){
  aovNm <- paste0('Figs/phylo_ANOVA_table_0m_',day,'.tex')
} else {
  aovNm <- paste0('Figs/phylo_ANOVA_table_hab_',day,'.tex')
}
print(aovX, file=aovNm, include.rownames=FALSE)
