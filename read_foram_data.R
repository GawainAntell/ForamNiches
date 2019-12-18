library(phytools)
library(paleoPhylo)
library(ape)

occ <- read.csv('Data/foram_uniq_occs_latlong_8ka_191218.csv', stringsAsFactors=FALSE)
bins <- unique(occ$bin)

# Limit analysis to species included in tree from Aze et al. 2011
trRaw <- read.csv('Data/Aze_et_al_2011_bifurcating_tree_data.csv', stringsAsFactors=FALSE)
trRaw$Species.name <- gsub('Globigerinoides sacculifer', 'Trilobatus sacculifer', trRaw$Species.name)
trRaw$Species.name <- gsub('Globigerinoides trilobus', 'Trilobatus trilobus', trRaw$Species.name)
tr_paleo <- with(trRaw, 
                 as.paleoPhylo(Species.code, Ancestor.code, Start.date, End.date)
)
tr <- buildApe(tr_paleo)
sppCodes <- tr$tip.label
rowOrdr <- match(sppCodes, trRaw$Species.code)
tr$tip.label <- trRaw$Species.name[rowOrdr]

# 10 species are not present in the phylogeny, mostly because microporiferate
sppAll <- unique(occ$species)
lostSpp <- setdiff(sppAll, tr$tip.label)

# Beella megastoma is arguably the same species as B. digitata,
# and Truncorotalia crassula is arguably senior synonym to crassaformis
# (Schiebel and Hemleben 2017). The depth ranges for both are unknown.
lostSpp <- c(lostSpp, 'Beella megastoma', 'Truncorotalia crassula')

spp <- setdiff(sppAll, lostSpp)
rows2toss <- ! occ$species %in% spp
occ <- occ[!rows2toss,]
row.names(occ) <- as.character(1:nrow(occ))
