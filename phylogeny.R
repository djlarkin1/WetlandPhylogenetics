library(ape)
library(dplyr)

# Tree derived from the Tank phylogeny
# Created in Dropbox\ExSituPhylogenetics\R-ExSituPhylogenetics\A.TankTree.R
load("wes.tre.RData")
wes.tre #contains 234 tips

# Complete species list from field data
# Contains 243 species (240 useable)
speciesList <- read.csv("Species_Table_USDA_Family.csv", header = T)
# which of these species were not cross-ref'd w/ TNRS (240 spp sent to TNRS)
tnrsSpecies <- read.csv("WesSppForTNRS.csv", header = T)
speciesList$USDA_Species_Name[speciesList$USDA_Species_Name %in% tnrsSpecies$Original == FALSE]
# only non-angiosperms or unplaceable

# Which of 240 species missing from the tree (234 tips)?
# only non-angios and Lamiaceae_sp.
tnrsTable <- read.csv("wesTnrsSimple.csv", header = T)
targetSpp <- gsub(" ", "_", as.character(tnrsTable$Accepted_name))
targetSpp[targetSpp %in% wes.tre$tip == F]

pdf(file = "phylogenyPlot.pdf", width = 7, height = 18)
plot(wes.tre, no.margin = T, cex = 0.5)
dev.off()
