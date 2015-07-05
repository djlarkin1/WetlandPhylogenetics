library(ape)
library(dplyr)
library(phytools)
library(phangorn)
library(adephylo)

# Tree derived from the Tank phylogeny
# Created in Dropbox\ExSituPhylogenetics\R-ExSituPhylogenetics\A.TankTree.R
load("wes.tre.RData")
wes.tre #contains 234 tips 

# Dichanthelium_boreale misplaced
wetland.tre <- wes.tre
wetland.tre <- drop.tip(wetland.tre, "Dichanthelium_boreale")
tip.node <- grep("Poaceae_sp.", wetland.tre$tip)
node.i <- Ancestors(wetland.tre, tip.node, "parent")
clade.i <- extract.clade(wetland.tre, node.i)
depth.i <- distRoot(clade.i, tip=wetland.tre$tip[tip.node])
wetland.tre <- bind.tip(wetland.tre, "Dichanthelium_boreale", where = node.i, 
                             edge.length = depth.i)

# Carex spp. divided
carex.vec <- c("Carex_atherodes", "Carex_lacustris", "Carex_pellita",
               "Carex_stipata_var._stipata", "Carex_tenera", 
               "Carex_utriculata", "Carex_debilis_var._rudgei", 
               "Carex_lasiocarpa_var._americana", 
               "Carex_scoparia_var._scoparia", "Carex_sp.")
wetland.tre <- drop.tip(wetland.tre, carex.vec)

for(j in 1:length(carex.vec)) {
  wetland.tre <- add.species.to.genus(wetland.tre, carex.vec[j])
  print(paste("species", j, "of", length(carex.vec)))
}

# missing Cirsium sp.
wetland.tre <- add.species.to.genus(wetland.tre, "Cirsium_sp.")

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
targetSpp[targetSpp %in% wetland.tre$tip == F]

pdf(file = "phylogenyPlot.pdf", width = 7, height = 18)
plot(wetland.tre, no.margin = T, cex = 0.5)
dev.off()

save(wetland.tre, file = "wetland.tre.RData")