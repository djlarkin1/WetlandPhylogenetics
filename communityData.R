library(dplyr)
library(reshape2)
library(vegan)

# Reconciling species names -----------------------------------------------

# Wes' original spp table
speciesList <- read.csv("Species_Table_USDA_Family.csv", header = T)
speciesList <- speciesList[,c(1:3,7,9,11,13:15)]
colnames(speciesList) <- c("origSpp", "sppCode", "usdaSpp",
                           "NativeWI", "Life_Cycle", "Habit", "C_Value",
                           "Wetland_Indicator_Status", "Invasive")

# names that Dan used to query TNRS
tnrsSpecies <- read.csv("WesSppForTNRS.csv", header = T)
tnrsSpecies$Original %in% speciesList$usdaSpp
speciesList <- full_join(speciesList, tnrsSpecies, by = c("usdaSpp" = "Original"))
colnames(speciesList)[colnames(speciesList)=="Rev"] <- "subToTnrsSpp"

# TNRS reconciled names
tnrsTable <- read.csv("wesTnrsSimple.csv", header = T)
tnrsTable <- tnrsTable[,c(1,6:8)]
speciesList <- full_join(speciesList, tnrsTable, by = c("subToTnrsSpp" = "Name_submitted"))
phyName <- gsub(" ", "_", as.character(speciesList$Accepted_name))
speciesList <- mutate(speciesList, phyName = phyName)


# Organizing community data -----------------------------------------------

# Three plot sizes: s/m/l
# Abundance for all
# Height for s & l

# Full veg data
vegData <- read.csv("Vegetation_Table_Complete.csv", header = T)
vegData <- mutate(vegData, Plot_Size = tolower(Plot_Size), 
                  Plot_ID = tolower(Plot_ID))
vegData <- select(vegData, -Plot, -Large_Plot, -Plot_Type, -Plot_ID2)
vegData <- mutate(vegData, Sample_ID = paste(SSU, Plot_ID, sep = "_"))

# community data frames
comm.abund <- select(vegData, Sample_ID, SSU, Plot_Size, Species, Abundance)
comm.abund <- dcast(comm.abund, Sample_ID + SSU + Plot_Size ~ Species, 
                    mean, fill = 0)

comm.height <- select(vegData, Sample_ID, SSU, Plot_Size, Species, Height)
comm.height <- dcast(comm.height, Sample_ID + SSU + Plot_Size ~ Species, 
                    mean, fill = 0)

# Site-level spp richness
siteRichness <- select(vegData, SSU, Species)
siteRichness <- dcast(siteRichness, SSU ~ Species)
siteRichness <- specnumber(siteRichness[,-1])




save(list = c("speciesList", "vegData", "comm.abund", 
              "comm.height", "siteRichness"), file = "communityData.RData")
