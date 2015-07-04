library(dplyr)

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

save(list = c("speciesList"), file = "communityData.RData")