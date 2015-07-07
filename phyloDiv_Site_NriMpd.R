library(plyr)
library(picante)
load("communityData.RData")


# Should community or phylo distance be transformed? ----------------------

# Comm data should be sqrt-trans
hist((sqrt(unlist(site.phy.abund$comm[site.phy.abund$comm > 0]))))
# No obvious improvement for phy dist
hist((cophenetic(site.phy.abund$phy)))


# MPD and NRI calculations ------------------------------------------------

# Functions
# Calculating NRI
nri.site.fun <- function(comm, phy, site, ab.wt) {
  site.vector <- site$SSU[match(rownames(comm), site$Sample_ID)]
  x <- ddply(comm, c(.(SSU = site.vector)), "colMeans")
  y <- ses.mpd(sqrt(x[,-1]), (cophenetic(phy)), null.model = "taxa.labels", 
               runs = 100, abundance.weighted = ab.wt)
  y$mpd.obs.z <- y$mpd.obs.z * -1
  y <- mutate(y, SSU = x$SSU)
  y <- mutate(y, Type = site$Type[match(y$SSU, site$SSU)])
  return(y)
}

# Summarizing NRI
nri.summ.fun <- function(x) {
  y <- summarise(group_by(x, Type), mean(ntaxa), 
                 mean(mpd.obs), sd(mpd.obs), 
                 mean(mpd.obs.z), sd(mpd.obs.z))
}

# Site-level results
abund.site.nri_aw <- nri.site.fun(site.phy.abund$comm, site.phy.abund$phy, 
                                  plotsSites, T)
abund.site.nri_pa <- nri.site.fun(site.phy.abund$comm, site.phy.abund$phy, 
                                  plotsSites, F)

unloadNamespace("reshape2")
unloadNamespace("plyr")
library(dplyr)

# Phylogenetically: Restored more clustered
# When abund weighted, natural veer overdispersed, restored clustered
# If just pa, both appear clustered, restored more so
(abund.site.nri_summ_pa <- nri.summ.fun(abund.site.nri_pa))
(abund.site.nri_summ_aw <- nri.summ.fun(abund.site.nri_aw))


# Global phylo structure test ---------------------------------------------

library(plyr)

# Summarizing data by site type
site.vector <- plotsSites$SSU[match(rownames(site.phy.abund$comm), plotsSites$Sample_ID)]
site.abund.comm <- ddply(site.phy.abund$comm, c(.(SSU = site.vector)), "colMeans")
rownames(site.abund.comm) <- site.abund.comm[, "SSU"]
site.abund.comm <- mutate(site.abund.comm, Type = Site_Table$Type[match(site.abund.comm$SSU, Site_Table$SSU)])
# Restorations
rest.abund.comm <- droplevels(filter(site.abund.comm, Type == "Res"))
rest.abund.comm <- select(rest.abund.comm, -SSU, -Type)
cols.del <- which(colSums(rest.abund.comm) == 0)
rest.abund.comm <- rest.abund.comm[, -cols.del]
rest.abund <- match.phylo.comm(phy = site.phy.abund$phy, comm = rest.abund.comm)
# Natural wetlands
nat.abund.comm <- droplevels(filter(site.abund.comm, Type == "Nat"))
nat.abund.comm <- select(nat.abund.comm, -SSU, -Type)
cols.del <- which(colSums(nat.abund.comm) == 0)
nat.abund.comm <- nat.abund.comm[, -cols.del]
nat.abund <- match.phylo.comm(phy = site.phy.abund$phy, comm = nat.abund.comm)


phylostruct.fun <- function(com, phy, metric, type, abund.wt){
  y <- phylostruct(sqrt(com), phy, metric = metric, 
              null.model = "frequency", runs = 1000)
  y <- data.frame(metric = y$metric, null = y$null.model, 
                  mean.obs = y$mean.obs, mean.null = y$mean.null,
                  null.025 = y$quantiles.null[1], 
                  null.975 = y$quantiles.null[2],
                  phylo.structure = y$phylo.structure, type = type,
                  abund.wt = abund.wt)
  return(y)
}

rest.site.ps_pa <- phylostruct.fun(rest.abund$comm, rest.abund$phy, 
                                   metric = "psv", "rest", "pa")
rest.site.ps_aw <- phylostruct.fun(rest.abund$comm, rest.abund$phy, 
                                   metric = "pse", "rest", "aw")
nat.site.ps_pa <- phylostruct.fun(nat.abund$comm, nat.abund$phy, 
                                   metric = "psv", "nat", "pa")
nat.site.ps_aw <- phylostruct.fun(nat.abund$comm, nat.abund$phy, 
                                   metric = "pse", "nat", "aw")
phylostruct.results <- rbind(nat.site.ps_pa, nat.site.ps_aw, 
                             rest.site.ps_pa, rest.site.ps_aw)




##################
# Next steps
# Global p-val/null model test (generate 1000 null comms, calculate mpd)
# Check out IAC?


save(list = c("abund.site.nri_aw", "abund.site.nri_pa", 
              "abund.site.nri_summ_aw", "abund.site.nri_summ_pa", 
              "site.abund.comm", "rest.abund.comm", "nat.abund.comm",
              "phylostruct.results"), 
     file = "phyloDiv_Site_NriMpd.RData")



# # Using pez ---------------------------------------------------
# library(pez)
# library(plyr)
# 
# site.vector <- plotsSites$SSU[match(rownames(site.phy.abund$comm), plotsSites$Sample_ID)]
# site.abund.comm <- ddply(site.phy.abund$comm, c(.(SSU = site.vector)), "colMeans")
# rownames(site.abund.comm) <- site.abund.comm[, "SSU"]
# site.abund.comm <- select(site.abund.comm, -SSU)
# rownames(Site_Table) <- Site_Table$SSU  
# data <- comparative.comm(phy = site.phy.abund$phy, 
#                          comm = as.matrix(site.abund.comm),
#                          env = Site_Table)
# # Dispersion
# dispersion(data, metric = "sesmpd", permute = 100)
# 
# # Null models for phylo div
# # Must be presence-absence data
# data <- data[,colSums(data$comm) > 0]
# data$comm[data$comm > 0] <- 1
# sims.phy <- ConDivSim(data, type = "phy")
