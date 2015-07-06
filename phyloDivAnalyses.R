library(plyr)
library(picante)
load("communityData.RData")


# MPD and NRI calculations ------------------------------------------------

# Calculating NRI
nri.site.fun <- function(comm, phy, site) {
  site.vector <- site$SSU[match(rownames(comm), site$Sample_ID)]
  x <- ddply(comm, c(.(SSU = site.vector)), "colMeans")
  y <- ses.mpd(x[,-1], cophenetic(phy), null.model = "taxa.labels", runs = 100,
               abundance.weighted = T)
  y$mpd.obs.z <- y$mpd.obs.z * -1
  y <- mutate(y, SSU = x$SSU)
  y <- mutate(y, Type = site$Type[match(y$SSU, site$SSU)])
  return(y)
}

s.abund.site.nri <- nri.site.fun(s.phy.abund$comm, s.phy.abund$phy, plotsSites)
m.abund.site.nri <- nri.site.fun(m.phy.abund$comm, m.phy.abund$phy, plotsSites)
l.abund.site.nri <- nri.site.fun(l.phy.abund$comm, l.phy.abund$phy, plotsSites)

# Summarizing NRI
unloadNamespace("reshape2")
unloadNamespace("plyr")
library(dplyr)

nri.summ.fun <- function(x) {
  y <- summarise(group_by(x, Type), mean(ntaxa), 
                 mean(mpd.obs), sd(mpd.obs), 
                 mean(mpd.obs.z), sd(mpd.obs.z))
}

# Phylogenetically: Natural wetlands overdispersed, Restored clustered
# Richness: Differences more pronounced at finer scales
(s.abund.site.nri_summ <- nri.summ.fun(s.abund.site.nri))
(m.abund.site.nri_summ <- nri.summ.fun(m.abund.site.nri))
(l.abund.site.nri_summ <- nri.summ.fun(l.abund.site.nri))