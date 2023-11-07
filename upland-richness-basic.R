#get data from database
#library(RODBC)
library(tidyverse)
library(ggplot2)
#connect to upland database
#conUpE <- odbcDriverConnect("driver=SQL Server;server=INPSCPNMS01\\Production, 50000; database=SCPN_UplandEvent;trusted_connection=yes;")
#pull in info for grassland parks
#spp_all <- sqlQuery(conUpE, "SELECT * FROM vtbl_NestedSpecies_All", as.is=FALSE, stringsAsFactors = FALSE)
#close connection
#close(conUpE)
#spp_all<-read.csv(spp_occur, "C:/Users/mswan/Documents/R/Projects/upland-veg-richness/allspecies.csv")
#filter dataframe to get only VALID or Correct and INPANEL and numeric and reduce to COlumns desired for CSP
spp_all<-spp_all%>%
  filter(UsableQuadratForNested %in% "1", !NestedSpeciesValidation %in% "Issue")
#do we have all the plots we expect?
plots<-spp_CSP_all%>%
  distinct(Plot)
#keep only species that occur somewhere in EcoSite
spp_occur<-spp_all%>%
  group_by(EcoSite, CurrentSpecies)%>%
  mutate(InEcosite = if_else(sum(SpeciesPresentInQuadratForNested)>0, 1, 0))%>%
  ungroup()%>%
  filter(InEcosite >0)
#took out filter for panel until we decide whether or not to exclude  not in panel plots , EventPanel %in% c("A", "B", "C")
#write this file to use in future to avoid long download times.
#write.csv(spp_occur, "C:/Users/mswan/Documents/R/Projects/upland-veg-richness/allspecies.csv")



#richness by quadrat
quadrichness<-spp_occur%>%
  group_by(EcoSite, EventYear, TransectQuadratID)%>%
  summarize(quadrich = sum(SpeciesPresentInQuadratForNested))
#ave quad richness by plot
avequadrich<-quadrichness%>%
  separate(TransectQuadratID, into=c(NA, "Plot", NA), c(5, 8))%>%
  group_by(EcoSite, EventYear, Plot)%>%
  summarize(avequadrich = mean(quadrich))
ggplot(avequadrich)+
  geom_jitter(aes(x=EventYear, y = avequadrich, color = EcoSite), alpha = .9, size = 5, shape = 1, stroke = 2)
#richness by plot
plotrichness<-spp_occur%>%
  group_by(EcoSite, EventYear, Plot)%>%
  summarize(plotrich = sum(SpeciesPresentInQuadratForNested))

ggplot(plotrichness)+
  geom_jitter(aes(x=EventYear, y = plotrich, color = EcoSite), alpha = .9, size = 5, shape = 20)
#ave plot richness for ecosite
ecorich<-plotrichness%>%
  group_by(EcoSite, EventYear)%>%
  summarize(ecorich = mean(plotrich))

ggplot(ecorich)+
  geom_jitter(aes(x=EventYear, y = ecorich, color = EcoSite), alpha = .9, size = 5, shape = 20)
