library(tidyverse)
library(vegan)
library(gridExtra)
library(grid)
library(gtable)
library(gtools)
# library(devtools)
library(codyn)

theme_set(theme_bw(12))

#For this one, questions:
#*EventPanels other than ABC
#*Events where we only went once

# Read in Data ------------------------------------------------------------


#Install packages
packages <- c("tidyverse", "here","lme4", "vegan","LMERConvenienceFunctions", "lmerTest", "emmeans", "multcomp")
lapply(packages, library, character.only = T)

#Load files
species_corrected <- read.csv("species_corrected.csv", header = T) %>% 
  mutate(EcoSite = str_remove(EcoSite, "\\w{4}_"), #Separate nested parameters
         Plot = as.factor(str_extract(Plot, "\\d{2}")),
         Quadrat = as.factor(Quadrat),
         Replicate = paste(Park, EcoSite,EventPanel, Plot, Transect, Quadrat, sep = "_"),
         Year = as.numeric(EventYear)) %>% 
  filter(EventPanel %in% c("A","B","C")) #Select only event panel a, b, and c



species_corrected_nozero <- species_corrected %>% #remove all the extra zeros to make it easier to work with
 filter(!is.na(ReportTaxon)& #remove unknown species
 !is.na(Cover_pct)&  #remove na in cover
 Cover_pct!=0)   #remove extra zeros
  


sampleyears = species_corrected_nozero %>% 
  dplyr::select(Replicate, Park, EcoSite, Plot, Transect, Quadrat, Year) %>% 
  distinct() %>% 
  group_by(Replicate) %>% 
  #mutate(mean_cover = mean(Cover_pct)) %>% #no values changed!
  tally() %>%   #add tally for the number of times a replicate was sampled
  filter(n>1)

species_corrected_nozero_sampleyears <- species_corrected_nozero %>% 
  filter(Replicate%in% sampleyears$Replicate)


# species_corrected_zero <- species_corrected %>% #remove all the extra zeros to make it easier to work with
#   filter(is.na(ReportTaxon)|
#            is.na(Cover_pct)|
#            Cover_pct==0)
#   
#verified that the rows add up properly


#create a list of park names
parks <- unique(species_corrected_nozero_sampleyears$Park)

#create a list of dataframes separated by park
species_park <- lapply(1:length(parks), function(x) {species_corrected_nozero_sampleyears %>% filter(Park==parks[x])})
#renames objects in list park names
names(species_park)<- parks #Change the names of the items in the list to the names of the data files



# Species turnover --------------------------------------------------------


###looking at temporal variability
species_turnover<-data.frame()




for (i in 1:length(parks)){
  
  com_rep <- unique(species_park[[parks[i]]]$Replicate)
  
  subset<-species_park[[parks[i]]]%>%
    filter(Replicate%in%com_rep)
  
  out <- turnover(df = subset, time.var = "Year", species.var = "CurrentSpecies", abundance.var = "Cover_pct", replicate.var = "Replicate")
  
  #out$id<-com_rep[]
  
  species_turnover<-rbind(species_turnover, out)  
}


View(species_turnover %>% 
  group_by(Replicate, Year) %>% 
  add_tally())

species_turnover$Replicate
species_turnvoer_avg<-species_turnover%>% 
  group_by(Replicate, Year)%>%
  summarise(total=mean(total, .groups = c(Replicate, Year)))%>%
  separate(Replicate, into=c("Park", "EcoSite", "EventPanel", "Plot", "Transect", "Quadrat"), sep="_", remove=F)%>%
  mutate(id3=paste(Park, EcoSite, Plot, sep="_"))%>%
  group_by(id3, Year)%>%
  summarize(turnover=mean(total))%>%
  ungroup()%>%
  separate(id3, into=c("Park", "EcoSite", "Plot"), sep="_")%>%
  group_by(Park, EcoSite, Plot)%>%
  summarize(turnover=mean(turnover))


names(species_turnvoer_avg)


species_turnover_plot <- species_turnvoer_avg
# all_turnover<-data.frame()
# com_rep <- unique(species_corrected$Replicate)


# for (i in 1:length(com_rep)){
#   
#   subset<-species_corrected%>%
#     filter(Replicate==com_rep[i])
#   
#   out <- turnover(df = subset, time.var = "Year", species.var = "CurrentSpecies", abundance.var = "Cover_pct", replicate.var = "Replicate")
#   
#   out$id<-com_rep[i]
#   
#   all_turnover<-rbind(all_turnover, out)  
# }
# 
# 
# all_turnvoer_ave<-all_turnover%>% 
#   group_by(id, Year)%>%
#   summarise(total=mean(total))%>%
#   separate(id, into=c("Park", "EcoSite", "Plot", "Transect", "Quadrat"), sep="_", remove=F)%>%
#   mutate(id3=paste(Park, EcoSite, Plot, sep="_"))%>%
#   group_by(id3, Year)%>%
#   summarize(turnover=mean(total))%>%
#   ungroup()%>%
#   separate(id3, into=c("Park", "EcoSite", "Plot"), sep="_")%>%
#   group_by(EcoSite)%>%
#   summarize(turnover=mean(turnover))
# 
# ### looking at spatial variability using Whittacker's beta diversity (gamma / average alpha)
# gammadiv<-BAND%>%
#   filter(abundance!=0)%>%
#   group_by(id3, time, rep, species)%>%
#   summarize(splist=mean(abundance))%>%
#   ungroup()%>%
#   group_by(id3, time, rep)%>%
#   summarize(gamma=length(species))
# 
# ##richness
# S<-function(x){
#   x1<-x[x!=0]
#   length(x1)
# }
# 
# rich <- group_by(BAND, id3, time, rep, site) %>% 
#   summarize(S=S(abundance))%>%
#   ungroup()%>%
#   group_by(id3, time, rep)%>%
#   summarize(alpha=mean(S))
# 
# beta_div<-gammadiv%>%
#   left_join(rich)%>%
#   mutate(wbeta=gamma/alpha)%>%
#   group_by(id3, time)%>%
#   summarize(wbeta=mean(wbeta))%>%
#   ungroup()%>%
#   separate(id3, into=c("alpha","even","comtype"), sep="_")%>%
#   group_by(comtype)%>%
#   summarize(betadiv=mean(wbeta))



# # Richness Evenness Metrics -----------------------------------------------
# com_rep<-unique(BAND$id)
# 
# BAND_div_evar<-data.frame()
# for (i in 1:length(com_rep)){
#   
#   subset<-BAND%>%
#     filter(id==com_rep[i])
#   
#   out <- community_structure(df = subset, time.var = "time", abundance.var = "abundance", replicate.var = "site", metric = "Evar")
#   out$id<-com_rep[i]
#   
#   BAND_div_evar<-rbind(BAND_div_evar, out)  
# }
# 
# BAND_diversity_mean<-BAND_div_evar%>%
#   separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
#   mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
#   group_by(id3, time, site)%>%#average over replicates
#   summarize(Sp=mean(richness),
#             Evar=mean(Evar, na.rm=T))%>%
#   ungroup()%>%
#   group_by(id3, time)%>%#average over sites
#   summarize(Sp=mean(Sp),
#             Evar=mean(Evar, na.rm=T))
# 
# 
# ###making figure for appendix 5
# BAND_div_tograph<-BAND_diversity_mean%>%
#   separate(id3, into=c("alpha","theta","scenario"), sep="_", remove=F)
# 
# ggplot(data = BAND_div_tograph, aes(x = Sp, y = Evar, color = scenario))+
#   geom_point(size = 2)+
#   scale_color_brewer(type = "div", name = "Community Type", labels = c("High Spatial, High Temporal", "Low Spatial, Low Temporal", "High Spatial, Low Temporal", "Low Spatial, High Temporal"))+
#   xlab("BANDulated Species Richness")+
#   ylab("BANDulated Evenness (Evar)")+
#   theme(panel.grid = element_blank())



# Change Metrics -------------------------------------------------------------
#need to remove samples that were collected less than 3 times
#tried to make it less than 2 but I get a similar error
lowsampleyear <-  species_park$BAND %>% 
  group_by(Replicate, Year) %>% 
  summarise(mean = mean(Cover_pct)) %>% 
  add_tally() %>% 
  filter(n<3) %>% 
  dplyr::select(Replicate) %>% 
  distinct()

BAND <- species_park$BAND %>% 
  filter(!(Replicate%in%lowsampleyear$Replicate))

dim(species_park$BAND)


### RAC Change
BAND_rac_change<-data.frame()

com_rep<-unique(BAND$Replicate)

for (i in 1:length(com_rep)){
  
  subset<-BAND%>%
    filter(Replicate==com_rep[i])
  
  out <- RAC_change(df = subset, time.var = "Year", species.var = "CurrentSpecies", abundance.var = "Cover_pct", replicate.var = "Replicate")
  
  out$id<-com_rep[i]
  
  BAND_rac_change<-rbind(BAND_rac_change, out)  
}



printList <- function(list) {
  
  for (item in 1:length(list)) {
    
    print(head(list[[item]]))
    
  }
}
printList(com_rep)






BAND_rac_change_mean<-BAND_rac_change%>%
  separate(id, into=c("Park", "EcoSite", "EventPanel","Plot", "Transect", "Quadrat"), sep="_", remove=F)%>%
  mutate(id3=paste("Park", "EcoSite", "Plot", sep="_"))%>%
  group_by(id3, Year, Year2, Replicate)%>%
  summarize(S=mean(richness_change),
            E=mean(evenness_change,na.rm=T),
            R=mean(rank_change),
            G=mean(gains),
            L=mean(losses))%>%
  ungroup()%>%
  group_by(id3, Year, Year2)%>%
  summarize(S=mean(S),
            E=mean(E,na.rm=T),
            R=mean(R),
            G=mean(G),
            L=mean(L))

##Multivariate Change
BAND_mult_change<-data.frame()

com_rep<-unique(BAND$id)

for (i in 1:length(com_rep)){
  
  subset<-BAND%>%
    filter(id==com_rep[i])
  
  out <- multivariate_change(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site")
  
  out$id<-com_rep[i]
  
  BAND_mult_change<-rbind(BAND_mult_change, out)  
}

BAND_multchange_mean<-BAND_mult_change%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time, time2)%>%
  summarize(composition_change=mean(composition_change, na.rm = T),
            dispersion_change=mean(dispersion_change, na.rm =T))

## Curve change

BAND_curve_change<-data.frame()

com_rep<-unique(BAND$id)

for (i in 1:length(com_rep)){
  
  subset<-BAND%>%
    filter(id==com_rep[i])
  
  out <- curve_change(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site")
  
  out$id<-com_rep[i]
  
  BAND_curve_change<-rbind(BAND_curve_change, out)  
}

BAND_cc_ave<-BAND_curve_change%>% 
  group_by(id, time, time2)%>%
  summarise(curve_change=mean(curve_change))%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time, time2)%>%
  summarize(curve_change=mean(curve_change))


# Difference Metrics ------------------------------------------------------

BAND_diff<-BAND%>%
  mutate(treatment = as.factor(ifelse(as.integer(site) < 5, "T1", "T2")))

#RAC diff
BAND_rac_diff<-data.frame()
com_rep<-unique(BAND_diff$id)

for (i in 1:length(com_rep)){
  
  subset<-BAND_diff%>%
    filter(id==com_rep[i])
  
  out <- RAC_difference(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site", pool = TRUE, treatment.var = "treatment")
  
  out$id<-com_rep[i]
  
  BAND_rac_diff<-rbind(BAND_rac_diff, out)  
}

BAND_rac_diff_mean<-BAND_rac_diff%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(S=mean(richness_diff),
            E=mean(evenness_diff,na.rm=T),
            R=mean(rank_diff),
            D=mean(species_diff))

##multivariate_diff
BAND_mult_diff<-data.frame()
com_rep<-unique(BAND_diff$id)

for (i in 1:length(com_rep)){
  
  subset<-BAND_diff%>%
    filter(id==com_rep[i])
  
  out <- multivariate_difference(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site", treatment.var = "treatment")
  
  out$id<-com_rep[i]
  BAND_mult_diff<-rbind(BAND_mult_diff, out)  
}

BAND_mult_diff_mean<-BAND_mult_diff%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(composition_diff=mean(composition_diff, na.rm = T),
            abs_dispersion_diff=mean(abs_dispersion_diff, na.rm=T))

###curve diff
BAND_curve_diff<-data.frame()
com_rep<-unique(BAND_diff$id)

for (i in 1:length(com_rep)){
  subset<-BAND_diff%>%
    filter(id==com_rep[i])
  
  out <- curve_difference(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site", pool = TRUE, treatment.var = "treatment")
  
  out$id<-com_rep[i]
  BAND_curve_diff<-rbind(BAND_curve_diff, out)  
}

BAND_cc_diff<-BAND_curve_diff%>% 
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(curve_diff=mean(curve_diff))

BAND_diff_div<-BAND_div_evar%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time, site)%>%#average over replicates
  summarize(Sp=mean(richness),
            Evar=mean(Evar, na.rm=T))%>%
  ungroup()%>%
  group_by(id3, time)%>%#average over sites
  summarize(Sp=mean(Sp),
            Evar=mean(Evar, na.rm=T))

# Merging all metrics to single datasets ----------------------------------

BAND_allmetrics<-BAND_multchange_mean%>%
  left_join(BAND_rac_change_mean)%>%
  left_join(BAND_cc_ave)%>%
  left_join(BAND_diversity_mean)%>%
  separate(id3, into=c("alpha","even","comtype"), sep="_")%>%
  mutate(comtype2 = as.factor(comtype))


BAND_diff_allmetrics<-BAND_rac_diff_mean%>%
  left_join(BAND_cc_diff)%>%
  left_join(BAND_mult_diff_mean)%>%
  left_join(BAND_diff_div)%>%
  separate(id3, into=c("alpha","even","comtype"), sep="_")%>%
  mutate(comtype2 = as.factor(comtype))

# Effect of rich and even on CHANGE and DIFFERNCE BAND data Table 2 ------------------------------

#How are CHANGE metrics affected by the richness and evenness of the community?
#Richness
#rich with delta rank
with(BAND_allmetrics,  cor.test(Sp, R))

##removing extra division by Sp to demonstate to why it is necessary to divide again by the size of the species pool.
BAND_allmetrics2<-BAND_allmetrics%>%
  mutate(MRS = R*Sp)
with(BAND_allmetrics2, cor.test(Sp, MRS))

ggplot(data=BAND_allmetrics2, aes(x=Sp, y=MRS, color = Evar))+
  geom_point()+
  xlab("BANDulated Community Richness")+
  ylab("Rank Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol =4)

#rich with gains or lossess same thing
with(BAND_allmetrics, cor.test(Sp, G))

#rich with composition
with(BAND_allmetrics, cor.test(Sp, composition_change))

#rich with dispersion
with(BAND_allmetrics, cor.test(Sp, dispersion_change))

#rich with curve
with(BAND_allmetrics, cor.test(Sp, curve_change))

#even with delta rank
with(BAND_allmetrics, cor.test(Evar, R))

#even with gains/losses
with(BAND_allmetrics, cor.test(Evar, G))

#even with composition
with(BAND_allmetrics, cor.test(Evar, composition_change))

#even with dispersion
with(BAND_allmetrics, cor.test(Evar, dispersion_change))

#even with curve
with(BAND_allmetrics, cor.test(Evar, curve_change))

###doing the correlations for the difference metrics

#rich with delta rank
with(BAND_diff_allmetrics, cor.test(Sp, R))

#rich with species_diff
with(BAND_diff_allmetrics,cor.test(Sp, D))

#rich with composition
with(BAND_diff_allmetrics, cor.test(Sp, composition_diff))

#rich with dispersion
with(BAND_diff_allmetrics, cor.test(Sp, abs_dispersion_diff))

#rich with curve
with(BAND_diff_allmetrics, cor.test(Sp, curve_diff))

#even with delta rank
with(BAND_diff_allmetrics, cor.test(Evar, R))

#even with species difference
with(BAND_diff_allmetrics, cor.test(Evar, D))

#even with composition
with(BAND_diff_allmetrics, cor.test(Evar, composition_diff))

#even with dispersion
with(BAND_diff_allmetrics, cor.test(Evar, abs_dispersion_diff))

#even with curve
with(BAND_diff_allmetrics, cor.test(Evar, curve_diff))


# Looking at figures of rich/even on metrics Appendix 6 figures and Table -----------


#Change metrics
BAND_tograph<-BAND_allmetrics%>%
  gather(metric, value, composition_change:curve_change)%>%
  filter(metric != "S")%>%
  filter(metric != "E")%>%
  filter(metric != "L")


labels_change <-c(composition_change = "Composition Change",
                  curve_change = "Curve Change",
                  dispersion_change = "Dispersion Change",
                  G = "Gains/Losses",
                  R = "Rank Change")

#change with richness
ggplot(data=BAND_tograph, aes(x=as.factor(Sp), y=value, color = comtype))+
  geom_boxplot()+
  scale_color_manual(name = "Community Type", labels = c("High Spatial, High Temporal","Low Spatial, Low Temporal","High Spatial, Low Temporal","Low Spatial, High Temporal"), values= c("gold","orange","darkorange4","red"))+
  xlab("BANDulated Species Richness")+
  ylab("Measure Value")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~metric, ncol = 3, scales = "free", labeller = labeller(metric = labels_change))+
  theme(strip.background = element_rect(fill = 0))


#change with evenness
ggplot(data=BAND_tograph, aes(x=as.factor(even), y=value, color = comtype))+
  geom_boxplot()+
  scale_color_manual(name = "Community Type", labels = c("High Spatial, High Temporal","Low Spatial, Low Temporal","High Spatial, Low Temporal","Low Spatial, High Temporal"), values= c("gold","orange","darkorange4","red"))+
  xlab("BANDulated  Evenness")+
  ylab("Measure Value")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~metric, ncol = 3, scales = "free", labeller = labeller(metric = labels_change))+
  theme(strip.background = element_rect(fill = 0))+
  scale_x_discrete(labels = c("Low","Mid","High"))

BAND_means<-BAND_allmetrics%>%
  mutate(abs_dc = abs(dispersion_change))%>%
  gather(metric, value, c(composition_change:curve_change, abs_dc))%>%
  group_by(comtype, metric)%>%
  summarize(mean = mean (value),
            se = sd(value)/sqrt(81))



###Difference Metrics

BAND_diff_tograph<-BAND_diff_allmetrics%>%
  gather(metric, value, S:abs_dispersion_diff)%>%
  filter(metric != "S")%>%
  filter(metric != "E")

labels_diff <-c(composition_diff = "Composition Difference",
                curve_diff = "Curve Difference",
                abs_dispersion_diff = "Dispersion Difference",
                D = "Species Differences",
                R = "Rank Difference")
#diff with richness
ggplot(data=BAND_diff_tograph, aes(x=as.factor(Sp), y=value, color = comtype))+
  geom_boxplot()+
  scale_color_manual(name = "Community Type", labels = c("High Spatial, High Temporal","Low Spatial, Low Temporal","High Spatial, Low Temporal","Low Spatial, High Temporal"), values= c("gold","orange","darkorange4","red"))+
  xlab("BANDulated Richness")+
  ylab("Measure Value")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~metric, ncol = 3, scales = "free", labeller = labeller(metric = labels_diff))+
  theme(strip.background = element_rect(fill = 0))


#diff with evenness
ggplot(data=BAND_diff_tograph, aes(x=as.factor(even), y=value, color = comtype))+
  geom_boxplot()+
  scale_color_manual(name = "Community Type", labels = c("High Spatial, High Temporal","Low Spatial, Low Temporal","High Spatial, Low Temporal","Low Spatial, High Temporal"), values= c("gold","orange","darkorange4","red"))+
  xlab("BANDulated Evenness")+
  ylab("Measure Value")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~metric, ncol = 3, scales = "free", labeller = labeller(metric = labels_diff))+
  theme(strip.background = element_rect(fill = 0))+
  scale_x_discrete(labels = c("Low","Mid","High"))

diff_means<-BAND_diff_tograph%>%
  group_by(comtype, metric)%>%
  summarize(mean = mean (value),
            se = sd(value)/sqrt(90))
