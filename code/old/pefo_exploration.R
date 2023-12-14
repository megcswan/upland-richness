##worked example with pplots

#Install packages
packages <- c("tidyverse","codyn", "vegan")
lapply(packages, library, character.only = T)

#library(gridExtra)
#library(grid)
#library(gtable)



theme_set(theme_bw(12))

# read in data ------------------------------------------------------------


dat<-read.csv("http://lter.konza.ksu.edu/sites/default/files/data/PPL011.csv")
df <- read.csv("species_filtered.csv")

#create a list of park names
parks <- unique(df$Park)

#create a list of dataframes separated by park
species_park <- lapply(1:length(parks), function(x) {df %>% filter(Park==parks[x])})
#renames objects in list park names
names(species_park)<- parks #Change the names of the items in the list to the names of the data files


##pplots
testyears<-dat%>%
  filter(RecYear==2002|RecYear==2011, Treatment=="N1P0"|Treatment=="N2P3"|Treatment=="N2P0")%>%
  mutate(GenusSpecies = paste(Genus, Species, sep = " "))%>%
  dplyr::select(-c (Genus, Species))%>%
  group_by(RecYear, PlotID)%>%
  mutate(rank=rank(-Abundance, ties.method="average"))%>%
  ungroup()%>%
  as_tibble()

names(dat)
pefo <- species_park$PEFO 

#Rank abundance

pefo_rank <- pefo %>% 
  group_by(Replicate, EventYear) %>% 
  mutate(rank = rank(-Cover_pct, ties.method = "average"))

sampleyears <- pefo_rank %>% 
  dplyr::select(Replicate, EventYear) %>% 
  unique()

trts<-testyears%>%
  dplyr::select(PlotID, Treatment)%>%
  unique()

# Make Figure 5 and 6 and 7 -----------------------------------------------------------

###NMDS in Figure 5

##step 1. do NMDS of pretreatment and last year of data
wide_dat<-testyears%>%
  dplyr::select(-rank)%>%
  dplyr::select(Treatment, RecYear, PlotID, GenusSpecies, Abundance)%>%
  spread(GenusSpecies, Abundance, fill=0)


pefo_wide<-pefo_rank%>%
  dplyr::select(-rank)%>%
  dplyr::select(EcoSite:Quadrat, Replicate, CurrentSpecies, Cover_pct)%>%
  spread(CurrentSpecies, Cover_pct, fill=0) %>% 
  mutate(year_factor = as.factor(EventYear))
glimpse(pefo_wide)
plots<-wide_dat[,1:3]
pefo_env <- pefo_wide %>% dplyr::select(EcoSite:Replicate)

mds<-metaMDS(wide_dat[,4:61], autotransform=FALSE, shrink=FALSE, trymax = 50)
mds

# names(pefo_wide)
set.seed(1234)
mds_pefo<-metaMDS(pefo_wide[,8:188], autotransform=TRUE, shrink=FALSE, trymax = 100)
mds_pefo

set.seed(1234)
mds_pefo_k3<-metaMDS(pefo_wide[,8:188], autotransform=TRUE, shrink=FALSE, trymax = 100, k = 3)
mds_pefo_k3

#try to add plot = TRUE to metaMDS

##test for differences in centroid means
adonis2(pefo_wide[,8:188]~year_factor*EcoSite, pefo_wide)

#test for differences dispersion
# dist<-vegdist(wide_dat[,4:61])
# betadisp<-betadisper(dist, wide_dat$Treatment,type="centroid")
# betadisp
# permutest(betadisp)

dist<-vegdist(pefo_wide[,8:188])
betadisp<-betadisper(dist, pefo_wide$year_factor,type="centroid")
betadisp
permutest(betadisp)


#isolate NMDS scores to great figures
scores <- data.frame(scores(mds_pefo_k3, display="sites"))  # Extracts NMDS scores for year "i" #
scores2<- cbind(pefo_env, scores) # binds the NMDS scores of year i to all years previously run
toplot<-scores2

##Make NMDS graphs for figure 5
#graph for controls
ggplot(subset(toplot), aes(x=NMDS1, y=NMDS2, color=as.factor(EcoSite)))+
  geom_point(size=1)+
  #scale_color_manual(name="", values=c("black","red"))+
  #scale_x_continuous(limits=c(-1,1.2))+
  #scale_y_continuous(limits=c(-1.1,1.3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")

# #graph of N+P
# ggplot(subset(toplot, Treatment=="N2P3"), aes(x=NMDS1, y=NMDS2, color=as.factor(RecYear)))+
#   geom_point(size=5)+
#   scale_color_manual(name="", values=c("black","red"))+
#   scale_x_continuous(limits=c(-1,1.2))+
#   scale_y_continuous(limits=c(-1.1,1.3))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
#
# #graph of N only
# ggplot(subset(toplot, Treatment=="N2P0"), aes(x=NMDS1, y=NMDS2, color=as.factor(RecYear)))+
#   geom_point(size=5)+
#   scale_color_manual(name="", values=c("black","red"))+
#   scale_x_continuous(limits=c(-1,1.2))+
#   scale_y_continuous(limits=c(-1.1,1.3))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
#

### Making Rank Abundance Curves for Figure 6

quantile(pefo$Cover_pct)

pefo %>% 
  filter(Cover_pct>50)
pefo_rank <- pefo %>% 
  group_by(Replicate, EventYear) %>% 
  mutate(rank = rank(-Cover_pct, ties.method = "average"))

#label top three species in control plots blue, top 3 in N or N+P purple and top 3 in NP red.
#top 3 control: ag, as, sn, top 3 n/np: os, sc, sm, do
#ag=blue, as=lightblue, do=darkred, os = pink, other = green3, sc = red, sn = cornflowerblue, sm = orange
pefo_c<-pefo_rank%>%
  filter(EcoSite=="C")


pefo_c_highcover <- pefo_c %>% 
  ungroup() %>% 
  filter(Cover_pct>25) %>% 
  dplyr::select(CurrentSpecies, ReportTaxon) %>% 
  distinct()

pefo_color <- pefo_c %>% 
  mutate(colorfill=ifelse(CurrentSpecies%in%pefo_c_highcover$CurrentSpecies, CurrentSpecies, "other"))

View(pefo_c %>% 
  filter(CurrentSpecies=="Sporobolus coromandelianus"))
####make RACs graphs for figure 6
#contorls
ggplot(data=subset(pefo_color, EventPanel=="A" & rank<=5), aes(x=rank, y=Cover_pct))+
  geom_line(aes(group=EventYear), linewidth=1)+
  geom_point(aes(color=colorfill), size=3)+#
  facet_wrap(~as.factor(EventYear)) +
  scale_color_manual(values = c("blue","lightblue","darkred", "pink","green3","orange","cornflowerblue", "purple","yellow","black","gray"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", axis.title.x = element_blank(), axis.title.y = element_blank())

pefo_species <- pefo_wide %>% 
  ungroup() %>% 
  dplyr::select(-c(EcoSite:Replicate)) %>% 
  as.data.frame()

pefo_env <- pefo_wide %>% 
  ungroup() %>% 
  dplyr::select(c(EcoSite:Replicate)) %>% 
  as.data.frame() %>% 
  mutate(year_factor=as.factor(EventYear))
data(dune.env)
# data(dune)
# names(pefo_c_color)
# d.specs<-  %>%
#   select(-(Park:Quadrat), -PlotTransQuad) %>% 
#   as.data.frame()

d.env  <-  d  %>%
  select((Park:Quadrat), PlotTransQuad) %>% 
  as.data.frame() %>% 
  mutate_all(.,as.factor)
glimpse(d.env)
library(BiodiversityR)

data(dune)
RankAbun.1 <-  rankabundance(pefo_species)
RankAbun.1
rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3))
rankabunplot(RankAbun.1, scale='logabun', addit=FALSE, specnames=c(1:30), 
             srt=45, ylim=c(1,10000))
rankabuncomp(dune, y=dune_env, factor=, 
             scale='proportion', legend=TRUE)
## CLICK IN THE GRAPH TO INDICATE WHERE THE LEGEND NEEDS TO BE PLACED
## IF YOU OPT FOR LEGEND=TRUE.

## Not run: 
# ggplot2 plotting method

# Only label the two most abundant species

RA.data <- rankabuncomp(d.specs, y=d.env, factor='EventPanel', 
                        return.data=TRUE, specnames=c(1:4), legend=FALSE)
names(RA.data)
dom_spec_band_p <- RA.data %>%
  filter(rank %in% c(1,2)) %>% 
  select( species) %>% 
  distinct()
library(ggplot2)
library("ggrepel")

# possibly need for extrafont::loadfonts(device="win") to have Arial
# as alternative, use library(ggThemeAssist)
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12, family="Arial"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

plotgg1 <- ggplot(data=RA.data, aes(x = rank, y = abundance)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), linewidth=1) +
  geom_point(aes(colour=Grouping, shape=Grouping), size=5, alpha=0.7) +
  geom_text_repel(data=subset(RA.data, labelit == TRUE), 
                  aes(colour=Grouping, label=species), 
                  nudge_x=1, nudge_y=1, show.legend=FALSE) +
  BioR.theme +
  scale_color_brewer(palette = "Set1") +
  labs(x = "rank", y = "abundance", colour = "Management", shape = "Management")

plotgg1

# use different facets
# now label first 10 species
RA.data <- rankabuncomp(dune, y=dune.env, factor='Management', 
                        return.data=TRUE, specnames=c(1:10), legend=FALSE)

plotgg2 <- ggplot(data=RA.data, aes(x = rank, y = abundance)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=1) +
  geom_point(aes(colour=Grouping), size=5, alpha=0.7) +
  geom_text_repel(data=subset(RA.data, labelit == TRUE), 
                  aes(label=species), 
                  angle=45, nudge_x=1, nudge_y=1, show.legend=FALSE) +
  BioR.theme +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~ Grouping) +
  labs(x = "rank", y = "abundance", colour = "Management")

plotgg2









#Code from paper
ggplot(data=subset(ractoplot, Treatment=="N1P0"&RecYear==2011), aes(x=rank, y=Abundance))+
  geom_line(aes(group=PlotID), color="darkgray", linewidth=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","darkred","pink","green3","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())

#N+P
ggplot(data=subset(ractoplot, Treatment=="N2P3"&RecYear==2002), aes(x=rank, y=Abundance))+
  geom_line(aes(group=PlotID), color="darkgray", linewidth=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","darkred", "green3","red","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
ggplot(data=subset(ractoplot, Treatment=="N2P3"&RecYear==2011), aes(x=rank, y=Abundance))+
  geom_line(aes(group=PlotID), color="darkgray", linewidth=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","darkred", "pink","green3","red","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())

#N only
ggplot(data=subset(ractoplot, Treatment=="n2p0"&RecYear==2002), aes(x=rank, y=Abundance))+
  geom_line(aes(group=PlotID), color="darkgray", size=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","darkred", "green3","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
ggplot(data=subset(ractoplot, Treatment=="n2p0"&RecYear==2011), aes(x=rank, y=Abundance))+
  geom_line(aes(group=PlotID), color="darkgray", size=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue",'darkred',"pink","green3","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())

##Make Curve Change graphs for figure 7.

###Getting relative rank and cumulative abundance to plot cumulative curves
ccplot<-ractoplot%>%
  mutate(SumAbund=sum(Abundance),
         Abund2=Abundance/SumAbund)%>%
  mutate(rank=rank(-Abund2, ties.method="average"),
         maxrank=max(rank),
         relrank=rank/maxrank)%>%
  arrange(relrank)%>%
  mutate(cumabund=cumsum(Abund2))%>%
  mutate(year=as.character(RecYear))

##graph of controls
ggplot(data=subset(ccplot, Treatment=="n1p0"), aes(x=relrank, y=cumabund, color=year))+
  geom_step(size=2)+
  scale_color_manual(values=c("black","red"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  ylab("Cumulative Relative Abundance")+
  xlab("Relative Rank")+
  scale_x_continuous(breaks=c(0.25, 0.75))+
  facet_wrap(~PlotID, ncol=3)+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

##graph of N+P
ggplot(data=subset(ccplot, Treatment=="n2p3"), aes(x=relrank, y=cumabund, color=year))+
  geom_step(size=2)+
  scale_color_manual(values=c("black","red"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  ylab("Cumulative Relative Abundance")+
  xlab("Relative Rank")+
  scale_x_continuous(breaks=c(0.25, 0.75))+
  facet_wrap(~PlotID, ncol=3)+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

##graph of N only
ggplot(data=subset(ccplot, Treatment=="n2p0"), aes(x=relrank, y=cumabund, color=year))+
  geom_step(size=2)+
  scale_color_manual(values=c("black","red"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  ylab("Cumulative Relative Abundance")+
  xlab("Relative Rank")+
  scale_x_continuous(breaks=c(0.25, 0.75))+
  facet_wrap(~PlotID, ncol=3)+
  theme(strip.background = element_blank(),strip.text.x = element_blank())


# doing RAC change and curve change anlayses for Figure 8 -------------------------------------------------------------

##get RAC changes
rac <- RAC_change(testyears, time.var = "RecYear", species.var = "GenusSpecies", abundance.var = "Abundance", replicate.var = "PlotID")

##doing curve change
cc <- curve_change(testyears, time.var = "RecYear", species.var = "GenusSpecies", abundance.var = "Abundance", replicate.var = "PlotID")

##merge for stats
allmetrics_full<-rac%>%
  left_join(cc)%>%
  left_join(trts)%>%
  gather(metric, value, richness_change:curve_change)

summary(aov(value~Treatment, data=subset(allmetrics_full, metric=="richness_change")))# p = 0.906
summary(aov(value~Treatment, data=subset(allmetrics_full, metric=="evenness_change")))  # p = 0.053
summary(aov(value~Treatment, data=subset(allmetrics_full, metric=="rank_change")))  # p < 0.001
TukeyHSD(aov(value~Treatment, data=subset(allmetrics_full, metric=="rank_change")))
summary(aov(value~Treatment, data=subset(allmetrics_full, metric=="gains"))) # p = 0.353
summary(aov(value~Treatment, data=subset(allmetrics_full, metric=="losses"))) # p = 0.174
summary(aov(value~Treatment, data=subset(allmetrics_full, metric=="curve_change"))) # p = 0.223



#merge together for figures
allmetrics<-rac%>%
  left_join(cc)%>%
  left_join(trts)%>%
  gather(metric, value, richness_change:curve_change)%>%
  group_by(Treatment, RecYear, RecYear2, metric)%>%
  summarize(vmean=mean(value),
            vn=length(PlotID),
            vsd=sd(value))%>%
  mutate(vse=vsd/sqrt(vn))

##graph all metrics for figure 6            
theme_set(theme_bw(12))
S<-ggplot(data=subset(allmetrics, metric=="richness_change"), aes(x=Treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N", "N+P"))+
  ylab("Richness Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
E<-ggplot(data=subset(allmetrics, metric=="evenness_change"), aes(x=Treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N","N+P"))+
  ylab("Evenness Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
R<-ggplot(data=subset(allmetrics, metric=="rank_change"), aes(x=Treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N", "N+P"))+
  ylab("Rank Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate('text', label="C", x=1, y = 0.175, size=5)+
  annotate('text', label="B", x=2, y = 0.24, size=5)+
  annotate('text', label="A", x=3, y = 0.28, size=5)
G<-ggplot(data=subset(allmetrics, metric=="gains"), aes(x=Treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N","N+P"))+
  ylab("Speceis Gains")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
L<-ggplot(data=subset(allmetrics, metric=="losses"), aes(x=Treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N","N+P"))+
  ylab("Species Losses")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
C<-ggplot(data=subset(allmetrics, metric=="curve_change"), aes(x=Treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N","N+P"))+
  ylab("Curve Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grid.arrange(S, E,R, G, L, C, ncol=3)

# Making Appendix figure 7 of change/difference metrics for all years -------------------------
allyears<-dat%>%
  filter(Treatment=="n1p0"|Treatment=="n2p3"|Treatment=="n2p0")%>%
  mutate(GenusSpecies = paste(Genus, Species, sep = " "))%>%
  dplyr::select(-Genus, -Species)

rac_allyears<-RAC_change(allyears, time.var = "RecYear", species.var = "GenusSpecies", abundance.var = "Abundance", replicate.var = "PlotID")
cc_allyears<- curve_change(allyears, time.var = "RecYear", species.var = "GenusSpecies", abundance.var = "Abundance", replicate.var = "PlotID")
mult_change_allyears<-multivariate_change(allyears, time.var = "RecYear", species.var = "GenusSpecies", abundance.var = "Abundance", replicate.var = "PlotID", treatment.var = "Treatment")

rac_cc_mean<-rac_allyears%>%
  left_join(cc_allyears)%>%
  left_join(trts)%>%
  gather(metric, value, richness_change:curve_change)%>%
  group_by(Treatment, RecYear, RecYear2, metric)%>%
  summarize(vmean=mean(value),
            vn=length(PlotID),
            vsd=sd(value))%>%
  mutate(vse=vsd/sqrt(vn))

theme_set(theme_bw(12))
bc<-
  ggplot(data=mult_change_allyears, aes(x=RecYear2, y=composition_change, group=Treatment, color=Treatment))+
  geom_point(size=3)+
  scale_color_manual(name="Treatment", values=c("black","red","purple"))+
  geom_line(size=1)+
  ylab("Compositional Change")+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank())
disp<-
  ggplot(data=mult_change_allyears, aes(x=RecYear2, y=dispersion_change, group=Treatment, color=Treatment))+
  geom_point(size=3)+
  scale_color_manual(name="Treatment", values=c("black","red","purple"))+
  geom_line(size=1)+
  ylab("Dispersion Change")+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank())


S<-
  ggplot(data=subset(rac_cc_mean,metric=="richness_change"), aes(x=RecYear2, y=vmean, color=Treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","red","purple"))+
  ylab("Richness Change")+
  geom_line(size=1, aes(group=Treatment))+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x=element_blank())
E<-
  ggplot(data=subset(rac_cc_mean,metric=="evenness_change"), aes(x=RecYear2, y=vmean, color=Treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","red","purple"))+
  ylab("Evenness Change")+
  geom_line(size=1, aes(group=Treatment))+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x=element_blank())
R<-
  ggplot(data=subset(rac_cc_mean,metric=="rank_change"), aes(x=RecYear2, y=vmean, color=Treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","red","purple"))+
  ylab("Rank Change")+
  geom_line(size=1, aes(group=Treatment))+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x=element_blank())
G<-
  ggplot(data=subset(rac_cc_mean,metric=="gains"), aes(x=RecYear2, y=vmean, color=Treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","red","purple"))+
  ylab("Species Gains")+
  geom_line(size=1, aes(group=Treatment))+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x=element_blank())
L<-
  ggplot(data=subset(rac_cc_mean,metric=="losses"), aes(x=RecYear2, y=vmean, color=Treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(name = "Treatment", label=c("Control", "N", "N+P"),values=c("black","red","purple"))+
  ylab("Species Losses")+
  geom_line(size=1, aes(group=Treatment))+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank())

legend=gtable_filter(ggplot_gtable(ggplot_build(L)), "guide-box") 
grid.draw(legend)


grid.arrange(arrangeGrob(bc+theme(legend.position="none"),
                         disp+theme(legend.position="none"),
                         S+theme(legend.position="none"),
                         E+theme(legend.position="none"),
                         R+theme(legend.position="none"),
                         G+theme(legend.position="none"),
                         L+theme(legend.position="none"),
                         ncol=2),legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width),nrow=1)



# Looking at differences to make Table  4-------------------------------------------------------------

rac_diff<-RAC_difference(df = testyears, time.var="RecYear", species.var = "GenusSpecies", abundance.var = "Abundance", treatment.var = "Treatment", replicate.var = "PlotID", pool = TRUE)

cc_diff<-curve_difference(df = testyears, time.var="RecYear", species.var = "GenusSpecies", abundance.var = "Abundance", treatment.var = "Treatment", replicate.var = "PlotID", pool = TRUE)



# Looking at differnces all years to make appendix figure 8 ---------------


##graph all differences for all years
rac_diff_allyears<-RAC_difference(allyears, time.var="RecYear", species.var = "GenusSpecies", abundance.var = "Abundance", treatment.var = "Treatment", replicate.var = "PlotID", pool = TRUE)%>%
  mutate(group1=paste(Treatment, Treatment2, sep="_"))

cc_diff_allyears<-curve_difference(allyears, time.var="RecYear", species.var = "GenusSpecies", abundance.var = "Abundance", treatment.var = "Treatment", replicate.var = "PlotID", pool = TRUE)%>%
  mutate(group1=paste(Treatment, Treatment2, sep="_"))

mult_diff_allyears<-multivariate_difference(allyears, time.var="RecYear", species.var = "GenusSpecies", abundance.var = "Abundance", treatment.var = "Treatment", replicate.var = "PlotID")%>%
  mutate(group1=paste(Treatment, Treatment2, sep="_"))

bc_d<-
  ggplot(data=mult_diff_allyears, aes(x=as.numeric(RecYear), y=composition_diff, color=group1, group=group1))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_manual(name = "Treatment\nComparision", label=c("Control - N", "Control - N+P", "N - N+P"),values=c("orange","darkgoldenrod4","darkred"))+
  ylab("Compositional Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
disp_d<-
  ggplot(data=mult_diff_allyears, aes(x=as.numeric(RecYear), y=abs_dispersion_diff, color=group1, group=group1))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_manual(name = "Treatment\nComparision", label=c("Control - N", "Control - N+P", "N - N+P"),values=c("orange","darkgoldenrod4","darkred"))+
  ylab("Absolute Disperion Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
s_d<-
  ggplot(data=rac_diff_allyears, aes(x=as.numeric(RecYear), y=richness_diff, color=group1, group=group1))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_manual(name = "Treatment\nComparision", label=c("Control - N", "Control - N+P", "N - N+P"),values=c("orange","darkgoldenrod4","darkred"))+
  ylab("Compositional Difference")+
  ylab("Richness Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
e_d<-
  ggplot(data=rac_diff_allyears, aes(x=as.numeric(RecYear), y=evenness_diff, color=group1, group=group1))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_manual(name = "Treatment\nComparision", label=c("Control - N", "Control - N+P", "N - N+P"),values=c("orange","darkgoldenrod4","darkred"))+
  ylab("Compositional Difference")+
  ylab("Evenness Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
r_d<-
  ggplot(data=rac_diff_allyears, aes(x=as.numeric(RecYear), y=rank_diff, color=group1, group=group1))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_manual(name = "Treatment\nComparision", label=c("Control - N", "Control - N+P", "N - N+P"),values=c("orange","darkgoldenrod4","darkred"))+
  ylab("Compositional Difference")+
  ylab("Rank Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
sp_d<-
  ggplot(data=rac_diff_allyears, aes(x=as.numeric(RecYear), y=species_diff, color=group1, group=group1))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_manual(name = "Treatment\nComparision", label=c("Control - N", "Control - N+P", "N - N+P"),values=c("orange","darkgoldenrod4","darkred"))+
  ylab("Compositional Difference")+
  ylab("Species Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())


legend=gtable_filter(ggplot_gtable(ggplot_build(sp_d)), "guide-box") 
grid.draw(legend)


grid.arrange(arrangeGrob(bc_d+theme(legend.position="none"),
                         disp_d+theme(legend.position="none"),
                         s_d+theme(legend.position="none"),
                         e_d+theme(legend.position="none"),
                         r_d+theme(legend.position="none"),
                         sp_d+theme(legend.position="none"),
                         ncol=2),legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width),nrow=1)