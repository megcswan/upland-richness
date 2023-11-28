#Questions
#*Difference between columns cover and presence absence, discuss which to use
#*unknown species
#*look at reps graph, why is is so weird?

#Install packages
packages <- c("tidyverse", "here","codyn", "BiodiversityR")
lapply(packages, library, character.only = T)
# #Load data
# 
allspecies = read_csv("allspecies.csv") %>% select(-"...1")

#organizing columns so I understand nested design and to make it easier to work with each level of the design
# TransectQuadratID=Park_EcoSitePlot_Date_Transect_Quadrat
allspecies <- allspecies %>%
  mutate(TransectQuadrat = str_remove(TransectQuadratID, pattern ="\\w{4}_\\w{1}\\d{2}_\\d{8}_")) %>%
  separate(TransectQuadrat, into = c("Transect","Quadrat"), sep = "_")#separate transect and quadrat into 2 columns


# Make correction to data frame -------------------------------------------


#need to add together 2 instances of Carex spp.
#See commented code below for the instances
#long dataframe with presence absence and cover corrected
species_corrected <- allspecies %>% 
  group_by(Park, EcoSite, EventYear, EventPanel, Plot, Transect, Quadrat, CurrentSpecies) %>% 
  mutate(PresentAbsent = sum(SpeciesPresentInQuadratForNested),
         Cover_pct = sum(CoverClassMidpoint_Quadrat_pct)) %>% 
  select(Park:EventPanel,Plot, Transect, Quadrat, ReportTaxon, Nativity, Lifeform, Duration,CurrentSpecies, PresentAbsent, Cover_pct) %>% 
  distinct() %>% #removed 2 taxa from record, add back later? 
  mutate(PresentAbsent = ifelse(PresentAbsent==2, 1, PresentAbsent))#,

#check again on what the columns mean, should these values be na not 0 for pct
View(allspecies %>% 
  filter(SpeciesPresentInQuadratForCover!=1 & SpeciesPresentInQuadratForNested==1))

#To see duplicates that were combined run commented script below
# View(allspecies %>%
#        group_by(Park, EcoSite, EventYear, EventPanel, Plot, Transect, Quadrat, CurrentSpecies) %>%
#        mutate(PresentAbsent = sum(SpeciesPresentInQuadratForNested),
#               Cover_pct = sum(CoverClassMidpoint_Quadrat_pct)) %>%
#   filter(CoverClassMidpoint_Quadrat_pct!=Cover_pct))


# #Pivot wide to calculate richness  --------------------------------------


#Values are Present/Absent
species_wide <- species_corrected %>% 
  select(-c(Cover_pct, ReportTaxon, Nativity, Lifeform, Duration)) %>% 
  pivot_wider(names_from = CurrentSpecies, values_from = PresentAbsent, values_fill = 0) %>% 
  ungroup()

#make dataframe that has only the species as columns
species_only=species_wide %>% 
  select(-c(Park:Quadrat)) 

#same table as above but as binary
species_binary <- species_only%>% 
  mutate_all(., function(x) as.integer(ifelse(x>0, 1, 0)))

#table with event details
species_env = species_wide %>% 
  select(c(Park:Quadrat)) %>% 
  mutate(rownumber = row_number())

#Calculate species richness
rich_calc = estimateR(species_binary) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(cols = starts_with("V"), names_to = "rownumber" ) %>% 
  pivot_wider(names_from = rowname) %>% 
  mutate(rownumber = as.integer(str_remove(rownumber, "V")))   

rich_env=species_env %>% 
  left_join(rich_calc)





# #Pivot wide to calculate diversity metrics ------------------------------


species_wide_div <- species_corrected %>% 
  select(-c(PresentAbsent, ReportTaxon, Nativity, Lifeform, Duration)) %>% 
  pivot_wider(names_from = CurrentSpecies, values_from = Cover_pct, values_fill = 0) %>% 
  ungroup()
#make dataframe that has only the species as columns
species_only_div=species_wide_div %>% 
  select(-c(Park:Quadrat)) 
#same table as above but as binary
species_binary_div <- species_only_div%>% 
  mutate_all(., function(x) as.integer(ifelse(x>0, 1, 0)))

rich_calc_div = estimateR(species_binary_div) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(cols = starts_with("V"), names_to = "rownumber" ) %>% 
  pivot_wider(names_from = rowname) %>% 
  mutate(rownumber = as.integer(str_remove(rownumber, "V"))) %>% 
  rename_all(paste0, "_cover") %>% 
  rename(rownumber=rownumber_cover)


#Calculate species diversity
div_calc = species_only %>% 
  mutate(simpson = diversity(., index = "simpson"),
         shannon = diversity(., index = "shannon"),
         invD = diversity(., index = "invsimp"),
         rownumber = row_number()) %>% 
  select(simpson:rownumber)



# Make dataframe with vegan calculated metrics ----------------------------
names(rich_calc_div)

community=rich_env %>% 
  left_join(div_calc) %>% 
  left_join(rich_calc_div)

View(community %>% 
  filter(S.obs!=S.obs_cover))


# #Files  written to .csv---------------------------------------------------

# 
# write.csv(species_corrected, file = "species_corrected.csv")
# write.csv(community, file = "community.csv")
# created a file that has metrics on it



#*Need to:
#*Edit this document have a master file for data wrangling
#*Calculate codyn measurements
#*Try running some models

# #Add diversity metrics from codyn ---------------------------------------
reps = species_corrected %>%
  ungroup() %>% 
  select(EventYear,Plot, EventPanel) %>% 
  distinct() %>% 
  separate(Plot, sep="_", into =c("Park","X")) %>% 
  mutate(EcoSite = str_extract(X, "\\w{1}"),
         Plot = str_extract(X, "\\d{2}")) %>% 
  select(-X) %>% 
  filter(Park == "BAND"&
           !EventPanel %in% c("S","Z"))

names(reps)


reps %>% 
  # group_by(EventYear, EventPanel, EcoSite) %>% 
  # tally() %>% 
  ggplot(aes(x = EventYear, fill = EcoSite)) +
  geom_bar(stat='bin', position = 'dodge') +
  facet_wrap(~EventPanel)



# #calculate metrics from codyn -------------------------------------------

#stopped here
#need to remove unknown species

names(species_corrected)


species_corrected<- species_corrected %>% 
  mutate(PlotTransQuad = paste(Plot, Transect, Quadrat, sep = "_"))

species_plot <- species_corrected %>% 
  group_by(Plot, EventYear, EventPanel, CurrentSpecies) %>% 
  summarise(PA =sum(PresentAbsent),
            Cover_avg = mean(Cover_pct, na.rm=TRUE))
  

adund_change <- abundance_change(df =
                                   species_plot,
                                 species.var = "CurrentSpecies", 
                                 abundance.var = "Cover_avg", 
                                 replicate.var = "Plot", 
                                 time.var = "EventYear")%>% 
  separate(Plot, into = c("Park", "Plot"), sep = c("_"))

comm_stab <- community_stability(df =
                                   species_plot,
                                 abundance.var = "Cover_avg", 
                                 replicate.var = "Plot", 
                                 time.var = "EventYear") %>% 
  separate(Plot, into = c("Park", "Plot"), sep = c("_"))

comm_stab_plot <- comm_stab %>% 
  ggplot(aes(x = Park, y = stability)) +
  geom_boxplot()

abund_stab_plot <- abund_stab %>% 
  ggplot(aes(x = Park, y = stability)) +
  geom_boxplot()
names(d)
data(dune.env)
data(dune)
d.specs<- d %>%
  select(-(Park:Quadrat), -PlotTransQuad) %>% 
  as.data.frame()

d.env  <-  d  %>%
  select((Park:Quadrat), PlotTransQuad) %>% 
  as.data.frame() %>% 
  mutate_all(.,as.factor)
glimpse(d.env)
RankAbun.1 <-  rankabundance(d.specs)
RankAbun.1
rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3))
rankabunplot(RankAbun.1, scale='logabun', addit=FALSE, specnames=c(1:30), 
             srt=45, ylim=c(1,100))
rankabuncomp(dune, y=dune.env, factor='Management', 
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
#stopped here 
#need to figure out how to get the top species
#maybe do a rank abundance curve for each transect.
c_dom <- c %>% 
  group_by(Plot, EventYear, CurrentSpecies) %>% 
  summarize(Cover_pct = mean(Cover_pct)) %>% 
  filter(Cover_pct>0)%>% 
  top_n(Cover_pct, 4)

c_dom %>% 
  ungroup() %>% 
  select(CurrentSpecies) %>% 
  distinct()

names(species_dup_div)
a <- species_dup_div %>% 
  filter(is.na(Cover_pct))

#need to figure out why there are NAs, probably has something do do with nested and other thing

b <- allspecies %>% 
  filter(Plot%in% a$Plot & 
           EventYear %in% a$EventYear &
           Transect %in% a$Transect &
           Quadrat %in% a$Quadrat &
           CurrentSpecies %in% a$CurrentSpecies)

data(collins08)






species_dup
data(pplots) 

# Without replicates 
df <- subset(pplots, plot == 25) 
abundance_change(df = df, species.var = "species", abundance.var = "relative_cover", time.var = "year")

#ave quad richness by plot
avequadrich<-quadrichness%>%
  separate(TransectQuadratID, into=c(NA, "Plot", NA), c(5, 8))%>%
  group_by(EcoSite, EventYear, Plot)%>%
  summarize(avequadrich = mean(quadrich))
ggplot(avequadrich)+
  geom_jitter(aes(x=EventYear, y = avequadrich, color = EcoSite), alpha = .9, size = 5, shape = 1, stroke = 2)

ggplot(plotrichness)+
  geom_jitter(aes(x=EventYear, y = plotrich, color = EcoSite), alpha = .9, size = 5, shape = 20)
#ave plot richness for ecosite
ecorich<-plotrichness%>%
  group_by(EcoSite, EventYear)%>%
  summarize(ecorich = mean(plotrich))

ggplot(ecorich)+
  geom_jitter(aes(x=EventYear, y = ecorich, color = EcoSite), alpha = .9, size = 5, shape = 20)
