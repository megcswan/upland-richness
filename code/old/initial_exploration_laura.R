#Install packages
packages <- c("tidyverse", "here", "vegan")
lapply(packages, library, character.only = T)

#Load data

allspecies = read_csv("allspecies.csv") %>% select(-"...1")

#organizing columns so I understand nested design
# TransectQuadratID=Park_EcoSitePlot_Date_Transect_Quadrat

allspecies <- allspecies %>% 
  mutate(Quadrat = str_remove(TransectQuadratID, pattern ="_\\d{8}"))



#
plotrichness<-allspecies%>%
  group_by(EcoSite, EventYear, Plot)%>%
  summarize(plotrich = sum(SpeciesPresentInQuadratForNested))

# Diversity by plot -------------------------------------------------------


species_wide_plot = allspecies %>% 
  group_by(Park, EcoSite, EventYear, Plot, CurrentSpecies) %>% 
  mutate(CoverClass_avg = mean(NestedQuadratSizeClass)) %>% 
  select(Park:Plot,CurrentSpecies, CoverClass_avg) %>% 
  distinct() %>% #removed 2 taxa? from record, this is probably not correct
  pivot_wider(names_from = CurrentSpecies, values_from = CoverClass_avg, values_fill = 0) %>% 
  ungroup()

#For duplicates from above?
#re-run species_wide stopping before the select
#run below
# View(species_wide %>% 
#   filter(CoverClassMidpoint_Quadrat_pct!=CoverMidpoint_pct))

species_wide_int_plot=species_wide_plot %>% 
  select(-c(Park:Plot)) %>% 
  mutate_all(., function(x) as.integer(ifelse(x>0, 1, 0)))

species_env_plot = species_wide_plot %>% 
  select(c(Park:Plot)) %>% 
  mutate(rownumber = row_number())

 
  

rich_calc_plot = estimateR(species_wide_int_plot) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(cols = starts_with("V"), names_to = "rownumber" ) %>% 
  pivot_wider(names_from = rowname) %>% 
  mutate(rownumber = as.integer(str_remove(rownumber, "V")))   

rich_plot=species_env_plot %>% 
  left_join(rich_calc_plot)


# Diversity  by quadrat ---------------------------------------------------


quadrichness<-allspecies%>%
  group_by(EcoSite, EventYear, TransectQuadratID)%>%
  summarize(quadrich = sum(SpeciesPresentInQuadratForNested))

species_wide = allspecies %>% 
  group_by(Park, EcoSite, EventYear, Plot, EventPanel,Quadrat, CurrentSpecies) %>% 
  mutate(CoverClass_avg = mean(NestedQuadratSizeClass)) %>% 
  select(Park:EventPanel, Quadrat,CurrentSpecies, CoverClass_avg) %>% 
  distinct() %>% #removed 2 taxa? from record, this is probably not correct
  pivot_wider(names_from = CurrentSpecies, values_from = CoverClass_avg, values_fill = 0) %>% 
  ungroup()


#For duplicates from above?
#re-run species_wide stopping before the select
#run below
# View(species_wide %>% 
#   filter(CoverClassMidpoint_Quadrat_pct!=CoverMidpoint_pct))

species_wide_int=species_wide %>% 
  select(-c(Park:Quadrat)) %>%
  mutate_all(., function(x) as.integer((as.character(x))))

species_wide_int=species_wide %>% 
  select(-c(Park:Quadrat)) %>% 
  mutate_all(., function(x) as.integer(ifelse(x>0, 1, 0)))

species_env = species_wide %>% 
  select(c(Park:Quadrat)) %>% 
  mutate(rownumber = row_number())




rich_calc = estimateR(species_wide_int) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(cols = starts_with("V"), names_to = "rownumber" ) %>% 
  pivot_wider(names_from = rowname) %>% 
  mutate(rownumber = as.integer(str_remove(rownumber, "V")))   

rich_quadrat=species_env %>% 
  left_join(rich_calc) 
rich_quadrat%>% 
  ggplot(aes(x=EventYear, y = S.obs, color = EcoSite))+#str_remove(EcoSite, "\\w{4}_")))+
  geom_boxplot()

#sit down and read this protocol, what are these numbers and letters, ugh!
names(rich_quadrat)

quadrats_notcollected=rich_quadrat %>% 
  group_by(Park, EcoSite, EventYear, Plot, EventPanel) %>% 
  tally() %>% 
  filter(n!=15)
quadrats_notcollected

rich_quadrat %>% 
  group_by(Park, EventYear) %>% 
  tally() 

rich_quadrat %>% 
  group_by(EcoSite, EventYear) %>% 
  tally() 


View(rich_quadrat %>% 
  group_by(Quadrat, EventYear) %>% 
  tally() )
View(rich_quadrat %>% 
       group_by(EventPanel, EventYear) %>% 
       tally() )

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
