packages <- c("tidyverse", "here","lme4","LMERConvenienceFunctions", "lmerTest", "emmeans", "multcomp", "MuMIn",'nortest', "vegan", "codyn", "data.table", "janitor") #list of packages to install
lapply(packages, library, character.only = T) #load packages

load(file=here("results","upload4"))
load(file=here("results","species_wide_df"))






species_wide_park_ecosite <- species_wide %>% 
  left_join(., park_ecosystem) %>%
  split(., species_wide_park_ecosite$Ecosystem)

species_ecosystem_ls <- names(species_wide_park_ecosite)



#thinking about redoing based on ecosystem with park_ecosite_eventpanel as the group variable
species_only <- species_wide_park_ecosite[[1]] %>%
  dplyr::select(-c(EventYear, Replicate, Park_Ecosite, EventPanel, Plot, Transect, Quadrat, Ecosystem, EcoSite, Park, EcoSiteName))

# groups <- species_wide_park[[2]] %>%
#   dplyr::select(EventYear:Quadrat, Replicate) %>% 
#   mutate(Park_Ecosite_EventPanel = paste(Park_Ecosite, EventPanel, sep ="_"))
# 
# group=species_wide_park[[2]] $Park_Ecosite_EventPanel
# 
# indval <- multipatt(species_only1, group, 
#                     control = how(nperm=999))
# a <-capture.output(summary(indval)) %>% 
#   as.data.frame() %>% 
#   rename(col = ".")
#   filter()


groups <- species_wide_park_ecosite[[1]] %>%
  dplyr::select(EventYear:Quadrat, Replicate, Ecosystem) %>% 
  mutate(Park_Ecosite_EventPanel = paste(Park_Ecosite, EventPanel, sep ="_"))

group <- groups$Park_Ecosite_EventPanel

t1 = Sys.time()
indval <- multipatt(species_only, group, 
                    control = how(nperm=999))
t2 = Sys.time()

a <-capture.output(summary(indval)) %>% 
  as.data.frame() %>% 
  rename(col = ".")
  filter()

  
species_wide_park_ecosite <- species_wide %>% 
  left_join(., park_ecosystem) 

species_wide_park_ecosite <- species_wide_park_ecosite %>%
  split(., species_wide_park_ecosite$Ecosystem)

species_ecosystem_ls <- names(species_wide_park_ecosite)


for (i in 1:length(species_ecosystem_ls)){ #run a loop over the dataframes in the list
  
#thinking about redoing based on ecosystem with park_ecosite_eventpanel as the group variable
species_only <- species_wide_park_ecosite[[i]] %>%
  dplyr::select(-c(EventYear, Replicate, Park_Ecosite, EventPanel, Plot, Transect, Quadrat, Ecosystem, EcoSite, Park, EcoSiteName))

groups <- species_wide_park_ecosite[[i]] %>%
  dplyr::select(EventYear:Quadrat, Replicate, Ecosystem) %>% 
  mutate(Park_Ecosite_EventPanel = paste(Park_Ecosite, EventPanel, sep ="_"))

group <- groups$Park_Ecosite_EventPanel

t1 = Sys.time()
indval <- multipatt(species_only, group, 
                    control = how(nperm=999))
t2 = Sys.time()

assign(here("results", "indicator_species", paste(species_ecosystem_ls[i], "indicator_summary")),
       
       capture.output(summary(indval)) %>% 
  as.data.frame() %>% 
  rename(col = "."))

}

#different attempt


species_wide_park_ecosite <- species_wide %>% 
  left_join(., park_ecosystem) 

species_wide_park_ecosite <- species_wide_park_ecosite %>%
  split(., species_wide_park_ecosite$Ecosystem)

species_ecosystem_ls <- names(species_wide_park_ecosite)

library(indicspecies)
for (i in 1:length(species_ecosystem_ls)){ #run a loop over the dataframes in the list
  
  #thinking about redoing based on ecosystem with park_ecosite_eventpanel as the group variable
  species_only <- species_wide_park_ecosite[[3]] %>%
    dplyr::select(-c(EventYear, Replicate, Park_Ecosite, EventPanel, Plot, Transect, Quadrat, Ecosystem, EcoSite, Park, EcoSiteName))
  
  groups <- species_wide_park_ecosite[[3]] %>%
    dplyr::select(EventYear:Quadrat, Replicate, Ecosystem) %>% 
    mutate(Park_Ecosite_EventPanel = paste(Park_Ecosite, EventPanel, sep ="_"))
  
  group <- groups$Park_Ecosite
  
  t1 = Sys.time()
  indval <- multipatt(species_only, group, 
                      control = how(nperm=999))
  t2 = Sys.time()
  
  assign(here("results", "indicator_species", paste(species_ecosystem_ls[i], "indicator_summary")),
         
         capture.output(summary(indval)) %>% 
           as.data.frame() %>% 
           rename(col = "."))
  
}

grasslands <- species_wide_park_ecosite[[1]] %>% 
  dplyr::select(Park_Ecosite, Ecosystem) %>% 
  distinct()

mixedconifer <- species_wide_park_ecosite[[2]] %>% 
  dplyr::select(Park_Ecosite, Ecosystem) %>% 
  distinct()

pj <- species_wide_park_ecosite[[3]] %>% 
  dplyr::select(Park_Ecosite, Ecosystem) %>% 
  distinct()

ponderosa <- species_wide_park_ecosite[[4]] %>% 
  dplyr::select(Park_Ecosite, Ecosystem) %>% 
  distinct()
