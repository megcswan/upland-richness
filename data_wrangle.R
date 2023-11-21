#Install packages
packages <- c("tidyverse", "here","codyn", "BiodiversityR")
lapply(packages, library, character.only = T)

# #Load data
allspecies = read_csv("allspecies.csv") %>% select(-"...1")

#organizing columns so I understand nested design and to make it easier to work with each level of the design
# TransectQuadratID=Park_EcoSitePlot_Date_Transect_Quadrat
allspecies <- allspecies %>%
  mutate(TransectQuadrat = str_remove(TransectQuadratID, pattern ="\\w{4}_\\w{1}\\d{2}_\\d{8}_")) %>%
  separate(TransectQuadrat, into = c("Transect","Quadrat"), sep = "_")


#Make correction to data frame, need to add together 2 instances of Carex spp.
#See commented code below for the instances
species_dup <- allspecies %>% 
  group_by(Park, EcoSite, EventYear, EventPanel, Plot, Transect, Quadrat, CurrentSpecies) %>% 
  mutate(PresentAbsent = sum(SpeciesPresentInQuadratForNested)) %>% 
  select(Park:EventPanel,Plot, Transect, Quadrat, CurrentSpecies, PresentAbsent) %>% 
  distinct() %>% #removed 2 taxa from record, add back later? 
  mutate(PresentAbsent = ifelse(PresentAbsent==2, 1, PresentAbsent))#,
#EcoSite = str_remove(EcoSite, "\\w{4}_"),
#Plot = str_remove(Plot, "\\w{4}_")) 

#Write file for presence absence data
# write.csv(species_dup, file = "species_PresentAbsent_long.csv")

#For Carex spp. that were combined
#run below
# View(species_dup %>%
#   filter(CoverClassMidpoint_Quadrat_pct!=CoverMidpoint_pct))

#Pivot wide to calculate diversity metrics
species_wide <- species_dup %>% 
  pivot_wider(names_from = CurrentSpecies, values_from = PresentAbsent, values_fill = 0) %>% 
  ungroup()
# write.csv(species_wide, file = "species_PresentAbsent_wide.csv")

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

#make diversity dataframe
species_dup_div <- allspecies %>% 
  group_by(Park, EcoSite, EventYear, Plot, Transect, Quadrat, CurrentSpecies) %>% 
  mutate(Cover_pct = sum(CoverClassMidpoint_Quadrat_pct)) %>% 
  select(Park:EventPanel,Plot, Transect, Quadrat, CurrentSpecies, Cover_pct) %>% 
  distinct()  #removed 2 taxa from record, add back later? 

#write file for cover data
# write.csv(species_dup_div, file = "species_Cover_long.csv")



#same Carex spp. combined as above
# View(species_dup_div %>%
#   filter(CoverClassMidpoint_Quadrat_pct!=CoverMidpoint_pct))

#Pivot wide to calculate diversity metrics
species_wide_div <- species_dup_div %>% 
  pivot_wider(names_from = CurrentSpecies, values_from = Cover_pct, values_fill = 0) %>% 
  ungroup()

#make dataframe that has only the species as columns
species_only=species_wide_div %>% 
  select(-c(Park:Quadrat)) 

# #same table as above but as binary
# species_binary <- species_only%>% 
#   mutate_all(., function(x) as.integer(ifelse(x>0, 1, 0)))

#table with event details
species_env = species_wide_div %>% 
  select(c(Park:Quadrat)) %>% 
  mutate(rownumber = row_number())


#Calculate diversity metrics
div_calc = species_only %>% 
  mutate(simpson = diversity(., index = "simpson"),
         shannon = diversity(., index = "shannon"),
         invD = diversity(., index = "invsimp"),
         rownumber = row_number()) %>% 
  select(simpson:rownumber)



community=rich_env %>% 
  left_join(div_calc) 

#Write file with vegan calculated diversity metrics
# write.csv(community, file = "community.csv")
# created a file that has metrics on it



