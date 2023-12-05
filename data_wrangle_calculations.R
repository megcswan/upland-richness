#Questions
#*Difference between columns cover and presence absence, discuss which to use
#*unknown species
#*look at reps graph, why is is so weird?
#*#EcoSite or ecosystem type
#*other variables (would the stuff on irma be useful?)

#Install packages
packages <- c("tidyverse", "here","codyn", "BiodiversityR")
lapply(packages, library, character.only = T)
# #Load data
# 
allspecies = read_csv("allspecies.csv") %>% dplyr::select(-"...1")

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
  dplyr::select(Park:EventPanel,Plot, Transect, Quadrat, ReportTaxon, Nativity, Lifeform, Duration,CurrentSpecies, PresentAbsent, Cover_pct) %>% 
  distinct() %>% #removed 2 taxa from record, add back later? 
  mutate(PresentAbsent = ifelse(PresentAbsent==2, 1, PresentAbsent))

species_corrected <- species_corrected%>% 
  mutate(EcoSite = str_remove(EcoSite, "\\w{4}_"), #Separate nested parameters
         Plot = as.factor(str_extract(Plot, "\\d{2}")),
         Quadrat = as.factor(Quadrat),
         Replicate = paste(Park, EcoSite,EventPanel, Plot, Transect, Quadrat, sep = "_"),
         Year = as.numeric(EventYear)) %>% 
  filter(EventPanel %in% c("A","B","C")) #Select only event panel a, b, and c


#To see duplicates that were combined run commented script below
# View(allspecies %>%
#        group_by(Park, EcoSite, EventYear, EventPanel, Plot, Transect, Quadrat, CurrentSpecies) %>%
#        mutate(PresentAbsent = sum(SpeciesPresentInQuadratForNested),
#               Cover_pct = sum(CoverClassMidpoint_Quadrat_pct)) %>%
#   filter(CoverClassMidpoint_Quadrat_pct!=Cover_pct))



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





# #Pivot wide to calculate richness  --------------------------------------


#Values are Present/Absent
species_wide <- species_corrected_nozero_sampleyears %>% 
  dplyr::select(-c(Cover_pct, ReportTaxon, Nativity, Lifeform, Duration)) %>% 
  pivot_wider(names_from = CurrentSpecies, values_from = PresentAbsent, values_fill = 0) %>% 
  ungroup()

#make dataframe that has only the species as columns
species_only=species_wide %>% 
  dplyr::select(-c(Park:Year)) 

#same table as above but as binary
species_binary <- species_only%>% 
  mutate_all(., function(x) as.integer(ifelse(x>0, 1, 0)))

#table with event details
species_env = species_wide %>% 
  dplyr::select(c(Park:Quadrat)) %>% 
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


species_wide_div <- species_corrected_nozero_sampleyears %>% 
  dplyr::select(-c(PresentAbsent, ReportTaxon, Nativity, Lifeform, Duration)) %>% 
  pivot_wider(names_from = CurrentSpecies, values_from = Cover_pct, values_fill = 0) %>% 
  ungroup()
#make dataframe that has only the species as columns
species_only_div=species_wide_div %>% 
  dplyr::select(-c(Park:Year))


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
  dplyr::select(simpson:rownumber)



# Make dataframe with vegan calculated metrics ----------------------------

community=rich_env %>% 
  left_join(div_calc) %>% 
  left_join(rich_calc_div)

View(community %>% 
  filter(S.obs!=S.obs_cover))


# #Files  written to .csv---------------------------------------------------

# 
# write.csv(species_corrected, file = "species_corrected.csv", row.names = F)
# write.csv(species_corrected_nozero_sampleyears, file = "species_filtered.csv", row.names = F)
# write.csv(community, file = "community.csv", row.names = F)
# write.csv(species_wide_div, file ="species_wide_cover.csv")
# created a file that has metrics on it



