#merge this into a more permanent file
#run richness on all of these first than do the community metrics
#run parks that have not been run yet (only 1 ecosite)




#Install packages
packages <- c("tidyverse", "here","lme4", "vegan","LMERConvenienceFunctions", "lmerTest", "emmeans", "multcomp") #list of packages to install
lapply(packages, library, character.only = T) #load packages



community <- read.csv("community.csv", header = T) %>% #read data file
  mutate(Plot = as.factor(Plot), #fix nested columns to better analyze data
         Quadrat = as.factor(Quadrat),
         SampleYear = EventYear -2007,
         Year_factor = as.factor(EventYear))

species_long <- read.csv("species_filtered.csv") %>%  #load long file of all species
  mutate(Replicate=paste( EventPanel,Plot, Transect, Quadrat, sep="_"),
                         Plot = as.factor(Plot), #fix nested columns to better analyze data
                         Quadrat = as.factor(Quadrat),
                         SampleYear = EventYear -2007,
                         Year_factor = as.factor(EventYear))


#Explore community driven data
community_ecosite <- community %>% 
  filter(!(Park%in%c("AZRU", "CHCU", "PETR", "WACA")))   #remove sites with only 1 ecosite
community_ecosite_park <-community_ecosite  %>%
  split(., community_ecosite$Park) #split example dataset by group factor

#create a list of parks_ecosite
parks_ecosite <- names(community_ecosite_park)

#create empty dataframes to load results into
results_df <- as.data.frame(matrix(ncol=6,nrow=3)) # make an empty dataframe for results
rich_anova_results <- data.frame() # make an empty dataframe for results

#Create richness models for sites with more than one EcoSite
#Model parameters:
#Y = Sqrt richness
#FixedEffects = Sample Year (EventYEar -2007) by EcoSite
#RandomEffects = 1 intercept, eventpanel *plot*Transect*Quadrat
for (i in 1:length(community_ecosite_park)){ #run a loop over the dataframes in the list
  rich_mod <- lmer(sqrt(S.obs)~SampleYear*EcoSite + (1|EventPanel:Plot:Transect:Quadrat) ,data= community_ecosite_park[[i]])
  results_df<-anova(rich_mod) %>%  #run anova
    as.data.frame() %>% #save results as dataframe
    rownames_to_column(var = "FixedEffect") %>%  #change rownames to a column, these are the fixed effects of the model (factors)
    rename(p.value="Pr(>F)") %>% 
    mutate(Park = parks_ecosite[i],
           Significance = case_when(p.value<0.0001 ~ "<0.0001",
                                    p.value<0.001 ~ "<0.001",
                                    p.value<0.01 ~ "<0.01",
                                    p.value<0.05 ~ "<0.05",
                                    p.value<0.1 ~ "<0.1",
                                    p.value >=0.1 ~ "NS"))
    assign(paste("rich_mod", parks_ecosite[i], sep ="_"), rich_mod)
    rich_anova_results = rbind(rich_anova_results, results_df)
}


#Explore model assumptions
par(mfrow = c(4, 3))
shapiro_community_ecosite <- data.frame()


for (i in 1:length(community_ecosite_park)){ #run a loop over the dataframes in the list
  qqnorm(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))), main = paste("Q-Q Plot", parks_ecosite[i], sep=" "))
  hist(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))), main = paste("Histogram", parks_ecosite[i], sep=" "))
  shapiro <- shapiro.test(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))))
  p_value = cbind(parks_ecosite[i], shapiro$p.value) 
  
  shapiro_community_ecosite <- rbind(shapiro_community_ecosite, p_value)

}


# Run models for parks that do not have more than one ecosite -------------

#Explore community driven data
community_noecosite <- community %>% 
  filter((Park%in%c("AZRU", "CHCU", "PETR", "WACA"))) #remove sites with only 1 ecosite


community_noecosite_park <-community_noecosite  %>% #create a list of dataframes separated by park
  split(., community_noecosite$Park) #split example dataset by group factor

#create a list of parks_noecosite
parks_noecosite <- names(community_noecosite_park)

#create empty dataframes to load results into
results_df <- as.data.frame(matrix(ncol=6,nrow=3)) # make an empty dataframe for results
rich_anova_results2 <- data.frame() # make an empty dataframe for results

#Create richness models for sites with more than one EcoSite
#Model parameters:
#Y = Sqrt richness
#FixedEffects = Sample Year (EventYEar -2007) by EcoSite
#RandomEffects = 1 intercept, eventpanel *plot*Transect*Quadrat
for (i in 1:length(community_noecosite_park)){ #run a loop over the dataframes in the list
  rich_mod <- lmer(sqrt(S.obs)~SampleYear + (1|EventPanel:Plot:Transect:Quadrat) ,data= community_noecosite_park[[i]])
  results_df<-anova(rich_mod) %>%  #run anova
    as.data.frame() %>% #save results as dataframe
    rownames_to_column(var = "FixedEffect") %>%  #change rownames to a column, these are the fixed effects of the model (factors)
    rename(p.value="Pr(>F)") %>% 
    mutate(Park = parks_noecosite[i],
           Significance = case_when(p.value<0.0001 ~ "<0.0001",
                                    p.value<0.001 ~ "<0.001",
                                    p.value<0.01 ~ "<0.01",
                                    p.value<0.05 ~ "<0.05",
                                    p.value<0.1 ~ "<0.1",
                                    p.value >=0.1 ~ "NS"))
  assign(paste("rich_mod", parks_noecosite[i], sep ="_"), rich_mod)
  rich_anova_results2 = rbind(rich_anova_results2, results_df)
}


#Explore model assumptions
par(mfrow = c(3, 3))
shapiro_community_noecosite <- data.frame()


for (i in 1:length(community_noecosite_park)){ #run a loop over the dataframes in the list
  qqnorm(resid(get(paste("rich_mod", parks_noecosite[i], sep="_"))), main = paste("Q-Q Plot", parks_noecosite[i], sep=" "))
  hist(resid(get(paste("rich_mod", parks_noecosite[i], sep="_"))), main = paste("Histogram", parks_noecosite[i], sep=" "))
  shapiro <- shapiro.test(resid(get(paste("rich_mod", parks_noecosite[i], sep="_"))))
  p_value = cbind(parks_noecosite[i], shapiro$p.value) 
  
  shapiro_community_noecosite <- rbind(shapiro_community_noecosite, p_value)
  
}

#Create richness models for sites with more than one EcoSite
#Model g
#Model parameters:
#Y = Sqrt richness
#FixedEffects = Sample Year (EventYEar -2007) by EcoSite
#RandomEffects = 1 intercept, eventpanel *plot*Transect*Quadrat


#This is taking forever and giving warnings skipping for now
# for (i in 1:length(community_ecosite_park)){ #run a loop over the dataframes in the list
#   rich_mod <- glmer(sqrt(S.obs)~SampleYear*EcoSite + (1|EventPanel:Plot:Transect:Quadrat) ,data= community_ecosite_park[[i]], family ="poisson")
#   results_df<-anova(rich_mod) %>%  #run anova
#     as.data.frame() %>% #save results as dataframe
#     rownames_to_column(var = "FixedEffect") %>%  #change rownames to a column, these are the fixed effects of the model (factors)
#     rename(p.value="Pr(>F)") %>% 
#     mutate(Park = parks_ecosite[i],
#            Significance = case_when(p.value<0.0001 ~ "<0.0001",
#                                     p.value<0.001 ~ "<0.001",
#                                     p.value<0.01 ~ "<0.01",
#                                     p.value<0.05 ~ "<0.05",
#                                     p.value<0.1 ~ "<0.1",
#                                     p.value >=0.1 ~ "NS"))
#   assign(paste("rich_mod", parks_ecosite[i], sep ="_"), rich_mod)
#   rich_anova_results = rbind(rich_anova_results, results_df)
# }

# 
# #Explore model assumptions
# par(mfrow = c(4, 3))
# shapiro_community_ecosite <- data.frame()
# 
# 
# for (i in 1:length(community_ecosite_park)){ #run a loop over the dataframes in the list
#   qqnorm(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))), main = paste("Q-Q Plot", parks_ecosite[i], sep=" "))
#   hist(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))), main = paste("Histogram", parks_ecosite[i], sep=" "))
#   shapiro <- shapiro.test(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))))
#   p_value = cbind(parks_ecosite[i], shapiro$p.value) 
#   
#   shapiro_community_ecosite <- rbind(shapiro_community_ecosite, p_value)
#   
# }
# Create plots for the parks_ecosite ----------------------------------------------

rich_plot_ecosite = community_ecosite %>% 
  ggplot(aes(x = SampleYear, y = S.obs)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = S.obs, color = EcoSite),method = "glm" )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park)+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
rich_plot_ecosite

rich_plot = community_noecosite %>% 
  ggplot(aes(x = SampleYear, y = S.obs)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = S.obs, color = EcoSite),method = "glm" )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park)+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
rich_plot


# #Calculate codyn metrics ------------------------------------------------

#Calculate community metrics from codyn

#Rank abundance

species_rank <- species_long %>% 
  group_by(Replicate, EventYear) %>% 
  mutate(Park_Ecosite = paste(Park, EcoSite, sep = "_")) %>% 
  mutate(rank = rank(-Cover_pct, ties.method = "average"))

sampleyears <- species_rank %>% 
  dplyr::select(Replicate, EventYear) %>% 
  unique()

#remove plots that were only sampled 1 year
#The following numbers are the data frames in species park that have plots that were only sampled 1 year: 
#3, 6, 8, 14, 16, 17

#create dataframe that tallys and filter plots/replicates sampled in only 1 year
species_multyear <- species_rank %>% ungroup() %>% 
  dplyr::select(Park_Ecosite,Replicate, Plot, EventYear) %>% 
  distinct()  %>% 
  group_by(Park_Ecosite,Replicate) %>%  
  add_tally() %>% 
  filter(n>1) %>% 
  mutate(rep = paste0(Park_Ecosite, "_", Replicate))

#remove records taht were only sampled in one year
species_rank <- species_rank %>% 
  mutate(rep = paste0(Park_Ecosite, "_", Replicate)) %>% 
  filter(rep%in%species_multyear$rep)

species_park <-species_rank  %>%
  split(., species_rank$Park_Ecosite) #split  dataset by Park_Ecosite

species_park_list=names(species_park) #create a list of the dataframe names
stability_df <- data.frame()#Create an empty dataframe for stability 
turnover_df <- data.frame()#Create an empty dataframe for turnover
rate_change_df <- data.frame()#Create an empty dataframe for rate_change
rate_change_interval_df <- data.frame()#Create an empty dataframe for rate_change_interval
synchrony_df <- data.frame()#Create an empty dataframe for synchrony
variance_ratio_df <- data.frame()#Create an empty dataframe for variance_ratio
rank_shift_df <- data.frame()#Create an empty dataframe for rank_shift
appearance_df <- data.frame()
disappearance_df <- data.frame()


for (i in 1:length(species_park_list)){ #run a loop over the dataframes in the list
  
  stability_park_df <- codyn::community_stability(species_park[[i]], time.var = "EventYear", abundance.var = "Cover_pct",  replicate.var = "Replicate") %>% 
    mutate(Park_Ecosite = species_park_list[i] ) %>% 
    separate(Park_Ecosite, sep = "_", into = c("Park", "EcoSite"), remove=FALSE)
  stability_df <- rbind(stability_df, stability_park_df)
  
  turnover_park_df <- codyn::turnover(species_park[[i]], time.var = "EventYear", abundance.var = "Cover_pct",  replicate.var = "Replicate", species.var = "CurrentSpecies", metric="total") %>%
    mutate(Park_Ecosite = species_park_list[i] )%>% 
    separate(Park_Ecosite, sep = "_", into = c("Park", "EcoSite"), remove=FALSE)
  turnover_df <- rbind(turnover_df, turnover_park_df)
  
  
  appearance_park_df <- codyn::turnover(species_park[[i]], time.var = "EventYear", abundance.var = "Cover_pct",  replicate.var = "Replicate", species.var = "CurrentSpecies", metric="appearance") %>%
    mutate(Park_Ecosite = species_park_list[i] )%>% 
    separate(Park_Ecosite, sep = "_", into = c("Park", "EcoSite"), remove=FALSE)
  appearance_df <- rbind(appearance_df, appearance_park_df)
  
  disappearance_park_df <- codyn::turnover(species_park[[i]], time.var = "EventYear", abundance.var = "Cover_pct",  replicate.var = "Replicate", species.var = "CurrentSpecies", metric="disappearance") %>%
    mutate(Park_Ecosite = species_park_list[i] )%>% 
    separate(Park_Ecosite, sep = "_", into = c("Park", "EcoSite"), remove=FALSE)
  disappearance_df <- rbind(disappearance_df, disappearance_park_df)
  
  turnover_park_df <- codyn::turnover(species_park[[i]], time.var = "EventYear", abundance.var = "Cover_pct",  replicate.var = "Replicate", species.var = "CurrentSpecies", metric="total") %>%
    mutate(Park_Ecosite = species_park_list[i] )%>% 
    separate(Park_Ecosite, sep = "_", into = c("Park", "EcoSite"), remove=FALSE)
  turnover_df <- rbind(turnover_df, turnover_park_df)
  
  rate_change_park_df <- codyn::rate_change(species_park[[i]], time.var = "EventYear", abundance.var = "Cover_pct",  replicate.var = "Replicate", species.var = "CurrentSpecies") %>%
    mutate(Park_Ecosite = species_park_list[i] )%>% 
    separate(Park_Ecosite, sep = "_", into = c("Park", "EcoSite"), remove=FALSE)
  rate_change_df <- rbind(rate_change_df, rate_change_park_df)
  
  rate_change_interval_park_df <- codyn::rate_change_interval(species_park[[i]], time.var = "EventYear", abundance.var = "Cover_pct",  replicate.var = "Replicate", species.var = "CurrentSpecies") %>%
    mutate(Park_Ecosite = species_park_list[i] )%>% 
    separate(Park_Ecosite, sep = "_", into = c("Park", "EcoSite"), remove=FALSE)
  rate_change_interval_df <- rbind(rate_change_interval_df, rate_change_interval_park_df)
  
  rank_shift_park_df <- codyn::rank_shift(species_park[[i]], time.var = "EventYear", abundance.var = "Cover_pct",  replicate.var = "Replicate", species.var = "CurrentSpecies") %>%
    mutate(Park_Ecosite = species_park_list[i] )%>% 
    separate(Park_Ecosite, sep = "_", into = c("Park", "EcoSite"), remove=FALSE)
  rank_shift_df <- rbind(rank_shift_df, rank_shift_park_df)
  
}



#Tried to get synchrony but it is not working, says needs more than 1 species per replicate
#Tried several ways to remove plots with only 1 species, but it still does not work
#Below, I tried to filter by species richness
#checked warnings,   One or more species has non-varying abundance within a subplot and has been omitted
#It is likely that the abundance does not change, so may not be able to use this metric
#Skip for now look into later

# community_multispecies<- community %>% 
#   mutate(Park_Ecosite = paste(Park, EcoSite, sep = "_")) %>% 
#   filter(S.obs>1) %>% 
#   mutate(Park_Ecosite_Year_Replicate = paste(Park, EcoSite, EventYear, Plot, Transect, Quadrat, sep = "_"))
# 
# species_multi <- species_rank%>% #remove plots that contain only 1 species, this is a requirement for these calculations
#   mutate(Park_Ecosite_Year_Replicate = paste(Park, EcoSite, EventYear, Plot, Transect, Quadrat, sep = "_")) %>% 
# filter(Park_Ecosite_Year_Replicate%in%community_multispecies$Park_Ecosite_Year_Replicate)
# 
# species_multi_park <-species_multi  %>%
#   split(., species_multi$Park_Ecosite) #
# 
# 
# 
# for (i in 1:length(species_park_list)){ #run a loop over the dataframes in the list
# 
#   
#   synchrony_park_df <- codyn::synchrony(species_multi_park[[i]], time.var = "EventYear", abundance.var = "Cover_pct",  replicate.var = "Replicate", species.var = "CurrentSpecies", metric="Gross") %>%
#     mutate(Park_Ecosite = species_park_list[i] )%>% 
#     separate(Park_Ecosite, sep = "_", into = c("Park", "EcoSite"))
#   synchrony_df <- rbind(synchrony_df, synchrony_park_df)
#   
#   }
# 
# species_multi_park[[2]] %>% 
#   group_by(Replicate, CurrentSpecies) %>% 
#   tally()
# synchrony_park_df <- codyn::synchrony(species_multi_park[[4]], time.var = "EventYear", abundance.var = "Cover_pct",  replicate.var = "Replicate", species.var = "CurrentSpecies", metric="Loreau") %>%
#   

#Same thing for variance ratio, it does work for some but not all parks
# variance_ratio_park_df <- codyn::variance_ratio(species_park[[6]], time.var = "EventYear", abundance.var = "Cover_pct",  replicate.var = "Replicate", species.var = "CurrentSpecies",  bootnumber = 10,average.replicates = T) %>%
#   mutate(Park_Ecosite = species_park_list[i] )%>% 
#   separate(Park_Ecosite, sep = "_", into = c("Park", "EcoSite"))
# variance_ratio_df <- rbind(variance_ratio_df, variance_ratio_park_df)

#Create visualizations
stability_df %>% 
  ggplot(aes(x = Park, y =stability)) +
  geom_boxplot() +
  ylim(0,100)

turnover_df %>% 
  ggplot(aes(y = total, x = EventYear, color = Park, linetype = EcoSite))+
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("black","red","purple", "goldenrod", "blue", "lightblue", "lavender","darkgreen", "orange","green"))+
  facet_wrap(~Park)


rate_change_df %>% 
  ggplot(aes(x = paste(Park, EcoSite, sep = "_"), y =rate_change)) +
  geom_boxplot() 

rate_change_interval_df %>% 
  ggplot(aes(x = interval, y = distance, color= Park, linetype=EcoSite)) +
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("black","red","purple", "goldenrod", "blue", "lightblue", "lavender","darkgreen", "orange","green"))+
  facet_wrap(~Park)

#need to figure out a better way to graph
rank_shift_df %>% 
  mutate(year =as.numeric(substr(year_pair, 6,9))) %>% 
  ggplot(aes(y = MRS, x = year, color = Park, linetype = EcoSite))+
  stat_smooth(method = "lm")+
  #scale_color_manual(values = c("black","red","purple", "goldenrod", "blue", "lightblue", "lavender","darkgreen", "orange","green"))+
  facet_wrap(~Park)



#run analyses on these factors
#create rank clock
# Create rank abundance curves by park ecosite ----------------------------


species_rank_park_avg <- species_long %>% 
  ungroup() %>% 
  mutate(Park_Ecosite = paste(Park, EcoSite, sep = "_")) %>% 
  group_by(Park_Ecosite,EventYear, EventPanel,CurrentSpecies) %>% 
  summarise(Cover_plot_avg= mean(Cover_pct)) %>% 
  mutate(rank = rank(-Cover_plot_avg, ties.method = "average")) %>% 
  ungroup() 

species_rank_park_avg <-species_rank_park_avg %>% 
  split(., species_rank_park_avg$Park_Ecosite)

species_park_list <- names(species_rank_park_avg)

####make RACs graphs for each ecosite

top_species_df = data.frame()
for (i in 1:length(species_park_list)){ #run a loop over the dataframes in the list
  top_species <- species_rank_park_avg[[i]] %>% 
    filter(rank<=3) %>% 
    dplyr::select(CurrentSpecies, rank) %>% 
    arrange(rank) %>% 
    mutate(Park_Ecosite =species_park_list[i] ) %>% 
    dplyr::select(-rank) %>% 
    distinct()
  
  top_species_df= rbind(top_species_df, top_species)
  

  assign(paste("RAC_plot", species_park_list[i], sep="_"),
         species_rank_park_avg[[i]] %>% 
           mutate(colorfill = ifelse(CurrentSpecies%in%top_species$CurrentSpecies, CurrentSpecies, "other"))%>% 
           ggplot(aes(x=rank, y=Cover_plot_avg))+
           geom_line(color="darkgray", linewidth=1)+
           geom_point(aes(color = colorfill), size=3)+#aes(color=colorfill), size=3)+
           facet_wrap(~paste(EventYear, EventPanel))+
           #scale_color_manual(values = c("blue","lightblue","darkred", "pink","green3","orange","cornflowerblue", "darkgreen","purple", "lavender"))+
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", axis.title.x = element_blank(), axis.title.y = element_blank()))
 
}


top_species_park <- top_species_df %>% 
  split(., top_species_df$Park_Ecosite)



RAC_plot_AZRU_A
RAC_plot_BAND_M
RAC_plot_BAND_P
RAC_plot_CHCU_S
RAC_plot_GLCA_B
RAC_plot_GLCA_H
RAC_plot_GRCA_P
RAC_plot_GRCA_M
RAC_plot_MEVE_L
RAC_plot_MEVE_S
RAC_plot_PEFO_C
RAC_plot_PEFO_S
RAC_plot_PETR_P
RAC_plot_WACA_P
RAC_plot_WUPA_L
RAC_plot_WUPA_S


  
for (i in 1:13){ #run a loop over the dataframes in the list
  
#need to fix colors, increase number of species presented, or figure out another way to label the top 5 species or something
    aggdat <- aggregate(Cover_plot_avg ~ CurrentSpecies * EventYear * Park_Ecosite*EventPanel, 
                    data = subset(species_rank_park_avg[[i]], 
                                  CurrentSpecies=="Bromus tectorum"),
                                  #CurrentSpecies%in%top_species_park[[3]]$CurrentSpecies), 
                    FUN = mean)
    
    assign(paste("BT_plot", species_park_list[i], sep="_"),
       aggdat %>% 
         ggplot( aes(EventYear, Cover_plot_avg, color = CurrentSpecies)) + 
         geom_line(size = 2) + 
         coord_polar() + 
         theme_bw() + 
         facet_wrap(~Park_Ecosite*EventPanel) +
         ggtitle(paste0("Bromus tectorum abundances \n", aggdat$Park_Ecosite)))
         #scale_color_manual(values = c("blue","lightblue","darkred", "pink","green3","orange","cornflowerblue", "darkgreen","purple", "lavender"))+
         #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", axis.title.x = element_blank(), axis.title.y = element_blank()))

 }

BT_plot_AZRU_A
BT_plot_BAND_M
BT_plot_BAND_P
BT_plot_CHCU_S
BT_plot_GLCA_B
BT_plot_GLCA_H
BT_plot_GRCA_M
BT_plot_MEVE_L
BT_plot_MEVE_S
BT_plot_PEFO_C


#Is this because of panels????


#*Things to do:
#*1. rank abundance curves by event panel?
#*2. what to do about event panel anyway? discuss with megan. were plots in event panels randomized?
#*3. calculate other metrics with codyn from Hallett paper
#*4. organize this script
#*5. look into pft and functional diversity?

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
  select(-Genus, -Species)

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