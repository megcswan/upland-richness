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

species_long <- read.csv("species_filtered.csv") #load long file of all species


#Explore community driven data
community_ecosite <- community %>% 
  filter(!(Park%in%c("AZRU", "CHCU", "PETR", "WACA"))) #remove sites with only 1 ecosite


community_ecosite_park <-community_ecosite  %>% #create a list of dataframes separated by park
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
par(mfrow = c(4, 3))
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

rich_plot = community_ecosite %>% 
  ggplot(aes(x = SampleYear, y = S.obs)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = S.obs, color = EcoSite),method = "glm" )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park)+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
rich_plot

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






MuMIn::r.squaredGLMM(get(paste("rich_mod", parks_ecosite[1], sep="_")))



mod_means_contr <- emmeans::emmeans(object = rich_lmer_cont,
                                    pairwise ~ EcoSite|SampleYear,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_contr$emmeans,
              Letters = letters)

rich_plot = pefo_comm %>% 
  ggplot(aes(x = SampleYear, y = S.obs)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = S.obs, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
rich_plot


### Sample year factor --------------------------------------------------


rich_lmer_factor<-lmer(sqrt(S.obs)~Year_factor*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= pefo_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(rich_lmer_factor)
plot(rich_lmer_factor)
qqnorm(resid(rich_lmer_factor))
hist(resid(rich_lmer_factor))
shapiro.test(resid(rich_lmer_factor))
MuMIn::r.squaredGLMM(rich_lmer_factor)



mod_means_factorr <- emmeans::emmeans(object = rich_lmer_factor,
                                    pairwise ~ EcoSite|Year_factor,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_factorr$emmeans,
              Letters = letters)

rich_plot = pefo_comm %>% 
  ggplot(aes(x = SampleYear, y = S.obs)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = S.obs, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
rich_plot



rich_plot_box = pefo_comm %>% 
  ggplot(aes(x = Year_factor, y = S.obs)) + 
  geom_jitter(aes(color = EcoSite), alpha = .1)+
  geom_boxplot(aes(x = Year_factor, y = S.obs, color = EcoSite))+
  facet_wrap(~EcoSite)+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
rich_plot_box


# Simpson -----------------------------------------------------------------

### Sample year continuous --------------------------------------------------

#model not meeting assumptions
simpson_lmer_cont<-lmer(exp(simpson)~SampleYear*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= pefo_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(simpson_lmer_cont)
plot(simpson_lmer_cont)
qqnorm(resid(simpson_lmer_cont))
hist(resid(simpson_lmer_cont))
shapiro.test(resid(simpson_lmer_cont))
MuMIn::r.squaredGLMM(simpson_lmer_cont)



mod_means_contr <- emmeans::emmeans(object = simpson_lmer_cont,
                                    pairwise ~ EcoSite|SampleYear,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_contr$emmeans,
              Letters = letters)

simpson_plot = pefo_comm %>% 
  ggplot(aes(x = SampleYear, y = simpson)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = simpson, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
simpson_plot


### Sample year factor --------------------------------------------------


simpson_lmer_factor<-lmer((simpson)~Year_factor*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= pefo_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(simpson_lmer_factor)
plot(simpson_lmer_factor)
qqnorm(resid(simpson_lmer_factor))
hist(resid(simpson_lmer_factor))
shapiro.test(resid(simpson_lmer_factor))
MuMIn::r.squaredGLMM(simpson_lmer_factor)



mod_means_factorr <- emmeans::emmeans(object = simpson_lmer_factor,
                                      pairwise ~ EcoSite|Year_factor,
                                      adjust = "tukey",
                                      pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_factorr$emmeans,
              Letters = letters)



simpson_plot_box = pefo_comm %>% 
  ggplot(aes(x = Year_factor, y = simpson)) + 
  geom_jitter(aes(color = EcoSite), alpha = .1)+
  geom_boxplot(aes(x = Year_factor, y = simpson, color = EcoSite))+
  facet_wrap(~EcoSite)+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
simpson_plot_box


# shannon -----------------------------------------------------------------

### Sample year continuous --------------------------------------------------

#model not meeting assumptions
shannon_lmer_cont<-lmer(exp(shannon)~SampleYear*EcoSite + (1|EventPanel:Plot:Transect:Quadrat) ,data= pefo_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(shannon_lmer_cont)
plot(shannon_lmer_cont)
qqnorm(resid(shannon_lmer_cont))
hist(resid(shannon_lmer_cont))
shapiro.test(resid(shannon_lmer_cont))
MuMIn::r.squaredGLMM(shannon_lmer_cont)



mod_means_contr <- emmeans::emmeans(object = shannon_lmer_cont,
                                    pairwise ~ EcoSite|SampleYear,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_contr$emmeans,
              Letters = letters)

shannon_plot = pefo_comm %>% 
  ggplot(aes(x = SampleYear, y = shannon)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = shannon, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
shannon_plot


### Sample year factor --------------------------------------------------


shannon_lmer_factor<-lmer(exp(shannon)~Year_factor*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= pefo_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(shannon_lmer_factor)
plot(shannon_lmer_factor)
qqnorm(resid(shannon_lmer_factor))
hist(resid(shannon_lmer_factor))
shapiro.test(resid(shannon_lmer_factor))
MuMIn::r.squaredGLMM(shannon_lmer_factor)



mod_means_factorr <- emmeans::emmeans(object = shannon_lmer_factor,
                                      pairwise ~ EcoSite|Year_factor,
                                      adjust = "tukey",
                                      pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_factorr$emmeans,
              Letters = letters)



shannon_plot_box = pefo_comm %>% 
  ggplot(aes(x = Year_factor, y = shannon)) + 
  geom_jitter(aes(color = EcoSite), alpha = .1)+
  geom_boxplot(aes(x = Year_factor, y = shannon, color = EcoSite))+
  facet_wrap(~EcoSite)+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
shannon_plot_box

# invD -----------------------------------------------------------------

### Sample year continuous --------------------------------------------------

#model not meeting assumptions
invD_lmer_cont<-lmer(sqrt(invD)~SampleYear*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= pefo_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(invD_lmer_cont)
plot(invD_lmer_cont)
qqnorm(resid(invD_lmer_cont))
hist(resid(invD_lmer_cont))
shapiro.test(resid(invD_lmer_cont))
MuMIn::r.squaredGLMM(invD_lmer_cont)



mod_means_contr <- emmeans::emmeans(object = invD_lmer_cont,
                                    pairwise ~ EcoSite|SampleYear,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_contr$emmeans,
              Letters = letters)

invD_plot = pefo_comm %>% 
  ggplot(aes(x = SampleYear, y = invD)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = invD, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
invD_plot


### Sample year factor --------------------------------------------------


invD_lmer_factor<-lmer(sqrt(invD)~Year_factor*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= pefo_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(invD_lmer_factor)
plot(invD_lmer_factor)
qqnorm(resid(invD_lmer_factor))
hist(resid(invD_lmer_factor))
shapiro.test(resid(invD_lmer_factor))
MuMIn::r.squaredGLMM(invD_lmer_factor)



mod_means_factorr <- emmeans::emmeans(object = invD_lmer_factor,
                                      pairwise ~ EcoSite|Year_factor,
                                      adjust = "tukey",
                                      pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_factorr$emmeans,
              Letters = letters)



invD_plot_box = pefo_comm %>% 
  ggplot(aes(x = Year_factor, y = invD)) + 
  geom_jitter(aes(color = EcoSite), alpha = .1)+
  geom_boxplot(aes(x = Year_factor, y = invD, color = EcoSite))+
  facet_wrap(~EcoSite)+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
invD_plot_box

# MEVE models --------------------------------------------

MEVE_comm <- community %>% 
  filter( Park=="MEVE")



# Richness models ---------------------------------------------------------


### Sample year continuous --------------------------------------------------


rich_lmer_cont<-lmer((S.obs)~SampleYear*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= MEVE_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(rich_lmer_cont)
plot(rich_lmer_cont)
qqnorm(resid(rich_lmer_cont))
hist(resid(rich_lmer_cont))
shapiro.test(resid(rich_lmer_cont))
MuMIn::r.squaredGLMM(rich_lmer_cont)



mod_means_contr <- emmeans::emmeans(object = rich_lmer_cont,
                                    pairwise ~ EcoSite|SampleYear,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_contr$emmeans,
              Letters = letters)

rich_plot = MEVE_comm %>% 
  ggplot(aes(x = SampleYear, y = S.obs)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = S.obs, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
rich_plot


### Sample year factor --------------------------------------------------

#removing S EcoSite here, too many missing years 
rich_lmer_factor<-lmer(sqrt(S.obs)~Year_factor + (1|EventPanel:Plot:Transect:Quadrat) ,data= subset(MEVE_comm, EcoSite=="L"))#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(rich_lmer_factor)
plot(rich_lmer_factor)
qqnorm(resid(rich_lmer_factor))
hist(resid(rich_lmer_factor))
shapiro.test(resid(rich_lmer_factor))
MuMIn::r.squaredGLMM(rich_lmer_factor)



mod_means_factorr <- emmeans::emmeans(object = rich_lmer_factor,
                                      pairwise ~ Year_factor,
                                      adjust = "tukey",
                                      pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_factorr$emmeans,
              Letters = letters)

rich_plot = MEVE_comm %>% 
  filter(EcoSite=="L") %>% 
  ggplot(aes(x = SampleYear, y = S.obs)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = S.obs, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
rich_plot



rich_plot_box = MEVE_comm %>% 
  filter(EcoSite=="L") %>% 
  ggplot(aes(x = Year_factor, y = S.obs)) + 
  geom_jitter(aes(color = EcoSite), alpha = .1)+
  geom_boxplot(aes(x = Year_factor, y = S.obs, color = EcoSite))+
  facet_wrap(~EcoSite)+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
rich_plot_box


# Simpson -----------------------------------------------------------------

### Sample year continuous --------------------------------------------------

#model not meeting assumptions
simpson_lmer_cont<-lmer(exp(simpson)~SampleYear*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= MEVE_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(simpson_lmer_cont)
plot(simpson_lmer_cont)
qqnorm(resid(simpson_lmer_cont))
hist(resid(simpson_lmer_cont))
shapiro.test(resid(simpson_lmer_cont))
MuMIn::r.squaredGLMM(simpson_lmer_cont)



mod_means_contr <- emmeans::emmeans(object = simpson_lmer_cont,
                                    pairwise ~ EcoSite|SampleYear,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_contr$emmeans,
              Letters = letters)

simpson_plot = MEVE_comm %>% 
  ggplot(aes(x = SampleYear, y = simpson)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = simpson, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
simpson_plot


### Sample year factor --------------------------------------------------

#subset by EcoSite L
#simpson is not working well with model assumptions

simpson_lmer_factor<-lmer(exp(simpson)~Year_factor + (1|EventPanel:Plot:Transect:Quadrat) ,data= subset(MEVE_comm, EcoSite =="L"))#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(simpson_lmer_factor)
plot(simpson_lmer_factor)
qqnorm(resid(simpson_lmer_factor))
hist(resid(simpson_lmer_factor))
shapiro.test(resid(simpson_lmer_factor))
MuMIn::r.squaredGLMM(simpson_lmer_factor)



mod_means_factorr <- emmeans::emmeans(object = simpson_lmer_factor,
                                      pairwise ~ Year_factor,
                                      adjust = "tukey",
                                      pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_factorr$emmeans,
              Letters = letters)



simpson_plot_box = MEVE_comm %>% 
  filter(EcoSite =="L") %>% 
  ggplot(aes(x = Year_factor, y = simpson)) + 
  geom_jitter(aes(color = EcoSite), alpha = .1)+
  geom_boxplot(aes(x = Year_factor, y = simpson, color = EcoSite))+
  facet_wrap(~EcoSite)+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
simpson_plot_box


# shannon -----------------------------------------------------------------

### Sample year continuous --------------------------------------------------

#model not meeting assumptions
shannon_lmer_cont<-lmer(exp(shannon)~SampleYear*EcoSite + (1|EventPanel:Plot:Transect:Quadrat) ,data= MEVE_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(shannon_lmer_cont)
plot(shannon_lmer_cont)
qqnorm(resid(shannon_lmer_cont))
hist(resid(shannon_lmer_cont))
shapiro.test(resid(shannon_lmer_cont))
MuMIn::r.squaredGLMM(shannon_lmer_cont)



mod_means_contr <- emmeans::emmeans(object = shannon_lmer_cont,
                                    pairwise ~ EcoSite|SampleYear,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_contr$emmeans,
              Letters = letters)

shannon_plot = MEVE_comm %>% 
  ggplot(aes(x = SampleYear, y = shannon)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = shannon, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
shannon_plot


### Sample year factor --------------------------------------------------


shannon_lmer_factor<-lmer(exp(shannon)~Year_factor + (1|EventPanel:Plot:Transect:Quadrat) ,data= subset(MEVE_comm, EcoSite =="L"))#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(shannon_lmer_factor)
plot(shannon_lmer_factor)
qqnorm(resid(shannon_lmer_factor))
hist(resid(shannon_lmer_factor))
shapiro.test(resid(shannon_lmer_factor))
MuMIn::r.squaredGLMM(shannon_lmer_factor)



mod_means_factorr <- emmeans::emmeans(object = shannon_lmer_factor,
                                      pairwise ~ Year_factor,
                                      adjust = "tukey",
                                      pbkrtest.limit = 27990)

meve_shannon_factor_cld <- multcomp::cld(object = mod_means_factorr$emmeans,
              Letters = letters) %>% 
  as.data.frame()

names(meve_shannon_factor_cld)


shannon_plot_box = MEVE_comm %>% 
  filter(EcoSite=="L" ) %>% 
  ggplot(aes(x = Year_factor, y = shannon)) + 
  geom_jitter(aes(color = EcoSite), alpha = .1)+
  geom_boxplot(aes(x = Year_factor, y = shannon, color = EcoSite))+
  facet_wrap(~EcoSite)+
  geom_text(data = meve_shannon_factor_cld,
            aes(x =Year_factor,
                y = 3,
                label = .group)) +  
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
shannon_plot_box

# invD -----------------------------------------------------------------

### Sample year continuous --------------------------------------------------

#model not meeting assumptions
invD_lmer_cont<-lmer(sqrt(invD)~SampleYear*EcoSite + (1|EventPanel:Plot:Transect:Quadrat) ,data= MEVE_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(invD_lmer_cont)
plot(invD_lmer_cont)
qqnorm(resid(invD_lmer_cont))
hist(resid(invD_lmer_cont))
shapiro.test(resid(invD_lmer_cont))
MuMIn::r.squaredGLMM(invD_lmer_cont)



mod_means_contr <- emmeans::emmeans(object = invD_lmer_cont,
                                    pairwise ~ EcoSite|SampleYear,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_contr$emmeans,
              Letters = letters)

invD_plot = MEVE_comm %>% 
  ggplot(aes(x = SampleYear, y = invD)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = invD, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
invD_plot


### Sample year factor --------------------------------------------------


invD_lmer_factor<-lmer(sqrt(invD)~Year_factor + (1|EcoSite:EventPanel:Plot:Transect:Quadrat) ,data= MEVE_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(invD_lmer_factor)
plot(invD_lmer_factor)
qqnorm(resid(invD_lmer_factor))
hist(resid(invD_lmer_factor))
shapiro.test(resid(invD_lmer_factor))
MuMIn::r.squaredGLMM(invD_lmer_factor)



mod_means_factorr <- emmeans::emmeans(object = invD_lmer_factor,
                                      pairwise ~ Year_factor,
                                      adjust = "tukey",
                                      pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_factorr$emmeans,
              Letters = letters)



invD_plot_box = MEVE_comm %>% 
  ggplot(aes(x = Year_factor, y = invD)) + 
  geom_jitter(aes(color = EcoSite), alpha = .1)+
  geom_boxplot(aes(x = Year_factor, y = invD))+
  #facet_wrap(~EcoSite)+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
invD_plot_box

parks= community %>% dplyr::select(Park) %>% distinct()
parks

# GRCA models --------------------------------------------

GRCA_comm <- community %>% 
  filter( Park=="GRCA")



# Richness models ---------------------------------------------------------


### Sample year continuous --------------------------------------------------

#not meeting model assumptions well
rich_lmer_cont<-lmer((S.obs)~SampleYear*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= GRCA_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(rich_lmer_cont)
plot(rich_lmer_cont)
qqnorm(resid(rich_lmer_cont))
hist(resid(rich_lmer_cont))
shapiro.test(resid(rich_lmer_cont))
MuMIn::r.squaredGLMM(rich_lmer_cont)



mod_means_contr <- emmeans::emmeans(object = rich_lmer_cont,
                                    pairwise ~ EcoSite|SampleYear,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_contr$emmeans,
              Letters = letters)

rich_plot = GRCA_comm %>% 
  ggplot(aes(x = SampleYear, y = S.obs)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = S.obs, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
rich_plot


### Sample year factor --------------------------------------------------


rich_lmer_factor<-lmer(sqrt(S.obs)~Year_factor*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= GRCA_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(rich_lmer_factor)
plot(rich_lmer_factor)
qqnorm(resid(rich_lmer_factor))
hist(resid(rich_lmer_factor))
shapiro.test(resid(rich_lmer_factor))
MuMIn::r.squaredGLMM(rich_lmer_factor)



mod_means_factorr <- emmeans::emmeans(object = rich_lmer_factor,
                                      pairwise ~ EcoSite|Year_factor,
                                      adjust = "tukey",
                                      pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_factorr$emmeans,
              Letters = letters)




rich_plot_box = GRCA_comm %>% 
  ggplot(aes(x = Year_factor, y = S.obs)) + 
  geom_jitter(aes(color = EcoSite), alpha = .1)+
  geom_boxplot(aes(x = Year_factor, y = S.obs, color = EcoSite))+
  facet_wrap(~EcoSite)+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
rich_plot_box


# Simpson -----------------------------------------------------------------

### Sample year continuous --------------------------------------------------

#model not meeting assumptions
simpson_lmer_cont<-lmer(exp(simpson)~SampleYear*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= GRCA_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(simpson_lmer_cont)
plot(simpson_lmer_cont)
qqnorm(resid(simpson_lmer_cont))
hist(resid(simpson_lmer_cont))
shapiro.test(resid(simpson_lmer_cont))
MuMIn::r.squaredGLMM(simpson_lmer_cont)



mod_means_contr <- emmeans::emmeans(object = simpson_lmer_cont,
                                    pairwise ~ EcoSite|SampleYear,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_contr$emmeans,
              Letters = letters)

simpson_plot = GRCA_comm %>% 
  ggplot(aes(x = SampleYear, y = simpson)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = simpson, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
simpson_plot


### Sample year factor --------------------------------------------------


simpson_lmer_factor<-lmer(exp(simpson)~Year_factor*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= GRCA_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(simpson_lmer_factor)
plot(simpson_lmer_factor)
qqnorm(resid(simpson_lmer_factor))
hist(resid(simpson_lmer_factor))
shapiro.test(resid(simpson_lmer_factor))
MuMIn::r.squaredGLMM(simpson_lmer_factor)



mod_means_factorr <- emmeans::emmeans(object = simpson_lmer_factor,
                                      pairwise ~ EcoSite|Year_factor,
                                      adjust = "tukey",
                                      pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_factorr$emmeans,
              Letters = letters)



simpson_plot_box = GRCA_comm %>% 
  ggplot(aes(x = Year_factor, y = simpson)) + 
  geom_jitter(aes(color = EcoSite), alpha = .1)+
  geom_boxplot(aes(x = Year_factor, y = simpson, color = EcoSite))+
  #facet_wrap(~EcoSite)+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
simpson_plot_box


# shannon -----------------------------------------------------------------

### Sample year continuous --------------------------------------------------

#model not meeting assumptions
shannon_lmer_cont<-lmer(exp(shannon)~SampleYear*EcoSite + (1|EventPanel:Plot:Transect:Quadrat) ,data= GRCA_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(shannon_lmer_cont)
plot(shannon_lmer_cont)
qqnorm(resid(shannon_lmer_cont))
hist(resid(shannon_lmer_cont))
shapiro.test(resid(shannon_lmer_cont))
MuMIn::r.squaredGLMM(shannon_lmer_cont)



mod_means_contr <- emmeans::emmeans(object = shannon_lmer_cont,
                                    pairwise ~ EcoSite|SampleYear,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_contr$emmeans,
              Letters = letters)

shannon_plot = GRCA_comm %>% 
  ggplot(aes(x = SampleYear, y = shannon)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = shannon, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
shannon_plot


### Sample year factor --------------------------------------------------


shannon_lmer_factor<-lmer(exp(shannon)~Year_factor*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= GRCA_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(shannon_lmer_factor)
plot(shannon_lmer_factor)
qqnorm(resid(shannon_lmer_factor))
hist(resid(shannon_lmer_factor))
shapiro.test(resid(shannon_lmer_factor))
MuMIn::r.squaredGLMM(shannon_lmer_factor)



mod_means_factorr <- emmeans::emmeans(object = shannon_lmer_factor,
                                      pairwise ~ EcoSite|Year_factor,
                                      adjust = "tukey",
                                      pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_factorr$emmeans,
              Letters = letters)



shannon_plot_box = GRCA_comm %>% 
  ggplot(aes(x = Year_factor, y = shannon)) + 
  geom_jitter(aes(color = EcoSite), alpha = .1)+
  geom_boxplot(aes(x = Year_factor, y = shannon, color = EcoSite))+
  #facet_wrap(~EcoSite)+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
shannon_plot_box

# invD -----------------------------------------------------------------

### Sample year continuous --------------------------------------------------
#says there is na in y but I'm not finding it
#model not meeting assumptions
invD_lmer_cont<-lmer((invD)~SampleYear*EcoSite + (1|EventPanel:Plot:Transect:Quadrat) ,data= GRCA_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(invD_lmer_cont)
plot(invD_lmer_cont)
qqnorm(resid(invD_lmer_cont))
hist(resid(invD_lmer_cont))
shapiro.test(resid(invD_lmer_cont))
MuMIn::r.squaredGLMM(invD_lmer_cont)



mod_means_contr <- emmeans::emmeans(object = invD_lmer_cont,
                                    pairwise ~ EcoSite|SampleYear,
                                    adjust = "tukey",
                                    pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_contr$emmeans,
              Letters = letters)

invD_plot = GRCA_comm %>% 
  ggplot(aes(x = SampleYear, y = invD)) + 
  geom_jitter(aes(color=EcoSite), alpha =0.1)+
  geom_smooth(aes(x = SampleYear, y = invD, color = EcoSite),method = "glm" )+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
invD_plot


### Sample year factor --------------------------------------------------


invD_lmer_factor<-lmer(sqrt(invD)~Year_factor*EcoSite + (EventPanel|Plot:Transect:Quadrat) ,data= GRCA_comm)#,REML = FALSE,
# control = lmerControl(optimizer ="Nelder_Mead"))

#isSingular message
#passes model assumptions
anova(invD_lmer_factor)
plot(invD_lmer_factor)
qqnorm(resid(invD_lmer_factor))
hist(resid(invD_lmer_factor))
shapiro.test(resid(invD_lmer_factor))
MuMIn::r.squaredGLMM(invD_lmer_factor)



mod_means_factorr <- emmeans::emmeans(object = invD_lmer_factor,
                                      pairwise ~ EcoSite|Year_factor,
                                      adjust = "tukey",
                                      pbkrtest.limit = 27990)

multcomp::cld(object = mod_means_factorr$emmeans,
              Letters = letters)



invD_plot_box = GRCA_comm %>% 
  ggplot(aes(x = Year_factor, y = invD)) + 
  geom_jitter(aes(color = EcoSite), alpha = .1)+
  geom_boxplot(aes(x = Year_factor, y = invD, color = EcoSite))+
  facet_wrap(~EcoSite)+
  scale_color_manual(values = c("blue", "orange"))+
  #scale_x_factorinuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_factorinuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "top")
invD_plot_box

