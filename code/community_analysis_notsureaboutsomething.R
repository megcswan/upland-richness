#analysis following last update to data_wrangle2_codyn
#Install packages
packages <- c("tidyverse", "here","lme4","LMERConvenienceFunctions", "lmerTest", "emmeans", "multcomp", "MuMIn",'nortest') #list of packages to install
lapply(packages, library, character.only = T) #load packages

load(file=here("results","calculated_dfs3"))



# Richness models ---------------------------------------------------------



rich_ecosite_anova_results <- data.frame() # make an empty dataframe for results


#Create richness models for sites with more than one EcoSite
#Model parameters:
#Y = Sqrt richness
#FixedEffects = Sample Year (EventYEar -4007) by EcoSite
#RandomEffects = 1 intercept, eventpanel *plot*Transect*Quadrat

for (i in 1:length(community_ecosite_park)){ #run a loop over the dataframes in the list
  rich_mod <- lmer(sqrt(richness)~SampleYear*EcoSite + (1|EventPanel:Plot:Transect:Quadrat) ,data= community_ecosite_park[[i]])
  results_df<-anova(rich_mod) %>%  #run anova
    as.data.frame() %>% #save results as dataframe
    rownames_to_column(var = "FixedEffect") %>%  #change rownames to a column, these are the fixed effects of the model (factors)
    rename(p.value="Pr(>F)") %>% 
    mutate(Park = parks_ecosite[i],
           Significance = case_when(p.value<0.0001 ~ "****",
                                    p.value<0.001 ~ "***",
                                    p.value<0.01 ~ "**",
                                    p.value<0.05 ~ "*",
                                    p.value<0.1 ~ " †",
                                    p.value >=0.1 ~ "NS"))
  assign(paste("rich_mod", parks_ecosite[i], sep ="_"), rich_mod)
  rich_ecosite_anova_results = rbind(rich_ecosite_anova_results, results_df)
}

rich_ecosite_anova_results_ls <- rich_ecosite_anova_results %>% 
  split(., rich_ecosite_anova_results$Park)
ecosite_rich_model_results=data.frame()

for (i in 1:length(parks_ecosite)){
  ecosite_rich_model_results_park<-rich_ecosite_anova_results %>% 
    filter(Park%in%parks_ecosite[i]) %>% 
    mutate(p_value =  case_when(p.value<0.0001 ~ "<0.0001",
                                p.value<0.001 ~ "<0.001",
                                p.value<0.01 ~ as.character(round(p.value, 4)),
                                p.value<0.05 ~ as.character(round(p.value, 4)),
                                p.value<0.1 ~ as.character(round(p.value, 4)),
                                p.value >=0.1 ~ as.character(round(p.value, 4)))) %>% 
    dplyr::select(Park, FixedEffect, `Sum Sq`, `Mean Sq`, NumDF, DenDF, `F value`,
                  p_value, Significance)
  assign(paste("rich_mod_results", parks_ecosite[i], sep ="_"),ecosite_rich_model_results_park)
  ecosite_rich_model_results <-  rbind(ecosite_rich_model_results,ecosite_rich_model_results_park)
}

# #Explore model assumptions
par(mfrow = c(1,2 ))

  qqnorm(resid(get(paste("rich_mod", parks_ecosite[1], sep="_"))), main = paste("Q-Q Plot", parks_ecosite[1], sep=" "))
  hist(resid(get(paste("rich_mod", parks_ecosite[1], sep="_"))), main = paste("Histogram", parks_ecosite[1], sep=" "))
  ad_test_park <- ad.test(resid(get(paste("rich_mod", parks_ecosite[1], sep="_"))))
  ad_test_park
  
  qqnorm(resid(get(paste("rich_mod", parks_ecosite[2], sep="_"))), main = paste("Q-Q Plot", parks_ecosite[2], sep=" "))
  hist(resid(get(paste("rich_mod", parks_ecosite[2], sep="_"))), main = paste("Histogram", parks_ecosite[2], sep=" "))
  ad_test_park <- ad.test(resid(get(paste("rich_mod", parks_ecosite[2], sep="_"))))
  ad_test_park
  
  qqnorm(resid(get(paste("rich_mod", parks_ecosite[3], sep="_"))), main = paste("Q-Q Plot", parks_ecosite[3], sep=" "))
  hist(resid(get(paste("rich_mod", parks_ecosite[3], sep="_"))), main = paste("Histogram", parks_ecosite[3], sep=" "))
  ad_test_park <- ad.test(resid(get(paste("rich_mod", parks_ecosite[3], sep="_"))))
  ad_test_park



#Below does not seem to be populating correct values so going to look at individuallly
par(mfrow = c(4, 3))
ad_test_community_ecosite <- data.frame()
#shapiro_community_ecosite <- data.frame()

# for (i in 1:length(community_ecosite_park)){ #run a loop over the dataframes in the list
#   qqnorm(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))), main = paste("Q-Q Plot", parks_ecosite[i], sep=" "))
#   hist(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))), main = paste("Histogram", parks_ecosite[i], sep=" "))
#   ad_test_park <- ad.test(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))))
#   park = parks_ecosite[i]
#   ad_test_value = ad_test$p.value
#   p_value_park = cbind(park, ad_test_value)
#   
# 
#   ad_test_community_ecosite <- rbind(ad_test_community_ecosite, p_value_park)
#   
#   # shapiro <- shapiro.test(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))))
#   # p_value = cbind(parks_ecosite[i], shapiro$p.value)
#   # 
#   # shapiro_community_ecosite <- rbind(shapiro_community_ecosite, p_value)
# 
# }
# ad_test_community_ecosite
#Create a loop to present the findings

# rich_mod_results <- ls(pattern="rich_mod_results")
# for (i in 1:length(rich_mod_results)){
#   print(knitr::kable(get(rich_mod_results[i])))
# }

#knitr::kable(ecosite_rich_model_results)

# Run models for parks that do not have more than one ecosite -------------


parks_noecosite <- names(community_noecosite_park)

#create empty dataframes to load results into
rich_noecosite_anova_results <- data.frame() # make an empty dataframe for results

#Create richness models for sites with more than one EcoSite
#Model parameters:
#Y = Sqrt richness
#FixedEffects = Sample Year (EventYEar -4007) by EcoSite
#RandomEffects = 1 intercept, eventpanel *plot*Transect*Quadrat
for (i in 1:length(community_noecosite_park)){ #run a loop over the dataframes in the list
  rich_mod <- lmer(sqrt(richness)~SampleYear + (1|EventPanel:Plot:Transect:Quadrat) ,data= community_noecosite_park[[i]])
  results_df<-anova(rich_mod) %>%  #run anova
    as.data.frame() %>% #save results as dataframe
    rownames_to_column(var = "FixedEffect") %>%  #change rownames to a column, these are the fixed effects of the model (factors)
    rename(p.value="Pr(>F)") %>% 
    mutate(Park = parks_noecosite[i],
           Significance = case_when(p.value<0.0001 ~ "****",
                                    p.value<0.001 ~ "***",
                                    p.value<0.01 ~ "**",
                                    p.value<0.05 ~ "*",
                                    p.value<0.1 ~ "†",
                                    p.value >=0.1 ~ "NS"))
  
  assign(paste("rich_mod", parks_noecosite[i], sep ="_"), rich_mod)
  rich_noecosite_anova_results <-  rbind(rich_noecosite_anova_results, results_df)
}

rich_noecosite_anova_results_ls <- rich_noecosite_anova_results %>% 
  split(., rich_noecosite_anova_results$Park)
noecosite_rich_model_results=data.frame()

for (i in 1:length(parks_noecosite)){
  noecosite_rich_model_results_park<-rich_noecosite_anova_results %>% 
    filter(Park%in%parks_noecosite[i]) %>% 
    mutate(p_value =  case_when(p.value<0.0001 ~ "<0.0001",
                                p.value<0.001 ~ "<0.001",
                                p.value<0.01 ~ as.character(round(p.value, 4)),
                                p.value<0.05 ~ as.character(round(p.value, 4)),
                                p.value<0.1 ~ as.character(round(p.value, 4)),
                                p.value >=0.1 ~ as.character(round(p.value, 4)))) %>% 
    dplyr::select(Park, FixedEffect, `Sum Sq`, `Mean Sq`, NumDF, DenDF, `F value`,
                  p_value, Significance)
  assign(paste("rich_mod_results", parks_noecosite[i], sep ="_"),noecosite_rich_model_results_park)
  noecosite_rich_model_results <-  rbind(noecosite_rich_model_results,noecosite_rich_model_results_park)
}


#Explore model assumptions
par(mfrow = c(4, 3))
ad_test_community_ecosite <- data.frame()


for (i in 1:length(community_ecosite_park)){ #run a loop over the dataframes in the list
  qqnorm(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))), main = paste("Q-Q Plot", parks_ecosite[i], sep=" "))
  hist(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))), main = paste("Histogram", parks_ecosite[i], sep=" "))
  ad_test <- ad.test(resid(get(paste("rich_mod", parks_ecosite[i], sep="_"))))
  p_value = cbind(parks_ecosite[i], ad_test$p.value)

  ad_test_community_ecosite <- rbind(ad_test_community_ecosite, p_value)

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



# RAC models ---------------------------------------------------------


rac_ecosite_anova_results <- data.frame() # make an empty dataframe for results

#Create richness models for sites with more than one EcoSite
#Model parameters:
#Y = Sqrt richness
#FixedEffects = Sample Year (EventYEar -4007) by EcoSite
#RandomEffects = 1 intercept, eventpanel *plot*Transect*Quadrat

for (i in 1:length(community_ecosite_park)){ #run a loop over the dataframes in the list
  assign(paste("rac_mod", parks_ecosite[i], sep="_"),
         lmer(MRS~SampleYear*EcoSite + (1|Plot/Transect/Quadrat) ,
              data= community_ecosite_park[[i]]))
  
  
  results_df<-anova(get(paste("rac_mod", parks_ecosite[i], sep="_"))) %>%  #run anova
    as.data.frame() %>% #save results as dataframe
    rownames_to_column(var = "FixedEffect") %>%  #change rownames to a column, these are the fixed effects of the model (factors)
    rename(p.value="Pr(>F)") %>% 
    mutate(Park = parks_ecosite[i],
           Significance = case_when(p.value<0.0001 ~ "****",
                                    p.value<0.001 ~ "***",
                                    p.value<0.01 ~ "**",
                                    p.value<0.05 ~ "*",
                                    p.value<0.1 ~ " †",
                                    p.value >=0.1 ~ "NS"))
  rac_ecosite_anova_results = rbind(rac_ecosite_anova_results, results_df)
}

rac_ecosite_anova_results_ls <- rac_ecosite_anova_results %>% 
  split(., rac_ecosite_anova_results$Park)


ecosite_rac_model_results=data.frame()

for (i in 1:length(parks_ecosite)){
  ecosite_rac_model_results_park<-rac_ecosite_anova_results %>% 
    filter(Park%in%parks_ecosite[i]) %>% 
    mutate(p_value =  case_when(p.value<0.0001 ~ "<0.0001",
                                p.value<0.001 ~ "<0.001",
                                p.value<0.01 ~ as.character(round(p.value, 4)),
                                p.value<0.05 ~ as.character(round(p.value, 4)),
                                p.value<0.1 ~ as.character(round(p.value, 4)),
                                p.value >=0.1 ~ as.character(round(p.value, 4)))) %>% 
    dplyr::select(Park, FixedEffect, `Sum Sq`, `Mean Sq`, NumDF, DenDF, `F value`,
                  p_value, Significance)
  assign(paste("rac_mod_results", parks_ecosite[i], sep ="_"),ecosite_rac_model_results_park)
  ecosite_rac_model_results <-  rbind(ecosite_rac_model_results,ecosite_rac_model_results_park)
}


#Explore model assumptions
par(mfrow = c(4, 3))
ad_test_community_ecosite <- data.frame()


for (i in 1:length(community_ecosite_park)){ #run a loop over the dataframes in the list
  qqnorm(resid(get(paste("rac_mod", parks_ecosite[i], sep="_"))), main = paste("Q-Q Plot", parks_ecosite[i], sep=" "))
  hist(resid(get(paste("rac_mod", parks_ecosite[i], sep="_"))), main = paste("Histogram", parks_ecosite[i], sep=" "))
  ad_test <- ad.test(resid(get(paste("rac_mod", parks_ecosite[i], sep="_"))))
  p_value = cbind(parks_ecosite[i], ad_test$p.value)
  
  ad_test_community_ecosite <- rbind(ad_test_community_ecosite, p_value)
  
}

ad.test(resid(get(paste("rac_mod", parks_ecosite[6], sep="_"))))
# Run models for parks that do not have more than one ecosite -------------


parks_noecosite <- names(community_noecosite_park)

#create empty dataframes to load results into
rac_noecosite_anova_results <- data.frame() # make an empty dataframe for results

#Create richness models for sites with more than one EcoSite
#Model parameters:
#Y = Sqrt richness
#FixedEffects = Sample Year (EventYEar -4007) by EcoSite
#RandomEffects = 1 intercept, eventpanel *plot*Transect*Quadrat
for (i in 1:length(community_noecosite_park)){ #run a loop over the dataframes in the list
  rac_mod <- lmer(log1p(MRS)~SampleYear + (1|Plot:Transect:Quadrat) ,data= community_noecosite_park[[i]])
  results_df<-anova(rac_mod) %>%  #run anova
    as.data.frame() %>% #save results as dataframe
    rownames_to_column(var = "FixedEffect") %>%  #change rownames to a column, these are the fixed effects of the model (factors)
    rename(p.value="Pr(>F)") %>% 
    mutate(Park = parks_noecosite[i],
           Significance = case_when(p.value<0.0001 ~ "****",
                                    p.value<0.001 ~ "***",
                                    p.value<0.01 ~ "**",
                                    p.value<0.05 ~ "*",
                                    p.value<0.1 ~ "†",
                                    p.value >=0.1 ~ "NS"))
  
  assign(paste("rac_mod", parks_noecosite[i], sep ="_"), rac_mod)
  rac_noecosite_anova_results <-  rbind(rac_noecosite_anova_results, results_df)
}

rac_noecosite_anova_results_ls <- rac_noecosite_anova_results %>% 
  split(., rac_noecosite_anova_results$Park)
noecosite_rac_model_results=data.frame()

for (i in 1:length(parks_noecosite)){
  noecosite_rac_model_results_park<-rac_noecosite_anova_results %>% 
    filter(Park%in%parks_noecosite[i]) %>% 
    mutate(p_value =  case_when(p.value<0.0001 ~ "<0.0001",
                                p.value<0.001 ~ "<0.001",
                                p.value<0.01 ~ as.character(round(p.value, 4)),
                                p.value<0.05 ~ as.character(round(p.value, 4)),
                                p.value<0.1 ~ as.character(round(p.value, 4)),
                                p.value >=0.1 ~ as.character(round(p.value, 4)))) %>% 
    dplyr::select(Park, FixedEffect, `Sum Sq`, `Mean Sq`, NumDF, DenDF, `F value`,
                  p_value, Significance)
  assign(paste("rac_mod_results", parks_noecosite[i], sep ="_"),noecosite_rac_model_results_park)
  noecosite_rac_model_results <-  rbind(noecosite_rac_model_results,noecosite_rac_model_results_park)
}


#Explore model assumptions
par(mfrow = c(3, 3))
shapiro_community_noecosite <- data.frame()


for (i in 1:length(community_noecosite_park)){ #run a loop over the dataframes in the list
  qqnorm(resid(get(paste("rac_mod", parks_noecosite[i], sep="_"))), main = paste("Q-Q Plot", parks_noecosite[i], sep=" "))
  hist(resid(get(paste("rac_mod", parks_noecosite[i], sep="_"))), main = paste("Histogram", parks_noecosite[i], sep=" "))
  shapiro <- shapiro.test(resid(get(paste("rac_mod", parks_noecosite[i], sep="_"))))
  p_value = cbind(parks_noecosite[i], shapiro$p.value)
  
  shapiro_community_noecosite <- rbind(shapiro_community_noecosite, p_value)
  
}







#Code to review is below
#Need to try to look at rank curves


# #Set graph parameters ---------------------------------------------------


theme_set(theme_bw(12))

#set colors and other settings for plots

#color ecosite park by park
scale_color_park <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c('skyblue3', 'darkolivegreen3', 'turquoise4', 'tan3','lightpink3','aquamarine4'), parks_ecosite), 
    ...
  )
}

#color ecosite by park and ecosite
names_park_ecosite <- sort(unique(community_ecosite$Park_Ecosite))
scale_color_park_ecosite <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c('skyblue3', 'darkolivegreen3', 
                        'turquoise3', 'tan2','lightpink2',
                        'aquamarine3', 'lightpink3',
                        'khaki3','skyblue',
                        'green4','goldenrod',
                        'turquoise4', 'lightpink4'), names_park_ecosite), 
    ...
  )
}

 
scale_color_park_noecosite <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c("#4D004B", "#A6BDDB" ,"#7FCDBB" ,"#2D004B"), parks_noecosite), 
    ...
  )
}


#tool for creating palettes


# Richness plots ----------------------------------------------------------

#couple options for richness plots
#need to think about 
#still need to think about how to display this, don't like the colors right now.
rich_plot_ecosite1 = community_ecosite %>% 
  ggplot(aes(x = SampleYear, y = sqrt(richness))) + 
  geom_jitter(aes(color=Park), alpha =0.01)+
  geom_smooth(aes( color = Park, linetype = EcoSite),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park)+
  scale_color_park()+
  scale_linetype_manual(values =c("dotted", "solid","dashed", "solid","dashed", "solid", "dashed"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "right")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  xlab("Sample Year")+
  ylab(expression(sqrt("Richness (S)")))+
  labs(color  = "Park", linetype = "EcoSite")


#rich_plot_ecosite1

# save results
# tiff(here("results","richness","fig_rich_park_ecosite_linetype.tiff"), units="in", width=11, height=7, res=300)
# rich_plot_ecosite1
# dev.off()

rich_ecosite_anova_results <- rich_ecosite_anova_results %>% 
  mutate(fac_abb = case_when(FixedEffect=="SampleYear"~"Y",
                             FixedEffect=="EcoSite"~"E",
                             FixedEffect=="SampleYear:EcoSite"~"Y_E"))

graphLabels_rich_ecosite <- rich_ecosite_anova_results %>% 
  dplyr::select(Park, Significance, fac_abb) %>% 
  pivot_wider(id_cols = Park, names_from = fac_abb, values_from = Significance)

rich_plot_ecosite2 = community_ecosite %>% 
  ggplot(aes(x = SampleYear, y = sqrt(richness))) + 
  geom_jitter(aes(color=Park_Ecosite), alpha =0.05)+
  geom_smooth(aes(color = Park_Ecosite),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
 facet_wrap(~Park)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  scale_color_park_ecosite()+
  theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  xlab("Sample Year")+
  ylab(expression(sqrt("Richness (S)")))+
  geom_text(data = graphLabels_rich_ecosite, aes(x = 0.5, y = 5.5, label = "Year"),inherit.aes = FALSE, size=3)+
  geom_text(data = graphLabels_rich_ecosite, aes(x = 0.5, y = 5, label = Y),inherit.aes = FALSE, size=3)+
  geom_text(data = graphLabels_rich_ecosite, aes(x = 5, y = 5.5, label = "Ecosite"),inherit.aes = FALSE, size=3)+
  geom_text(data = graphLabels_rich_ecosite, aes(x = 5, y = 5, label = E),inherit.aes = FALSE, size=3)+
  geom_text(data = graphLabels_rich_ecosite, aes(x = 12, y = 5.5, label = "Year:Ecosite"),inherit.aes = FALSE, size=3)+
  geom_text(data = graphLabels_rich_ecosite, aes(x = 12, y = 5, label = Y_E),inherit.aes = FALSE, size=3)


rich_plot_ecosite2

# cvdPlot(rich_plot_ecosite2)
#save results
# tiff(here("results","richness","fig_rich_park_ecosite_color.tiff"), units="in", width=11, height=7, res=300)
# rich_plot_ecosite2
# dev.off()

graphLabels_rich_noecosite <- rich_noecosite_anova_results %>% 
  dplyr::select(Park, Significance, FixedEffect) %>% 
  pivot_wider(id_cols = Park, names_from = FixedEffect, values_from = Significance)

rich_plot_noecosite = community_noecosite %>% 
  ggplot(aes(x = SampleYear, y = sqrt(richness))) + 
  geom_jitter(aes(color=Park), alpha =0.05)+
  geom_smooth(aes( color = Park),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park)+
  scale_color_park_noecosite()+
  theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  xlab("Sample Year")+
  ylab(expression(sqrt("Richness (S)")))+
  geom_text(data = graphLabels_rich_noecosite, aes(x = 0, y = 5.5, label = "Year"),inherit.aes = FALSE, size=3)+
  geom_text(data = graphLabels_rich_noecosite, aes(x = 0, y = 5, label = SampleYear),inherit.aes = FALSE, size=3)

rich_plot_noecosite

# cvdPlot(rich_plot_noecosite)

# save results
# tiff(here("results","richness","fig_rich_park_noecosite.tiff"), units="in", width=11, height=7, res=300)
# rich_plot_noecosite
# dev.off()



#Richness by stability
rich_stability_plot_ecosite = community_structure_df %>% 
  filter(!(Park%in%c("AZRU", "CHCU", "PETR", "WACA")))  %>%  #remove sites with only 1 ecosite

  ggplot(aes(y = stability, x = sqrt(richness))) + 
  #geom_jitter(aes(color=Park_Ecosite), alpha =0.05)+
  geom_smooth(aes(color = Park_Ecosite),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  scale_color_park_ecosite()+
  theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  ylab("Stability")+
  xlab(expression(sqrt("Richness (S)")))


rich_stability_plot_ecosite

names(community_structure_df)

rich_plot_ecosystem = community_structure_df %>% 
  ggplot(aes(x = SampleYear, y = sqrt(richness))) + 
  geom_jitter(aes(color=Ecosystem), alpha =0.01)+
  geom_smooth(aes(color = Ecosystem),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
 # facet_wrap(~Park)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  #scale_color_park_ecosite()+
  scale_color_manual(values = c("blue", "forestgreen", "red", "black")) +
  theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  xlab("Sample Year")+
  ylab(expression(sqrt("Richness (S)")))

rich_plot_ecosystem

#Create visualizations
community_structure_df %>% 
  ggplot(aes(x = EcoSite, y =stability)) +
  geom_boxplot() +
  ylim(0,100)+
  facet_wrap(~Park)

stability_df %>% 
  ggplot(aes(x = Park, y =stability)) +
  geom_boxplot() +
  ylim(0,100)

community_structure_df %>% 
  ggplot(aes(y = turnover_total, x = EventYear, color = Park, linetype = EcoSite))+
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

community_structure_df %>% 
  ggplot(aes(x = EventYear, y = richness, color= Park, linetype=EcoSite)) +
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("black","red","purple", "goldenrod", "blue", "lightblue", "lavender","darkgreen", "orange","green"))+
  facet_wrap(~Park)+
  ylim(0,20)

community_structure_df %>% 
  ggplot(aes(x = EventYear, y = SimpsonEvenness, color= Park, linetype=EcoSite)) +
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("black","red","purple", "goldenrod", "blue", "lightblue", "lavender","darkgreen", "orange","green"))+
  facet_wrap(~Park)+
  ylim(0,1)

#need to figure out a better way to graph
rank_shift_df %>% 
  mutate(year =as.numeric(substr(year_pair, 6,9))) %>% 
  ggplot(aes(y = MRS, x = year, color = Park, linetype = EcoSite))+
  stat_smooth(method = "lm")+
  #scale_color_manual(values = c("black","red","purple", "goldenrod", "blue", "lightblue", "lavender","darkgreen", "orange","green"))+
  facet_wrap(~Park)

rich_turnover_plot = community_structure_df %>% 
  ggplot(aes(y = turnover_total, x = sqrt(richness))) + 
  #geom_jitter(aes(color=Park_Ecosite), alpha =0.05)+
  geom_smooth(aes(color = Park_Ecosite),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park_Ecosite)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  #scale_color_park_ecosite()+
  theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  ylab("Turnover")+
  xlab(expression(sqrt("Richness (S)")))


rich_turnover_plot

names(rate_change_interval_df)
rate_change_plot = rate_change_interval_df %>% 
  ggplot(aes(y = distance, x = interval)) + 
  #geom_jitter(aes(color=Park_Ecosite), alpha =0.05)+
  geom_smooth(aes(color = Park_Ecosite),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park_Ecosite)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  #scale_color_park_ecosite()+
  #theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  ylab("Abundance Rate Change")


rate_change_plot

rate_change_plot_pa = rate_change_interval_pa_df %>% 
  ggplot(aes(y = distance_pa, x = interval_pa)) + 
  geom_point(aes(color=Park_Ecosite), alpha= 0.05)+
  geom_smooth(aes(color = Park_Ecosite),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park_Ecosite)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  #scale_color_park_ecosite()+
  #theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  ylab("Species Rate Change")


rate_change_plot_pa

View(community_structure_df %>% 
  filter(is.na(turnover_total_pa)))

View(community_structure_df %>% 
  ungroup() %>% 
  filter(is.na(turnover_total_pa)) %>% 
  dplyr::select(Park_Ecosite, EventYear, Plot, turnover_total_pa) %>% 
  distinct())

View(community_structure_df %>% 
  ungroup() %>% 
  filter(Park=="AZRU"))

pft_long <- community_structure_df %>% 
  dplyr::select(EventYear:Quadrat, forb_NA:vine_annual) %>% 
  pivot_longer(cols = forb_NA:vine_annual , names_to = "pft", values_to = "pft_pct")

library(RColorBrewer)
n <- length(unique(pft_long$pft))
colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
col_vec2 = c("#DFC27D")
col <- sample(col_vec, n)
area <- rep(1,n)
pie(area, col = col)

library("colorBlindness")

names(pft_long)
pft_stacked <- pft_long %>%
  mutate(Park_Ecosite = paste(Park, EcoSite, sep = "")) %>% 
  group_by(Park_Ecosite, EventYear, pft) %>%
  summarise(mean_pft = mean (pft_pct)) %>% 
  ggplot(aes(x =EventYear, fill = pft, y = mean_pft))+ 
  scale_fill_manual(values = col)+
  geom_bar(stat = "identity")+
  facet_wrap(~ Park_Ecosite)

pft_stacked
#ran pft stacked try to show as lines

#run analyses on these factors
#create rank clock
# Create rank abundance curves by park ecosite ----------------------------





species_park_list <- names(species_rank_park_avg)

####make RACs graphs for each ecosite
species_rank_park_avg <- species_long %>% 
  ungroup() %>% 
  mutate(Park_Ecosite = paste(Park, EcoSite, sep = "_")) %>% 
  group_by(Park_Ecosite,EventYear,CurrentSpecies) %>% 
  summarise(Cover_plot_avg= mean(Cover_pct)) %>% 
  mutate(rank = rank(-Cover_plot_avg, ties.method = "average")) %>% 
  ungroup() 

species_rank_park_avg <-species_rank_park_avg %>% 
  split(., species_rank_park_avg$Park_Ecosite)

species_park_list <- names(species_rank_park_avg)

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
           ggplot(aes(x=rank, y=Cover_plot_avg, label = ifelse(rank<=3,colorfill,'')))+
           geom_line(color="darkgray", linewidth=1)+
           ggrepel::geom_text_repel(hjust = -.5)+
           geom_point(aes(color = colorfill), size=3)+#aes(color=colorfill), size=3)+
           facet_wrap(~EventYear)+
           #scale_color_manual(values = c("blue","lightblue","darkred", "pink","green3","orange","cornflowerblue", "darkgreen","purple", "lavender"))+
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right", axis.title.x = element_blank(), axis.title.y = element_blank()))
 
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






RAC_plot_BAND_M_manual_A <- 
       species_rank_park_avg$BAND_M %>% 
         mutate(colorfill = ifelse(CurrentSpecies%in%top_species$CurrentSpecies, CurrentSpecies, "other"))%>% 
         ggplot(aes(x=rank, y=Cover_plot_avg, label = ifelse(rank<=3,colorfill,'')))+
         geom_point(aes(color=colorfill), size=3)+
         facet_wrap(~paste(EventYear, EventPanel))+
          geom_text(hjust = 0)+
         #scale_color_manual(values = c("blue","lightblue","darkred", "pink","green3","orange","cornflowerblue", "darkgreen","purple", "lavender"))+
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", axis.title.x = element_blank(), axis.title.y = element_blank())

RAC_plot_BAND_M_manual_A


  
for (i in 1:13){ #run a loop over the dataframes in the list
  
#need to fix colors, increase number of species presented, or figure out another way to label the top 5 species or something
    aggdat <- aggregate(Cover_plot_avg ~ CurrentSpecies * EventYear * Park_Ecosite, 
                    data = subset(species_rank_park_avg[[i]], 
                                  CurrentSpecies=="Bromus tectorum"),
                                  #CurrentSpecies%in%top_species_park[[3]]$CurrentSpecies), 
                    FUN = mean)
    
    assign(paste("BT_plot", species_park_list[i], sep="_"),
       aggdat %>% 
         ggplot( aes(EventYear, Cover_plot_avg, color = CurrentSpecies)) + 
         geom_line(linewidth = 2) + 
         coord_polar() + 
         theme_bw() + 
         facet_wrap(~Park_Ecosite) +
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


names(community_structure_df)



#Is this because of panels????


#*Things to do:
#*1. rank abundance curves by event panel?
#*2. what to do about event panel anyway? discuss with megan. were plots in event panels randomized?
#*3. calculate other metrics with codyn from Hallett paper
#*4. organize this script
#*5. look into pft and functional diversity?
#*6. native and exotic
#*7. environmental variables (MuMin)

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