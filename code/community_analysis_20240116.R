#analysis following last update to data_wrangle2_codyn
#Install packages
packages <- c("tidyverse", "here","lme4","LMERConvenienceFunctions", "lmerTest", "emmeans", "multcomp", "MuMIn",'nortest', "vegan", "codyn", "data.table", "janitor", "indicspecies") #list of packages to install
lapply(packages, library, character.only = T) #load packages

load(file=here("results","upload4"))
load(file=here("results","species_wide_df"))


species_list <- species_long %>% 
  dplyr::select(CurrentSpecies, Nativity, Lifeform, Duration) %>% 
  distinct()

species_na <- species_list %>% 
  filter_at(vars(Nativity, Lifeform, Duration), any_vars(is.na(.)))

# write.csv(species_na, file = "species_list_na.csv")

# Richness models ---------------------------------------------------------

community_structure_df <- community_structure_df %>% 
  mutate(Year_f = as.factor(SampleYear),
         Rep = paste(Park, EcoSite, EventPanel, Plot, Transect)) 


names(community_structure_df)

Shannon_mod <- lmer((Shannon)~SampleYear*Ecosystem +  (1|Year_f:Rep) , data= community_structure_df)
anova(Shannon_mod)
qqnorm(resid(Shannon_mod))
hist(resid(Shannon_mod))
plot(resid(Shannon_mod))



shannon_plot = community_structure_df %>% 
  ggplot(aes(x = SampleYear, y = Shannon)) + 
  #geom_jitter(aes(color=Park), alpha =0.01)+
  geom_smooth(aes( color = Ecosystem),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  #facet_wrap(~Park)+
  #scale_color_park()+
  #scale_linetype_manual(values =c("dotted", "solid","dashed", "solid","dashed", "solid", "dashed"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "right")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  xlab("Sample Year")+
  ylab("Shannon")+
  labs(color  = "Park", linetype = "EcoSite")

shannon_plot


SimpsonEvenness_mod <- lmer(sqrt(SimpsonEvenness)~SampleYear*Ecosystem +  (1|Year_f:Rep) , data= community_structure_df)
anova(SimpsonEvenness_mod)
qqnorm(resid(SimpsonEvenness_mod))
hist(resid(SimpsonEvenness_mod))
plot(resid(SimpsonEvenness_mod))


SimpsonEvenness_plot = community_structure_df %>% 
  ggplot(aes(x = SampleYear, y = SimpsonEvenness)) + 
  #geom_jitter(aes(color=Park), alpha =0.01)+
  geom_smooth(aes( color = Ecosystem),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  #facet_wrap(~Park)+
  #scale_color_park()+
  #scale_linetype_manual(values =c("dotted", "solid","dashed", "solid","dashed", "solid", "dashed"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "right")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  xlab("Sample Year")+
  ylab(expression(sqrt("SimpsonEvenness")))+
  labs(color  = "Park", linetype = "EcoSite")

SimpsonEvenness_plot


InverseSimpson_mod <- lmer(log1p(InverseSimpson)~SampleYear*Ecosystem +  (1|Year_f:Rep) , data= community_structure_df)
anova(InverseSimpson_mod)
qqnorm(resid(InverseSimpson_mod))
hist(resid(InverseSimpson_mod))
plot(resid(InverseSimpson_mod))



InverseSimpson_plot = community_structure_df %>% 
  ggplot(aes(x = richness, y = InverseSimpson)) + 
  #geom_jitter(aes(color=Park), alpha =0.01)+
  geom_smooth(aes( color = Ecosystem),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  #facet_wrap(~Park)+
  #scale_color_park()+
  #scale_linetype_manual(values =c("dotted", "solid","dashed", "solid","dashed", "solid", "dashed"))+
  #scale_x_continuous(name = "Experiment Year (years)", breaks = seq(1:length(unique(bugs_div$Year))))+
  #scale_y_continuous(name = "Log Macroinvertebrate Density \n (individuals m-2)", breaks = c(0,10))+
  theme(legend.position = "right")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  xlab("Richness")+
  ylab(expression(log1p("InverseSimpson")))+
  labs(color  = "Park", linetype = "EcoSite")

InverseSimpson_plot






###NMDS in Figure 5

# Create nmds by park ecosite ----------------------------
#currently running with relative abundance
species_wide_grass <- species_wide %>% 
  left_join(park_ecosystem) %>% 
  filter(Ecosystem =="Grassland and shrubland")

names(species_wide_grass)
#remove species columns
species_grass_env <- species_wide_grass %>% 
  dplyr::select(EventYear:Replicate, EcoSite:EcoSiteName) 

#select species columns and standardize data
#d
# species_rel_abund_total <- species_wide %>% 
#   ungroup()%>% 
#   dplyr::select(-(Park_Ecosite),-(Park:Replicate), -year_factor) %>% 
#   decostand(., method="total") %>% 
#   cbind(species_env,.)
# 
# species_rel_abund_hellinger <- species_wide %>% 
#   ungroup()%>% 
#   dplyr::select(-(Park_Ecosite),-(Park:Replicate), -year_factor) %>% 
#   decostand(., method="hellinger") %>% 
#   cbind(species_env,.)



set.seed(86753092)


names(species_only)
species_only <-  species_wide_grass %>% 
  dplyr::select(-c(EventYear:Replicate, EcoSite:EcoSiteName))%>% 
         vegdist(., method= "bray")

mds1 <-  vegan::metaMDS(species_only, autotransform=FALSE, shrink=FALSE, k =2)#, trymax = 50, k=3))


par(mfrow = c(1, 1))

goodness(mds1)
stressplot(mds1)


##test for differences in centroid means
# adonis2(pefo_wide[,8:188]~year_factor*EcoSite, pefo_wide)
#adonis - permanova

adonis2(species_only~EventYear, species_env_AZRU)
#test for differences dispersion
# dist<-vegdist(wide_dat[,4:61])
# betadisp<-betadisper(dist, wide_dat$Treatment,type="centroid")
# betadisp
# permutest(betadisp)

# dist<-vegdist(get(paste("species_only",names_parks[1], sep ="_")))
dist <- get(paste("species_only", names_parks[1], sep ="_"))
betadisp<-betadisper(dist, species_env_AZRU$EventYear,type="centroid")
betadisp
permutest(betadisp)


#isolate NMDS scores to great figures
scores <- data.frame(scores(get(paste("mds",names_parks[1], sep ="_")), display="sites"))  # Extracts NMDS scores for year "i" #
scores2<- cbind(get(paste("species_env",names_parks[1], sep ="_")), scores) # binds the NMDS scores of year i to all years previously run

toplot<-scores2 %>% 
  group_by(year_factor)

centroid <- toplot %>% 
  summarize(NMDS1= mean(NMDS1), NMDS2 = mean(NMDS2))

##Make NMDS graphs for figure 5data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAkCAYAAAD7PHgWAAADFklEQVR42u3Y20+SYRwHcGvddOn/0VX32tE5y2WrvLCVtjbbLGuZm4Iu8VBoaSwttRTUFwQ5iUoK4llzLkvJA6h4yNzMY8KLTrxq354X2sIDIqhJG7/tO+CG98PzPu/zex4CAvzlL3/5a8diiT8G+iwuTdpalKnsAPOaLtWd8jlgalUzRpY2wO8eB/M+VaYL9ikgW9IE/dw6ZENmvOuaQApBJoubI3wCx+FwjrPEOvT9cAClg2aUdE+DLSFIUWP0v0G0t59giesDXYUtadkErBpcQWnPDFjiFiRS2vsHBkmgGrITKA1NAuckCrVIJrfRVdLkbdBvAUoGViD4NIskgn8s1CQdCPBReQNtWl6HkUx4A8nwoiNDizYMLtgwMG8jc81mn2/9TAiKGTnmc9f02jYgk7LeOfIjWkG+m7NvYLxATQ/MraHpmxXaSSs0E1Y0jFtRP07jg4mGeoxG3SiNmlELVCMWVBstUBosUBjMkA+bdwSKv66g/PMCQbbhoUCdsy9gXKn6l3HJBp0bYO1uwCEHsMoJWElCfVlEsrQDcfy6V14D772vhXF5A41T7oE1HgIr9T9B9S2BJetCbEltgVfAu8U124ETW4Bjf4EqF0CpC6DIjlxGiqIH5FqFkTz5SY+AMYXVHgOZUXQGytwAhf2OPFX1IrpI+cQj4K0ChQugdVegkgAVBCjfA1BEcFT/MthkFG/nKx94BIx6LfsDXD00oJDc4kRJJ27my94yXcgjYCRPeqjACvIkx1MtuMGT5Xv1kFzLk+wONLkBDjsBBzcDmbUwlq8FuUae18tMxEvxZuCkKyDtAI7sDSjonUdMsRpXcyXcfS3U4TkiskivHiiwpGcWUQUqhL8Q7b/VXeYKiy5xhfAmGoLfCizuniHzTY4wrjDpSPeEoc8r7BsH5378pnMKV3IlCH1GxR/5pvViVvkmIK/VhLBsEUIyy+/4xK76fEYZur+v2YHZjQaEZFXgXIbgus+cSc6mC+ztLl2tJ7AyEv4Fnzo0BXH4RWcIMpjDRxCn7LRvHo7tbQvH/H9h+Ot/r98AkmxTZNPwmwAAAABJRU5ErkJggg==
#graph for controls
ggplot(subset(toplot), aes(x=NMDS1, y=NMDS2, color=year_factor))+
  geom_point(size=2)+
  stat_ellipse(level = 0.50) +
  geom_point(data = centroid, size = 5, shape=21, aes(fill=year_factor), color = "black")+
  scale_color_manual(values = c("red", "black","green"))+
  scale_fill_manual(values = c("red", "black","green"))+
  #scale_color_manual(name="", values=c("black","red"))+
  #scale_x_continuous(limits=c(-1,1.2))+
  #scale_y_continuous(limits=c(-1.1,1.3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top")

#BAND
set.seed(86753092)

assign(paste("species_env",names_parks[3], sep="_"), 
       species[[3]] %>% dplyr::select(Park_Ecosite:year_factor))
assign(paste("species_only",names_parks[3], sep="_"), 
       species[[3]] %>% 
         dplyr::select(-(Park_Ecosite:year_factor)) %>% 
         vegdist(., method= "bray")
)
assign(paste("mds", names_parks[3], sep="_"),
       vegan::metaMDS(get(paste("species_only",names_parks[3], sep ="_")), autotransform=FALSE, shrink=FALSE, k =2))#, trymax = 50, k=3))


par(mfrow = c(1, 1))

goodness(get(paste("mds", names_parks[3], sep="_")))
stressplot(get(paste("mds", names_parks[3], sep="_")))

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
par(mfrow = c(3, 3))
shapiro_community_noecosite <- data.frame()


for (i in 1:length(community_noecosite_park)){ #run a loop over the dataframes in the list
  qqnorm(resid(get(paste("rich_mod", parks_noecosite[i], sep="_"))), main = paste("Q-Q Plot", parks_noecosite[i], sep=" "))
  hist(resid(get(paste("rich_mod", parks_noecosite[i], sep="_"))), main = paste("Histogram", parks_noecosite[i], sep=" "))
  shapiro <- shapiro.test(resid(get(paste("rich_mod", parks_noecosite[i], sep="_"))))
  p_value = cbind(parks_noecosite[i], shapiro$p.value)

  shapiro_community_noecosite <- rbind(shapiro_community_noecosite, p_value)

}

shapiro_community_noecosite



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

names_park_ecosite_all <- unique(community_structure_df$Park_Ecosite)
scale_color_park_ecosite_all <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c('skyblue3', 'darkolivegreen3', 
                        'turquoise3', 'tan2','lightpink2',
                        'aquamarine3', 'lightpink3',
                        'khaki3','skyblue',
                        'green4','goldenrod',
                        'turquoise4', 'lightpink4',"#4D004B", "#A6BDDB" ,"#7FCDBB" ,"#2D004B"), names_park_ecosite_all), 
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

#Stopped here, need to figure out how to how each of these different matrix are calculated, fix data_wrangle to rearrange these several dataframes
#Turnover

turnover_all <- turnover_df %>% 
  left_join(appearance_df) %>% 
  left_join(disappearance_df)

turnover_total_plot_ecosite2 = turnover_all %>% 
  ggplot(aes(x = EventYear)) + 
  #geom_jitter(aes(color=Park_Ecosite), alpha =0.05)+
  geom_smooth(aes(color = Park_Ecosite,y = (turnover_total)),method = "glm", se=F, linewidth=1.25 )+
  #geom_smooth(aes(color = Park_Ecosite,y = (appearance+turnover_total)),method = "glm", se=F, linewidth=1.25 , linetype="dashed")+
  #geom_smooth(aes(color = Park_Ecosite,y = (disappearance)),method = "glm", se=F, linewidth=1.25 , linetype="dotted")+
  
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  #scale_color_park_ecosite()+
  theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  xlab("Sample Year")
turnover_total_plot_ecosite2





rich <- community_structure_df %>% 
  dplyr::select(EventYear, Replicate, richness) 


turnover_rich_plot = turnover_df %>% 
  left_join(rich) %>% 
  ggplot(aes(x = sqrt(richness), y = (turnover_total))) + 
  #geom_jitter(aes(color=Park_Ecosite), alpha =0.05)+
  geom_smooth(aes(color = Park_Ecosite),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  #scale_color_park_ecosite_all()+
  theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  xlab("Richness")+
  ylim(0,1)
turnover_rich_plot






turnover_rich_plot = turnover_df %>% 
  left_join(rich) %>% 
  left_join(park_ecosystem) %>% 
  ggplot(aes(x = sqrt(richness), y = (turnover_total))) + 
  #geom_jitter(aes(color=Park_Ecosite), alpha =0.05)+
  geom_smooth(method = "glm", se=T, linewidth=1.25, aes(linetype=Ecosystem) )+
 # geom_point(aes(color=Park_Ecosite))+
  #scale_color_manual(values = c("blue", "orange"))+
  #facet_wrap(~Ecosystem)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  #scale_color_park_ecosite_all()+
  theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  xlab("Richness")+
  ylim(0,1)
turnover_rich_plot






stability_rich_plot = stability_df %>% 
left_join(rich) %>% 
  ggplot(aes(x = richness, y = (stability))) + 
  #geom_jitter(aes(color=Park_Ecosite), alpha =0.05)+
  geom_smooth(aes(color = Park_Ecosite),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  scale_color_park_ecosite_all()+
  theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  xlab("Richness")
stability_rich_plot

#Rate change interval

rate_change_interval_pa_plot = rate_change_interval_pa_df %>%
  mutate(EventPanel = str_remove(Replicate,"\\w{4}_\\w{1}_")) %>% 
  separate(EventPanel, into = c("EventPanel", "Plot", "Transect", "Quadrat"), sep = "_") %>% 
  # filter(EventPanel=="C") %>% 
  #filter(Park == "CHCU") %>%  
  ggplot(aes(x = interval, y = (distance))) + 
  #geom_jitter(aes(color=Park_Ecosite), alpha =0.05)+
  geom_smooth(aes(color = Park_Ecosite),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  #scale_color_park_ecosite()+
  theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  ylab("PA distance")
rate_change_interval_pa_plot

rate_change_interval_plot = rate_change_interval_df %>%
  mutate(EventPanel = str_remove(Replicate,"\\w{4}_\\w{1}_")) %>% 
  separate(EventPanel, into = c("EventPanel", "Plot", "Transect", "Quadrat"), sep = "_") %>% 
  # filter(EventPanel=="C") %>% 
  #filter(Park == "CHCU") %>% 
  ggplot(aes(x = interval, y = (distance))) +
  #geom_jitter(aes(color=Park_Ecosite), alpha =0.05)+
  geom_smooth(aes(color = Park),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  facet_wrap(~Park)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  #scale_color_park_ecosite()+
  theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  ylab("Species distance")
rate_change_interval_plot








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


names(community_structure_df)

ecosystem <- community_structure_df %>% 
  dplyr::select(EventYear:Quadrat, Ecosystem) %>% 
  distinct()

turnover_rich_plot = turnover_df %>% 
  left_join(rich) %>%
  left_join(ecosystem) %>% 
  ggplot(aes(x = richness, y = (turnover_total))) + 
  #geom_jitter(aes(color=Park_Ecosite), alpha =0.05)+
  geom_smooth(aes(color = Ecosystem),method = "glm", se=F, linewidth=1.25 )+
  #scale_color_manual(values = c("blue", "orange"))+
  #facet_wrap(~Park)+
  #geom_smooth(method="glm", linewidth=1.25, color="darkred")+
  #scale_color_park_ecosite_all()+
  theme(legend.position = "top")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.ticks.length = unit(.3, "cm"), text = element_text(size=15),
        axis.text = element_text(size=12))+
  xlab("Richness")
turnover_rich_plot

###NMDS in Figure 5

# Create nmds by park ecosite ----------------------------
#currently running with relative abundance
species_rank <- species_long %>%
  group_by(Park, EcoSite, Replicate, EventYear) %>% #replicate EventPanel, Plot, Transect, Quadrat
  mutate(rank = rank(-Cover_pct, ties.method = "average")) #ties are averaged, add rank to the file


#create wide format
species_wide<-species_rank%>%
  dplyr::select(-rank)%>%
  dplyr::select(Park_Ecosite, Park, EcoSite:Quadrat, Replicate, CurrentSpecies, Cover_pct)%>%
  spread(CurrentSpecies, Cover_pct, fill=0) %>% 
  mutate(year_factor = as.factor(EventYear)) %>% 
  ungroup()

#remove species columns
species_env <- species_wide %>% 
  dplyr::select(Park_Ecosite,Park:Replicate, year_factor) 

#select species columns and standardize data
#d
# species_rel_abund_total <- species_wide %>% 
#   ungroup()%>% 
#   dplyr::select(-(Park_Ecosite),-(Park:Replicate), -year_factor) %>% 
#   decostand(., method="total") %>% 
#   cbind(species_env,.)
# 
# species_rel_abund_hellinger <- species_wide %>% 
#   ungroup()%>% 
#   dplyr::select(-(Park_Ecosite),-(Park:Replicate), -year_factor) %>% 
#   decostand(., method="hellinger") %>% 
#   cbind(species_env,.)

species_adund <- species_wide %>% 
  ungroup()%>% 
  dplyr::select(-(Park_Ecosite),-(Park:Replicate), -year_factor) %>% 
  cbind(species_env,.)


#create list of dataframes to be used in nmds
species <- species_adund %>% #using regular abundance data for now
  split(., species_adund$Park)

names_parks <- names(species)

#AZRU
set.seed(86753092)

assign(paste("species_env",names_parks[1], sep="_"), 
       species[[1]] %>% dplyr::select(Park_Ecosite:year_factor))
assign(paste("species_only",names_parks[1], sep="_"), 
       species[[1]] %>% 
         dplyr::select(-(Park_Ecosite:year_factor)) %>% 
         vegdist(., method= "bray")
         )
assign(paste("mds", names_parks[1], sep="_"),
    vegan::metaMDS(get(paste("species_only",names_parks[1], sep ="_")), autotransform=FALSE, shrink=FALSE, k =2))#, trymax = 50, k=3))


par(mfrow = c(1, 1))

goodness(get(paste("mds", names_parks[1], sep="_")))
stressplot(get(paste("mds", names_parks[1], sep="_")))


##test for differences in centroid means
# adonis2(pefo_wide[,8:188]~year_factor*EcoSite, pefo_wide)
#adonis - permanova

adonis2(get(paste("species_only", names_parks[1], sep ="_"))~EventYear, species_env_AZRU)
#test for differences dispersion
# dist<-vegdist(wide_dat[,4:61])
# betadisp<-betadisper(dist, wide_dat$Treatment,type="centroid")
# betadisp
# permutest(betadisp)

# dist<-vegdist(get(paste("species_only",names_parks[1], sep ="_")))
dist <- get(paste("species_only", names_parks[1], sep ="_"))
betadisp<-betadisper(dist, species_env_AZRU$EventYear,type="centroid")
betadisp
permutest(betadisp)


#isolate NMDS scores to great figures
scores <- data.frame(scores(get(paste("mds",names_parks[1], sep ="_")), display="sites"))  # Extracts NMDS scores for year "i" #
scores2<- cbind(get(paste("species_env",names_parks[1], sep ="_")), scores) # binds the NMDS scores of year i to all years previously run

toplot<-scores2 %>% 
  group_by(year_factor)

centroid <- toplot %>% 
  summarize(NMDS1= mean(NMDS1), NMDS2 = mean(NMDS2))

##Make NMDS graphs for figure 5data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAkCAYAAAD7PHgWAAADFklEQVR42u3Y20+SYRwHcGvddOn/0VX32tE5y2WrvLCVtjbbLGuZm4Iu8VBoaSwttRTUFwQ5iUoK4llzLkvJA6h4yNzMY8KLTrxq354X2sIDIqhJG7/tO+CG98PzPu/zex4CAvzlL3/5a8diiT8G+iwuTdpalKnsAPOaLtWd8jlgalUzRpY2wO8eB/M+VaYL9ikgW9IE/dw6ZENmvOuaQApBJoubI3wCx+FwjrPEOvT9cAClg2aUdE+DLSFIUWP0v0G0t59giesDXYUtadkErBpcQWnPDFjiFiRS2vsHBkmgGrITKA1NAuckCrVIJrfRVdLkbdBvAUoGViD4NIskgn8s1CQdCPBReQNtWl6HkUx4A8nwoiNDizYMLtgwMG8jc81mn2/9TAiKGTnmc9f02jYgk7LeOfIjWkG+m7NvYLxATQ/MraHpmxXaSSs0E1Y0jFtRP07jg4mGeoxG3SiNmlELVCMWVBstUBosUBjMkA+bdwSKv66g/PMCQbbhoUCdsy9gXKn6l3HJBp0bYO1uwCEHsMoJWElCfVlEsrQDcfy6V14D772vhXF5A41T7oE1HgIr9T9B9S2BJetCbEltgVfAu8U124ETW4Bjf4EqF0CpC6DIjlxGiqIH5FqFkTz5SY+AMYXVHgOZUXQGytwAhf2OPFX1IrpI+cQj4K0ChQugdVegkgAVBCjfA1BEcFT/MthkFG/nKx94BIx6LfsDXD00oJDc4kRJJ27my94yXcgjYCRPeqjACvIkx1MtuMGT5Xv1kFzLk+wONLkBDjsBBzcDmbUwlq8FuUae18tMxEvxZuCkKyDtAI7sDSjonUdMsRpXcyXcfS3U4TkiskivHiiwpGcWUQUqhL8Q7b/VXeYKiy5xhfAmGoLfCizuniHzTY4wrjDpSPeEoc8r7BsH5378pnMKV3IlCH1GxR/5pvViVvkmIK/VhLBsEUIyy+/4xK76fEYZur+v2YHZjQaEZFXgXIbgus+cSc6mC+ztLl2tJ7AyEv4Fnzo0BXH4RWcIMpjDRxCn7LRvHo7tbQvH/H9h+Ot/r98AkmxTZNPwmwAAAABJRU5ErkJggg==
#graph for controls
ggplot(subset(toplot), aes(x=NMDS1, y=NMDS2, color=year_factor))+
  geom_point(size=2)+
  stat_ellipse(level = 0.50) +
  geom_point(data = centroid, size = 5, shape=21, aes(fill=year_factor), color = "black")+
  scale_color_manual(values = c("red", "black","green"))+
  scale_fill_manual(values = c("red", "black","green"))+
  #scale_color_manual(name="", values=c("black","red"))+
  #scale_x_continuous(limits=c(-1,1.2))+
  #scale_y_continuous(limits=c(-1.1,1.3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top")

#BAND
set.seed(86753092)

assign(paste("species_env",names_parks[3], sep="_"), 
       species[[3]] %>% dplyr::select(Park_Ecosite:year_factor))
assign(paste("species_only",names_parks[3], sep="_"), 
       species[[3]] %>% 
         dplyr::select(-(Park_Ecosite:year_factor)) %>% 
         vegdist(., method= "bray")
)
assign(paste("mds", names_parks[3], sep="_"),
       vegan::metaMDS(get(paste("species_only",names_parks[3], sep ="_")), autotransform=FALSE, shrink=FALSE, k =2))#, trymax = 50, k=3))


par(mfrow = c(1, 1))

goodness(get(paste("mds", names_parks[3], sep="_")))
stressplot(get(paste("mds", names_parks[3], sep="_")))


##test for differences in centroid means
# adonis2(pefo_wide[,8:188]~year_factor*EcoSite, pefo_wide)
#adonis - permanova

adonis2(get(paste("species_only", names_parks[3], sep ="_"))~EventYear, get(paste("species_env", names_parks[3], sep ="_")))
#test for differences dispersion
# dist<-vegdist(wide_dat[,4:61])
# betadisp<-betadisper(dist, wide_dat$Treatment,type="centroid")
# betadisp
# permutest(betadisp)

# dist<-vegdist(get(paste("species_only",names_parks[3], sep ="_")))
dist <- get(paste("species_only", names_parks[3], sep ="_"))
betadisp<-betadisper(dist, species_env_AZRU$EventYear,type="centroid")
betadisp
permutest(betadisp)


#isolate NMDS scores to great figures
scores <- data.frame(scores(get(paste("mds",names_parks[3], sep ="_")), display="sites"))  # Extracts NMDS scores for year "i" #
scores2<- cbind(get(paste("species_env",names_parks[3], sep ="_")), scores) # binds the NMDS scores of year i to all years previously run

toplot<-scores2 %>% 
  group_by(year_factor)

centroid <- toplot %>% 
  summarize(NMDS1= mean(NMDS1), NMDS2 = mean(NMDS2))
#use color brewer thing to redo plot below
#reorder years so they are in order

#start running the other nmds

##Make NMDS graphs for figure 5data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAkCAYAAAD7PHgWAAADFklEQVR42u3Y20+SYRwHcGvddOn/0VX32tE5y2WrvLCVtjbbLGuZm4Iu8VBoaSwttRTUFwQ5iUoK4llzLkvJA6h4yNzMY8KLTrxq354X2sIDIqhJG7/tO+CG98PzPu/zex4CAvzlL3/5a8diiT8G+iwuTdpalKnsAPOaLtWd8jlgalUzRpY2wO8eB/M+VaYL9ikgW9IE/dw6ZENmvOuaQApBJoubI3wCx+FwjrPEOvT9cAClg2aUdE+DLSFIUWP0v0G0t59giesDXYUtadkErBpcQWnPDFjiFiRS2vsHBkmgGrITKA1NAuckCrVIJrfRVdLkbdBvAUoGViD4NIskgn8s1CQdCPBReQNtWl6HkUx4A8nwoiNDizYMLtgwMG8jc81mn2/9TAiKGTnmc9f02jYgk7LeOfIjWkG+m7NvYLxATQ/MraHpmxXaSSs0E1Y0jFtRP07jg4mGeoxG3SiNmlELVCMWVBstUBosUBjMkA+bdwSKv66g/PMCQbbhoUCdsy9gXKn6l3HJBp0bYO1uwCEHsMoJWElCfVlEsrQDcfy6V14D772vhXF5A41T7oE1HgIr9T9B9S2BJetCbEltgVfAu8U124ETW4Bjf4EqF0CpC6DIjlxGiqIH5FqFkTz5SY+AMYXVHgOZUXQGytwAhf2OPFX1IrpI+cQj4K0ChQugdVegkgAVBCjfA1BEcFT/MthkFG/nKx94BIx6LfsDXD00oJDc4kRJJ27my94yXcgjYCRPeqjACvIkx1MtuMGT5Xv1kFzLk+wONLkBDjsBBzcDmbUwlq8FuUae18tMxEvxZuCkKyDtAI7sDSjonUdMsRpXcyXcfS3U4TkiskivHiiwpGcWUQUqhL8Q7b/VXeYKiy5xhfAmGoLfCizuniHzTY4wrjDpSPeEoc8r7BsH5378pnMKV3IlCH1GxR/5pvViVvkmIK/VhLBsEUIyy+/4xK76fEYZur+v2YHZjQaEZFXgXIbgus+cSc6mC+ztLl2tJ7AyEv4Fnzo0BXH4RWcIMpjDRxCn7LRvHo7tbQvH/H9h+Ot/r98AkmxTZNPwmwAAAABJRU5ErkJggg==
#graph for controls
ggplot(subset(toplot), aes(x=NMDS1, y=NMDS2, color=year_factor))+
  geom_point(size=2)+
  stat_ellipse(level = 0.50) +
  geom_point(data = centroid, size = 5, shape=21, aes(fill=year_factor), color = "black")+
  scale_color_manual(values = viridis::viridis(10))+
  scale_fill_manual(values =  viridis::viridis(10))+
  #scale_color_manual(name="", values=c("black","red"))+
  #scale_x_continuous(limits=c(-1,1.2))+
  #scale_y_continuous(limits=c(-1.1,1.3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top")




#older stuff
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