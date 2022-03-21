library(ggridges)
library(tidyverse)
library(MaizePal)
library(RColorBrewer)

H_perm_all<-read.table("~/Downloads/H_Perm_ALL.txt",header=T)

total_h<-select(H_perm_all, c(Her_ALL, Trait)) %>% rename(total_h=Her_ALL) #added by Sam

h<-rename(H_perm_all,MOA=Her_K1,Background=Her_K2,Genome=Her_K3) %>%
  gather(key="kin",value="h2",c(MOA,Background,Genome,Her_ALL)) %>% 
  select(kin,h2,Trait) 

hrank<-merge(h,total_h,by="Trait") %>% 
  mutate(trait_cat=ifelse(Trait %in% c("tassel_length_BLUP;Brown_2011","spike_length_BLUP;Brown_2011","branch_zone_BLUP;Brown_2011","branch_number_transformed_BLUP;Brown_2011","tassprimbranchno_BLUP;Hung_2012_1","tasslength_BLUP;Hung_2012_1"),"tassel architecture",
                          ifelse(Trait %in% c("cob_length_BLUP;Brown_2011","cob_diameter_BLUP;Brown_2011","ear_row_number_transformed_BLUP;Brown_2011","cobdiam_BLUP;Hung_2012_1","coblength_BLUP;Hung_2012_1","earrowno_BLUP;Hung_2012_1","kernelnoperrow_BLUP;Hung_2012_1","earmass_BLUP;Hung_2012_1","cobmass_BLUP;Hung_2012_1","totalkernelweight_BLUP;Hung_2012_1","weight20kernels_BLUP;Hung_2012_1","totalkernelno_BLUP;Hung_2012_1"),"yield components",
                                 ifelse(Trait %in% c("WI_rind_penetrometer_resistance_BLUE;Peiffer_2013","WI_rind_penetrometer_resistance_BLUP;Peiffer_2013","ALL_rind_penetrometer_resistance_BLUE;Peiffer_2013","ALL_rind_penetrometer_resistance_BLUP;Peiffer_2013","MO_rind_penetrometer_resistance_BLUE;Peiffer_2013","MO_rind_penetrometer_resistance_BLUP;Peiffer_2013","NY_rind_penetrometer_resistance_BLUE;Peiffer_2013","NY_rind_penetrometer_resistance_BLUP;Peiffer_2013"),"stalk strength",
                                        ifelse(Trait %in% c("southern_leaf_blight_BLUP;Bian_2014_Kump_2011","lesion_severity_BLUP;Olukolu_2014","sqrt_diseased_leaf_area1_BLUP;Poland_2011","sqrt_diseased_leaf_area2_BLUP;Poland_2011","sqrt_diseased_leaf_area3_BLUP;Poland_2011","diseased_leaf_area1_BLUP;Poland_2011","diseased_leaf_area2_BLUP;Poland_2011","diseased_leaf_area3_BLUP;Poland_2011","northern_leaf_blight_index_BLUP;Poland_2011"),"disease",
                                               ifelse(Trait %in% c("height_ratio_BLUP;Olukolu_2014","stalk_width_ratio_BLUP;Olukolu_2014","leaf_length_BLUP;Tian_2011","leaf_width_BLUP;Tian_2011","upper_leaf_angle_BLUP;Tian_2011","last_leaf_with_epicuticular_wax;Foerster_2015","ALL_ear_height_BLUE;Peiffer_2013","ALL_ear_height_BLUP;Peiffer_2013","plantheight_BLUP;Hung_2012_1","earheight_BLUP;Hung_2012_1","leaflength_BLUP;Hung_2012_1","leafwidth_BLUP;Hung_2012_1","upperleafangle_BLUP;Hung_2012_1","nodenumberbelowear_BLUP;Hung_2012_1","nodenumberaboveear_BLUP;Hung_2012_1","numbraceroots_BLUP;Hung_2012_1","plant_height_BLUP;Peiffer_2014","ear_height_BLUP;Peiffer_2014","plant_height_minus_ear_height_BLUP;Peiffer_2014","ear_height_div_plant_height_BLUP;Peiffer_2014","plant_height_div_DTA_BLUP;Peiffer_2014","leaf_angle_boxcox_transformed_BLUP;Tian_2011"),"plant architecture",
                                                      ifelse(Trait %in% c("DTA_ratio_BLUP;Olukolu_2014","ALL_DTA_BLUE;Peiffer_2013","ALL_DTA_BLUP;Peiffer_2013","DTS_BLUP;Hung_2012_1","DTA_BLUP;Hung_2012_1","asi_BLUP;Hung_2012_1","DTA_BLUP;Buckler_2009","DTS_BLUP;Buckler_2009","DTA_BLUP;Wallace_2014","ASI_BLUP;Buckler_2009","GDD_DTS_BLUP;Peiffer_2014","GDD_DTA_BLUP;Peiffer_2014","GDD_ASI_BLUP;Peiffer_2014","DTS_BLUP;Peiffer_2014","DTA_BLUP;Peiffer_2014","ASI_BLUP;Peiffer_2014","GDD_DTA_long_BLUP;Hung_2012","GDD_DTA_short_BLUP;Hung_2012","GDD_DTA_photo_resp_BLUP;Hung_2012","GDD_DTS_long_BLUP;Hung_2012","GDD_DTS_short_BLUP;Hung_2012","GDD_DTS_photo_resp_BLUP;Hung_2012","GDD_DTA_NC06_BLUP;Hung_2012","GDD_DTS_NC06_BLUP;Hung_2012","GDD_DTA_MO06_BLUP;Hung_2012","GDD_DTS_MO06_BLUP;Hung_2012","GDD_DTA_NY06_BLUP;Hung_2012","GDD_DTS_NY06_BLUP;Hung_2012","GDD_DTA_IL06_BLUP;Hung_2012","GDD_DTS_IL06_BLUP;Hung_2012","GDD_DTA_FL06_BLUP;Hung_2012","GDD_DTS_FL06_BLUP;Hung_2012","GDD_DTA_PR06_BLUP;Hung_2012","GDD_DTS_PR06_BLUP;Hung_2012","GDD_DTA_NC07_BLUP;Hung_2012","GDD_DTS_NC07_BLUP;Hung_2012","GDD_DTA_MO07_BLUP;Hung_2012","GDD_DTS_MO07_BLUP;Hung_2012","GDD_DTA_NY07_BLUP;Hung_2012","GDD_DTS_NY07_BLUP;Hung_2012","GDD_DTA_IL07_BLUP;Hung_2012","GDD_DTS_IL07_BLUP;Hung_2012","GDD_DTA_FL07_BLUP;Hung_2012","GDD_DTS_FL07_BLUP;Hung_2012"),"flowering time",
                                                             ifelse(Trait %in% c("glsblup;Benson_2015","flecking_ls_mean;Olukolu_2016","flecking_poisson_transformation_BLUP;Olukolu_2016","starch_kernel_BLUP;Cook_2012","protein_kernel_BLUP;Cook_2012","oil_kernel_BLUP;Cook_2012"),"miscellaneous",
                                                                    ifelse(Trait %in% c("alphaT_transformed_BLUE;Diepenbrock_2017","deltaT_transformed_BLUE;Diepenbrock_2017","gammaT_transformed_BLUE;Diepenbrock_2017","TotalT_transformed_BLUE;Diepenbrock_2017","alphaT3_transformed_BLUE;Diepenbrock_2017","deltaT3_transformed_BLUE;Diepenbrock_2017","gammaT3_transformed_BLUE;Diepenbrock_2017","TotalT3_transformed_BLUE;Diepenbrock_2017","TotalT_plus_T3_transformed_BLUE;Diepenbrock_2017","plastochromanol8_transformed_BLUE;Diepenbrock_2017","alphaT_div_gammaT_transformed_BLUE;Diepenbrock_2017","alphaT3_div_gammaT3_transformed_BLUE;Diepenbrock_2017","gammaT_div_gammaTPlusalphaT_transformed_BLUE;Diepenbrock_2017","gammaT3_div_gammaT3PlusalphaT3_transformed_BLUE;Diepenbrock_2017","deltaT_div_gammaTPlusalphaT_transformed_BLUE;Diepenbrock_2017","deltaT_div_alphaT_transformed_BLUE;Diepenbrock_2017","deltaT_div_gammaT_transformed_BLUE;Diepenbrock_2017","deltaT3_div_gammaT3PlusalphaT3_transformed_BLUE;Diepenbrock_2017","deltaT3_div_alphaT3_transformed_BLUE;Diepenbrock_2017","deltaT3_div_gammaT3_transformed_BLUE;Diepenbrock_2017","TotalT_div_totalT3_transformed_BLUE;Diepenbrock_2017") ,"vitamin E",
                                                                           ifelse(Trait %in% c("chlorophylla_BLUP;Wallace_2014","chlorophyllb_BLUP;Wallace_2014","malate_BLUP;Wallace_2014","fumarate_BLUP;Wallace_2014","fumarate2_BLUP;Wallace_2014","glutamate_BLUP;Wallace_2014","aminoacids_BLUP;Wallace_2014","protein_BLUP;Wallace_2014","nitrate_BLUP;Wallace_2014","starch_BLUP;Wallace_2014","sucrose_BLUP;Wallace_2014","glucose_BLUP;Wallace_2014","fructose_BLUP;Wallace_2014","principal_component_metabolies1;Wallace_2014","principal_component_metabolies2;Wallace_2014"),"metabolites","ERROR")))))))))) %>%
  filter(kin!="Her_ALL")

trait_cat_colors<-data.frame(trait_cat=levels(factor(hrank$trait_cat)),trait_cols=c(maize_pal("OaxacaGreen")[c(1,2,4)],maize_pal("MaizMorado")[1:3],maize_pal("Painted")[1:3]))
traits<-group_by(hrank,Trait,trait_cat) %>% summarize(meanh=mean(total_h)) %>% select(Trait,meanh,trait_cat)
trait_colors<-merge(trait_cat_colors,traits,by="trait_cat")

ggplot(hrank, aes(y=reorder(Trait,total_h),x=h2 )) + 
  stat_density_ridges(aes(fill=kin),alpha=0.5,size=0,from=0,to=1,bandwidth = 0.01)+
  ylab("Trait")+ xlab(expression(h^2)) +
  scale_fill_manual(name = "Kinship Matrix",values=maize_pal("HighlandMAGIC")[c(1,3,6)]) +
  scale_y_discrete(labels=rep("•",143))+
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(),
        axis.text.y=element_text(size = rel(3),margin = margin(-30,-30,-30,-30),color=trait_colors[order(trait_colors$meanh),]$trait_cols),
        axis.title.y=element_text(size = rel(1.5),margin = margin(40,40,40,40)),
        axis.title.x=element_text(size = rel(1.5)),
        axis.text.x=element_text(size=rel(1.5)),
        legend.title=element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1.25)),
        legend.position = c(0.75, 0.2)) +
  geom_point(data=trait_colors,aes(x=meanh,y=Trait,color=trait_cat,shape=NA))+ #this gets the legend to match the colors
  scale_color_manual(name="Trait",values=c(maize_pal("OaxacaGreen")[c(1,2,4)],maize_pal("MaizMorado")[1:3],maize_pal("Painted")[1:3]))

##Sam trying things
hrank<-inner_join(x=hrank,y=trait_cat_colors, by="trait_cat")
traits<-group_by(hrank,Trait,trait_cat) %>% 
  arrange(desc(total_h)) %>%
  summarize(meanh=mean(total_h), maxh=max(total_h)) %>% 
  select(Trait,meanh,maxh,trait_cat)

ggplot(hrank, aes(y=reorder(Trait,total_h),x=h2 )) + 
  stat_density_ridges(aes(fill=kin),alpha=0.5,size=0,from=0,to=1,bandwidth = 0.01)+
  ylab("Trait")+ xlab(expression(h^2)) +
  scale_fill_manual(name = "Kinship Matrix",values=brewer.pal(3,"Dark2")[c(2,3,1)]) +
  scale_y_discrete(labels=rep("•",143))+
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(),
        #axis.text.y=element_text(size = rel(3),margin = margin(-30,-30,-30,-30)),
        axis.text.y=element_blank(),
        axis.title.y=element_text(size = rel(1.5),margin = margin(40,40,40,40)),
        axis.title.x=element_text(size = rel(1.5)),
        axis.text.x=element_text(size=rel(1.5)),
        legend.title=element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1.25)))+#,
        #legend.position = c(0.75, 0.2)) +
  geom_point(data=traits,aes(x=-0.05,y=Trait,color=trait_cat))+ 
  scale_color_manual(name="Trait",
                     values=c(maize_pal("OaxacaGreen")[c(1,2)],
                              maize_pal("MaizAzul")[3],
                              maize_pal("MaizMorado")[c(1,3,4)],
                              maize_pal("Painted")[1:3]))
ggsave("Test_Results_Full/Allh2.Ridgeplot.traitcatlabeled.png",dpi = 300, device = "png",height = 8)
