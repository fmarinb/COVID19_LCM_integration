## VarSelLCM

## Package loading
library(tidyverse)
library(VarSelLCM)
library(forcats)
library(ggpubr)


## Data loading ###

div<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/diversity.csv",
              header=T, sep=',')[,-c(1, c(3:8))]
clon<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/clonality.csv",
               header=T, sep=',')[,-c(1, c(3:8))]
vus<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/v_usage_prop_w.csv",
              header=T, sep=',')[,-c(1, c(3:8))]%>%replace(is.na(.), 0)

jus<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/j_usage_prop_w.csv",
              header=T, sep=',')[,-c(1, c(3:8))]
cdr3<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/t_cdr3_table3_freq_mira.csv",
               header=T, sep=',')[,-c(1, c(3:8))]
genot<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/genotypes.csv",
                header=T, sep=';')%>%
  mutate(across(-c(1,2), ~recode(., 'A/A' = 1, 'A/G'= 2, 'G/G'=3)))%>%
  mutate(across(2, ~recode(., 'C/T'=1, 'T/T'=2)))%>%
  mutate(across(-1, ~factor(.))) # Recode genotypes to factor


cyt<-read_delim("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/cytometry.csv", 
                delim = ";", escape_double = FALSE, trim_ws = TRUE)
symp<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/sympthptoms.csv",
               header=T, sep=';')%>%
  mutate(across(-1, ~ recode(., 'NO' = 0, 'SI' = 1)))%>%
  replace(is.na(.), 0)%>% 
  mutate(across(-1, ~factor(.))) # Recode sympthoms to factor


meta<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/metadata.csv",
               header = T, sep=';')[, c(1:12)]


######## Dataframe building ############

##### Joining datasets #####

dfs <- list(meta, div, clon, vus, jus, cdr3)
int_df <- reduce(dfs, full_join, by = "sample")

dfs_2<- list(cyt, genot, symp)
int_df2<- reduce(dfs_2, inner_join, by= "code")

dfs_3<- list(int_df, int_df2)
final_df<- reduce(dfs_3, inner_join, by= "code")


##### Scaling #####

final_scale_df <- final_df %>%
  mutate_if(is.numeric, function(x) (x - median(x)) / IQR(x))%>%
  dplyr::select(-SPT.AGE_ASV)

write.csv(final_scale_df, 'final_scale_df.csv', quote=F)

##### Dataframe of predictors #####

sev<-final_scale_df[,'Severity'] # severity as known status
sev_55<-final_scale_df[,'sev_55'] # sev_55 as known status
x<-final_scale_df[, c(13:length(names(final_scale_df)))]
rownames(x)<-final_scale_df$sample


######### LCM #############

# Clustering with variable selection with BIC
set.seed(1234)
ari_list_bic <- c()
res_var_sel_bic_list <- list()

for (i in 1:100) {
res_var_sel_bic<-VarSelCluster(x, gvals = 1:3, vbleSelec = TRUE, nbcores=4,
                               crit.varsel = 'BIC')
res_var_sel_bic_list[[paste0("res_var_sel_bic_", i)]] <- res_var_sel_bic

bic<-BIC(res_var_sel_bic)
ari<-ARI(sev_55, fitted(res_var_sel_bic))
ari_list_bic <- c(ari_list_bic, ari)
}

ari_list_bic


match(max(ari_list_bic), ari_list_bic)
res.bic<-res_var_sel_bic_list$res_var_sel_bic_4


# Estimated probabilities for each observation in each new class


head(fitted(res.bic, type="probability"))
table(fitted(res.bic))
summary(res.bic)
saveRDS(res.bic, file='lcm_bic.RData')

classes<-fitted(res.bic)


# Selecting the BIC criteria to futher study

VarSelShiny(res.bic)


selected_var<-res.bic@model@names.relevant

relevance<-res.bic@criteria@discrim

names <- c("S.GGE_FIY", "SYGGE_4_22", "SY.GE_GST", "TRBV11_3", 
             "SL.SYE_DGNST", "SVG.E_DN", "CMCD36y", "TRBV6_3", 
             "S.GNE_AGILV", "S.GGYE_FLY", "SLG.YE_AGSTV", "TRBV9...", 
             "SLGL.YE_GHNRS", "G.YE_DEGNS", "TRBV12_3", "TRBV16",
             "normalizedShannonWienerIndex_mean", "TRBV11_1", "GS", 
             "TRBV18", "TRBJ1_6", "Hyperexpanded...", "D50", "TRBV15", 
             "SPG.E_DHNTWY", "SP.YE_GHNRST", "GINI...", "Large", "TRBV7_8", 
             "NCMCD163nCD206y", "TRBV14", "TRBJ2_5", "TRBV6_4", "TRBV7_4",
             "TRBJ2_2", "TRBV6_9", "Rare", "TRBJ2_3", "CD45n", "TMCD36y",
             "NCMCD163nCD206n", "CMCD163nCD206y", "NCMCD163yCD206y", 
             "NCMCD11by", "CMCD38y", "TRBV11_2", "SS.YE_AGST", "YSSGE_4_22...", 
             "TRBV4_3", "SP.TGE_DEGQRS", "pneu", "NCMCD38y", "TRBV5_5",
             "TRBV7_6", "TRBV12_5", "TRBV5_1", "TRBV3_1", "TRBV6_6", 
             "TRBV19", "dysn", "flu_vacc")
values <- c(175.25258321, 175.25258321, 154.59740775, 90.22304046,
            81.67897365, 74.20960749, 65.23287517, 59.51646927, 
            51.00811431, 49.63872927, 47.41944927, 38.31325068,
            26.16796974, 23.96635888, 22.29562192, 19.77537224, 
            19.74042621, 19.24019181, 18.87035398, 17.41901175, 
            16.53030677, 15.87899733, 13.39049136, 12.68132789, 
            12.25880085, 12.19661838, 12.17505821, 11.21600880, 
            10.86304548, 10.20544448, 8.51565932, 7.84323888, 
            7.75432978, 7.22086165, 6.88904671, 6.22033821, 
            5.93469779, 5.15426347, 4.75065456, 4.56387465, 
            4.32209174, 4.17438524, 4.05612663, 4.04542590, 
            3.71546151, 3.70239523, 3.52362209, 3.04247149, 
            2.46840949, 2.33444933, 2.25900602, 2.10039451, 
            2.06755518, 1.86993570, 1.41721023, 0.79569348, 
            0.74112735, 0.71920641, 0.46842775, 0.19445982, 0.13803519)
df_relevance <- data.frame(Names = names, Values = values)
df_relevance2<-df_relevance%>%
  filter(., Names %in% selected_var)%>%
  mutate(., perc=(Values/sum(Values)*100))


write.csv(df_relevance2, 'df_relevance2.csv', quote = F)
df_relevance3<-read.delim('C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/df_relevance2.csv', 
                          header = T, stringsAsFactors = F, sep = ';')
df_relevance3$Names <- factor(df_relevance3$Names) %>%
  fct_reorder(df_relevance3$perc)

top_25_relevance<-ggplot(data=df_relevance3[c(1:25), ], aes(x = Names, y=perc, fill=Class))+ 
  geom_col(color="black")+
  coord_flip()+
  ylab('% of relevance')+
  xlab('Variable')+
  theme_minimal()+
  theme(legend.text = element_text(size=16))+
  theme(axis.text.y = element_text(size=12, colour = 'black'))+
  theme(axis.title.y = element_text(size=22))+
  theme(axis.text.x = element_text(size=16, colour = 'black'))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.title.x = element_text(size=22))+
  theme(legend.title = element_blank())


#### Build the new dataframe of selected variables with each cluster in each observation

selected_df3<-final_scale_df%>%
  select(sample, Severity, Morethan55, sev_55, all_of(selected_var))%>% 
  mutate(cluster = classes ) %>% 
  relocate(cluster, .after = 1)
write.csv(selected_df3, 'selected_df3.csv', quote = F)

# Changes in varuable type between pre- and post LCM datasets
library('forcats')
var_dist<-read.delim('C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/var_dist.csv',
                   header = T, stringsAsFactors = F, sep = ';')

var_dist_plot<-ggplot(data = var_dist, aes(fct_infreq(Status), fill=Class))+
  geom_bar(color='black', position="fill")+
  scale_y_continuous(labels = scales::percent)+
  ylab('% of relevance')+
  theme_minimal()+
  theme(legend.title = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(legend.text = element_text(size=16))+
  theme(axis.text.y = element_text(size=16, colour = 'black'))+
  theme(axis.title.y = element_text(size=22))+
  theme(axis.text.x = element_text(size=20, colour = 'black'))


#### Severity distribution among clusters

table(selected_df3$Severity, selected_df3$cluster)
table(selected_df3$sev_55, selected_df3$cluster)

chisq2sev<-chisq.test(table(selected_df3$Severity, selected_df3$cluster),
           simulate.p.value = T, B = 1000)
chisq2sev55<-chisq.test(table(selected_df3$sev_55, selected_df3$cluster), 
           simulate.p.value = T, B = 1000)

library(rcompanion)

cramerV(table(selected_df3$Severity, selected_df3$cluster), ci = T, 
        bias.correct = T)
cramerV(table(selected_df3$sev_55, selected_df3$cluster), ci = T, 
        bias.correct = T)

tb<-selected_df3%>%
  select(., Severity, sev_55, cluster)%>%
  count(Severity, sev_55, cluster)

library(ggalluvial)

ggplot(data = tb,
       aes(axis1 = Severity, axis2 = cluster,
           y = n)) +
  scale_x_discrete(limits = c("Severity", "Cluster"), expand = c(.2, .05)) +
  ylab('Number of patients')+
  geom_alluvium(aes(fill=sev_55)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=6) +
  theme_minimal()+
  theme(legend.position="top") +
  theme(legend.text = element_text(size=16))+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size=16, colour = 'black'))+
  theme(axis.title.y = element_text(size=22))+
  theme(axis.text.x = element_text(size=22, colour = 'black'))


############## FAMD with filtered variables and clusters ##############

library(FactoMineR)
library(factoextra)
library(gridExtra)
library(grid)
library(plotly)


selected_df3<-selected_df3%>%
  rename(., 'SW'= "normalizedShannonWienerIndex_mean")%>%
  rename(., 'GINI'='GINI...')%>%
  rename(., 'Hyperexpanded'='Hyperexpanded...')%>%
  rename(., 'TRBV9'='TRBV9...')%>%
  rename(., 'TRBJ2_7'='TRBJ2_7...')
  

x_filt<-selected_df3[, c(6:length(names(selected_df3)))]
rownames(x_filt)<-selected_df3$sample
res.famd_new<-FAMD(x_filt, graph = FALSE, ncp = 5)


  #### PD contrbutions ####

  eig.val <- get_eigenvalue(res.famd_new)
  head(eig.val)
  
  a<-fviz_eig(res.famd_new, choice = 'eigenvalue', 
              geom = c('bar','line'), barfill = '#3D3576', barcolor = '#3D3576',
              main="")
  b<-fviz_eig(res.famd_new, choice = 'variance', 
              geom = c('bar','line'), barfill = '#00BC7F', barcolor = '#00BC7F',
              main="")
  grid.arrange(a, b, ncol=2)
  
  ### Visualize plot
  
  val_df <- as.data.frame(res.famd_new$ind)
  df<- cbind(selected_df3[,c(1:6)], val_df[1:3])
  
  ##### Plot ####
 

  
  FAMD_post_sev<-fviz_mfa_ind(res.famd_new,repel=T, addEllipses = T,
                  ellipse.type = "convex",
                  geom = 'point', 
                  pointsize = 5, 
                  alpha.ind = 0.4, 
                  habillage =as.factor(selected_df3$Severity), 
                  show.clust.cent = FALSE,
                  palette = c('#3E6B97', '#E13E43'))
  
  FAMD_post_sev + geom_point(data =  FAMD_post_sev$data, aes(x = x, y = y), 
                 color = "black", 
                 size=5, alpha=0.2)+
    theme(axis.title = element_text(size = 20),  
          axis.text = element_text(size = 16),
          title = element_blank(),
          legend.text = element_text(size=14))

  FAMD_post_sev55<-fviz_mfa_ind(res.famd_new,repel=T, addEllipses = T,
                              ellipse.type = "convex",
                              geom = 'point', 
                              pointsize = 5, 
                              alpha.ind = 0.4, 
                              habillage =as.factor(selected_df3$sev_55), 
                              show.clust.cent = FALSE,
                              palette = c('#2E5A87', '#86CAE1', '#FC7662', '#B4193A'))
  
  FAMD_post_sev55 + geom_point(data =  FAMD_post_sev55$data, aes(x = x, y = y), 
                             color = "black", 
                             size=5, alpha=0.2)+
    theme(axis.title = element_text(size = 20),  
          axis.text = element_text(size = 16),
          title = element_blank(),
          legend.text = element_text(size=14))
  

  
  #### Variables ####
  
  # Contribution barplots
  
  b<-fviz_contrib(res.famd_new, "var", axes = 1, 
                  top = 20, xtickslab.rt = 60, fill=  '#3D3576', color = '#3D3576')
  c<-fviz_contrib(res.famd_new, "var", axes = 2, 
                  top=20, xtickslab.rt = 60, fill = '#00BC7F', color='#00BC7F')
  grid.arrange(b, c, ncol=2)
  
  
  
  cos2_post<-fviz_famd_var(res.famd_new, "quanti.var", 
                col.var = 'cos2', 
                select.var = list(contrib=20),
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                shape.var = 17, 
                repel = TRUE,
                labelsize = 5)
  cos2_post<-cos2_post+ggtitle("Quantitative variables: FMAD post LCM")+
    theme(axis.title = element_text(size = 20),  
          axis.text = element_text(size = 16),
          title = element_blank(),
          legend.text = element_text(size=14))
  
  
  p<-fviz_mfa_ind(res.famd_new,repel=T, addEllipses = T,
              ellipse.type = "convex",
               geom = 'point', 
               pointsize = 5, 
               alpha.ind = 0.4, 
               habillage =as.factor(selected_df3$cluster), 
               show.clust.cent = FALSE,
               palette = c("#4B0055", "#009F94", "#FDE333"))
  
  p + geom_point(data = p$data, aes(x = x, y = y), 
                 color = "black", 
                 size=5, alpha=0.2)+
                 theme(axis.title = element_text(size = 20),  
                       axis.text = element_text(size = 16),
                       title = element_blank(),
                       legend.text = element_text(size=14))
                 labs(colour = "Clusters")
  
#################### PERMANOVA between LCM clusters ###################
 
library(vegan)
library(RVAideMemoire)
  
# Quantitative variables
  
quant_selected_df3<-selected_df3%>%select(., -pneu)
perm_eu<-adonis2(quant_selected_df3[, -c(1:5)]~quant_selected_df3$cluster,
          permutations=1000, parallel = 4, 
          method = 'euclidean')
RVAideMemoire::pairwise.perm.manova(dist(quant_selected_df3[, -c(1:5)], "euclidean"),
                       quant_selected_df3$cluster, 
                       nperm=1000, p.method = 'BH', progress=F)

quant_selected_df3_sev<-quant_selected_df3[, -c(1,2,4,5)]
perm_eu_sev<-adonis2(quant_selected_df3_sev%>%select(., -Severity)~quant_selected_df3_sev$Severity,
                 permutations=1000, parallel = 4, 
                 method = 'euclidean')

quant_selected_df3_sev55<-quant_selected_df3[, -c(1,2,3,4)]
perm_eu_sev55<-adonis2(quant_selected_df3_sev55%>%select(., -sev_55)~quant_selected_df3_sev55$sev_55,
                     permutations=1000, parallel = 4, 
                     method = 'euclidean')


################# Distribution of top 25 LCM variables between groups ##########

### Dataframe of un-scaled variables and clusters ####

final_df<-final_df%>%
  rename(., 'SW'= "normalizedShannonWienerIndex_mean")%>%
  rename(., 'GINI'='GINI...')%>%
  rename(., 'Hyperexpanded'='Hyperexpanded...')%>%
  rename(., 'TRBV9'='TRBV9...')%>%
  rename(., 'TRBJ2_7'='TRBJ2_7...')%>%
  rename(., 'YSSGE_4_22' = 'YSSGE_4_22...')

final_df_25<-final_df%>%select(., sample, all_of(df_relevance3$Names[c(1:25)]))

final_df_25<-final_df_25%>%
  full_join(., selected_df3%>%select(., sample, cluster), 'sample')

final_df_25<-final_df_25%>%relocate(cluster, .after = 1)

quant_vars_div<-colnames(final_df_25%>%select(., SW, GS, D50, Hyperexpanded))
quant_vars_v <- colnames(final_df_25 %>% select(starts_with("TRBV")))
quant_vars_j <- colnames(final_df_25 %>% select(starts_with("TRBJ")))
quant_vars_motif<- colnames(final_df_25[, -c(1,2)] %>% select(,-starts_with("TRBJ"),
                                                   -starts_with("TRBV"),
                                                   -SW, -GS, -D50, -Hyperexpanded, -CMCD36y))
  


library(rstatix)
library(patchwork)


# Diversity and clonality

top25_plots_div <- list()
for (varname in quant_vars_div){
  formula<-as.formula(paste(varname, "~ cluster"))
  stat.test<-dunn_test(formula,
                       data = final_df_25,
                       p.adjust.method = 'BH')
  stat.test$p.adj<-round(stat.test$p.adj, 8)
  stat.test <- stat.test %>% add_xy_position(x = "cluster")%>%filter(., p.adj<0.05)
  
  
  p<-ggboxplot(final_df_25, x = "cluster", y = varname, 
            color = "black", fill = "cluster", alpha=0.3,  palette = c("#4B0055", "#009F94", "#FDE333"),
            add = "jitter", shape = 21, add.params = list(color = "black")) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank())+
    theme(axis.title.x = element_blank())+
    theme(axis.text.y = element_text(size=16, colour = 'black'))+
    theme(axis.title.y = element_text(size=22))+
    stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, size=6,
                       step.increase = c(0.05, 0.075))
  top25_plots_div[[varname]] <- p
  
}


div_results<-top25_plots_div[[1]]
for (i in 2:length(top25_plots_div)){
  div_results <- div_results | top25_plots_div[[i]]
}

print(div_results)




# V alleles

top25_plots_v <- list()
for (varname in quant_vars_v){
  formula<-as.formula(paste(varname, "~ cluster"))
  stat.test<-dunn_test(formula,
                       data = final_df_25,
                       p.adjust.method = 'BH')
  stat.test$p.adj<-round(stat.test$p.adj, 4)
  stat.test <- stat.test %>% add_xy_position(x = "cluster")%>%filter(., p.adj<0.05)
  
  p<-ggboxplot(final_df_25, x = "cluster", y = varname, 
               color = "black", fill = "cluster", alpha=0.3,  palette = c("#4B0055", "#009F94", "#FDE333"),
               add = "jitter", shape = 21, add.params = list(color = "black")) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank())+
    theme(axis.title.x = element_blank())+
    theme(axis.text.y = element_text(size=16, colour = 'black'))+
    theme(axis.title.y = element_text(size=22))+
    scale_y_log10()+
    theme(axis.text.y = element_text(size=16, colour = 'black'))+
    stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, size=6,
                       hide.ns = T, step.increase = 0.2)
  
  top25_plots_v[[varname]] <- p
  
}

v_results<-top25_plots_v[[1]]
for (i in 2:length(top25_plots_v)){
  v_results <- v_results | top25_plots_v[[i]]
}

print(v_results)


# J alelles

top25_plots_j <- list()
for (varname in quant_vars_j){
  formula<-as.formula(paste(varname, "~ cluster"))
  stat.test<-dunn_test(formula,
                       data = final_df_25,
                       p.adjust.method = 'BH')
  stat.test$p.adj<-round(stat.test$p.adj, 4)
  stat.test <- stat.test %>% add_xy_position(x = "cluster")%>%filter(., p.adj<0.05)
  
  p<-ggboxplot(final_df_25, x = "cluster", y = varname, 
               color = "black", fill = "cluster", alpha=0.3,  palette = c("#4B0055", "#009F94", "#FDE333"),
               add = "jitter", shape = 21, add.params = list(color = "black")) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank())+
    theme(axis.title.x = element_blank())+
    scale_y_log10()+
    theme(axis.text.y = element_text(size=16, colour = 'black'))+
    stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, size=4,
                       hide.ns = T, step.increase = 0.2)
  
  top25_plots_j[[varname]] <- p
  
}

j_results<-top25_plots_j[[1]]
for (i in 2:length(top25_plots_j)){
  j_results <- j_results | top25_plots_j[[i]]
}

print(j_results)



# GLIPH2 Motifs

top25_plots_motif <- list()
for (varname in quant_vars_motif){
  formula<-as.formula(paste(varname, "~ cluster"))
  stat.test<-dunn_test(formula,
                       data = final_df_25,
                       p.adjust.method = 'BH')
  stat.test$p.adj<-round(stat.test$p.adj, 4)
  stat.test <- stat.test %>% add_xy_position(x = "cluster")%>%filter(., p.adj<0.05)
  
  p<-ggboxplot(final_df_25, x = "cluster", y = varname, 
               color = "black", fill = "cluster", alpha=0.3,  palette = c("#4B0055", "#009F94", "#FDE333"),
               add = "jitter", shape = 21, add.params = list(color = "black")) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank())+
    theme(axis.title.x = element_blank())+
    scale_y_log10()+
    theme(axis.text.y = element_text(size=16, colour = 'black'))+
    theme(axis.title.y = element_text(size=22))+
    theme(axis.text.y = element_text(size=16, colour = 'black'))+
    stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, size=6,
                       hide.ns = T, step.increase = 0.2)
  
  top25_plots_motif[[varname]] <- p
  
}

motif_results<-top25_plots_motif[[1]]
for (i in 2:length(top25_plots_motif)){
  motif_results <- motif_results | top25_plots_motif[[i]]
}

motif_results + plot_layout(nrow = 2)

signif_plots<-list(top25_plots_motif$S.GGE_FIY, top25_plots_motif$SYGGE_4_22, 
                   top25_plots_motif$SY.GE_GST, top25_plots_motif$SL.SYE_DGNST, 
                   top25_plots_motif$S.GNE_AGILV, top25_plots_motif$S.GGYE_FLY)

motif_sig<-signif_plots[[1]]
for (i in 2:length(signif_plots)){
  motif_sig <- motif_sig | signif_plots[[i]]
}


signif_plots_v<-list(top25_plots_v$TRBV6_3, top25_plots_v$TRBV12_3, 
                     top25_plots_v$TRBV15)


v_sig<-signif_plots_v[[1]]
for (i in 2:length(signif_plots_v)){
  v_sig <- v_sig | signif_plots_v[[i]]
}


# Cytof

formula<-as.formula(paste('CMCD36y', "~ cluster"))
stat.test<-dunn_test(formula,
                     data = final_df_25,
                     p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "cluster")%>%filter(., p.adj<0.05)

p<-ggboxplot(final_df_25, x = "cluster", y = 'CMCD36y', 
             color = "black", fill = "cluster", alpha=0.3,  palette = c("#4B0055", "#009F94", "#FDE333"),
             add = "jitter", shape = 21, add.params = list(color = "black")) +
  theme(legend.position = "none") +
  ylab('CM-CD36+')+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.text.y = element_text(size=16, colour = 'black'))+
  theme(axis.title.y = element_text(size=22))+
  theme(axis.text.y = element_text(size=16, colour = 'black'))+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, size=6,
                     hide.ns = T, step.increase = 0.1, y.position = 101)




replace_underscore <- function(x) {
  gsub("_", "-", x)
}

replace_yn <- function(x) {
  gsub("y", "+", gsub("n", "-", gsub("M", "M-", x)))
}



div_70_lcm<-final_df_lcm[, c(2:9, 72)]%>%
  rename(., cluster=cluster.x)%>%
  mutate(cluster = as.factor(cluster)) %>%
  relocate(., cluster, .before = 1)%>%
  mutate_if(., .predicate = is.numeric, scale)%>%
  group_by(., cluster)%>%
  summarize_all(., mean)%>%
  



v_70_lcm<-final_df_lcm[, c(10:37, 72)]%>%
rename(., cluster=cluster.x)%>%
  mutate(cluster = as.factor(cluster)) %>%
  relocate(., cluster, .before = 1)%>%
  mutate_if(., .predicate = is.numeric, scale)%>%
  group_by(., cluster)%>%
  summarize_all(., mean)%>% rename_all(replace_underscore)


j_70_lcm<-final_df_lcm[, c(38:43, 72)]%>%
rename(., cluster=cluster.x)%>%
  mutate(cluster = as.factor(cluster)) %>%
  relocate(., cluster, .before = 1)%>%
  mutate_if(., .predicate = is.numeric, scale)%>%
  group_by(., cluster)%>%
  summarize_all(., mean)%>% rename_all(replace_underscore)



motif_70_lcm<-final_df_lcm[, c(44:59, 72)]%>%
rename(., cluster=cluster.x)%>%
  mutate(cluster = as.factor(cluster)) %>%
  relocate(., cluster, .before = 1)%>%
  mutate_if(., .predicate = is.numeric, scale)%>%
  group_by(., cluster)%>%
  summarize_all(., mean)%>%
  rename_all(replace_underscore)%>%
  rename_all(~ str_replace(., "\\.", "%"))
  
  
  
prueba<-motif_70_lcm%>%rename_all(replace_mot)


cytof_70_lcm<-final_df_lcm[, c(60:70, 72)]%>%
  rename(., cluster=cluster.x)%>%
  mutate(cluster = as.factor(cluster)) %>%
  relocate(., cluster, .before = 1)%>%
  mutate_if(., .predicate = is.numeric, scale)%>%
  group_by(., cluster)%>%
  summarize_all(., mean)%>%rename_all(replace_yn)  



radar_divclon<-ggradar(
  div_70_lcm, 
  grid.min = min(as.matrix(div_70_lcm[, -1])), grid.mid = 0, grid.max = max(as.matrix(div_70_lcm[, -1])),
  # Polygons
  values.radar = c( grid.min=round(min(as.matrix(div_70_lcm[, -1])),2), grid.mid = 0, grid.max=round(max(as.matrix(div_70_lcm[, -1])),2)),
  group.line.width = 1.5, 
  group.point.size = 3,
  group.colours = c("#4B0055", "#009F94", "#FDE333"),
  # Background and grid lines
  gridline.mid.colour = "grey",
  legend.position = "top", axis.label.size = 6
)




radar_v<-ggradar(
  v_70_lcm, 
  grid.min = min(as.matrix(div_70_lcm[, -1])), grid.mid = 0, grid.max = max(as.matrix(v_70_lcm[, -1])),
  # Polygons
  values.radar = c( grid.min=round(min(as.matrix(v_70_lcm[, -1])),2), grid.mid = 0, grid.max=round(max(as.matrix(v_70_lcm[, -1])),2)),
  group.line.width = 1.5, 
  group.point.size = 3,
  group.colours = c("#4B0055", "#009F94", "#FDE333"),
  # Background and grid lines
  gridline.mid.colour = "grey",
  legend.position = "bottom",
  axis.label.size = 4, plot.legend = FALSE
)



radar_j<-ggradar(
  j_70_lcm, 
  grid.min = min(as.matrix(j_70_lcm[, -1])), grid.mid = 0, grid.max = max(as.matrix(j_70_lcm[, -1])),
  # Polygons
  values.radar = c( grid.min=round(min(as.matrix(j_70_lcm[, -1])),2), grid.mid = 0, grid.max=round(max(as.matrix(j_70_lcm[, -1])),2)),
  group.line.width = 1.5, 
  group.point.size = 3,
  group.colours = c("#4B0055", "#009F94", "#FDE333"),
  # Background and grid lines
  gridline.mid.colour = "grey",
  legend.position = "bottom", plot.legend = FALSE, axis.label.size = 4
)


radar_motif<-ggradar(
  motif_70_lcm, 
  grid.min = min(as.matrix(motif_70_lcm[, -1])), grid.mid = 0, grid.max = max(as.matrix(motif_70_lcm[, -1])),
  # Polygons
  values.radar = c( grid.min=round(min(as.matrix(motif_70_lcm[, -1])),2), grid.mid = 0, grid.max=round(max(as.matrix(motif_70_lcm[, -1])),2)),
  group.line.width = 1.5, 
  group.point.size = 3,
  group.colours = c("#4B0055", "#009F94", "#FDE333"),
  # Background and grid lines
  gridline.mid.colour = "grey",
  legend.position = "bottom",
  plot.legend = FALSE,
  axis.label.size = 4
)




radar_cytof<-ggradar(
  cytof_70_lcm, 
  grid.min = min(as.matrix( cytof_70_lcm[, -1])), grid.mid = 0, grid.max = max(as.matrix( cytof_70_lcm[, -1])),
  values.radar = c( grid.min=round(min(as.matrix( cytof_70_lcm[, -1])),2), grid.mid = 0, grid.max=round(max(as.matrix( cytof_70_lcm[, -1])),2)),
  group.line.width = 1.5, 
  group.point.size = 3,
  group.colours = c("#4B0055", "#009F94", "#FDE333"),
  gridline.mid.colour = "grey",
  legend.position = "bottom",
  plot.legend = FALSE,
  axis.label.size = 4
)

