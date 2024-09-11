# FAMD analysis
# Package loading
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(grid)
library(plotly)

##### Data Loading ######

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
  mutate(across(-1, ~factor(.))) 


cyt<-read_delim("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/cytometry.csv", 
                delim = ";", escape_double = FALSE, trim_ws = TRUE)
symp<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/sympthptoms.csv",
               header=T, sep=';')%>%
  replace(is.na(.), 'NO')%>% 
  mutate(across(-1, ~factor(.))) 


meta<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/metadata.csv",
               header = T, sep=';')[, c(1:12)]

####### Joining datasets ########

dfs <- list(meta, div, clon, vus, jus, cdr3)
int_df <- reduce(dfs, full_join, by = "sample")

dfs_2<- list(cyt, genot, symp)
int_df2<- reduce(dfs_2, inner_join, by= "code")

dfs_3<- list(int_df, int_df2)
final_df<- reduce(dfs_3, inner_join, by= "code")

######## Fatctorize response variable ######

final_df<-final_df%>%
  mutate(Severity = fct_recode(Severity, "1" = "Mild", "2" = "Severe"))%>%
  mutate(sev_55 = fct_recode(sev_55, "1" = "Mild;<55", "2" = "Mild;>55",
                             "3" = "Severe;<55", "4" = "Severe;>55"))


######### Scaling #########
# Due to the number of non-normally variables, we apply a robust transformation
# SPT.AGE_ASV has been removed due to the presence of inf an Nan values after scaling

final_scale_df <- final_df %>%
  mutate_if(is.numeric, function(x) (x - median(x)) / IQR(x))%>%
  dplyr::select(-SPT.AGE_ASV)

############ FAMD analysis  ##################


final_scale_df<-final_scale_df%>%
  rename(., 'SW'= "normalizedShannonWienerIndex_mean")%>%
  rename(., 'GINI'='GINI...')%>%
  rename(., 'Hyperexpanded'='Hyperexpanded...')%>%
  rename(., 'TRBV9'='TRBV9...')%>%
  rename(., 'TRBJ2_7'='TRBJ2_7...')


sev<-final_scale_df[,'Severity'] # severity as known status
sev_55<-final_scale_df[,'sev_55'] # sev_55 as known status
x<-final_scale_df[, c(13:length(names(final_scale_df)))]
rownames(x)<-final_scale_df$sample

res.famd<-FAMD(x, graph = FALSE, ncp = 5)

### Visualize PD contrbutions

eig.val <- get_eigenvalue(res.famd)
head(eig.val)

eig_pre<-fviz_eig(res.famd, choice = 'eigenvalue', 
         geom = c('bar','line'), barfill = '#C51517', barcolor = '#C51517', 
         main="")
var_pre<-fviz_eig(res.famd, choice = 'variance', 
         geom = c('bar','line'), barfill = '#1C5A99',
         main="")
grid.arrange(eig_pre, var_pre, ncol=2)

val_df <- as.data.frame(res.famd$ind)
df<- cbind(final_scale_df[,c(5,8,9,12)], val_df[1:3])

## Plot

FAMD_pre_sev<-fviz_mfa_ind(res.famd,repel=T, addEllipses = T,
                            ellipse.type = "convex",
                            geom = 'point', 
                            pointsize = 5, 
                            alpha.ind = 0.4, 
                            habillage =as.factor(selected_df3$Severity), 
                            show.clust.cent = FALSE,
                            palette = c('#3E6B97', '#E13E43'))

FAMD_pre_sev + geom_point(data = FAMD_pre_sev$data, aes(x = x, y = y), 
                           color = "black", 
                           size=5, alpha=0.2)+
  theme(axis.title = element_text(size = 20),  
        axis.text = element_text(size = 16),
        title = element_blank(),
        legend.text = element_text(size=14))

FAMD_pre_sev55<-fviz_mfa_ind(res.famd,repel=T, addEllipses = T,
                              ellipse.type = "convex",
                              geom = 'point', 
                              pointsize = 5, 
                              alpha.ind = 0.4, 
                              habillage =as.factor(selected_df3$sev_55), 
                              show.clust.cent = FALSE,
                              palette = c('#2E5A87', '#86CAE1', '#FC7662', '#B4193A'))

FAMD_pre_sev55 + geom_point(data = FAMD_pre_sev55$data, aes(x = x, y = y), 
                             color = "black", 
                             size=5, alpha=0.2)+
  theme(axis.title = element_text(size = 20),  
        axis.text = element_text(size = 16),
        title = element_blank(),
        legend.text = element_text(size=14))


# Plot of variables

b_pre<-fviz_contrib(res.famd, "var", axes = 1, 
                top = 20, xtickslab.rt = 45, 
                fill = '#C51517', color = '#C51517')

c_pre<-fviz_contrib(res.famd, "var", axes = 2, 
                top=20, xtickslab.rt = 45, 
                fill = '#1C5A99', color =  '#1C5A99')

grid.arrange(b_pre, c_pre, ncol=2)


cos2_pre<-fviz_famd_var(res.famd, "quanti.var", col.var = 'cos2', select.var = list(cos2=15),
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)
cos2_pre<-cos2_pre+ggtitle("Quantitative variables: FMAD pre LCM")+
  theme(axis.title = element_text(size = 20),  
        axis.text = element_text(size = 16),
        title = element_blank(),
        legend.text = element_text(size=14))



cos2_qual_pre<-fviz_famd_var(res.famd, "var", col.var ='cos2',  select.var = list(name=c('ACE2_rs2285666',                    
                                              'MX1_rs469390',                      
                                              'TMPRSS2_rs2070788' ,                
                                              'aff_dem'  ,                         
                                              'anosm' ,                            
                                              'ageusia',                          
                                              'mial',                              
                                              'cef' ,                              
                                              'fever',                             
                                              'dysn' ,                             
                                              'aesth' ,                            
                                              'flu_vacc' ,                         
                                              'long_covid',                        
                                              'pneu')),
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)




fmad_ace2<-fviz_mfa_ind(res.famd, habillage = as.factor(final_scale_df$ACE2_rs2285666)
             ,repel=T, 
             addEllipses = T, 
             ellipse.type = "convex",
             geom = 'point', 
             pointsize = 3, 
             alpha.ind = 0.4, title='ACE2_rs2285666')+scale_color_discrete(labels=c('C/T', 'T/T'))

fmad_tmprss2<-fviz_mfa_ind(res.famd, habillage =as.factor(final_scale_df$TMPRSS2_rs2070788)
             ,repel=T, 
             addEllipses = T, 
             geom = 'point',
             ellipse.type = "convex",
             pointsize = 3, 
             alpha.ind = 0.4, title='TMPRSS2_rs2070788')+scale_color_discrete(labels = c('A/A', 'A/G', 'G/G'))

fmad_mx1<-fviz_mfa_ind(res.famd, habillage = as.factor(final_scale_df$MX1_rs469390)
             ,repel=T, 
             addEllipses = T, 
             ellipse.type = "convex",
             geom = 'point', 
             pointsize = 3, 
             alpha.ind = 0.4, title='MX1_rs469390')+scale_color_discrete(labels = c('A/A', 'A/G', 'G/G'))


grid.arrange(fmad_ace2, fmad_tmprss2,fmad_mx1, nrow=1)

#### PERMANOVA TESTS OVER GENOTYPES
library(vegan)
library(RVAideMemoire)

quant_ace<- final_scale_df[,-c(142:length(colnames(final_scale_df)))]
perm_ace2<-adonis2(quant_ace[, -c(1:12)]%>%select(., -ACE2_rs2285666)~quant_ace$ACE2_rs2285666,
                 permutations=1000, parallel = 4, 
                 method = 'euclidean')


quant_mx1<- final_scale_df[,-c(143:length(colnames(final_scale_df)))]%>%select(., -ACE2_rs2285666)
perm_mx1<-adonis2(quant_mx1[, -c(1:12)]%>%select(., -MX1_rs469390)~quant_mx1$MX1_rs469390,
                   permutations=1000, parallel = 4, 
                   method = 'euclidean')

quant_tmprss2<- final_scale_df[,-c(144:length(colnames(final_scale_df)))]%>%select(., -ACE2_rs2285666, -MX1_rs469390)
perm_tmprss2<-adonis2(quant_tmprss2[, -c(1:12)]%>%select(., -TMPRSS2_rs2070788)~quant_tmprss2$TMPRSS2_rs2070788,
                  permutations=1000, parallel = 4, 
                  method = 'euclidean')

quant_tmprss2<- final_scale_df[,-c(144:length(colnames(final_scale_df)))]%>%select(., -ACE2_rs2285666, -MX1_rs469390)
perm_tmprss2<-adonis2(quant_tmprss2[, -c(1:12)]%>%select(., -TMPRSS2_rs2070788)~quant_tmprss2$TMPRSS2_rs2070788,
                      permutations=1000, parallel = 4, 
                      method = 'euclidean')

quant_sev<-final_scale_df[,-c(141:length(colnames(final_scale_df)))]
quant_sev<-quant_sev[, -c(1:8, 10:12)]
perm_sev<-adonis2(quant_sev%>%select(., -Severity)~quant_sev$Severity,
                      permutations=1000, parallel = 4, 
                      method = 'euclidean' )

quant_sev55<-final_scale_df[,-c(141:length(colnames(final_scale_df)))]
quant_sev55<-quant_sev55[, -c(1:11)]
perm_sev55<-adonis2(quant_sev55%>%select(., -sev_55)~quant_sev55$sev_55,
                  permutations=1000, parallel = 4, 
                  method = 'euclidean' )



