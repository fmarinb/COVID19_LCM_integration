####### MAIN DATASET CREATION #######

####### Packages ########
library(tidyverse)
library(DescTools)
library(MASS)
library(mice)

####### Data loading #####

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

####### Joining datasets ########

dfs <- list(meta, div, clon, vus, jus, cdr3)
int_df <- reduce(dfs, full_join, by = "sample")

dfs_2<- list(cyt, genot, symp)
int_df2<- reduce(dfs_2, inner_join, by= "code")

dfs_3<- list(int_df, int_df2)
final_df<- reduce(dfs_3, inner_join, by= "code")

is.numeric(final_df$pneu)

###### Exploratoy analyses and transformation ######

library(gtsummary)

final_df$Sex[final_df$Sex=='Man']<-'Male'
final_df$Sex[final_df$Sex=='Woman']<-'Female'

summary_df<-final_df%>%select(., Severity, Morethan55, Age, Sex)%>%
  tbl_summary(by=Severity)%>%
  as_gt() %>%
  gtsave("summary.docx")

summary_df<-final_df%>%select(., Severity, Morethan55, Age, Sex)%>%
  tbl_summary(by=Severity)%>%
  as_gt() %>%
  gtsave("summary.docx")


# Variable type and dimensions
str(final_df)

# Normality test
  # Lilliefors Test
num_final_df<-final_df%>%select_if(is.numeric)
normality_test<-sapply(num_final_df, LillieTest)


plot_data<-data.frame(variable=rep(names(num_final_df),
                                   each=nrow(num_final_df)),
                      value=unlist(num_final_df))

  # Gaussian distributions
ggplot(plot_data[c(1:4758),], aes(x = value)) +
  geom_density(fill = "blue", alpha = 0.5) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

ggplot(plot_data[c(4758:nrow(plot_data)),], aes(x = value)) +
  geom_density(fill = "blue", alpha = 0.5) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

  # QQplot
ggplot(plot_data[c(1:4758),], aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

ggplot(plot_data[c(4758:nrow(plot_data)),], aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

  # 98 variables with non-nomal distribution

var_nonnorm<-as.data.frame((t(normality_test)))%>%
  filter(p.value<0.05)%>%
  rownames_to_column(., var = 'Variable')
View(var_nonnorm)

no_norm_df<-as.numeric(no_norm_df)


  # Box-Cox transformation

      # Create a copy of final_df selecting the non-normal variables
transforming_df <- final_df[, var_nonnorm$Variable]
summary(transforming_df)
      # Apply boxcox transformation to each variable
      # In orden to trying again, re-run the 'join tables' section

for (var in names(transforming_df)){
  print(var)
  x<-transforming_df[[var]]
  x[x == 0] <- 1e-15 # Error fix replacing zeros for almost sero values
  bc<-boxcox(lm(x ~ 1)) # MASS::boxcox object
  lambda <-bc$x[which.max(bc$y)] #Estimating lambda
  if(lambda!=0){ #Applying transformation depending on lambda
    x_trans<-(x ^ lambda - 1) / lambda
  }else{
    x_trans<-log(x)
  }
  final_df[[var]]<-x_trans # Replace with new values
  print('Success')
}

transf_df<-final_df


  # 47 Variables are still following non-normal distribution after box-cox transf

num_transf_df<-transf_df%>%select_if(is.numeric)
normality_test_transf<-sapply(num_transf_df, LillieTest)
var_nonnorm_transf<-as.data.frame((t(normality_test_transf)))%>%
  filter(p.value<0.05)%>%
  rownames_to_column(., var = 'Variable')
View(var_nonnorm_transf)


plot_data2<-plot_data%>%
  filter(., variable %in% var_nonnorm_transf$Variable)
ggplot(plot_data2, aes(x = value)) +
  geom_density(fill = "blue", alpha = 0.5) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

  # Save transformed and original datasets

write_tsv(final_df, 'final_df.txt')
write_tsv(transf_df, 'transf_df.txt')


# 3 NA in TRBV6_3 and 6 in TRBV4_3. Two samples with NA in both 
##As it is the result of a frequency calculation, 
##the appearance of NA in allelic use is due to the absence of that allele in 
##that sample, so it is replaced by 0 in data loading section.


# Scaling
# Due to the number of non-normally variables, we apply a robust transformation
# SPT.AGE_ASV has been removed due to the presence of inf an Nan values after scaling

final_scale_df <- final_df %>%
  mutate_if(is.numeric, function(x) (x - median(x)) / IQR(x))%>%
  dplyr::select(-SPT.AGE_ASV)

is.numeric(final_scale_df$MX1_rs469390)
is.numeric(final_scale_df$fever)


transf_scale_df<- transf_df %>%
  mutate_if(is.numeric, function(x) (x - median(x)) / IQR(x))%>%
  dplyr::select(-SPT.AGE_ASV)

  # Save scaled dataset

write_tsv(final_scale_df, 'final_scale_df.txt')
write_tsv(transf_scale_df, 'transf_scale_df.txt')

 



