library(tidyverse)
library(readxl)

total_pred_ct1_1v2<-read.delim("C:/Users/Usuario/TCRMatch-1.1.2/data/ct1_1v2.tsv", 
                               header=F)
total_pred_ct2_1v2<-read.delim("C:/Users/Usuario/TCRMatch-1.1.2/data/ct2_1v2.tsv", 
                               header=F)
total_pred_ct1_1v3<-read.delim("C:/Users/Usuario/TCRMatch-1.1.2/data/ct1_1v3.tsv", 
                               header=F)
total_pred_ct3_1v3<-read.delim("C:/Users/Usuario/TCRMatch-1.1.2/data/ct3_1v3.tsv", 
                               header=F)

total_pred_ct1_1v2<-unique(total_pred_ct1_1v2$V1)
total_pred_ct2_1v2<-unique(total_pred_ct2_1v2$V1)
total_pred_ct1_1v3<-unique(total_pred_ct1_1v3$V1)
total_pred_ct3_1v3<-unique(total_pred_ct3_1v3$V1)


total_pred<-list(total_pred_ct1_1v2$V1,
                 total_pred_ct2_1v2$V1,
                 total_pred_ct1_1v3$V1,
                 total_pred_ct3_1v3$V1)

iedb_pred_ct1_1v2<-read.delim("C:/Users/Usuario/TCRMatch-1.1.2/output_file_ct1_1v2.txt", 
                               header=T)
iedb_pred_ct2_1v2<-read.delim("C:/Users/Usuario/TCRMatch-1.1.2/output_file_ct2_1v2.txt", 
                              header=T)
iedb_pred_ct3_1v3<-read.delim("C:/Users/Usuario/TCRMatch-1.1.2/output_file_ct3_1v3.txt", 
                              header=T)
iedb_pred_ct1_1v3<-read.delim("C:/Users/Usuario/TCRMatch-1.1.2/output_file_ct1_1v3.txt", 
                              header=T)


iedb_pred_ct1_1v2<-unique(iedb_pred_ct1_1v2$trimmed_input_sequence)
iedb_pred_ct2_1v2<-unique(iedb_pred_ct2_1v2$trimmed_input_sequence) 
iedb_pred_ct1_1v3<-unique(iedb_pred_ct1_1v3$trimmed_input_sequence)
iedb_pred_ct3_1v3<-unique(iedb_pred_ct3_1v3$trimmed_input_sequence)

iedb_pred<-list(iedb_pred_ct1_1v2$trimmed_input_sequence, 
             iedb_pred_ct2_1v2$trimmed_input_sequence, 
             iedb_pred_ct1_1v3$trimmed_input_sequence,
             iedb_pred_ct3_1v3$trimmed_input_sequence)


covid_pred_ct1_1v2<-read_excel("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/results/tables/TCRmatch_097_IEDB_1v2.xlsx",
                               sheet = 3)
covid_pred_ct2_1v2<-read_excel("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/results/tables/TCRmatch_097_IEDB_1v2.xlsx",
                               sheet = 4)
covid_pred_ct1_1v3<-read_excel("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/results/tables/TCRmatch_097_IEDB_1v3.xlsx",
                               sheet = 4)
covid_pred_ct3_1v3<-read_excel("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/results/tables/TCRmatch_097_IEDB_1v3.xlsx",
                               sheet = 3)


covid_pred_ct1_1v2<-unique(covid_pred_ct1_1v2$trimmed_input_sequence)
covid_pred_ct2_1v2<-unique(covid_pred_ct2_1v2$trimmed_input_sequence)
covid_pred_ct1_1v3<-unique(covid_pred_ct1_1v3$trimmed_input_sequence)
covid_pred_ct3_1v3<-unique(covid_pred_ct3_1v3$trimmed_input_sequence)


covid_pred<-list(Cluster.1_1v2=covid_pred_ct1_1v2, 
              Cluster.2_1v2=covid_pred_ct2_1v2,
              Cluster.1_1v3=covid_pred_ct1_1v3,
              Cluster.3_1v3=covid_pred_ct3_1v3)

 

library(scales) 
  
`%notin%` <- Negate(`%in%`)

seq_df_1_1v2 <- data.frame(CDR3b_sequences = total_pred_ct1_1v2) %>%
  mutate(Class = case_when(
    CDR3b_sequences %in% iedb_pred_ct1_1v2 ~ 'SARS-CoV2 match',
    CDR3b_sequences %notin% iedb_pred_ct1_1v2 & CDR3b_sequences %in% iedb_pred_ct1_1v2~ 'iedb_positive',
    CDR3b_sequences %notin% c(iedb_pred_ct1_1v2, iedb_pred_ct1_1v2) ~ 'Unmatched sequences'
  )) %>%
  group_by(Class) %>%
  summarise(Count = n())

seq_df_1_1v2$fraction = seq_df_1_1v2$Count / sum(seq_df_1_1v2$Count)
seq_df_1_1v2$ymax = cumsum(seq_df_1_1v2$fraction)
seq_df_1_1v2$ymin = c(0, head(seq_df_1_1v2$ymax, n=-1))
seq_df_1_1v2$labelPosition <- (seq_df_1_1v2$ymax + seq_df_1_1v2$ymin) / 2
seq_df_1_1v2$label <- paste0(seq_df_1_1v2$Class, "\n Value: ", seq_df_1_1v2$Count)
ct112<-ggplot(seq_df_1_1v2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Class)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  coord_polar(theta="y") + 
  xlim(c(2, 4))+
  theme_void() +
  theme(legend.position = "none")+
  scale_fill_manual(values = c('#C8A066', '#ACE6DE'))


seq_df_2_1v2 <- data.frame(CDR3b_sequences = total_pred_ct2_1v2) %>%
  mutate(Class = case_when(
    CDR3b_sequences %in% iedb_pred_ct2_1v2 ~ 'SARS-CoV2 match',
    CDR3b_sequences %notin% iedb_pred_ct2_1v2 & CDR3b_sequences %in% iedb_pred_ct2_1v2~ 'iedb_positive',
    CDR3b_sequences %notin% c(iedb_pred_ct2_1v2, iedb_pred_ct2_1v2) ~ 'Unmatched sequences'
  )) %>%
  group_by(Class) %>%
  summarise(Count = n())


seq_df_2_1v2$fraction = seq_df_2_1v2$Count / sum(seq_df_2_1v2$Count)
seq_df_2_1v2$ymax = cumsum(seq_df_2_1v2$fraction)
seq_df_2_1v2$ymin = c(0, head(seq_df_2_1v2$ymax, n=-1))
seq_df_2_1v2$labelPosition <- (seq_df_2_1v2$ymax + seq_df_2_1v2$ymin) / 2
seq_df_2_1v2$label <- paste0(seq_df_2_1v2$Class, "\n Value: ", seq_df_2_1v2$Count)
ct212<-ggplot(seq_df_2_1v2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Class)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  coord_polar(theta="y") + 
  xlim(c(2, 4))+
  theme_void() +
  theme(legend.position = "none")+
  scale_fill_manual(values = c('#C8A066', '#ACE6DE'))







seq_df_1_1v3 <- data.frame(CDR3b_sequences = total_pred_ct1_1v3) %>%
  mutate(Class = case_when(
    CDR3b_sequences %in% iedb_pred_ct1_1v3 ~ 'SARS-CoV2 match',
    CDR3b_sequences %notin% iedb_pred_ct1_1v3 & CDR3b_sequences %in% iedb_pred_ct1_1v3~ 'iedb_positive',
    CDR3b_sequences %notin% c(iedb_pred_ct1_1v3, iedb_pred_ct1_1v3) ~ 'Unmatched sequences'
  )) %>%
  group_by(Class) %>%
  summarise(Count = n())

seq_df_1_1v3$fraction = seq_df_1_1v3$Count / sum(seq_df_1_1v3$Count)
seq_df_1_1v3$ymax = cumsum(seq_df_1_1v3$fraction)
seq_df_1_1v3$ymin = c(0, head(seq_df_1_1v3$ymax, n=-1))
seq_df_1_1v3$labelPosition <- (seq_df_1_1v3$ymax + seq_df_1_1v3$ymin) / 2
seq_df_1_1v3$label <- paste0(seq_df_1_1v3$Class, "\n Value: ", seq_df_1_1v3$Count)
ct113<-ggplot(seq_df_1_1v3, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Class)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  coord_polar(theta="y") + 
  xlim(c(2, 4))+
  theme_void() +
  theme(legend.position = "none")+
  scale_fill_manual(values = c('#C8A066', '#ACE6DE'))



seq_df_3_1v3 <- data.frame(CDR3b_sequences = total_pred_ct3_1v3) %>%
  mutate(Class = case_when(
    CDR3b_sequences %in% iedb_pred_ct3_1v3 ~ 'SARS-CoV2 match',
    CDR3b_sequences %notin% iedb_pred_ct3_1v3 & CDR3b_sequences %in% iedb_pred_ct3_1v3~ 'iedb_positive',
    CDR3b_sequences %notin% c(iedb_pred_ct3_1v3, iedb_pred_ct3_1v3) ~ 'Unmatched sequences'
  )) %>%
  group_by(Class) %>%
  summarise(Count = n())

seq_df_3_1v3$fraction = seq_df_1_1v2$Count / sum(seq_df_3_1v3$Count)
seq_df_3_1v3$ymax = cumsum(seq_df_3_1v3$fraction)
seq_df_3_1v3$ymin = c(0, head(seq_df_3_1v3$ymax, n=-1))
seq_df_3_1v3$labelPosition <- (seq_df_3_1v3$ymax + seq_df_3_1v3$ymin) / 2
seq_df_3_1v3$label <- paste0(seq_df_3_1v3$Class, "\n Value: ", seq_df_3_1v3$Count)
ct313<-ggplot(seq_df_3_1v3, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Class)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  coord_polar(theta="y") + 
  xlim(c(2, 4))+
  theme_void() +
  theme(legend.position = "none")+
  scale_fill_manual(values = c('#C8A066', '#ACE6DE'))



library(UpSetR)
library(ComplexHeatmap)



covid_pred<-list(Cluster.1_1v2=covid_pred_ct1_1v2, 
                 Cluster.2=covid_pred_ct2_1v2,
                 Cluster.1_1v3=covid_pred_ct1_1v3,
                 Cluster.3=covid_pred_ct3_1v3)

total_pred<-list(Cluster.1_1v2=total_pred_ct1_1v2,
                 Cluster.2=total_pred_ct2_1v2,
                 Cluster.1_1v3=total_pred_ct1_1v3,
                 Cluster.3=total_pred_ct3_1v3)


m_covid<-make_comb_mat(covid_pred, mode='intersect')
m_pred<-make_comb_mat(total_pred, mode='intersect')



UpSet(
  m_pred,
  top_annotation = upset_top_annotation(m_pred, add_numbers = TRUE),
  set_order = c("Cluster.1_1v2", "Cluster.2", "Cluster.1_1v3", "Cluster.3"),
  comb_order = order(comb_size(m_pred)),
  pt_size = unit(6, "mm"),
  lwd = 3,
  right_annotation = rowAnnotation(
    "Set size" = anno_barplot(set_size(m_pred),
                              border = FALSE,
                              gp = gpar(fill = "black"),
                              width = unit(2, "cm")),
    comparisons = c('1v2', '1v2', '1v3', '1v3')))



UpSet(
  m_covid,
  top_annotation = upset_top_annotation(m_covid, add_numbers = TRUE),
  set_order = c("Cluster.1_1v2", "Cluster.2", "Cluster.1_1v3", "Cluster.3"),
  comb_order = order(comb_size(m_covid)),
  pt_size = unit(6, "mm"),
  lwd = 3,
  right_annotation = rowAnnotation(
    "Set size" = anno_barplot(set_size(m_covid),
                              border = FALSE,
                              gp = gpar(fill = "black"),
                              width = unit(2, "cm")),
  comparisons = c('1v2', '1v2', '1v3', '1v3')))



common_ct1_covid_sequences<-unique(intersect(covid_pred_ct1_1v2,
                                      covid_pred_ct1_1v3))

covid_ct1_1v2_df<-read_excel("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/results/tables/TCRmatch_097_IEDB_1v2.xlsx",
                                                     sheet = 3)%>%
  select(trimmed_input_sequence, antigen, score)%>%
  separate_rows(antigen, sep = ",")%>%
  group_by(antigen) %>%
  summarise(Count = n())%>% 
  mutate(cluster = 'Cluster 1') %>%
  mutate(comparison = '1v2') %>% mutate(., perc=(Count/sum(Count)*100))

covid_ct2_1v2_df<-read_excel("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/results/tables/TCRmatch_097_IEDB_1v2.xlsx",
                       sheet = 4)%>%
  select(trimmed_input_sequence, antigen, score)%>%
  separate_rows(antigen, sep = ",")%>%
  group_by(antigen) %>%
  summarise(Count = n()) %>% 
  mutate(cluster = 'Cluster 2') %>%
  mutate(comparison = '1v2') %>% mutate(., perc=(Count/sum(Count)*100))

covid_ct1_1v3_df<-read_excel("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/results/tables/TCRmatch_097_IEDB_1v3.xlsx",
                                 sheet = 4)%>%
  select(trimmed_input_sequence, antigen, score)%>%
  separate_rows(antigen, sep = ",")%>%
  group_by(antigen) %>%
  summarise(Count = n()) %>% 
  mutate(cluster = 'Cluster 1') %>%
  mutate(comparison = '1v3') %>% mutate(., perc=(Count/sum(Count)*100))

covid_ct3_1v3_df<-read_excel("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/results/tables/TCRmatch_097_IEDB_1v3.xlsx",
                                 sheet = 3)%>%
  select(trimmed_input_sequence, antigen, score)%>%
  separate_rows(antigen, sep = ",")%>%
  group_by(antigen) %>%
  summarise(Count = n()) %>% 
  mutate(cluster = 'Cluster 3') %>%
  mutate(comparison = '1v3')%>% mutate(., perc=(Count/sum(Count)*100))


antigens_df<-bind_rows(covid_ct1_1v2_df, 
                   covid_ct2_1v2_df, 
                   covid_ct1_1v3_df,
                   covid_ct3_1v3_df)

antigens_df$antigen[antigens_df$antigen=='ORF7b']<-'ORF7b protein'
antigens_df$antigen[antigens_df$antigen=='orf1ab polyprotein']<-'ORF1AB polyprotein'
antigens_df$antigen[antigens_df$antigen=='non-structural protein NS4b']<-'nonstructural protein NS4b'
antigens_df$antigen[antigens_df$antigen=='Matrix protein 1']<-'matrix protein 1'
antigens_df$antigen[antigens_df$antigen=='Spike glycoprotein precursor']<-'spike glycoprotein precursor'
ggplot(data=antigens_df, aes(x = antigen, y=perc, fill=cluster))+ 
  geom_col(color="black")+
  coord_flip()+
  ylab('Perc (%)')+
  xlab('Antigen')+
  theme_minimal()+
  facet_grid(cols=vars(comparison))+
  theme(legend.text = element_text(size=16))+
  theme(axis.text.y = element_text(size=12, colour = 'black'))+
  theme(axis.title.y = element_text(size=22))+
  theme(axis.text.x = element_text(size=16, colour = 'black'))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.title.x = element_text(size=22))+
  theme(strip.text = element_text(size = 16, face = "bold"))+
  theme(legend.title = element_blank())





