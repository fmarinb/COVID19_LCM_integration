###### Functions to check the normality and homocedasticity####################
library(nortest)
library(easystats)
library(effectsize)
library(rstatix)

#### KRUSKALL-WALLIS TEST ######################################################

kw_test<-function(df, test_groups, test_metrics){
  df_total<-df
  
  groups<-df_total%>%
    select(., all_of(test_groups))%>%
    colnames()
  
  div_names<-df_total%>%
    select(., all_of(test_metrics))%>%
    colnames()
  
  p_val<-c()
  ks_val<-c()
  res_v<-c()
  res_ci_l<-c()
  res_ci_h<-c()
  res_int<-c()
  post_hoc<-c()
  group<-c()
  div<-c()
  dt<-c()
  
  for (j in 1:length(groups)){ 
    for (i in 1:length(div_names)){
      print(paste('Kruskall-Wallis test between' ,
                  groups[j], 'and', div_names[i]))
      
      # Create the table of values
      test_table<-df_total%>%
        select(groups[j], div_names[i])%>%
        as.data.frame()%>%
        na.omit()
      
      # Agregate them to show the median value by group  
      test_table_param<-aggregate(test_table[,2], 
                                  by = list(group=test_table[,1]), 
                                  median)
      
      colnames(test_table_param)<-c('Group', 'Median')
      print(test_table_param)
      
      # Perform the Kruskall-Wallis test
      ks<-kruskal.test(test_table[,2]~test_table[,1])
      p_val<-append(p_val, ks$p.value)
      ks_val<-append(ks_val, ks$statistic)
      
      # Calculate size effect
      res<-rank_epsilon_squared(test_table[,2]~test_table[,1], data = test_table)
      res_v<-append(res_v, res$rank_epsilon_squared)
      res_ci_l<-append(res_ci_l, res$CI_low)
      res_ci_h<-append(res_ci_h, res$CI_high)
      
      # The interpretation is done following the criteria for it
      res_int<-append(res_int, 
                      interpret_epsilon_squared(res$rank_epsilon_squared))
      
      # Post hoc analysis with Dunn`s test and Bonferroni-Holm correction
      print("Dunn's test for post-hoc analysis")
      values<-test_table[,2]
      gr<-test_table[,1]
      df_temp<-data.frame(values=values, gr=gr)
      pt<-dunn_test(data=df_temp, values ~ gr,  p.adjust.method = 'BH')
      pt<-pt%>%mutate(., Metric=rep(div_names[i], 3), .before=group1)
      post_hoc<-bind_rows(post_hoc, pt, .id = NULL)
      
      
      group<-append(group, groups[j])
      div<-append(div, div_names[i])
      
    }
  }
  results_kw<-tibble(Group=group,
                     Index=div,
                     Kruskall_Wallis_p.val=p_val,
                     Statistic=ks_val,
                     `Rank_epsilon_squared(RES)`=res_v,
                     RES_95_low_CI=res_ci_l,
                     RES_d_95_high_CI=res_ci_h,
                     RES_meaning=res_int)
  
  results_posthoc<-compact(post_hoc)
  
  results<-list(Kruskall_Wallis=results_kw, Post_hoc=results_posthoc)
  
  return(results)
} 

selected_var<- sub("\\.\\.\\.", "", selected_var)
selected_var<- gsub("normalizedShannonWienerIndex_mean", "SW", selected_var)
selected_var<- gsub("normalizedShannonWienerIndex_mean", "SW", selected_var)
names(final_df)[names(final_df) == 'YSSGE_4_22...'] <- 'YSSGE_4_22'


final_df_lcm<-final_df%>%select(., sample, all_of(selected_var))
final_df_lcm<-final_df_lcm%>%
  full_join(., selected_df3%>%select(., sample, cluster, Severity, sev_55), 'sample')

selected_lcm70_test<-kw_test(df = final_df_lcm, test_groups = 'cluster', test_metrics = colnames(final_df_lcm[, -c(1, 71:74)]))

write.csv(selected_lcm70_test$Kruskall_Wallis, 'sel_70lcm_kw.csv', quote = F)
write.csv(selected_lcm70_test$Post_hoc, 'sel_70lcm_ph.csv', quote=F)



