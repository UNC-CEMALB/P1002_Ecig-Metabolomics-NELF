rm(list=ls())

#load packages
library(pheatmap)
library(RColorBrewer)
library(missForest)
library(tidyverse)
library(reshape)
library(janitor)
library(factoextra)
library(gridExtra)
library(ggplotify)
library(imputeLCMD)
library(xlsx)


setwd("/Cytokine_Imputation")


#load cytokine data
cyt_original <- read_csv("input/CytokineData_61Subjects_080421.csv")

#remove unnecessary metadata cols
cyt <- cyt_original %>% select(!c(SampleID, Gender, ExposureGroup)) %>%
  filter(is.na(ID)==FALSE) %>%
  column_to_rownames("ID")

#There are two IL8 cytokines specifying two ranges of detection. We only want to consider the upper range, thus the lower range is dropped
cyt <- as.data.frame(t(cyt))
cyt <- cyt %>% rownames_to_column("cytokine")
cyt[9,1] <- "IL8_drop"
cyt[23,1] <- "IL8"
cyt <- cyt %>% filter(cytokine!="IL8_drop") %>% column_to_rownames("cytokine")

#List of cytokines that should be considered for analysis
cyt_list <- rownames(cyt)

#Convert cytokine to binary (if the value is missing 0, if the value is present 1)
cyt_missing <- cyt %>% mutate_all(~replace(.,!is.na(.), 1)) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate_all(function(x) as.numeric(as.character(x)))


pheatmap(cyt_missing,cluster_rows=FALSE, cluster_cols = FALSE,
         angle_col = 45, fontsize_col = 7, legend = FALSE, cellheight = 10, cellwidth = 10, color = c("grey44", "indianred"),
         breaks = c(0, .5, 1), main="Missing Cytokines",
         filename = "figures/missing_cytokines_hm.png",width = 12,height = 8)


#Identify cytokines that are detected in at least 50% of samples
feature_to_keep <- cyt_missing %>% adorn_totals(where = "col",...=colnames(cyt_missing)) %>% 
  filter(Total>=ncol(cyt)*.5)

#Filter the cytokine data for the cytokines that should be retained based on previous identification step
cyt_filt <- cyt %>% rownames_to_column("cytokine") %>%
  filter(cytokine %in% rownames(feature_to_keep)) %>%
  column_to_rownames("cytokine")


#Boxplot of non-missing raw data for the molecular features of interest
df_long_cyt <- cyt_filt %>% 
  pivot_longer(everything(), names_to = "ID", values_to="values") %>% 
  mutate(group=if_else(grepl("NS",ID),"NS","ECig"))


ggplot(df_long_cyt, aes(y=log10(values),x=ID, col=group, order=group)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Raw")+
  ylab("Concentration  (log10)")

ggsave("figures/full_dist_boxplots/raw_data_boxplot.png", width = 12, height = 8)


#Prepare data for imputation per feature (cytokine)
cyt_filt_t <- t(cyt_filt)
cyt_filt_t <- apply(cyt_filt_t, 2, as.numeric)

set.seed(17)

################################################### Random Forest Imputation ####################################################
#impute
missForest_imp <- missForest(cyt_filt_t)

#get out data
missForest_res <- missForest_imp$ximp
rownames(missForest_res) <- colnames(cyt_filt)
missForest_df <- as.data.frame(missForest_res)

#write out data
write.csv(missForest_df, "output/RF_imputation_cytokine_data.csv")

#make another boxplot for the imputed data
bp_imp_missForest <- as.data.frame(t(missForest_df))

df_long_imp_missForest <- bp_imp_missForest %>% 
  pivot_longer(everything(), names_to = "ID", values_to="values") %>% 
  mutate(group=if_else(grepl("NS",ID),"NS","ECig"))


ggplot(df_long_imp_missForest, aes(y=log10(values), x=ID, col=group, order=group)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("RF Imputed")+
  ylab("Concentration (log10)")

ggsave("figures/full_dist_boxplots/RF_imputation_boxplot.png", width = 12, height = 8)


################################################### LOD/sqrt2 Global Imputation ####################################################

#prep dataframe
LODsqrt2_glob_imp <- as.data.frame(cyt_filt_t)
rownames(LODsqrt2_glob_imp) <- colnames(cyt_filt) 

glob_min <- min(LODsqrt2_glob_imp, na.rm = TRUE)

#impute
LODsqrt2_glob_df <- LODsqrt2_glob_imp %>% mutate_if(is.numeric, function(x) ifelse(is.na(x), glob_min/sqrt(2),x))

#write out data
write.csv(LODsqrt2_glob_df, "output/global_LODsqrt2_imputation_cytokine_data.csv")

#make another boxplot for the imputed data
bp_imp_LODsqrt2_glob <- as.data.frame(t(LODsqrt2_glob_df))

df_long_imp_LODsqrt2_glob <- bp_imp_LODsqrt2_glob %>% 
  pivot_longer(everything(), names_to = "ID", values_to="values") %>% 
  mutate(group=if_else(grepl("NS",ID),"NS","ECig"))


ggplot(df_long_imp_LODsqrt2_glob, aes(y=log10(values), x=ID, col=group, order=group)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("LOD/sqrt(2) Globally Imputed")+
  ylab("Concentration (log10)")

ggsave("figures/full_dist_boxplots/global_LODsqrt2_imputation_boxplot.png", width = 12, height = 8)


################################################### LOD/sqrt2 Per-Feature Imputation ####################################################

#prep dataframe
LODsqrt2_feat_imp <- as.data.frame(cyt_filt_t)
rownames(LODsqrt2_feat_imp) <- colnames(cyt_filt) 

#impute
LODsqrt2_feat_df <- LODsqrt2_feat_imp %>% mutate_if(is.numeric, function(x) ifelse(is.na(x), min(x, na.rm = T)/sqrt(2),x))

#write out data
write.csv(LODsqrt2_feat_df, "output/per_featureLODsqrt2_imputation_cytokine_data.csv")

#make another boxplot for the imputed data
bp_imp_LODsqrt2_feat <- as.data.frame(t(LODsqrt2_feat_df))

df_long_imp_LODsqrt2_feat <- bp_imp_LODsqrt2_feat %>% 
  pivot_longer(everything(), names_to = "ID", values_to="values") %>% 
  mutate(group=if_else(grepl("NS",ID),"NS","ECig"))

ggplot(df_long_imp_LODsqrt2_feat, aes(y=log10(values), x=ID, col=group, order=group)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("LOD/sqrt(2) Per-Feature Imputed")+
  ylab("Concentration (log10)")

ggsave("figures/full_dist_boxplots/per_feature_LODsqrt2_imputation_boxplot.png", width = 12, height = 8)



################################################### QRILC Imputation ####################################################

#prep dataframe
QRILC_prep <- as.data.frame(cyt_filt_t)
rownames(QRILC_prep) <- colnames(cyt_filt) 
QRILC_prep <- log10(QRILC_prep)

#Make sure that the distribution is roughly normal
viz <- QRILC_prep %>% rownames_to_column("cyt") %>% pivot_longer(!cyt, names_to = "samp", values_to = "conc" ) %>% drop_na()
ggplot(viz, aes(conc))+geom_histogram(bins=50)


#Impute
QRILC_imp <- impute.QRILC(QRILC_prep, tune.sigma = 1)
QRILC_df <- as.data.frame(QRILC_imp[1])

#Write out the imputations
write.csv(QRILC_df, "output/QRILC_imputation_cytokine_data.csv")

#make another boxplot for the imputed data
bp_imp_QRILC <- as.data.frame(t(QRILC_df))

df_long_imp_QRILC <- bp_imp_QRILC %>% 
  pivot_longer(everything(), names_to = "ID", values_to="values") %>% 
  mutate(group=if_else(grepl("NS",ID),"NS","ECig"))

ggplot(df_long_imp_QRILC, aes(y=values, x=ID, col=group, order=group)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(limits = c(-2.5,5))+
  ggtitle("QRILC Imputed")+
  ylab("Concentration (log10)")

ggsave("figures/full_dist_boxplots/QRILC_imputation_boxplot.png", width = 12, height = 8)


#Convert the now imputed data back to the original scale
QRILC_back_df <- QRILC_df %>% mutate_all(., function(x) 10^x)

#Write out the back converted imputations
write.csv(QRILC_back_df, "output/back_convert_QRILC_imputation_cytokine_data.csv")

#make another boxplot for the back converted imputed data
bp_imp_QRILC_back <- as.data.frame(t(QRILC_back_df))

df_long_imp_QRILC_back <- bp_imp_QRILC_back %>% 
  pivot_longer(everything(), names_to = "ID", values_to="values") %>% 
  mutate(group=if_else(grepl("NS",ID),"NS","ECig"))

ggplot(df_long_imp_QRILC_back, aes(y=log10(values), x=ID, col=group, order=group)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("QRILC Imputed- Back Converted")+
  ylab("Concentration (log10)")

ggsave("figures/full_dist_boxplots/back_convert_QRILC_imputation_boxplot.png", width = 12, height = 8)


###########################################

# Not all features were missing values. These are the cytokines that were missing values that passed necessary filters to be kept
imputed_features <- c("GMCSF","IL17","IL5","Eotaxin3","MCP4")

# Reference data frame of original values of just features that were imputed
comp_df <- as.data.frame(t(cyt_filt))
comp_df <- comp_df %>% select(all_of(imputed_features))

# Format RF results
missForrest_filt <- missForest_df %>%
  select(all_of(imputed_features)) %>%
  rownames_to_column("sample") %>%
  pivot_longer(!sample, names_to = "cyt",values_to = "RF") %>% 
  mutate(merge_col=paste0(sample,"_",cyt)) %>% 
  select(!c("sample","cyt"))

# Format Global LOD/sqrt(2) results
glob_LODsqrt2_filt <- LODsqrt2_glob_df %>%
  select(all_of(imputed_features)) %>%
  rownames_to_column("sample") %>%
  pivot_longer(!sample, names_to = "cyt",values_to = "Global_LOD_sqrt2") %>% 
  mutate(merge_col=paste0(sample,"_",cyt)) %>% 
  select(!c("sample","cyt"))

# Format Per-Feature LOD/sqrt(2) results
per_feat_LODsqrt2_filt <- LODsqrt2_feat_df %>%
  select(all_of(imputed_features)) %>%
  rownames_to_column("sample") %>%
  pivot_longer(!sample, names_to = "cyt",values_to = "Per_feature_LOD_sqrt2") %>% 
  mutate(merge_col=paste0(sample,"_",cyt)) %>% 
  select(!c("sample","cyt"))

# Format QRILC results
QRILC_filt <- QRILC_df %>%
  select(all_of(imputed_features)) %>%
  rownames_to_column("sample") %>%
  pivot_longer(!sample, names_to = "cyt",values_to = "QRILC") %>% 
  mutate(merge_col=paste0(sample,"_",cyt)) %>% 
  select(!c("sample","cyt"))

# Format back-convert QRILC results
QRILC_back_filt <- QRILC_back_df %>%
  select(all_of(imputed_features)) %>%
  rownames_to_column("sample") %>%
  pivot_longer(!sample, names_to = "cyt",values_to = "QRILC_back") %>% 
  mutate(merge_col=paste0(sample,"_",cyt)) %>% 
  select(!c("sample","cyt"))

# Merge all results and format
merged_results <- merge(missForrest_filt, glob_LODsqrt2_filt, by="merge_col")
merged_results <- merge(merged_results, per_feat_LODsqrt2_filt, by="merge_col")
merged_results <- merge(merged_results, QRILC_filt, by="merge_col")
merged_results <- merge(merged_results, QRILC_back_filt, by="merge_col")

merged_results <- merged_results %>% separate(merge_col, into = c("sample", "cyt"),sep = "_(?=[^_]+$)") %>% pivot_longer(!c("sample","cyt"), names_to = "imp_meth",values_to = "conc")

for(i in 1:length(imputed_features)){
  f <- imputed_features[i]
  na_samps <- comp_df %>% select(all_of(f)) %>% filter(is.na(.)) %>% rownames()
  temp <- merged_results %>% filter(cyt==f & sample %in% all_of(na_samps))
  five_num_sum <- temp %>% group_by(imp_meth) %>% summarise(min=min(conc), q1=quantile(conc,0.25), med=median(conc), q3=quantile(conc, 0.75), max=max(conc)) %>% column_to_rownames("imp_meth")
  write.xlsx(five_num_sum, file="output/Imputation_Summary_Statistics.xlsx", sheetName=f, append=TRUE)
  bp <- ggplot(temp, aes(y=conc, x=imp_meth, fill=imp_meth)) + geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1)) + 
    ggtitle(paste0("Distribution of Imputed Values for ",f," (n=",length(na_samps),")"))+
    ylab("imputed concentration")+
    xlab("imputation method")
  ggsave(filename = paste0("figures/feature_imputation_boxplots/",f,"_imputation_boxplot.png"), plot = bp, width=8, height = 6, units = "in")
}
