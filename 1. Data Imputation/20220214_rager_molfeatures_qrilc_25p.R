
{
#load
{
  #
  rm(list=ls())
  ## load packages, data
  #Set seed for reproducability
  set.seed(42)
  library(tidyverse)
  
  library(vsn)
  library(missForest)
  library(imputeLCMD)
  
  library(pheatmap)
  library(GGally)
  library(ggh4x)
  
  
  dir_file <- "C:/Users/Nyssa/Documents/unc/rager/"
  rec_date <- "20220214"
  date <- "20220214"
}


#filter missing id & filter
{
  
  df <- as_tibble(read.csv(paste0(dir_file,"20210803_molecular_feature_ECIG_study_NoGapFillNoScaling_Dmp3872.csv")))
  
  
  ## name rows, remove ctrls (NA_), transpose, and shape as tibble
  mf <- df %>% column_to_rownames("MF_Number") %>% 
    select(!contains("NA_")) %>% t %>% as_tibble()
  
  ## shape & vis heatmap missing data
  ## gets total count (sums of 0 or 1s) for coverage percent (next step)
  mf_missing <- mf %>% mutate_all(~replace(.,!is.na(.), 1)) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    mutate_all(function(x) as.numeric(as.character(x)))%>% mutate(across(.cols=everything(),sum)) %>% unique()
  
  #get a count of the total rows
  cou <- nrow(mf)
  mf_m3 <- mf_missing
  
  ## divide sum of each row (0-59) by total row count to get coverage %
      for (i in 1:ncol(mf_missing)){
        mf_m3[i] <- mf_missing[i]/cou
      }
  
  ## reshape df, setting mol feature name to 'molf'
  mf4 <- mf_m3%>% pivot_longer(names_to = "molf",cols=everything())
  ## set drop val (25% (+1%))
  mf_drop_val <- 0.26
  
  #hist(mf4$value,col="aliceblue")
  #abline(v = mf_drop_val,col = "mediumvioletred",lwd = 3,lty=2)
  
  
  ## drop values under filter
  mf5 <- mf4 %>% filter(value >= mf_drop_val)
  
  mf6 <- mf %>% select(cols=mf5$molf)
  colnames(mf6) <- mf5$molf
  
  ## announce dropped count
  n_dropped<- nrow(mf4)-nrow(mf5)
  print(paste0("At ",mf_drop_val*100,"+% missing, ",n_dropped," molecular features were dropped from further analysis"))
  
  
  ## select only those retained mol features, drop NA, reshape (cant use df composed above bc we dropped abundance values for 0 or 1s)
  mf7 <- df[df$MF_Number %in% mf5$molf, ]
  mf8 <- mf7%>% select(!contains("NA_"))
  mf9 <- mf8 %>% t %>% as_tibble()
  mf9$ID <- colnames(mf8)
  colnames(mf9)<-mf9[1,]
  mf10 <- mf9[-1,]
  names(mf10)[names(mf10) == 'MF_Number'] <- 'ID'
  
  ## split ID into useful units
  mf11 <- mf10 %>% 
    mutate(group=if_else(grepl("NS",ID),"0","1")) %>% separate(col = 'ID',into=c("treat","sample",sep="_"))
  
  ## mutate column types to become more useful
  mf_meanb <- mf11 %>% mutate_if(is.character,as.numeric)
  mf_meanc <- mf_meanb 
  mf_meand <- mf_meanc %>% select(-c("sample","group","treat","_")) %>% as.matrix() %>% justvsn()
  mf_before <- as_tibble( mf_meand) 
  mf_before$sample <- mf_meanc$sample
  mf_before$group <- mf_meanc$group
  #mf_before$treat <- mf_meanc$treat
  mf_before <- mf_before %>% select('sample','group',everything()) %>% 
    pivot_longer(names_to = "MF",values_to = "val",cols = 3:ncol(mf_before))
  
  mf_vsn_t <- t(mf_meand)
  mf_vsn_t <- apply(mf_vsn_t, 2, as.numeric)
}

  ################################################### QRILC Imputation ####################################################
  
{
  ## assign correct starting point (both justvsn'd, transpose transposed)
  mf_filt <- mf_meand
  
  mf_filt_t <-mf_vsn_t
  
  #prep dataframe
  QRILC_prep <- as.data.frame(mf_filt_t)
  rownames(QRILC_prep) <- colnames(mf_filt) 
  QRILC_prep <- log10(QRILC_prep)
  
  #Make sure that the distribution is roughly normal
  viz <- QRILC_prep %>% rownames_to_column("cyt") %>% pivot_longer(!cyt, names_to = "samp", values_to = "conc" ) %>% drop_na()

  #Impute
  QRILC_imp <- impute.QRILC(QRILC_prep, tune.sigma = 1)
  QRILC_df <- as.data.frame(QRILC_imp[1])
  
  #Write out the imputations
  #write.csv(QRILC_df, paste0(dir_name,"output/",date,"_QRILC_imputation_cytokine_data.csv"))
  
  #make another boxplot for the imputed data
  bp_imp_QRILC <- as.data.frame(t(QRILC_df))
  bp_imp_QRILC$group <- mf10$ID
  
  df_long_imp_QRILC <- bp_imp_QRILC %>% 
    pivot_longer(-c("group"), names_to = "ID", values_to="values") %>% 
    mutate(group=if_else(grepl("NS",group),"NS","ECig"))
  
  QRILC_back_df <- QRILC_df %>% mutate_all(., function(x) 10^x)
  
  #Write out the back converted imputations
  #write.csv(QRILC_back_df, paste0(dir_name,"output/",date,"_back_convert_QRILC_imputation_cytokine_data.csv"))
  
  #make another boxplot for the back converted imputed data
  bp_imp_QRILC_back <- as.data.frame(t(QRILC_back_df))
  bp_imp_QRILC_back$sample <- mf10$ID
}
  {
    
    ## reshape IDs again
  df_long_imp_QRILC_back <- bp_imp_QRILC_back %>% 
    pivot_longer(-c('sample'), names_to = "ID", values_to="values") %>% 
    mutate(group=if_else(grepl("NS",sample),"NS","ECig"))
  
  ## select rows and keep non-scaled only
  df2 <- df_long_imp_QRILC_back %>% select('sample','ID','values','group') %>% filter(values >= 10)
  df2$group <- as.character(df2$group)
  df2["group"][df2["group"] == "0"] <- "NS"
  df2["group"][df2["group"] == "1"] <- "ECig"
  df3 <- df2 %>% pivot_wider(names_from = 'ID',values_from = 'values')
  
  write_tsv(df3,file = paste0("C:/Users/Nyssa/Documents/unc/rager/",date,"_mf_qrilc_df_square2.tsv"))
  
}
print(done)
}

