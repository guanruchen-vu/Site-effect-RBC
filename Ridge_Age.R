
library(tidyverse)
library(vtable)
library(knitr) 
library(glmnet)
library(locfit)
library(table1)
library(boot)
library(ggplot2)
library(GGally)
library(ranger)
library(tuneRanger)
library(mlr)
library(table1)
library(readr)     
library(dplyr)      
library(Hmisc)      
library(kableExtra) 

# Setting the seed
set.seed(123)
# Loading in data
rbc <- read_csv("/media/disk2/RBC_version0.1/regional_GMV_all_sites_harmonized_covariates.csv")

# Reading in the structural data
func <- read_csv("/media/disk2/RBC_version0.1/FC_network17_harmonized_covariates.csv")

# Setting each dataset to only include patients that overlap 
rbc = rbc[rbc$participant_id %in% func$participant_id,]
func = func[func$participant_id %in% rbc$participant_id,]

# Sanity check
sum(rbc$participant_id == func$participant_id) == nrow(func)

# Merging datasets
cols_to_add <- setdiff(names(func), names(rbc))
both <- dplyr::inner_join(rbc, func[, c("participant_id", cols_to_add)], by = "participant_id")

test <- dplyr::inner_join(rbc, func, by = "participant_id")

# Loading in functions
# Loading in functions
source("/media/disk2/multivariateBWAS/code_for_review/build_funcs.R")
# Sourcing influence and correlation functions
source("/media/disk2/multivariateBWAS/code_for_review/influence_funcs.R")
source("/media/disk2/multivariateBWAS/code_for_review/corr_funcs.R")

# Loading in color palette
cbbPalette <- c("#000000","#E69F00","#56B4E9","#009E73",
                "#F0E442","#0072B2","#D55E00","#CC79A7")

Single_Split <- function(data_clean, target_site, site_col, outcome, typeDat, folds){
  target_data = data_clean %>% filter(.data[[site_col]] == target_site)
  other_data  = data_clean %>% filter(.data[[site_col]] != target_site)
  
  target_data$fold <- sample(rep(1:folds, length.out = nrow(target_data)))
  
  fold_local_os = c(); fold_down_os = c(); fold_full_os = c()
  fold_local_pearson = c(); fold_down_pearson = c(); fold_full_pearson = c()
  
  inf_local_all = c(); inf_down_all = c(); inf_full_all = c()
  preds_full_all  = c()
  y_real_full_all = c()
  subj_id_all     = c()
  
  prep_matrix <- function(train_df, test_df) {
    if(typeDat != "Baseline"){
      if(outcome == "age"){
        x    = data.matrix(train_df[, grepl(typeDat, colnames(train_df)) | colnames(train_df) %in% c("sex", "euler", "meanFD")])
        newx = data.matrix(test_df[,  grepl(typeDat, colnames(test_df))  | colnames(test_df)  %in% c("sex", "euler", "meanFD")])
      } else {
        x    = data.matrix(train_df[, grepl(typeDat, colnames(train_df)) | colnames(train_df) %in% c("age", "sex", "euler", "meanFD")])
        newx = data.matrix(test_df[,  grepl(typeDat, colnames(test_df))  | colnames(test_df)  %in% c("age", "sex", "euler", "meanFD")])
      }
    } else {
      if(outcome == "age"){
        x    = data.matrix(train_df[, colnames(train_df) %in% c("sex", "euler", "meanFD")])
        newx = data.matrix(test_df[,  colnames(test_df)  %in% c("sex", "euler", "meanFD")])
      } else {
        x    = data.matrix(train_df[, colnames(train_df) %in% c("age", "sex", "euler", "meanFD")])
        newx = data.matrix(test_df[,  colnames(test_df)  %in% c("age", "sex", "euler", "meanFD")])
      }
    }
    y      = data.matrix(train_df[, outcome])
    y_real = as.numeric(test_df[[outcome]])
    return(list(x = x, y = y, newx = newx, y_real = y_real))
  }
  
  # Fitting a ridge regression model
  lambda_seq <- c(
    seq(10, 1, -1),        # From 10 to 1 (step of 1)
    seq(0.9, 0.1, -0.1),   # From 0.9 to 0.1 (step of 0.1)
    seq(0.09, 0.01, -0.01),# From 0.09 to 0.01 (step of 0.01)
    seq(0.009, 0.001, -0.001), # From 0.009 to 0.001 (step of 0.001)
    seq(0.0009, 0.0001, -0.0001), # From 0.0009 to 0.0001 (step of 0.0001)
    seq(0.00009, 0.00001, -0.00001), # From 0.00009 to 0.00001
    seq(0.000009, 0.000001, -0.000001), # From 0.000009 to 0.000001
    seq(0.0000009, 0.0000001, -0.0000001), # From 0.0000009 to 0.0000001
    0
  )
  
  for(i in 1:folds){
    test        = target_data %>% filter(fold == i)
    train_local = target_data %>% filter(fold != i)
    N_train     = nrow(train_local)
    
    train_multi_full = bind_rows(train_local, other_data)
    
    site_counts <- train_multi_full %>%
      count(.data[[site_col]]) %>%
      mutate(prop = n / sum(n), target_n = floor(prop * N_train))
    remainder <- N_train - sum(site_counts$target_n)
    if(remainder > 0){
      add_idx <- sample(1:nrow(site_counts), remainder)
      site_counts$target_n[add_idx] <- site_counts$target_n[add_idx] + 1
    }
    
    train_multi_down <- train_multi_full %>%
      left_join(site_counts %>% select(all_of(site_col), target_n), by = site_col) %>%
      group_by(.data[[site_col]]) %>%
      filter(target_n[1] > 0) %>%
      slice(sample.int(n(), size = target_n[1], replace = FALSE)) %>%
      ungroup() %>%
      select(-target_n)
    
    run_model_get_OS <- function(train_df, test_df) {
      mats = prep_matrix(train_df, test_df)
      safe_nfolds = min(10, max(3, nrow(mats$x)))
      
      ridge <- tryCatch({
        cv.glmnet(mats$x, mats$y, alpha = 0, lambda = lambda_seq, nfolds = safe_nfolds)
      }, error = function(e) NULL)
      
      if(is.null(ridge)) return(list(
        os      = NA,
        pearson = NA,
        inf     = rep(NA, nrow(test_df)),
        preds   = rep(NA, nrow(test_df)),
        y_real  = mats$y_real
      ))
      
      preds       = as.vector(predict(ridge, newx = mats$newx, s = "lambda.min"))
      os_val      = cor_os(preds, mats$y_real)
      # 新增 Pearson estimator
      pearson_val = suppressWarnings(cor(preds, mats$y_real, method = "pearson")) 
      inf_vec     = if_logit_sq(preds, mats$y_real)
      
      return(list(os = os_val, pearson = pearson_val, inf = inf_vec, preds = preds, y_real = mats$y_real))
    }
    
    res_local = run_model_get_OS(train_local,     test)
    res_down  = run_model_get_OS(train_multi_down, test)
    res_full  = run_model_get_OS(train_multi_full, test)
    
    fold_local_os = c(fold_local_os, res_local$os)
    fold_down_os  = c(fold_down_os,  res_down$os)
    fold_full_os  = c(fold_full_os,  res_full$os)

    fold_local_pearson = c(fold_local_pearson, res_local$pearson)
    fold_down_pearson  = c(fold_down_pearson,  res_down$pearson)
    fold_full_pearson  = c(fold_full_pearson,  res_full$pearson)
    
    inf_local_all = c(inf_local_all, res_local$inf)
    inf_down_all  = c(inf_down_all,  res_down$inf)
    inf_full_all  = c(inf_full_all,  res_full$inf)
    
    preds_full_all  = c(preds_full_all,  res_full$preds)
    y_real_full_all = c(y_real_full_all, res_full$y_real)
    subj_id_all     = c(subj_id_all,     test$participant_id)
  }
  
  N_total   = nrow(target_data)
  var_local = var(inf_local_all, na.rm = TRUE) / N_total
  var_down  = var(inf_down_all,  na.rm = TRUE) / N_total
  var_full  = var(inf_full_all,  na.rm = TRUE) / N_total
  
  preds_df = data.frame(
    subj_id   = subj_id_all,
    y_real    = y_real_full_all,
    pred_full = preds_full_all,
    site      = target_site
  )
  
  return(list(
    metrics  = data.frame(
      Local_OS         = mean(fold_local_os, na.rm = TRUE),
      Multi_Down_OS    = mean(fold_down_os,  na.rm = TRUE),
      Multi_Full_OS    = mean(fold_full_os,  na.rm = TRUE),
      Local_Pearson    = mean(fold_local_pearson, na.rm = TRUE),
      Multi_Down_Pearson= mean(fold_down_pearson, na.rm = TRUE),
      Multi_Full_Pearson= mean(fold_full_pearson, na.rm = TRUE),
      Local_Var        = var_local,
      Multi_Down_Var   = var_down,
      Multi_Full_Var   = var_full
    ),
    preds_df = preds_df
  ))
}

Run_Site_Ablation_Splits <- function(data, target_site, site_col, outcome, typeDat, folds = 5, splits = 100){
  data_clean = data %>% drop_na(all_of(outcome), all_of(site_col))
  print(paste("=================================================="))
  print(paste("Starting Ablation for Site:", target_site))
  print(paste("Total Splits:", splits, "| Target N:", 
              nrow(data_clean %>% filter(.data[[site_col]] == target_site))))
  print(paste("=================================================="))
  
  all_splits_res = data.frame()
  all_preds_list = list()
  
  for(s in 1:splits){
    if(s %% 10 == 0) print(paste("  ... Processing split", s, "/", splits))
    
    split_out = Single_Split(data_clean, target_site, site_col, outcome, typeDat, folds)
    
    # metrics
    split_res           = split_out$metrics
    split_res$split_num = s
    all_splits_res      = rbind(all_splits_res, split_res)
    
    preds_temp           = split_out$preds_df
    preds_temp$split_num = s
    all_preds_list[[s]]  = preds_temp
  }
  
  all_preds_long = bind_rows(all_preds_list)
  
  median_preds = all_preds_long %>%
    group_by(subj_id, y_real, site) %>%
    summarise(pred_full_median = median(pred_full, na.rm = TRUE), .groups = "drop")
  
  median_res = all_splits_res %>%
    summarise(
      Target_Site               = target_site,
      Median_Local_OS           = median(Local_OS,         na.rm = TRUE),
      Median_Multi_Down_OS      = median(Multi_Down_OS,    na.rm = TRUE),
      Median_Multi_Full_OS      = median(Multi_Full_OS,    na.rm = TRUE),
      IQR_Local_OS              = IQR(Local_OS,            na.rm = TRUE),
      IQR_Multi_Down_OS         = IQR(Multi_Down_OS,       na.rm = TRUE),
      
      Median_Local_Pearson      = median(Local_Pearson,      na.rm = TRUE),
      Median_Multi_Down_Pearson = median(Multi_Down_Pearson, na.rm = TRUE),
      Median_Multi_Full_Pearson = median(Multi_Full_Pearson, na.rm = TRUE),
      IQR_Local_Pearson         = IQR(Local_Pearson,         na.rm = TRUE),
      IQR_Multi_Down_Pearson    = IQR(Multi_Down_Pearson,    na.rm = TRUE),
      IQR_Multi_Full_Pearson    = IQR(Multi_Full_Pearson,    na.rm = TRUE),
      
      Median_Local_Var          = median(Local_Var,         na.rm = TRUE),
      Median_Multi_Down_Var     = median(Multi_Down_Var,    na.rm = TRUE),
      Median_Multi_Full_Var     = median(Multi_Full_Var,    na.rm = TRUE)
    )
  
  return(list(
    summary      = median_res,
    all_splits   = all_splits_res,
    all_preds    = all_preds_long,  
    median_preds = median_preds     
  ))
}

SITE_COLUMN_NAME <- "study_site"

valid_sites <- both %>%
  count(.data[[SITE_COLUMN_NAME]]) %>%
  filter(n >= 10) %>%
  pull(.data[[SITE_COLUMN_NAME]])

final_summary_table = data.frame()
all_raw_splits      = data.frame()
all_median_preds    = data.frame()   
all_preds_raw       = data.frame()   

struct_string <- "X17"
func_string   <- "_to_"

for(st in valid_sites){
  site_result = Run_Site_Ablation_Splits(
    data        = both,
    target_site = st,
    site_col    = SITE_COLUMN_NAME,
    outcome     = "age",
    typeDat     = paste0(struct_string, "|", func_string),
    folds       = 5,
    splits      = 100
  )
  
  final_summary_table = rbind(final_summary_table, site_result$summary)
  
  raw_temp      = site_result$all_splits
  raw_temp$Site = st
  all_raw_splits = rbind(all_raw_splits, raw_temp)
  
  all_median_preds = rbind(all_median_preds, site_result$median_preds)
  all_preds_raw    = rbind(all_preds_raw,    site_result$all_preds)
}

saveRDS(final_summary_table, "final_summary_table_age_ridge.rds")
saveRDS(all_median_preds,    "all_median_preds_age_ridge.rds")
saveRDS(all_preds_raw,       "all_preds_raw_age_ridge.rds")
saveRDS(all_raw_splits,      "all_raw_splits_age_ridge.rds")
