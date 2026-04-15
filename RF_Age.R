library(dplyr)       
library(tidyr)       
library(mlr)        
library(tuneRanger)  
library(ranger) 
library(tidyverse)
library(readr)
# Setting the seed
set.seed(123)
logit <- function(p) { log(p / (1 - p)) }
expit <- function(x) { exp(x) / (1 + exp(x)) }
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
source("/media/disk2/multivariateBWAS/code_for_review/build_funcs.R")
source("/media/disk2/multivariateBWAS/code_for_review/influence_funcs.R")
source("/media/disk2/multivariateBWAS/code_for_review/corr_funcs.R")

# Loading in color palette
cbbPalette <- c("#000000","#E69F00","#56B4E9","#009E73",
                "#F0E442","#0072B2","#D55E00","#CC79A7")

# tuning function
helper_ranger = function(trainDat, outcome, testDat, task){
  if (sd(predict(task$model, newdata = data.frame(testDat))$data$response) > 1e-5) {
    return(task$model)
  }
  tune_table <- task$results[order(task$results$mse), ]
  for (i in seq_len(nrow(tune_table))) {
    params <- tune_table[i, ]
    lrn <- mlr::makeLearner(
      "regr.ranger",
      num.trees       = 250,
      mtry            = params$mtry,
      min.node.size   = params$min.node.size,
      sample.fraction = params$sample.fraction,
      num.threads     = 1,
      respect.unordered.factors = "order"
    )
    tr_task <- mlr::makeRegrTask(id = "rf_task_retry", data = trainDat, target = outcome)
    wm <- mlr::train(lrn, tr_task)
    
    if (sd(predict(wm, newdata = data.frame(testDat))$data$response) > 1e-5) {
      return(wm)
    }
  }
  return(task$model)
}

###############################################################################
# Modified Single_Split_RF: records both Pearson and One-Step
###############################################################################
Single_Split_RF <- function(data_clean, target_site, site_col, outcome, typeDat, folds){
  target_data = data_clean %>% filter(.data[[site_col]] == target_site)
  other_data = data_clean %>% filter(.data[[site_col]] != target_site)
  
  target_data$fold <- sample(rep(1:folds, length.out = nrow(target_data)))
  
  # Pearson correlation per fold
  fold_local_pearson = c(); fold_down_pearson = c(); fold_full_pearson = c()
  # One-Step per fold
  fold_local_os = c(); fold_down_os = c(); fold_full_os = c()
  # Influence vectors (for variance, based on One-Step)
  inf_local_all = c(); inf_down_all = c(); inf_full_all = c()
  
  for(i in 1:folds){
    test = target_data %>% filter(fold == i)
    train_local = target_data %>% filter(fold != i)
    N_train = nrow(train_local) 
    
    train_multi_full = bind_rows(train_local, other_data)
    site_counts <- train_multi_full %>%
      count(.data[[site_col]]) %>%
      mutate(
        prop = n / sum(n),
        target_n = floor(prop * N_train) 
      )
    remainder <- N_train - sum(site_counts$target_n)
    if (remainder > 0) {
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
    
    #---------------------------------------------------------------------------
    # Inner function: returns Pearson, One-Step, and influence vector
    #---------------------------------------------------------------------------
    run_model_get_metrics_RF <- function(train_df, test_df) {
      if(typeDat != "Baseline"){
        if(outcome == "age"){
          target_cols = c(colnames(train_df)[grepl(typeDat, colnames(train_df))], "sex", "euler", "meanFD", outcome)
        } else {
          target_cols = c(colnames(train_df)[grepl(typeDat, colnames(train_df))], "age", "sex", "euler", "meanFD", outcome)
        }
      } else {
        if(outcome == "age"){
          target_cols = c("sex", "euler", "meanFD", outcome)
        } else {
          target_cols = c("age", "sex", "euler", "meanFD", outcome)
        }
      }
      
      train_sub = data.frame(train_df[, colnames(train_df) %in% target_cols])
      test_sub  = data.frame(test_df[, colnames(test_df) %in% setdiff(target_cols, outcome)])
      y_real    = as.numeric(test_df[[outcome]])
      
      rf_res <- tryCatch({
        if (nrow(train_sub) < 15) {
          fallback_mod <- ranger::ranger(
            x = train_sub[, setdiff(colnames(train_sub), outcome), drop = FALSE],
            y = train_sub[[outcome]],
            num.trees = 250,
            min.node.size = 2  
          )
          preds = predict(fallback_mod, data = test_sub)$predictions
          list(preds = preds)
        } else {
          task = makeRegrTask(id = "rf_task", data = train_sub, target = outcome)
          rfTune = tuneRanger(task, num.trees = 250, iters = 25, iters.warmup = 10, show.info = FALSE)
          final_model = helper_ranger(train_sub, outcome, test_sub, rfTune)
          preds = predict(final_model, newdata = test_sub)$data$response
          list(preds = preds)
        }
      }, error = function(e) {
        print(paste("error:", e$message)) 
        return(NULL)
      })
      
      if(is.null(rf_res)) return(list(pearson = NA, os = NA, inf = rep(NA, nrow(test_df))))
      
      preds = as.vector(rf_res$preds)
      
      # Pearson correlation
      pearson_val = cor(preds, y_real, use = "complete.obs")
      # One-Step correlation
      os_val = cor_os(preds, y_real)
      # Influence function (based on One-Step)
      inf_vec = if_logit_sq(preds, y_real) 
      
      return(list(pearson = pearson_val, os = os_val, inf = inf_vec))
    }
    
    res_local = run_model_get_metrics_RF(train_local, test)
    res_down  = run_model_get_metrics_RF(train_multi_down, test)
    res_full  = run_model_get_metrics_RF(train_multi_full, test)
    
    # Pearson
    fold_local_pearson = c(fold_local_pearson, res_local$pearson)
    fold_down_pearson  = c(fold_down_pearson, res_down$pearson)
    fold_full_pearson  = c(fold_full_pearson, res_full$pearson)
    
    # One-Step
    fold_local_os = c(fold_local_os, res_local$os)
    fold_down_os  = c(fold_down_os, res_down$os)
    fold_full_os  = c(fold_full_os, res_full$os)
    
    # Influence (for variance calculation)
    inf_local_all = c(inf_local_all, res_local$inf)
    inf_down_all  = c(inf_down_all, res_down$inf)
    inf_full_all  = c(inf_full_all, res_full$inf)
    
  } 
  
  N_total = nrow(target_data)
  
  var_local = var(inf_local_all, na.rm = TRUE) / N_total
  var_down  = var(inf_down_all, na.rm = TRUE) / N_total
  var_full  = var(inf_full_all, na.rm = TRUE) / N_total
  
  return(data.frame(
    # Pearson
    Local_Pearson      = mean(fold_local_pearson, na.rm = TRUE),
    Multi_Down_Pearson = mean(fold_down_pearson, na.rm = TRUE),
    Multi_Full_Pearson = mean(fold_full_pearson, na.rm = TRUE),
    # One-Step
    Local_OS      = mean(fold_local_os, na.rm = TRUE),
    Multi_Down_OS = mean(fold_down_os, na.rm = TRUE),
    Multi_Full_OS = mean(fold_full_os, na.rm = TRUE),
    # Variance (based on influence function / One-Step)
    Local_Var      = var_local,
    Multi_Down_Var = var_down,
    Multi_Full_Var = var_full
  ))
}


###############################################################################
# Modified Run_Site_Ablation_Splits_RF: summary includes both Pearson & OS
###############################################################################
Run_Site_Ablation_Splits_RF <- function(data, target_site, site_col, outcome, typeDat, folds = 5, splits = 100){
  data_clean = data %>% drop_na(all_of(outcome), all_of(site_col))
  print(paste("=================================================="))
  print(paste("Starting RF Ablation for Site:", target_site))
  print(paste("Total Splits:", splits, "| Target N:", nrow(data_clean %>% filter(.data[[site_col]] == target_site))))
  print(paste("=================================================="))
  all_splits_res <- data.frame()
  
  for(s in 1:splits){
    if(s %% 10 == 0) print(paste("  ... Processing split", s, "/", splits))
    split_res <- Single_Split_RF(data_clean, target_site, site_col, outcome, typeDat, folds)
    split_res$split_num <- s
    all_splits_res <- rbind(all_splits_res, split_res)
  }
  
  median_res <- all_splits_res %>%
    summarise(
      Target_Site = target_site,
      # Pearson summary
      Median_Local_Pearson      = median(Local_Pearson, na.rm = TRUE),
      Median_Multi_Down_Pearson = median(Multi_Down_Pearson, na.rm = TRUE),
      Median_Multi_Full_Pearson = median(Multi_Full_Pearson, na.rm = TRUE),
      IQR_Local_Pearson         = IQR(Local_Pearson, na.rm = TRUE),
      IQR_Multi_Down_Pearson    = IQR(Multi_Down_Pearson, na.rm = TRUE),
      IQR_Multi_Full_Pearson    = IQR(Multi_Full_Pearson, na.rm = TRUE),
      # One-Step summary
      Median_Local_OS      = median(Local_OS, na.rm = TRUE),
      Median_Multi_Down_OS = median(Multi_Down_OS, na.rm = TRUE),
      Median_Multi_Full_OS = median(Multi_Full_OS, na.rm = TRUE),
      IQR_Local_OS         = IQR(Local_OS, na.rm = TRUE),
      IQR_Multi_Down_OS    = IQR(Multi_Down_OS, na.rm = TRUE),
      IQR_Multi_Full_OS    = IQR(Multi_Full_OS, na.rm = TRUE),
      # Variance summary (based on influence / One-Step)
      Median_Local_Var      = median(Local_Var, na.rm = TRUE),
      Median_Multi_Down_Var = median(Multi_Down_Var, na.rm = TRUE),
      Median_Multi_Full_Var = median(Multi_Full_Var, na.rm = TRUE)
    )
  
  return(list(summary = median_res, all_splits = all_splits_res))
}

library(dplyr)
library(tidyr)
both$sex <- as.factor(both$sex)

SITE_COLUMN_NAME <- "study_site" 
struct_string <- "X17"
func_string <- "_to_"

valid_sites <- both %>% 
  count(.data[[SITE_COLUMN_NAME]]) %>% 
  filter(n >= 10) %>% 
  pull(.data[[SITE_COLUMN_NAME]])


sites_to_run <- valid_sites[1:9]

final_summary_table_rf <- data.frame()
all_raw_splits_rf <- data.frame() 

for(st in sites_to_run){
  
  site_result <- Run_Site_Ablation_Splits_RF(
    data = both, 
    target_site = st, 
    site_col = SITE_COLUMN_NAME, 
    outcome = "age", 
    typeDat = paste0(struct_string, "|", func_string),
    folds = 5,    
    splits = 100   
  )
  final_summary_table_rf <- rbind(final_summary_table_rf, site_result$summary)
  
  raw_temp <- site_result$all_splits
  raw_temp$Site <- st
  all_raw_splits_rf <- rbind(all_raw_splits_rf, raw_temp)
}


print("========== RF Final Summary Table (Pearson + One-Step) ==========")
print(final_summary_table_rf)

saveRDS(final_summary_table_rf, "/home/guanru/meta/final_summary_table_age_RF.RDS")
saveRDS(all_raw_splits_rf, "/home/guanru/meta/all_raw_splits_age_RF.RDS")