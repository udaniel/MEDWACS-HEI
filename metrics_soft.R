ROC_boot <- function(splits) {
    
    # validation_pred <- analysis(all_boots_test$splits[1][[1]])
    validation_pred <- analysis(splits)
    
    total_results <- 
        validation_pred %>% yardstick::roc_auc(truth = binarize_diab, hei_system_noHEI_pred) %>% 
        mutate(.model = "MEDWACS_no_HEI") %>% 
        bind_rows(
            validation_pred %>% yardstick::roc_auc(truth = binarize_diab, hei_system_HEI_ALL_pred) %>% 
                mutate(.model = "inc_HEI_All")
        ) %>% 
        bind_rows(
            validation_pred %>% yardstick::roc_auc(truth = binarize_diab, boruta_HEI_model_pred) %>% 
                mutate(.model = "boruta_model")
        ) %>% 
        # bind_rows(
        #     validation_pred %>% yardstick::roc_auc(truth = binarize_diab, hei_system_individual_inter_pred) %>% 
        #         mutate(.model = "final_model_inter")
        # ) %>% 
        bind_rows(
            validation_pred %>% yardstick::roc_auc(truth = binarize_diab, hei_system_HEI_only_pred) %>% 
                mutate(.model = "HEI_Only")
        )
    
    
    total_results %>% 
        rename(estimate = .estimate,
               term = .model)
        # tibble(model = total_results$.model,
    #        term = paste0(total_results$.metric, "_perc"),
    #        estimate = total_results$.estimate * 100) %>% 
    #     mutate(term = paste0(model, ".", term))
}






multi_boot_HEI <- function(splits) {
    
    validation_pred <- analysis(splits)
    
    total_results <- 
        validation_pred %>% yardstick::roc_auc(truth = binarize_diab, predict_noHEI) %>% 
        bind_rows(
            validation_pred %>% yardstick::pr_auc(truth = binarize_diab, predict_noHEI),
            validation_pred %>% yardstick::brier_class(truth = binarize_diab, predict_noHEI) %>% mutate(.estimate = 1 - .estimate)
        ) %>% mutate(.model = "No HEI") %>% 
        bind_rows(bind_rows(
            validation_pred %>% yardstick::roc_auc(truth = binarize_diab, predict_HEI_ALL),
            validation_pred %>% yardstick::pr_auc(truth = binarize_diab, predict_HEI_ALL),
            validation_pred %>% yardstick::brier_class(truth = binarize_diab, predict_HEI_ALL) %>% mutate(.estimate = 1 - .estimate)
        ) %>% mutate(.model = "HEI All")) %>% 
        bind_rows(bind_rows(
            validation_pred %>% yardstick::roc_auc(truth = binarize_diab, predict_HEI_boruta),
            validation_pred %>% yardstick::pr_auc(truth = binarize_diab, predict_HEI_boruta),
            validation_pred %>% yardstick::brier_class(truth = binarize_diab, predict_HEI_boruta) %>% mutate(.estimate = 1 - .estimate)
        ) %>% mutate(.model = "HEI Boruta")) %>% 
        # bind_rows(bind_rows(
        #     validation_pred %>% yardstick::roc_auc(truth = binarize_diab, predict_individual_inter),
        #     validation_pred %>% yardstick::pr_auc(truth = binarize_diab, predict_individual_inter),
        #     validation_pred %>% yardstick::brier_class(truth = binarize_diab, predict_individual_inter) %>% mutate(.estimate = 1 - .estimate)
        # ) %>% mutate(.model = "HEI interaction")) %>% 
        bind_rows(bind_rows(
            validation_pred %>% yardstick::roc_auc(truth = binarize_diab, predict_HEI_only),
            validation_pred %>% yardstick::pr_auc(truth = binarize_diab, predict_HEI_only),
            validation_pred %>% yardstick::brier_class(truth = binarize_diab, predict_HEI_only) %>% mutate(.estimate = 1 - .estimate)
        ) %>% mutate(.model = "HEI only"))
        
    
    tibble(model = total_results$.model,
           term = paste0(total_results$.metric, "_perc"),
           estimate = total_results$.estimate * 100) %>% 
        mutate(term = paste0(model, ".", term))
}


multi_boot_final_HEI <- function(splits) {
    
    validation_pred <- analysis(splits)
    
    total_results <- 
        validation_pred %>% yardstick::roc_auc(truth = binarize_diab, predict_test_final_HEI) %>% 
        bind_rows(validation_pred %>% yardstick::pr_auc(truth = binarize_diab, predict_test_final_HEI),
                  validation_pred %>% yardstick::brier_class(truth = binarize_diab, predict_test_final_HEI) %>% mutate(.estimate = 1 - .estimate)) %>% 
        mutate(.model = "HEI Final",
               .estimate = .estimate * 100,
               term = paste0(.metric, "_perc")) %>% 
        rename(model = .model,
               metric = .metric,
               estimate = .estimate) %>% 
        select(-.estimator)
    
}
