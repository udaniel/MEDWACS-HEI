# ==============================================================================
# Script Name: Figure_Reproduction.R
# Purpose: Reproduces Figure 1 (Model Performance & Calibration) and 
#          Figure 2 (SHAP Interpretability) for the manuscript.
# Models: 
#   1. Base MEDWACS (7 parameters)
#   2. Base + HEI Total Score
#   3. Base + Boruta-Selected HEI Components
#   4. HEI Components Only
# ==============================================================================

#### 1. Load Required Libraries ####
library(tidyverse)
library(parallel)
library(doParallel)
library(yardstick)
library(rsample)
library(Boruta)
library(patchwork)
library(kernelshap)
library(shapviz)
library(caret)


#### 2. Load Pre-trained Models and Data ####
# Load ML models
warning_system_all <- read_rds("warning_system_all.rds")
boruta_result <- read_rds("boruta_result.rds")
hei_system_noHEI <- read_rds("hei_system_noHEI.rds")
hei_system_HEI_ALL <- read_rds("hei_system_HEI_ALL.rds")
boruta_HEI_model <- read_rds("boruta_HEI_model.rds")
hei_system_HEI_only <- read_rds("hei_system_HEI_only.rds")

# Define readable labels for HEI components for plots
hei_labels <- c(
    HEI_TOTALFRT     = "HEI Total Fruits",
    HEI_FRT          = "HEI Whole Fruits",
    HEI_VEG          = "HEI Total Vegetables",
    HEI_GREENNBEAN   = "HEI Greens and Beans",
    HEI_TOTALPRO     = "HEI Total Protein Foods",
    HEI_SEAPLANTPRO  = "HEI Seafood and Plant Proteins",
    HEI_WHOLEGRAIN   = "HEI Whole Grains",
    HEI_DAIRY        = "HEI Dairy",
    HEI_FATTYACID    = "HEI Fatty Acids",
    HEI_REFINEDGRAIN = "HEI Refined Grains",
    HEI_SODIUM       = "HEI Sodium",
    HEI_ADDEDSUGAR   = "HEI Added Sugars",
    HEI_SATFAT       = "HEI Saturated Fats"
)

# Load helper functions and raw datasets
source("metrics_soft.R")
merged_nhanes_alltimeline_adult <- read_rds("merged_nhanes_alltimeline_adult.rds")
response_clean <- read_rds("response_clean.rds")
dictionary_nhanes <- read_rds("dictionary_nhanes.rds")
HEI_2020_one_day <- read.csv("HEI_2020_all_cycles_first_day_recall.csv") %>% as_tibble()
HEI_2020_two_day <- read.csv("HEI_2020_no_pre_2003_first_and_second_day_recall.csv") %>% as_tibble()


#### 3. Data Harmonization and Preprocessing ####
# Merge 1-day and 2-day HEI-2020 recall data. 
# Prioritize 2-day average recall (2003-2018), and fall back to 1-day recall for earlier cycles (1999-2002).
HEI_all_clean <- HEI_2020_one_day %>% 
    anti_join(HEI_2020_two_day, by = c("SEQN", "SDDSRVYR")) %>% 
    bind_rows(HEI_2020_two_day) %>% 
    relocate(SDDSRVYR, .after = SEQN) %>% 
    rename_with(~ gsub("2020", "", .x))

# Filter NHANES III (SDDSRVYR == -1), define the composite outcome (prediabetes/diabetes), 
# and join target variables.
merged_nhanes_alltimeline_adult %>% 
    filter(SDDSRVYR != -1) %>% 
    mutate(plasma_fasting = response_clean$LBXGLU[match(SEQN_new, response_clean$SEQN_new)]) %>% 
    rename(glycohemoglobin = LBXGH) %>% 
    select(SEQN, SEQN_new, SDDSRVYR, 
           names(warning_system_all$trainingData)[-1], 
           plasma_fasting, 
           glycohemoglobin) %>% 
    filter((!is.na(plasma_fasting) & !is.na(glycohemoglobin))) %>%
    mutate(binarize_diab = ifelse(plasma_fasting >= 100 | glycohemoglobin >= 5.7, "X1", "X0"),
           binarize_diab = factor(binarize_diab, levels = c("X1", "X0")),
           RIAGENDR = factor(RIAGENDR, levels = 1:2)) -> full_data_small_features

# Join HEI data and drop missing values
full_data_small_features %>% 
    left_join(HEI_all_clean %>% select(-RIDAGEYR), by = join_by("SEQN", "SDDSRVYR")) %>% 
    relocate(HEI_ALL:HEI_SATFAT, .before = plasma_fasting) %>% 
    filter(if_all(RIDAGEYR:HEI_SATFAT, ~!is.na(.x))) -> full_data_with_HEI


#### 4. Data Splitting ####
# Stratified 70/30 Train/Test split
set.seed(123)
ind <- createDataPartition(full_data_with_HEI$binarize_diab, p = 0.7, list = F)
hei_train_data <- full_data_with_HEI[ind, ]
hei_test_data <- full_data_with_HEI[-ind, ]

# Setup cross-validation framework for model evaluation
set.seed(123)
model_ctrl_warning_all <- trainControl(method = "adaptive_cv",
                                       repeats = 3,
                                       number = 10,
                                       allowParallel = T,
                                       savePredictions = "final",
                                       adaptive = list(
                                           min = 5,
                                           alpha = 0.05,
                                           method = "gls",
                                           complete = T
                                       ),
                                       classProbs = T,
                                       summaryFunction = twoClassSummary,
                                       search = "random",
                                       index = createMultiFolds(hei_train_data$binarize_diab, k = 10, times = 3))


#### 5. Assessment (Internal Validation on Test Set); Figure 1 ####

# Generate predictions for all 4 models on the held-out test set
hei_test_data_wPred <- 
    hei_test_data %>%
    mutate(predict_noHEI = predict(hei_system_noHEI, ., type = "prob")$X1,
           predict_HEI_ALL = predict(hei_system_HEI_ALL, ., type = "prob")$X1,
           predict_HEI_boruta = predict(boruta_HEI_model, ., type = "prob")$X1,
           predict_HEI_only = predict(hei_system_HEI_only, ., type = "prob")$X1)

# Perform 1,000 bootstrap resamples to calculate 95% Confidence Intervals
set.seed(27)
boots_system_test_hei <- bootstraps(hei_test_data_wPred, times = 1000, apparent = TRUE)
result_system_test_hei <- boots_system_test_hei %>% mutate(result = map(splits, multi_boot_HEI))
result_system_test_multi_hei <- int_pctl(result_system_test_hei, result)

# Clean and format bootstrap results
result_system_test_multi_hei %>% 
    separate_wider_delim(term, ".", names = c("model", "metric")) %>% 
    arrange(factor(metric, levels = c(
        "roc_auc_perc", "pr_auc_perc", "brier_class_perc"
    )),
    factor(model, levels = c(
        "No HEI", "HEI All", "HEI Boruta", "HEI only"
    ))) %>% 
    mutate(all_together = paste0(formatC(signif(.estimate, 3), digits = 3, format = "fg", flag = "#"), " (",
                                 formatC(signif(.lower, 3), digits = 3, format = "fg", flag = "#"), "-",
                                 formatC(signif(.upper, 3), digits = 3, format = "fg", flag = "#"), ")")) -> result_system_test_multi_hei_clean

print(result_system_test_multi_hei_clean)

# Format detailed results for Figure 1A (Violin/Boxplots)
result_detail_hei_clean <- 
    bind_rows(result_system_test_hei$result[1:1000]) %>% 
    select(-model) %>% 
    separate_wider_delim(term, ".", names = c("model", "metric")) %>% 
    dplyr::rename(performance = estimate) %>% 
    mutate(metric = factor(metric, levels = c("roc_auc_perc", "pr_auc_perc", "brier_class_perc"),
                           labels = c("100 X ROCAUC", "100 X PRAUC", "1 - Brier Score (%)")),
           model = factor(model, levels = c("No HEI", "HEI All", "HEI Boruta", "HEI only"),
                          labels = c("Base MEDWACS (7 params)",
                                     "Base + HEI Total Score",
                                     "Base + Boruta-Selected HEI Components", 
                                     "HEI Components Only")))

# Plot Figure 1A: Discrimination and Calibration Metrics
result_detail_hei_clean %>%
    ggplot(aes(x = metric, y = performance, fill = model)) +
    geom_violin(position = position_dodge(0.3), alpha = 0.5, show.legend = F, width = 2) +
    geom_boxplot(position = position_dodge(0.3), width = 0.1) +
    xlab("Metric") + ylab("Performance (%)") +
    scale_y_continuous(breaks = seq(0, 100, 5), limits = c(50, 100)) +
    scale_colour_viridis_d(option = "turbo") +
    scale_fill_viridis_d(option = "turbo") +
    labs(fill = "Model") + 
    theme_bw() + 
    theme(
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = 0, , face = "bold"),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 15),
        legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = T)) -> final_results_test_hei



#### 6. Calibration Plots; Figure 1B ####

# Generate calibration data using caret
calPlotData_all <- calibration(binarize_diab ~ predict_noHEI + predict_HEI_ALL + predict_HEI_boruta + predict_HEI_only, hei_test_data_wPred)
calPlotData_all <- calPlotData_all$data %>% as_tibble()

# Format calibration data for plotting
calPlotData_all %>% 
    mutate(across(c(Percent:Upper, midpoint), ~ .x / 100),
           calibModelVar = factor(calibModelVar, 
                                  levels = c("predict_noHEI", "predict_HEI_ALL", "predict_HEI_boruta", "predict_HEI_only"),
                                  labels = c("Base MEDWACS (7 params)", "Base + HEI Total Score",
                                             "Base + Boruta-Selected HEI Components", "HEI Components Only"))) -> calPlotData_all

# Plot Figure 1B: Calibration Curves
ggplot(calPlotData_all %>% na.omit(), aes(x = midpoint, y = Percent, color = calibModelVar)) +
    geom_line(linewidth = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    theme_bw() +
    labs(color = "Models") + 
    scale_x_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, 0.2)) +
    scale_colour_viridis_d(option = "turbo") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 15),
          legend.position = "bottom",
          strip.text = element_text(size = 15, color = "white", face = "bold"),
          strip.background = element_rect(color = "black", fill = "black", linetype = "solid"),
          plot.title = element_text(size = 20, face = "bold"),
          panel.spacing = unit(2, "lines"),
          plot.margin = margin(t = 30, unit = "pt")) +
    guides(color = guide_legend(nrow = 2, byrow = T)) + 
    xlab("Predicted probabilities") + ylab("Observed probabilities") -> CALIB_curve_all

# Combine Panel A and B using Patchwork
final_results_test_hei / CALIB_curve_all +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 30, face = 'bold')) -> performance_total_figure

# ggsave(filename = "Figure1_performance_total.png", plot = performance_total_figure, width = 10, height = 15)


#### 7. SHAP Values Analysis; Figure 2 ####

# Import pre-computed SHAP models
hei_shap_result_system_noHEI <- read_rds("hei_shap_result_system_noHEI.rds")
hei_shap_result_system_HEI_ALL <- read_rds("hei_shap_result_system_HEI_ALL.rds")
hei_shap_result_system_HEI_boruta <- read_rds("hei_shap_result_system_HEI_boruta.rds")
hei_shap_result_system_HEI_only <- read_rds("hei_shap_result_system_HEI_only.rds")


# ---------------------------------------------------------
# Panel A: Base MEDWACS (7 parameters)
# ---------------------------------------------------------
hei_sv_system_noHEI <- shapviz(hei_shap_result_system_noHEI)
names(hei_sv_system_noHEI) <- c("prediabetes_diabetes", "Normal")
all_old_names_system_noHEI <- colnames(hei_sv_system_noHEI$Normal)

# Map raw variable codes to human-readable names
dictionary_nhanes %>% 
    filter(variable_codename_use %in% all_old_names_system_noHEI) %>%
    distinct(variable_codename_use, .keep_all = T) %>% 
    mutate(variable_description_use = ifelse(variable_codename_use == "RIDAGEYR", "Age",
                                             ifelse(variable_codename_use == "RIAGENDR", "Gender",
                                                    ifelse(variable_codename_use == "BMXBMI", "Body Mass Index (kg/m²)", variable_description_use)))) %>% 
    arrange(match(variable_codename_use, all_old_names_system_noHEI)) %>% 
    pull(variable_description_use) -> all_new_names_system_noHEI


colnames(hei_sv_system_noHEI$prediabetes_diabetes) <- all_new_names_system_noHEI

# Drop normal class to plot risk directionality only
hei_sv_system_noHEI$Normal <- NULL
all_sv_imp_system_noHEI <- sv_importance(hei_sv_system_noHEI, kind = "bee", max_display = 20, viridis_args = list(option = "H"))

all_sv_imp_system_noHEI +
    theme_classic() +
    xlab("SHAP value\n(prediction probability)") +
    ggtitle("") +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15, face = "bold", vjust = -1),
          legend.title = element_text(size = 15, face = "bold"),
          plot.title = element_text(size = 15, face = "bold"),
          plot.margin = margin(b = 1, unit = "cm")) +
    guides(
        color = guide_colourbar(
            barheight = unit(13, "cm"))
    ) -> p_shap_system_noHEI


# ---------------------------------------------------------
# Panel B: Base + HEI Total Score
# ---------------------------------------------------------
hei_sv_system_HEI_ALL <- shapviz(hei_shap_result_system_HEI_ALL)
names(hei_sv_system_HEI_ALL) <- c("prediabetes_diabetes", "Normal")
all_old_names_system_HEI_ALL <- colnames(hei_sv_system_HEI_ALL$Normal)

# Rename features
dictionary_nhanes %>% 
    filter(variable_codename_use %in% all_old_names_system_HEI_ALL) %>%
    distinct(variable_codename_use, .keep_all = T) %>% 
    mutate(variable_description_use = ifelse(variable_codename_use == "RIDAGEYR", "Age",
                                             ifelse(variable_codename_use == "RIAGENDR", "Gender",
                                                    ifelse(variable_codename_use == "BMXBMI", "Body Mass Index (kg/m²)", variable_description_use)))) %>% 
    arrange(match(variable_codename_use, all_old_names_system_HEI_ALL)) %>% 
    pull(variable_description_use) -> all_new_names_system_HEI_ALL

all_new_names_system_HEI_ALL <- c(all_new_names_system_HEI_ALL, all_old_names_system_HEI_ALL[grepl(pattern = "HEI", x = all_old_names_system_HEI_ALL)])
names(all_new_names_system_HEI_ALL) <- all_old_names_system_HEI_ALL
all_new_names_system_HEI_ALL <- gsub(replacement = "HEI Total Score", pattern = "HEI_ALL", x = all_new_names_system_HEI_ALL)

colnames(hei_sv_system_HEI_ALL$prediabetes_diabetes) <- all_new_names_system_HEI_ALL[colnames(hei_sv_system_HEI_ALL$prediabetes_diabetes)]

hei_sv_system_HEI_ALL$Normal <- NULL
all_sv_imp_system_HEI_ALL <- sv_importance(hei_sv_system_HEI_ALL, kind = "bee", max_display = 20, viridis_args = list(option = "H"))

all_sv_imp_system_HEI_ALL +
    theme_classic() +
    xlab("SHAP value\n(prediction probability)") +
    ggtitle("") +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15, face = "bold", vjust = -1),
          legend.title = element_text(size = 15, face = "bold"),
          plot.title = element_text(size = 15, face = "bold"),
          plot.margin = margin(b = 1, unit = "cm")) -> p_shap_system_HEI_ALL


# ---------------------------------------------------------
# Panel C: Base + Boruta-Selected HEI Components
# ---------------------------------------------------------
hei_sv_system_HEI_boruta <- shapviz(hei_shap_result_system_HEI_boruta)
names(hei_sv_system_HEI_boruta) <- c("prediabetes_diabetes", "Normal")
all_old_names_system_HEI_boruta <- colnames(hei_sv_system_HEI_boruta$Normal)

# Rename features
dictionary_nhanes %>% 
    filter(variable_codename_use %in% all_old_names_system_HEI_boruta) %>%
    distinct(variable_codename_use, .keep_all = T) %>% 
    mutate(variable_description_use = ifelse(variable_codename_use == "RIDAGEYR", "Age",
                                             ifelse(variable_codename_use == "RIAGENDR", "Gender",
                                                    ifelse(variable_codename_use == "BMXBMI", "Body Mass Index (kg/m²)", variable_description_use)))) %>% 
    arrange(match(variable_codename_use, all_old_names_system_HEI_boruta)) %>% 
    pull(variable_description_use) -> all_new_names_system_HEI_boruta

all_new_names_system_HEI_boruta <- c(all_new_names_system_HEI_boruta, all_old_names_system_HEI_boruta[grepl(pattern = "HEI", x = all_old_names_system_HEI_boruta)])
names(all_new_names_system_HEI_boruta) <- all_old_names_system_HEI_boruta

idx_HEI_boruta <- names(all_new_names_system_HEI_boruta) %in% names(hei_labels)
all_new_names_system_HEI_boruta[idx_HEI_boruta] <- hei_labels[names(all_new_names_system_HEI_boruta)[idx_HEI_boruta]]

colnames(hei_sv_system_HEI_boruta$prediabetes_diabetes) <- all_new_names_system_HEI_boruta[colnames(hei_sv_system_HEI_boruta$prediabetes_diabetes)]

hei_sv_system_HEI_boruta$Normal <- NULL
all_sv_imp_system_HEI_boruta <- sv_importance(hei_sv_system_HEI_boruta, kind = "bee", max_display = 20, viridis_args = list(option = "H"))

all_sv_imp_system_HEI_boruta +
    theme_classic() +
    xlab("SHAP value\n(prediction probability)") +
    ggtitle("") +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15, face = "bold", vjust = -1),
          legend.title = element_text(size = 15, face = "bold"),
          plot.title = element_text(size = 15, face = "bold")) -> p_shap_system_HEI_boruta


# ---------------------------------------------------------
# Panel D: HEI Components Only
# ---------------------------------------------------------
hei_sv_system_HEI_only <- shapviz(hei_shap_result_system_HEI_only)
names(hei_sv_system_HEI_only) <- c("prediabetes_diabetes", "Normal")
all_old_names_system_HEI_only <- colnames(hei_sv_system_HEI_only$Normal)

# Apply specific HEI labels defined at the top of the script
colnames(hei_sv_system_HEI_only$prediabetes_diabetes) <- hei_labels[all_old_names_system_HEI_only]

hei_sv_system_HEI_only$Normal <- NULL
all_sv_imp_system_HEI_only <- sv_importance(hei_sv_system_HEI_only, kind = "bee", max_display = 20, viridis_args = list(option = "H"))

all_sv_imp_system_HEI_only +
    theme_classic() +
    xlab("SHAP value\n(prediction probability)") +
    ggtitle("") +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15, face = "bold", vjust = -1),
          legend.title = element_text(size = 15, face = "bold"),
          plot.title = element_text(size = 15, face = "bold")) -> p_shap_system_HEI_only


# ---------------------------------------------------------
# Assemble Figure 2 using Patchwork
# ---------------------------------------------------------
p_shap_system_noHEI + 
    p_shap_system_HEI_ALL +
    p_shap_system_HEI_boruta +
    p_shap_system_HEI_only + 
    plot_annotation(tag_levels = 'A') +
    plot_layout(byrow = T, guides = "collect") & 
    guides(
        color = guide_colourbar(
            barheight = unit(10, "cm"),
            barwidth = unit(1.5, "cm"),
            ticks = F,
        )
    ) &
    theme(plot.tag = element_text(size = 30, face = "bold"),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20, face = "bold", vjust = -1),
          legend.title = element_text(size = 20, face = "bold", margin = margin(b = 20)),
          legend.text = element_text(size = 20),
          plot.title = element_text(size = 20, face = "bold")) -> p_shap_all


# ggsave("Figure2_shap_all.png", plot = p_shap_all, width = 25, height = 20)
