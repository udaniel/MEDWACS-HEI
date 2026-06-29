# Limited utility of diet quality scores in anthropometric prediabetes and diabetes screening models

# **Overview**

This repository contains the R code, machine learning models, and preprocessed datasets required to reproduce the analyses and figures for the study: *"Limited utility of diet quality scores in anthropometric prediabetes and diabetes screening models"*.

This study evaluates the incremental predictive utility of integrating detailed diet quality metrics—specifically the Healthy Eating Index-2020 (HEI-2020)—into the **Machineborne Early Diabetic Warning And Control System (MEDWACS)**.

The original MEDWACS screening tool is available as a web application at <https://dtu-quantitative-sustainability-assessment.shinyapps.io/MEDWACS/>.

# **System requirements**

**Hardware requirements:** This code does not require specialized hardware and can be run on a standard computer.

**R package dependencies:** To reproduce the figures and analyses on your local machine, you will need to install R and the following package dependencies. The original code was written using **R version 4.3.2** and RStudio.

Required Packages and Versions: \* Boruta (v8.0.0) \* caret (v6.0-94) \* doParallel (v1.0.17) \* dplyr (v1.1.4) \* kernelshap (v0.7.0) \* parallel (v4.3.2) \* patchwork (v1.3.0) \* rsample (v1.2.1) \* shapviz (v0.9.4) \* tidyverse (v2.0.0) \* yardstick (v1.3.1)

# **Instructions guide**

1.  Install R and RStudio to run the code.
2.  Install all required packages listed above. If R and RStudio are already installed, package installation should take less than 10 minutes.
3.  Download or clone this repository to your local machine.

### Repository Files:

**Datasets:** \* `merged_nhanes_alltimeline_adult.rds` - The harmonized adult NHANES dataset. \* `response_clean.rds` - Cleaned response/outcome data (FPG and HbA1c). \* `dictionary_nhanes.rds` - NHANES parameter dictionary. \* `HEI_2020_all_cycles_first_day_recall.csv` - HEI-2020 scores based on 1-day dietary recall. \* `HEI_2020_no_pre_2003_first_and_second_day_recall.csv` - HEI-2020 scores based on the average of 2-day dietary recalls.

**Machine Learning Models:** \* `warning_system_all.rds` - The original MEDWACS model. \* `hei_system_noHEI.rds` - Base MEDWACS model trained on the HEI cohort. \* `hei_system_HEI_ALL.rds` - Model including Base parameters + HEI Total Score. \* `boruta_HEI_model.rds` - Model including Base parameters + Boruta-Selected HEI Components. \* `hei_system_HEI_only.rds` - Model including HEI Components Only. \* `boruta_result.rds` - Feature selection results from the Boruta algorithm.

**SHAP (Interpretability) Objects:** \* `hei_shap_result_system_noHEI.rds` \* `hei_shap_result_system_HEI_ALL.rds` \* `hei_shap_result_system_HEI_boruta.rds` \* `hei_shap_result_system_HEI_only.rds`

**Scripts:** \* `Figure_Reproduction.R` - The main script to reproduce Figure 1 (Discrimination and Calibration) and Figure 2 (SHAP Interpretability). \* `metrics_soft.R` - Helper functions containing metrics used for bootstrap resampling.

### Reproducing the Results:

Open the `Figure_Reproduction.R` file in RStudio, set your working directory to the downloaded repository folder, and run the script step-by-step to generate the final manuscript figures.

# **Citation**

If you use this code or data, please cite the corresponding paper: *(Citation details will be updated upon publication)*

To reference the original MEDWACS framework, please cite: Yoo D., Maggiore U., & Jolliet O. (2026). *Enhancing prediabetes and diabetes detection through a machine learning-enabled self-assessment approach*. Journal of Clinical Epidemiology. \# Limited utility of diet quality scores in anthropometric prediabetes and diabetes screening models

# **Overview**

This repository contains the R code, machine learning models, and preprocessed datasets required to reproduce the analyses and figures for the study: *"Limited utility of diet quality scores in anthropometric prediabetes and diabetes screening models"*.

This study evaluates the incremental predictive utility of integrating detailed diet quality metrics—specifically the Healthy Eating Index-2020 (HEI-2020)—into the **Machineborne Early Diabetic Warning And Control System (MEDWACS)**.

The original MEDWACS screening tool is available as a web application at <https://dtu-quantitative-sustainability-assessment.shinyapps.io/MEDWACS/>.

# **System requirements**

**Hardware requirements:** This code does not require specialized hardware and can be run on a standard computer.

**R package dependencies:** To reproduce the figures and analyses on your local machine, you will need to install R and the following package dependencies. The original code was written using **R version 4.3.2** and RStudio.

Required Packages and Versions: \* Boruta (v8.0.0) \* caret (v6.0-94) \* doParallel (v1.0.17) \* dplyr (v1.1.4) \* kernelshap (v0.7.0) \* parallel (v4.3.2) \* patchwork (v1.3.0) \* rsample (v1.2.1) \* shapviz (v0.9.4) \* tidyverse (v2.0.0) \* yardstick (v1.3.1)

# **Instructions guide**

1.  Install R and RStudio to run the code.
2.  Install all required packages listed above. If R and RStudio are already installed, package installation should take less than 10 minutes.
3.  Download or clone this repository to your local machine.

### Repository Files:

**Datasets:** \* `merged_nhanes_alltimeline_adult.rds` - The harmonized adult NHANES dataset. \* `response_clean.rds` - Cleaned response/outcome data (FPG and HbA1c). \* `dictionary_nhanes.rds` - NHANES parameter dictionary. \* `HEI_2020_all_cycles_first_day_recall.csv` - HEI-2020 scores based on 1-day dietary recall. \* `HEI_2020_no_pre_2003_first_and_second_day_recall.csv` - HEI-2020 scores based on the average of 2-day dietary recalls.

**Machine Learning Models:** \* `warning_system_all.rds` - The original MEDWACS model. \* `hei_system_noHEI.rds` - Base MEDWACS model trained on the HEI cohort. \* `hei_system_HEI_ALL.rds` - Model including Base parameters + HEI Total Score. \* `boruta_HEI_model.rds` - Model including Base parameters + Boruta-Selected HEI Components. \* `hei_system_HEI_only.rds` - Model including HEI Components Only. \* `boruta_result.rds` - Feature selection results from the Boruta algorithm.

**SHAP (Interpretability) Objects:** \* `hei_shap_result_system_noHEI.rds` \* `hei_shap_result_system_HEI_ALL.rds` \* `hei_shap_result_system_HEI_boruta.rds` \* `hei_shap_result_system_HEI_only.rds`

**Scripts:** \* `Figure_Reproduction.R` - The main script to reproduce Figure 1 (Discrimination and Calibration) and Figure 2 (SHAP Interpretability). \* `metrics_soft.R` - Helper functions containing metrics used for bootstrap resampling.

### Reproducing the Results:

Open the `Figure_Reproduction.R` file in RStudio, set your working directory to the downloaded repository folder, and run the script step-by-step to generate the final manuscript figures.

# **Citation**

If you use this code or data, please cite the corresponding paper: *(Citation details will be updated upon publication)*

To reference the original MEDWACS framework, please cite: Yoo D., Maggiore U., & Jolliet O. (2026). *Enhancing prediabetes and diabetes detection through a machine learning-enabled self-assessment approach*. Journal of Clinical Epidemiology. <https://doi.org/10.1016/j.jclinepi.2026.112266>.
