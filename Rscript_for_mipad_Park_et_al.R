#Project title: Unsupervised integration of peripheral blood-based multi-omics datasets for identifying new subtypes of AD 
#Authors for this R script: Jong-Chan Park (University College London), Natalia B. Torres (University College London)
#E-mail addresses of the authors: jong-chan.park@ucl.ac.uk, natalia.torres.18@ucl.ac.uk
#All commands are based on the MOFA+ vignette from this link: https://github.com/bioFAM/MOFA2
#Last modified: 1st May, 2021 (by Jong-Chan Park., Ph.D.)

#Let's Start!

########## 01. Basic settings ##########

library(MOFA2)
library(data.table)
library(reticulate)
use_python("/home/dell/anaconda3/envs/r-reticulate/bin", required=TRUE)
library(MultiAssayExperiment)
library(MOFAdata)
library(ggplot2)
library(tidyverse)
utils::data("CLL_data")
install.packages("data.table")
install.packages("tidyverse")

########## 02. Read raw datasets ##########

#making AD_metadata
AD_metadata <- fread("/home/dell/newinput/common_AD_metadata_new.txt")

#Read raw dataset files 
read.table("/home/dell/newinput/ipad1_tgwas_new01c.txt")
ipad1_tgwas <- read.table("/home/dell/newinput/ipad1_tgwas_new01c.txt")
read.table("/home/dell/newinput/ipad2_miRNA_new.txt")
ipad2_miRNA <- read.table("/home/dell/newinput/ipad2_miRNA_new.txt")
read.table("/home/dell/newinput/ipad3_proteomics_new.txt")
ipad3_proteomics <- read.table("/home/dell/newinput/ipad3_proteomics_new.txt")
read.table("/home/dell/newinput/ipad4_blood_biomarkers_new.txt")
ipad4_blood_biomarkers <- read.table("/home/dell/newinput/ipad4_blood_biomarkers_new.txt")

#Load datasets to MOFA+
Multiomic_dataset <- list(ip_tgwas=ipad1_tgwas, ip_miRNA=ipad2_miRNA, ip_proteomics=ipad3_proteomics, 
                          ip_biomarkers=ipad4_blood_biomarkers)
read.table("/home/dell/newinput/common_AD_covariates_new.txt")
AD_covariates2 <- read.table("/home/dell/newinput/common_AD_covariates_new.txt")
head(AD_covariates2)
mae_AD <- MultiAssayExperiment(
  experiments = Multiomic_dataset, 
  colData = AD_covariates2
)

########## 03. Training and Running ##########

MOFAobject <- create_mofa(mae_AD)
plot_data_overview(MOFAobject)
data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 19
get_default_training_options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42
MOFAobject <- prepare_mofa(
  MOFAobject, 
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#Running
MOFAobject <- run_mofa(MOFAobject) 

########## 04. Overall validation of model ##########

#Sanity check
stopifnot(all(sort(AD_metadata$sample)==sort(unlist(samples_names(MOFAobject)))))

#Add sample metadata to the model
samples_metadata(MOFAobject) <- AD_metadata

#correlation between factors
plot_factor_cor(MOFAobject)

#Plot variance decomposition
plot_variance_explained(MOFAobject, max_r2 = 5)

#Total variance explained per view
plot_variance_explained(MOFAobject, plot_total = T)[[2]]

#Roughly check the clustering
library(ggplot2)
p <- plot_factors(MOFAobject, factors = c(1,2,3,4), color_by = "CDR_0_0.5_1", dot_size = 3,                   
                  show_missing = T)
p <- p + geom_hline(yintercept=-1.6, linetype="dashed") +   
  geom_vline(xintercept=(1), linetype = "dashed")
print(p)

#Check correlation with covariates
correlate_factors_with_covariates(MOFAobject, 
                                  covariates = 
                                    c("CDR_0_0.5_1","Sex_1_2","Age","Dx_0_1_2","Educat",
                                      "APOE_PN_1_2","bl_cogscore_MMSEz","FDG_4rois_PN"), plot ="r")

########## 05. Imputation ##########

#Get imputation
MOFAobject_im <- impute(MOFAobject)
get_imputed_data(MOFAobject_im)
#Test
MOFAobject_im@data$ip_proteomics[[1]][1:5,160:170]
MOFAobject_im@imputed_data$ip_proteomics[[1]][["mean"]][1:5,160:170]
get_imputed_data(MOFAobject_im)

#Export the imputed data to xlsx file
library("openxlsx")
imputed_ipad2_miRNA <- MOFAobject_im@imputed_data$ip_miRNA[[1]][["mean"]]
imputed_ipad2_miRNA
write.xlsx(imputed_ipad2_miRNA, sheetName="sheet1", file="imputed_ipad2_miRNA_new.xlsx")
imputed_ipad3_proteomics <- MOFAobject_im@imputed_data$ip_proteomics[[1]][["mean"]]
imputed_ipad3_proteomics
write.xlsx(imputed_ipad3_proteomics, sheetName="sheet1", file="imputed_ipad3_proteomics_new.xlsx")
imputed_ipad4_blood_biomarkers <- MOFAobject_im@imputed_data$ip_biomarkers[[1]][["mean"]]
imputed_ipad4_blood_biomarkers
write.xlsx(imputed_ipad4_blood_biomarkers, sheetName="sheet1", file="imputed_ipad4_blood_biomarkers_new.xlsx")

########## 06. Clustering ##########

#Get factor-values  
get_factors(MOFAobject, factors = 1)
get_factors(MOFAobject, factors = 2)
get_factors(MOFAobject, factors = 3)
get_factors(MOFAobject, factors = 4)

#make txt file using factor-values: factor1vs2 (1vs2 can be replaced by 1vs3, 1vs4, 2vs3, 2vs4, 3vs4)
factor1vs2 <- read.table("/home/dell/newinput/ipad_factor1vs2.txt")

#k-medoid clustering: factor1vs2 (1vs2 can be replaced by 1vs3, 1vs4, 2vs3, 2vs4, 3vs4)
library(cluster)
pam_k3 <- pam(factor1vs2, k = 3) #possible k range: 2<=k<=5
pam_k3
cluster_km3 <- pam_k3$cluster
cluster_km3
factor1vs2_km3 <- mutate(factor1vs2, cluster = cluster_km3)
factor1vs2_km3
ggplot(factor1vs2_km3, aes(x = Factor3, y = Factor4, color = factor(cluster), fill = factor(cluster))) +
  geom_point() + stat_ellipse(type = "t",geom = "polygon",alpha = 0.4) + 
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))

#Elbow Plot: factor1vs2 (1vs2 can be replaced by 1vs3, 1vs4, 2vs3, 2vs4, 3vs4)
library(purrr)
tot_withinss <-map_dbl(1:10,  function(k){
  model <- kmeans(x = factor1vs2, centers = k)
  model$tot.withinss
})
elbow_df <- data.frame(
  k = 1:10,
  tot_withinss = tot_withinss
)

#Plot the elbow plot: factor1vs2 (1vs2 can be replaced by 1vs3, 1vs4, 2vs3, 2vs4, 3vs4)
ggplot(elbow_df, aes(x = k, y = tot_withinss)) +
  geom_line() +
  scale_x_continuous(breaks = 1:10) + 
  scale_y_continuous(limits = c(min(elbow_df[,"tot_withinss"]), max(elbow_df[,"tot_withinss"])+0.1)) + 
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))
dev.off()

#Silhouette analysis
library(cluster)
plot(silhouette(pam_k3)) #possible k range: 2<=k<=5

#Silhouette plot curve: using average width: factor1vs2 (1vs2 can be replaced by 1vs3, 1vs4, 2vs3, 2vs4, 3vs4)
sil_width <- map_dbl(2:10,  function(k){
  model <- pam(x = factor1vs2, k = k)
  model$silinfo$avg.width
})

#Generate a data frame containing both k and sil_width
sil_df <- data.frame(
  k = 2:10,
  sil_width = sil_width
)

#Plot the relationship between k and sil_width
ggplot(sil_df, aes(x = k, y = sil_width)) +
  geom_line() +
  scale_x_continuous(breaks = 2:10) + 
  scale_y_continuous(limits = c(min(sil_df[,"sil_width"]), max(sil_df[,"sil_width"])+0.1)) + 
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))

########## 07. Down-stream analysis (in MOFA+) ##########

#Get factor-value dot graph
plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Factor1" 
)
plot_factor(MOFAobject, 
            factors = 2, 
            color_by = "Factor2" 
)

#Plot feature weights
plot_top_weights(MOFAobject_im, view = "ip_tgwas", 
                 factor = 1, nfeatures = 30, scale =T)
plot_top_weights(MOFAobject_im, view = "ip_miRNA", 
                 factor = 1, nfeatures = 10, scale =T)
plot_top_weights(MOFAobject_im, view = "ip_proteomics", 
                 factor = 1, nfeatures = 10, scale =T)
plot_top_weights(MOFAobject_im, view = "ip_biomarkers", 
                 factor = 1, nfeatures = 10, scale =T)
plot_top_weights(MOFAobject_im, view = "ip_tgwas", 
                 factor = 2, nfeatures = 30, scale =T)
plot_top_weights(MOFAobject_im, view = "ip_miRNA", 
                 factor = 2, nfeatures = 10, scale =T)
plot_top_weights(MOFAobject_im, view = "ip_proteomics", 
                 factor = 2, nfeatures = 10, scale =T)
plot_top_weights(MOFAobject_im, view = "ip_biomarkers", 
                 factor = 2, nfeatures = 10, scale =T)

#Get weight-values from each dataset 
get_weights(MOFAobject_im, view = "ip_tgwas", 
            factor = 1, scale =T)
get_weights(MOFAobject_im, view = "ip_miRNA", 
            factor = 1, scale =T)
get_weights(MOFAobject_im, view = "ip_proteomics", 
            factor = 1, scale =T)
get_weights(MOFAobject_im, view = "ip_biomarkers", 
            factor = 1, scale =T)
get_weights(MOFAobject_im, view = "ip_tgwas", 
            factor = 2, scale =T)
get_weights(MOFAobject_im, view = "ip_miRNA", 
            factor = 2, scale =T)
get_weights(MOFAobject_im, view = "ip_proteomics", 
            factor = 2, scale =T)
get_weights(MOFAobject_im, view = "ip_biomarkers", 
            factor = 2, scale =T)

