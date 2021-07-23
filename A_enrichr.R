library(readxl)
library('biomaRt')
library(dplyr)
library(enrichR)
library(RDAVIDWebService)
library(openxlsx)
library(ggplot2)
library(extrafont)
source('helper_functions.R')

miRNA_table <- read_excel("109miRNA_traget_list_from_miRDBre_rev.xlsx")
miRNA_list <- colnames(miRNA_table)

############################# Fold_change_All_targets_수정.xlsx
org_table_ori <- read_excel("Fold_changes_All_targets_수정.xlsx", sheet = 2) # 1:M_TPAD, 2: M_IPAD

org_table <- org_table_ori[,c(1,2,10:12)]

for(mirna in miRNA_list){
  idx <- which(org_table$DB_numbers == mirna)
  mirna_gene <- na.omit(miRNA_table[,mirna])[[1]]
  temp_values <- - as.numeric(org_table[idx,3:5])
  temp <- cbind(mirna_gene, mirna, temp_values[1], temp_values[2], temp_values[3])
  temp2 <- data.frame(temp)
  names(temp2) <- names(org_table)
  org_table <- rbind(org_table, temp2)
  org_table <- org_table[-idx,]
}

unique_table <- list()
unique_list <- unique(org_table$targets)
for(target in unique_list){
  temp_uniq <- subset(org_table, org_table$targets == target)
  temp_uniq2 <- cbind(temp_uniq[1,1:2],mean(as.numeric(temp_uniq$Clu1_nor...10)),
                      mean(as.numeric(temp_uniq$Clu2_nor...11)),
                      mean(as.numeric(temp_uniq$Clu3_nor...12)))
  unique_table <- rbind(unique_table, temp_uniq2)
}

names(unique_table) <- names(org_table)

th_value <- 0.5 # the value of threshod for expression

clu1_up <- unique_table[unique_table$Clu1_nor...10 > th_value,]$targets
clu1_down <- unique_table[unique_table$Clu1_nor...10 < -th_value,]$targets

clu2_up <- unique_table[unique_table$Clu2_nor...11 > th_value,]$targets
clu2_down <- unique_table[unique_table$Clu2_nor...11 < -th_value,]$targets

clu3_up <- unique_table[unique_table$Clu3_nor...12 > th_value,]$targets
clu3_down <- unique_table[unique_table$Clu3_nor...12 < -th_value,]$targets


############################# Enriched pathway analysis
dbs_list <- listEnrichrDbs()

idx <- 1 ### 1: adjusted p-value, 2: p-value

dbs <- c("KEGG_2019_Human","GO_Molecular_Function_2018","GO_Cellular_Component_2018",
         "GO_Biological_Process_2018","BioCarta_2016","Reactome_2016")

cl1_up_list <- get_pathway(clu1_up, dbs, 0.05, idx)
cl1_down_list <- get_pathway(clu1_down, dbs, 0.05, idx)

plot_bar(cl1_up_list, cl1_down_list, 30, "T_Factors")

cl2_up_list <- get_pathway(clu2_up, dbs, 0.05, idx)
cl2_down_list <- get_pathway(clu2_down, dbs, 0.05, idx)
plot_bar(cl2_up_list, cl2_down_list, 30, "T_Factors")

cl3_up_list <- get_pathway(clu3_up, dbs, 0.05, idx)
cl3_down_list <- get_pathway(clu3_down, dbs, 0.05, idx)
plot_bar(cl3_up_list, cl3_down_list, 30, "T_Factors")


wb_T <- createWorkbook()
addWorksheet(wb_T, "Cl1_up")
addWorksheet(wb_T, "Cl1_down")
addWorksheet(wb_T, "Cl2_up")
addWorksheet(wb_T, "Cl2_down")
addWorksheet(wb_T, "Cl3_up")
addWorksheet(wb_T, "Cl3_down")

writeData(wb_T, "Cl1_up", cl1_up_list, rowNames = FALSE, colNames = TRUE)
writeData(wb_T, "Cl1_down", cl1_down_list, rowNames = FALSE, colNames = TRUE)
writeData(wb_T, "Cl2_up", cl2_up_list, rowNames = FALSE, colNames = TRUE)
writeData(wb_T, "Cl2_down", cl2_down_list, rowNames = FALSE, colNames = TRUE)
writeData(wb_T, "Cl3_up", cl3_up_list, rowNames = FALSE, colNames = TRUE)
writeData(wb_T, "Cl3_down", cl3_down_list, rowNames = FALSE, colNames = TRUE)

saveWorkbook(wb_T, file = "results_T.xlsx")



wb_I <- createWorkbook()
addWorksheet(wb_I, "Cl1_up")
addWorksheet(wb_I, "Cl1_down")
addWorksheet(wb_I, "Cl2_up")
addWorksheet(wb_I, "Cl2_down")
addWorksheet(wb_I, "Cl3_up")
addWorksheet(wb_I, "Cl3_down")

writeData(wb_I, "Cl1_up", cl1_up_list, rowNames = FALSE, colNames = TRUE)
writeData(wb_I, "Cl1_down", cl1_down_list, rowNames = FALSE, colNames = TRUE)
writeData(wb_I, "Cl2_up", cl2_up_list, rowNames = FALSE, colNames = TRUE)
writeData(wb_I, "Cl2_down", cl2_down_list, rowNames = FALSE, colNames = TRUE)
writeData(wb_I, "Cl3_up", cl3_up_list, rowNames = FALSE, colNames = TRUE)
writeData(wb_I, "Cl3_down", cl3_down_list, rowNames = FALSE, colNames = TRUE)

saveWorkbook(wb_I, file = "results_I.xlsx")




intersect(cl1_up_list$Term, cl1_down_list$Term)
intersect(cl2_up_list$Term, cl2_down_list$Term)
intersect(cl3_up_list$Term, cl3_down_list$Term)

intersect(cl1_down_list$Term, cl2_up_list$Term)
intersect(cl1_down_list$Term, cl3_up_list$Term)

intersect(cl2_up_list$Term, cl3_up_list$Term)
intersect(cl2_up_list$Term, cl3_down_list$Term)


############################# Plotting and figure export
loadfonts(device= "win")
windowsFonts()

data <- read_xlsx("TPAD_jcpark_filtered..xlsx")
data <- read_xlsx("IPAD_jcpark_filtered.xlsx")


ggsave("TPAD_reesults.png", dpi = 300, height=10, width=12, units = "in")
ggsave("IPAD_reesults.png", dpi = 300, height=10, width=12, units = "in")

S1 <- ggplot(data, aes(x= Subtype, y=Term, size=Overlap, color = Exp, group=Subtype)) +
  geom_point() + scale_color_manual(breaks = c("Up","Down"), values=c("red","blue")) +
  cowplot::theme_cowplot() +
  theme(axis.line = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.ticks = element_blank()) +
  ylab('') +
  theme(text = element_text(size = 13, family = "Arial")) +
  scale_x_discrete(labels = c("Clu 1", "Clu 2", "Clu 3"))
  
S1
dev.off()


