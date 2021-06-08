#Load packages
library(dplyr)  
library(psych)
library(tidyverse)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)

#Set working directory 
setwd("~/Desktop/Final_data/TMZ")

cat("===Uploading files===\n")
Response_TMZ <- read.delim("Data/TMZ_response.txt", header=TRUE)
RNA_seq_norm_portein_coding <- read.delim("Data/RNAseq.txt", header=TRUE, row.names = 1) 

cat("===Correlation analysis between TMZ-response (IC50) and RNA-seq data===\n")

#Select columns for correlation analysis
mat <- RNA_seq_norm_portein_coding[,]
survival_TMZ_IC50 <- as.data.frame(Response_TMZ[, 2])
colnames(survival_TMZ_IC50) <- "survival"

cat("Run Spearman correlation analysis without multiple testing correction\n")
mat_x <- as.data.frame(mat)
output_none_TMZ_IC50 <- corr.test(mat_x, survival_TMZ_IC50, method = "spearman", adjust = "none")

#Extract necessary data from correlation results of IC50 NOT corrected for multiple testing
output_none_TMZ_matrix_IC50 <- cbind(output_none_TMZ_IC50$r, output_none_TMZ_IC50$p)
colnames(output_none_TMZ_matrix_IC50) <- c("Correlation_Coefficient_r", "p_value")
output_none_TMZ_matrix_IC50 <- as.data.frame(output_none_TMZ_matrix_IC50)
output_none_TMZ_matrix_IC50$Gene_name <- rownames(output_none_TMZ_matrix_IC50)
rownames(output_none_TMZ_matrix_IC50) <- NULL
output_none_TMZ_matrix_IC50 <- output_none_TMZ_matrix_IC50[,c("Gene_name", "Correlation_Coefficient_r", "p_value")]

cat("Run Spearman correlation analysis corrected for multiple testing using Benjamini Hochberg\n")
mat_x <- as.data.frame(mat)
output_TMZ_IC50 <- corr.test(mat_x, survival_TMZ_IC50, method = "spearman", adjust = "BH")

#Extract necessary data from correlation results of IC50 corrected for multiple testing using Benjamini Hochberg
output_TMZ_matrix_IC50 <- cbind(output_TMZ_IC50$r, output_TMZ_IC50$p.adj,output_TMZ_IC50$p )
colnames(output_TMZ_matrix_IC50) <- c("Correlation_Coefficient_r", "p_value", "adjusted_p_value")
output_TMZ_matrix_IC50 <- as.data.frame(output_TMZ_matrix_IC50)
output_TMZ_matrix_IC50$Gene_name <- rownames(output_TMZ_matrix_IC50)
rownames(output_TMZ_matrix_IC50) <- NULL
output_TMZ_matrix_IC50 <- output_TMZ_matrix_IC50[,c("Gene_name", "Correlation_Coefficient_r", "p_value","adjusted_p_value")]

#Find significant protein_coding genes not corrected for multiple testing
significant_genes_TMZ_0.05_IC50 <- filter(output_none_TMZ_matrix_IC50, p_value <= 0.05)
significant_genes_TMZ_0.005_IC50 <- filter(output_none_TMZ_matrix_IC50, p_value <= 0.005)

#Dotchart of significant protein_coding genes not corrected for multiple testing
dotchart(significant_genes_TMZ_0.05_IC50$Correlation_Coefficient_r, xlab= "Correlation coefficient r")

#Correlation plot over all protein_coding genes not corrected for multiple testing
Correlation_plot_none_TMZ_IC50 <- ggplot(output_none_TMZ_matrix_IC50, aes(x=Correlation_Coefficient_r, y=p_value)) + 
  geom_point() + labs(x="Correlation Coefficient r", y="p-value") + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.5) +
  theme_bw()
ggsave("Results/IC50/Correlation_plot_none_TMZ_IC50.jpg")

#Correlation plot over all protein_coding genes corrected for multiple testing using Benjamini Hochberg
Correlation_plot_IC50 <- ggplot(output_TMZ_matrix_IC50, aes(x=Correlation_Coefficient_r, y=adjusted_p_value)) + 
  geom_point() + labs(x="Correlation Coefficient r", y="Adjusted p-value") + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.5) +
  theme_bw()
ggsave("Results/IC50/Correlation_plot_IC50_TMZ.jpg")

#Counts of protein_coding genes not corrected for multiple testing
All_genes_TMZ_IC50 <- nrow(output_none_TMZ_matrix_IC50)
Total_significant_genes_0.05_TMZ_IC50 <- nrow(significant_genes_TMZ_0.05_IC50)
Total_significant_genes_0.005_TMZ_IC50 <- nrow(significant_genes_TMZ_0.005_IC50)

cat("Total genes: ",All_genes_TMZ_IC50, sep="\n")
cat("Total significant genes (p<0.05): ",Total_significant_genes_0.05_TMZ_IC50, sep="\n")
cat("Total significant genes (p<0.005): ",Total_significant_genes_0.005_TMZ_IC50, sep="\n")

#Create count dataframe
Count_genes_TMZ_IC50 <- data.frame (Group  = c("all genes", "sign (p<0.05)", "sign (p<0.005)"),
                                    Count = c(All_genes_TMZ_IC50, Total_significant_genes_0.05_TMZ_IC50, Total_significant_genes_0.005_TMZ_IC50),
                                    Color = c("#00E100", "#E1BBBB", "#E10000")
)

#Create barplot
jpeg("Results/IC50/Barplot_TMZ_IC50.jpg")
A <- barplot(Count_genes_TMZ_IC50$Count, names.arg = Count_genes_TMZ_IC50$Group, col = Count_genes_TMZ_IC50$Color, ylab = "Amount of genes", main = "Significant genes treatment TMZ (IC50)")
text(x = A, y = Count_genes_TMZ_IC50$Count, label = Count_genes_TMZ_IC50$Count, pos = 3, cex = 1, col = "black")
dev.off()

cat("===Correlation analysis between TMZ-response (AUC) and RNA-seq data===\n")

#Select columns for correlation analysis
mat <- RNA_seq_norm_portein_coding[,]
survival_TMZ_AUC <- as.data.frame(Response_TMZ[, 3])
colnames(survival_TMZ_AUC) <- "survival"

cat("Run Spearman correlation analysis without multiple testing correction\n")
mat_x <- as.data.frame(mat)
output_none_TMZ_AUC <- corr.test(mat_x, survival_TMZ_AUC, method = "spearman", adjust = "none")

#Extract necessary data from correlation results of AUC NOT corrected for multiple testing
output_none_TMZ_matrix_AUC <- cbind(output_none_TMZ_AUC$r, output_none_TMZ_AUC$p)
colnames(output_none_TMZ_matrix_AUC) <- c("Correlation_Coefficient_r", "p_value")
output_none_TMZ_matrix_AUC <- as.data.frame(output_none_TMZ_matrix_AUC)
output_none_TMZ_matrix_AUC$Gene_name <- rownames(output_none_TMZ_matrix_AUC)
rownames(output_none_TMZ_matrix_AUC) <- NULL
output_none_TMZ_matrix_AUC <- output_none_TMZ_matrix_AUC[,c("Gene_name", "Correlation_Coefficient_r", "p_value")]

cat("Run Spearman correlation analysis corrected for multiple testing using Benjamini Hochberg\n")
mat_x <- as.data.frame(mat)
output_TMZ_AUC <- corr.test(mat_x, survival_TMZ_AUC, method = "spearman", adjust = "BH")

#Extract necessary data from correlation results of AUC corrected for multiple testing using Benjamini Hochberg
output_TMZ_matrix_AUC <- cbind(output_TMZ_AUC$r, output_TMZ_AUC$p.adj,output_TMZ_AUC$p )
colnames(output_TMZ_matrix_AUC) <- c("Correlation_Coefficient_r", "p_value", "adjusted_p_value")
output_TMZ_matrix_AUC <- as.data.frame(output_TMZ_matrix_AUC)
output_TMZ_matrix_AUC$Gene_name <- rownames(output_TMZ_matrix_AUC)
rownames(output_TMZ_matrix_AUC) <- NULL
output_TMZ_matrix_AUC <- output_TMZ_matrix_AUC[,c("Gene_name", "Correlation_Coefficient_r", "p_value", "adjusted_p_value")]

#Find significant protein_coding genes not corrected for multiple testing
significant_genes_TMZ_0.05_AUC <- filter(output_none_TMZ_matrix_AUC, p_value <= 0.05)
significant_genes_TMZ_0.005_AUC <- filter(output_none_TMZ_matrix_AUC, p_value <= 0.005)

#Dotchart of significant protein_coding genes not corrected for multiple testing
dotchart(significant_genes_TMZ_0.05_AUC$Correlation_Coefficient_r, xlab= "Correlation coefficient r")

#Correlation plot over all protein_coding genes not corrected for multiple testing
Correlation_plot_none_TMZ_AUC <- ggplot(output_none_TMZ_matrix_AUC, aes(x=Correlation_Coefficient_r, y=p_value)) + 
  geom_point() + labs(x="Correlation Coefficient r", y="p-value") + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.5) +
  theme_bw()
ggsave("Results/AUC/Correlation_plot_none_TMZ_AUC.jpg")

#Correlation plot over all protein_coding genes corrected for multiple testing using Benjamini Hochberg
Correlation_plot_AUC <- ggplot(output_TMZ_matrix_AUC, aes(x=Correlation_Coefficient_r, y=adjusted_p_value)) + 
  geom_point() + labs(x="Correlation Coefficient r", y="Adjusted p-value") + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.5) +
  theme_bw()
ggsave("Results/AUC/Correlation_plot_AUC_TMZ.jpg")

#Counts of protein_coding genes not corrected for multiple testing
All_genes_TMZ_AUC <- nrow(output_none_TMZ_matrix_AUC)
Total_significant_genes_0.05_TMZ_AUC <- nrow(significant_genes_TMZ_0.05_AUC)
Total_significant_genes_0.005_TMZ_AUC <- nrow(significant_genes_TMZ_0.005_AUC)

cat("Total genes:",All_genes_TMZ_AUC, sep="\n")
cat("Total significant genes (p<0.05):",Total_significant_genes_0.05_TMZ_AUC, sep="\n")
cat("Total significant genes (p<0.005):",Total_significant_genes_0.005_TMZ_AUC, sep="\n")

#Create count dataframe
Count_genes_TMZ_AUC <- data.frame (Group  = c("all genes", "sign (p<0.05)", "sign (p<0.005)"),
                                   Count = c(All_genes_TMZ_AUC, Total_significant_genes_0.05_TMZ_AUC, Total_significant_genes_0.005_TMZ_AUC),
                                   Color = c("#00E100", "#E1BBBB", "#E10000")
)

#Create barplot
jpeg("Results/AUC/Barplot_TMZ_AUC.jpg")
B <- barplot(Count_genes_TMZ_AUC$Count, names.arg = Count_genes_TMZ_AUC$Group, col = Count_genes_TMZ_AUC$Color, ylab = "Amount of genes", main = "Significant genes treatment TMZ (AUC)")
text(x = B, y = Count_genes_TMZ_AUC$Count, label = Count_genes_TMZ_AUC$Count, pos = 3, cex = 1, col = "black")
dev.off()

cat("===Correlation analysis between TMZ-response (100 uM) and RNA-seq data===\n")

#Select columns for correlation analysis
mat <- RNA_seq_norm_portein_coding[,]
survival_TMZ_100uM <- as.data.frame(Response_TMZ[, 4])
colnames(survival_TMZ_100uM) <- "survival"

cat("Run Spearman correlation analysis without multiple testing correction\n")
mat_x <- as.data.frame(mat)
output_none_TMZ_100uM <- corr.test(mat_x, survival_TMZ_100uM, method = "spearman", adjust = "none")

#Extract necessary data from correlation results of 100uM NOT corrected for multiple testing
output_none_TMZ_matrix_100uM <- cbind(output_none_TMZ_100uM$r, output_none_TMZ_100uM$p)
colnames(output_none_TMZ_matrix_100uM) <- c("Correlation_Coefficient_r", "p_value")
output_none_TMZ_matrix_100uM <- as.data.frame(output_none_TMZ_matrix_100uM)
output_none_TMZ_matrix_100uM$Gene_name <- rownames(output_none_TMZ_matrix_100uM)
rownames(output_none_TMZ_matrix_100uM) <- NULL
output_none_TMZ_matrix_100uM <- output_none_TMZ_matrix_100uM[,c("Gene_name", "Correlation_Coefficient_r", "p_value")]

cat("Run Spearman correlation analysis corrected for multiple testing using Benjamini Hochberg\n")
mat_x <- as.data.frame(mat)
output_TMZ_100uM <- corr.test(mat_x, survival_TMZ_100uM, method = "spearman", adjust = "BH")

#Extract necessary data from correlation results of 100uM corrected for multiple testing using Benjamini Hochberg
output_TMZ_matrix_100uM <- cbind(output_TMZ_100uM$r, output_TMZ_100uM$p.adj,output_TMZ_100uM$p )
colnames(output_TMZ_matrix_100uM) <- c("Correlation_Coefficient_r", "p_value", "adjusted_p_value")
output_TMZ_matrix_100uM <- as.data.frame(output_TMZ_matrix_100uM)
output_TMZ_matrix_100uM$Gene_name <- rownames(output_TMZ_matrix_100uM)
rownames(output_TMZ_matrix_100uM) <- NULL
output_TMZ_matrix_100uM <- output_TMZ_matrix_100uM[,c("Gene_name", "Correlation_Coefficient_r", "p_value", "adjusted_p_value")]

#Find significant protein_coding genes not corrected for multiple testing
significant_genes_TMZ_0.05_100uM <- filter(output_none_TMZ_matrix_100uM, p_value <= 0.05)
significant_genes_TMZ_0.005_100uM <- filter(output_none_TMZ_matrix_100uM, p_value <= 0.005)

#Dotchart of significant protein_coding genes not corrected for multiple testing
dotchart(significant_genes_TMZ_0.05_100uM$Correlation_Coefficient_r, xlab= "Correlation coefficient r")

#Correlation plot over all protein_coding genes not corrected for multiple testing
Correlation_plot_none_TMZ_100uM <- ggplot(output_none_TMZ_matrix_100uM, aes(x=Correlation_Coefficient_r, y=p_value)) + 
  geom_point() + labs(x="Correlation Coefficient r", y="p-value") + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.5) +
  theme_bw()
ggsave("Results/100uM/Correlation_plot_none_TMZ_100uM.jpg")

#Correlation plot over all protein_coding genes corrected for multiple testing using Benjamini Hochberg
Correlation_plot_100uM <- ggplot(output_TMZ_matrix_100uM, aes(x=Correlation_Coefficient_r, y=adjusted_p_value)) + 
  geom_point() + labs(x="Correlation Coefficient r", y="Adjusted p-value") + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.5) +
  theme_bw()
ggsave("Results/100uM/Correlation_plot_100uM_TMZ.jpg")

#Counts of protein_coding genes not corrected for multiple testing
All_genes_TMZ_100uM <- nrow(output_none_TMZ_matrix_100uM)
Total_significant_genes_0.05_TMZ_100uM <- nrow(significant_genes_TMZ_0.05_100uM)
Total_significant_genes_0.005_TMZ_100uM <- nrow(significant_genes_TMZ_0.005_100uM)

cat("Total genes:",All_genes_TMZ_100uM, sep="\n")
cat("Total significant genes (p<0.05):",Total_significant_genes_0.05_TMZ_100uM, sep="\n")
cat("Total significant genes (p<0.005):",Total_significant_genes_0.005_TMZ_100uM, sep="\n")

#Create count dataframe
Count_genes_TMZ_100uM <- data.frame (Group  = c("all genes", "sign (p<0.05)", "sign (p<0.005)"),
                                     Count = c(All_genes_TMZ_100uM, Total_significant_genes_0.05_TMZ_100uM, Total_significant_genes_0.005_TMZ_100uM),
                                     Color = c("#00E100", "#E1BBBB", "#E10000")
)

#Create barplot
jpeg("Results/100uM/Barplot_TMZ_100uM.jpg")
C <- barplot(Count_genes_TMZ_100uM$Count, names.arg = Count_genes_TMZ_100uM$Group, col = Count_genes_TMZ_100uM$Color, ylab = "Amount of genes", main = "Significant genes treatment TMZ (100uM)")
text(x = C, y = Count_genes_TMZ_100uM$Count, label = Count_genes_TMZ_100uM$Count, pos = 3, cex = 1, col = "black")
dev.off()

cat("===Creating a Venn Diagram===\n")

# Prepare a palette of 3 colors with R colorbrewer:
myCol <- brewer.pal(3, "Pastel2")

#Make a Venndiagram
venn.diagram(
  x = list(significant_genes_TMZ_0.05_IC50$Gene_name, significant_genes_TMZ_0.05_AUC$Gene_name, significant_genes_TMZ_0.05_100uM$Gene_name),
  category.names = c("p<0.05 IC50" , "p<0.05 AUC", "p<0.05 100uM"),
  filename = 'Results/Venndiagram_TMZ.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

#Save results 
write.csv(significant_genes_TMZ_0.005_IC50,"Results/IC50\\significant_genes_TMZ_0.005_IC50.csv", row.names = FALSE)
write.csv(significant_genes_TMZ_0.05_IC50,"Results/IC50\\significant_genes_TMZ_0.05_IC50.csv", row.names = FALSE)
write.csv(significant_genes_TMZ_0.005_AUC,"Results/AUC\\significant_genes_TMZ_0.005_AUC.csv", row.names = FALSE)
write.csv(significant_genes_TMZ_0.05_AUC,"Results/AUC\\significant_genes_TMZ_0.05_AUC.csv", row.names = FALSE)
write.csv(significant_genes_TMZ_0.005_100uM,"Results/100uM\\significant_genes_TMZ_0.005_100uM.csv", row.names = FALSE)
write.csv(significant_genes_TMZ_0.05_100uM,"Results/100uM\\significant_genes_TMZ_0.05_100uM.csv", row.names = FALSE)

cat("===Correlation analysis is finished. Thank you for waiting===\n")