library(edgeR)
library(pheatmap)
library(Rtsne)
library(ggplot2)
library(ggrepel)

# road data
data_file <- "/Users/song-yuseog/Desktop/R_script/K562_DCI/Chromatin_Contact.Normal-Philadelphia.merged.csv"

# CSV --> edgeR data object
data <- read.csv(data_file, header = TRUE, row.names = 1)

sample_groups <- c("Normal", "Normal", "Philadelphia", "Philadelphia")

# Create edgeR data object
dge <- DGEList(counts = data, group = sample_groups)

# Normalization 
dge <- calcNormFactors(dge)

# estimate dispersion
dge <- estimateCommonDisp(dge)
# Differential gene expression
dge <- estimateTagwiseDisp(dge)

# 1.COMMON DISPERSION APPROACH
# Fisher's exact test --> significantly different genes

dge.com <- exactTest(dge)
CDtop10<-topTags(dge.com, n = 10L)$table
#CDtop20<-topTags(dge.com, n = 20L)$table
#CDtop30<-topTags(dge.com, n = 30L)$table
#CDtop40<-topTags(dge.com, n = 40L)$table

ann_col <- data.frame(type = as.factor(sample_groups)) # sample의 type  나타내기

# Heatmap
heat <- pheatmap(dge, 
                 border_color = NA, 
                 cluster_rows = TRUE,
                 cluster_cols = TRUE, 
                 main = "Differential Chromatin Interaction: Normal vs Philadelphia chromosome",
                 fontsize_row = 10, 
                 color = colorRampPalette(c("forestgreen", "black", "red"))(100),
                 width = 1, angle_col = 315)


png(filename = "/Users/song-yuseog/Desktop/R_script/K562_DCI/DCI_N_vs_P.heatmap.png", width = 600, height = 960) # 파일 이름과 크기 지정
print(heat) 
dev.off() 

# Extract log2 fold changes and -log10(p-values) from the exact test results
logFC <- CDtop10$logFC
logCPM <- CDtop10$logCPM
p_values <- CDtop10$PValue
adjusted_p_values <- CDtop10$FDR
gene_names <- rownames(CDtop10)

# Create a data frame for the volcano plot
volcano_data <- data.frame(Gene = gene_names, logFC = logFC, logCPM = logCPM, PValue = p_values)
volcano_data$FDR <- p.adjust(volcano_data$PValue, method = "BH")  # adjust FDR
volcano_data$significant <- ifelse(volcano_data$FDR < 0.05 & abs(volcano_data$logFC) > 4, "yes", "no")

# Volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(FDR), color = significant)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("yes" = "red", "no" = "black")) +
  labs(title = "Volcano Plot",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted P-value)") +
  theme_minimal() +
  geom_text_repel(
    data = subset(volcano_data, significant == "yes"), 
    aes(label = Gene),
    box.padding = 0.5,  # Adjust the padding around labels
    point.padding = 1,  # Adjust the distance from points
    segment.color = "grey",  # Color of the line segments
    segment.size = 0.2,  # Size of the line segments
    hjust = 1, vjust = -0.5  # Adjust the horizontal and vertical justification
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

# Save the volcano plot as a PNG file
png(filename = "/Users/song-yuseog/Desktop/R_script/K562_DCI/DCI_N_vs_P.volcano.png", width = 600, height = 400)
print(volcano_plot)
dev.off()
