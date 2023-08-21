# 0. load packages --------------------------------------------------------
library(tidyverse)
library(readxl)
library(writexl)
library(ggsci)
library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

# 0. lazy load ------------------------------------------------------------
# load("00_variables/filtered.abun.RData")
# load("00_variables/sample.list.RData")
load("00_variables/transmission.species.RData")
species.in.disease <- as.data.frame(read_xlsx("00_variables/Transmission_species.xlsx",sheet = 1))


# 1. prepare directory and load data --------------------------------------
dir.create("04_Figure4/", showWarnings = F)
# dir.create("04_Figure4/01_input", showWarnings = F)


# 2. Plot -----------------------------------------------------------------
rownames(species.in.disease) <- species.in.disease$Species
species.in.disease <- species.in.disease[,-1]

species.in.disease.stats <- rowSums(species.in.disease!=0,na.rm =T) %>% as.data.frame()
colnames(species.in.disease.stats) <- "Freq"
species.in.disease.stats <- species.in.disease.stats %>% arrange(desc(Freq))

disease.stats <- colSums(species.in.disease!=0,na.rm =T) %>% as.data.frame()
colnames(disease.stats) <- "Freq"
disease.stats <- disease.stats %>% arrange(desc(Freq))

species.in.disease <- species.in.disease[rownames(species.in.disease.stats),]
species.in.disease <- species.in.disease[,rownames(disease.stats)]

bar1 <- HeatmapAnnotation(
  sum1 = anno_barplot(
    colSums(species.in.disease!=0,na.rm =T),
    bar_width = 0.9,
    gp = gpar(col = "white", fill = "#8988A3"),
    border = F,
    # axis_param = list(at = c(0,1.25e5,2.5e5),
    #                   labels = c("0","125k","250k")),
    height = unit(2, "cm")
  ), show_annotation_name = F
)

bar2 <- rowAnnotation(
  sum2 = anno_barplot(
    rowSums(species.in.disease!=0,na.rm =T),
    bar_width = 0.9,
    gp = gpar(col = "white", fill = "#8988A3"),
    border = F,
    # axis_param = list(at = c(0,2.5e5,5e5),
    #                   labels = c("0","250k","500k")),
    width = unit(2, "cm")
  ), show_annotation_name = F
)

col_fun = colorRamp2(c(0,2,4,6), c("#FAEBA8", "#EF906A","#BB415F","#741C5B"))

pdf("04_Figure4/Figure4.pdf", width = 12, height = 8)
Heatmap(as.matrix(species.in.disease),
        col = col_fun,
        na_col = "white",
        row_names_side = "left",
        cluster_rows = F,
        cluster_columns = F,
        rect_gp = gpar(col = "white", lwd = 2),
        name = "Nr. of the projects in which species enriched in diseases",
        top_annotation = bar1,
        right_annotation =bar2,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        width = ncol(species.in.disease) * unit(5, "mm"),
        height = nrow(species.in.disease) * unit(5, "mm")
)
dev.off()


temp <- intersect(rownames(species.in.disease),species.base.intervent.diff$Taxa[species.base.intervent.diff$group=="H2RA"])
temp <- species.in.disease[temp,]

temp1 <- as.data.frame(read_xlsx("00_variables/temp.xlsx",col_names = F))
colnames(temp1) <- c("species","project")

temp2 <- temp1[temp1$species %in% rownames(temp),]
unique(temp2$project) %>% length()


















