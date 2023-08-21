# 0. load packages --------------------------------------------------------
library(tidyverse)
library(readxl)
library(writexl)
library(ggsci)
library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(vegan)
library(ape)
library(phyloseq)

# 0. lazy load ------------------------------------------------------------
load("00_variables/sample.list.RData")
load("00_variables/filtered.abun.RData")

# 1. prepare directory and load data --------------------------------------
dir.create("02_Figure2", showWarnings = F)


# 2. Our result -----------------------------------------------------------
{ ## 2.1 Diff taxa ----
  
  { ### functions ----
    
    paired.test <- function(data, group, timepoint, level, type) {
      # 从数据中筛选出感兴趣的组别和时间点数据
      sample <- sample.list[(sample.list$group == group) & (sample.list$timepoint %in% timepoint) & (sample.list$type == type), c("sample", "id", "timepoint")]
      data <- data[[type]][, c("Taxa", colnames(data[[type]])[colnames(data[[type]]) %in% sample$sample])]
      
      # 把Taxa分割开来,按照想要的水平进行丰度的合并
      data <- separate(data, "Taxa", into = c("k", "p", "c", "o", "f", "g", "s", "t"), sep = "\\|")
      data <- aggregate(data[, -c(1:8)], by = list(data[, level]), sum)
      colnames(data)[1] <- "Taxa"
      taxa.list <- unique(data$Taxa)
      
      # 把data转置后与sample合并得到每个样本对应的个体id
      rownames(data) <- data$Taxa
      data <- data[, -1]
      data <- as.data.frame(t(data))
      data <- merge(data, sample, by.x = "row.names", by.y = "sample")
      data <- data[, c(1, ncol(data), ncol(data) - 1, 2:(ncol(data) - 2))]
      colnames(data)[1] <- "sample"
      
      # 一个菌一个菌来进行配对检验
      test.result <- data.frame()
      for (i in taxa.list) {
        subset <- data[, c("id", "timepoint", i)]
        colnames(subset)[3] <- "taxa"
        subset <- spread(subset, timepoint, taxa)
        colnames(subset) <- c("id", "time1", "time2")
        subset$diff <- subset$time1 - subset$time2 # 计算两组差值
        subset <- na.omit(subset)
        
        # 计算效应量 cohen's delta
        cohen <- effsize::cliff.delta(subset$time1, subset$time2)
        cohen.d <- cohen$estimate
        cohen.ci.lower <- cohen$conf.int[1]
        cohen.ci.upper <- cohen$conf.int[2]
        
        # 如果对于某一个菌而言，在所有样本和所有时间节点中都不存在时，跳过当前循环
        if (all(colSums(subset[, -1]) == 0)) {
          next
        }
        
        # 密度图 直观看是否符合正态分布
        # ggdensity(subset$diff,
        #           main = "Density plot of difference",
        #           xlab = "Difference")
        
        normal.test <- shapiro.test(subset$diff)
        
        if (normal.test$p.value < 0.05) {
          # 如果差值分布不符合正态分布，就使用配对wilcox检验
          res <- wilcox.test(subset$time1, subset$time2, paired = TRUE)
          test.result <- rbind(test.result, c(i, "wilcox.test", res$p.value, cohen.d, cohen.ci.lower, cohen.ci.upper))
        } else {
          # 如果差值分布符合正态分布，就使用配对t检验
          res <- t.test(subset$time1, subset$time2, paired = TRUE)
          test.result <- rbind(test.result, c(i, "t.test", res$p.value, cohen.d, cohen.ci.lower, cohen.ci.upper))
        }
      }
      colnames(test.result) <- c("Taxa", "test.method", "p.value", "cohen.d", "cohen.ci.lower", "cohen.ci.upper")
      
      test.result$p.value <- as.numeric(test.result$p.value)
      test.result$adj.p <- p.adjust(test.result$p.value, method = "BH")
      
      signif.diff.taxa <- test.result[test.result$p.value < 0.05, ]
      signif.diff.taxa <- signif.diff.taxa %>% arrange(p.value)
      signif.diff.taxa$cohen.d <- as.numeric(signif.diff.taxa$cohen.d)
      
      return(signif.diff.taxa)
    }
    lefse.analysis <- function(data, group, timepoint, level, type, lefse.plot.name = NULL) {
      # 从数据中筛选出感兴趣的组别和时间点数据
      sample <- sample.list[(sample.list$group == group) & (sample.list$timepoint %in% timepoint) & (sample.list$type == type), c("sample", "id", "timepoint")]
      data <- data[[type]][, c("Taxa", colnames(data[[type]])[colnames(data[[type]]) %in% sample$sample])]
      # 把Taxa分割开来,按照想要的水平进行丰度的合并
      data <- separate(data, "Taxa", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "t"), sep = "\\|")
      data <- unite(data, "Taxa", c(colnames(data)[1:which(colnames(data)[1:8] == level)]), sep = ";", remove = T)
      data <- data[, c("Taxa", sample$sample)]
      data <- aggregate(data[, -1], by = list(data$Taxa), sum)
      colnames(data)[1] <- "Taxa"
      rownames(data) <- data$Taxa
      data <- data[, -1]
      
      tax <- data.frame(Taxa = rownames(data))
      rownames(tax) <- tax$Taxa
      tax <- separate(tax, col = Taxa, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
      tax <- tax[, c(1:which(colnames(tax)[1:7] == level))]
      tax <- as.matrix(tax)
      
      rownames(sample) <- sample$sample
      
      ## -- calculate --
      OTU <- otu_table(data, taxa_are_rows = TRUE)
      TAX <- tax_table(tax)
      GROUP <- sample_data(sample[, "timepoint", drop = F])
      
      phylodat <- phyloseq(OTU, TAX, GROUP)
      try(markers <- microbiomeMarker::run_lefse(phylodat, group = "timepoint"), silent = T)
      
      if (exists("markers")) {
        marker <- markers@marker_table
        marker <- as.data.frame(matrix(unlist(marker), ncol = 5))
        colnames(marker) <- names(markers@marker_table)
        
        pattern <- tolower(str_sub(level, 1, 1))
        pattern <- paste0(pattern, "__")
        marker <- marker[str_detect(marker$feature, pattern), ]
        marker$pvalue <- as.numeric(marker$pvalue)
        marker <- marker %>% arrange(pvalue)
        marker <- separate(marker, "feature", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "t"), sep = "\\|")
        marker <- marker[, c(level, "enrich_group", "ef_lda", "pvalue")]
        marker$ef_lda <- as.numeric(marker$ef_lda)
        colnames(marker)[1] <- "Taxa"
        
        marker$ef_lda <- as.numeric(marker$ef_lda)
        marker$ef_lda[marker$enrich_group == unique(marker$enrich_group)[1]] <- marker$ef_lda[marker$enrich_group == unique(marker$enrich_group)[1]] * (-1)
        marker <- marker %>% arrange(ef_lda)
        marker$Taxa <- factor(marker$Taxa, levels = unique(marker$Taxa))
        p <- ggplot(marker, aes(ef_lda, Taxa, fill = enrich_group)) +
          geom_col(width = 0.7) +
          labs(x = "", y = NULL, fill = "Enriched group") +
          scale_x_continuous(expand = c(0, 0)) +
          theme_bw() +
          scale_fill_manual(values = brewer.pal(3, "Set2")[1:2])
        p
        if (!is.null(lefse.plot.name)) {
          assign(lefse.plot.name, p, envir = .GlobalEnv)
        }
        return(marker)
      } else {
        return(print("No markers have been found!"))
      }
    }
  }
  
  { ### species level ----
    
    ppi.base.intervent.statistical.species <- paired.test(data = filtered.abun, group = "PPI", timepoint = c("Baseline", "Intervention"), level = "s", type = "gut")

    ppi.base.intervent.lefse.species <- lefse.analysis(data = filtered.abun, group = "PPI", timepoint = c("Baseline", "Intervention"), level = "Species", type = "gut")

    # select intersection of the wilcox and lefse result
    ppi.base.intervent.species.intersect <- ppi.base.intervent.statistical.species[ppi.base.intervent.statistical.species$Taxa %in% ppi.base.intervent.lefse.species$Taxa, ]
    ppi.base.intervent.species.intersect$group <- "PPI"
    ppi.base.intervent.species.intersect$level <- "Species"

    species.base.intervent.diff <- rbind(
      ppi.base.intervent.species.intersect
    )
    species.base.intervent.diff$Taxa <- str_sub(species.base.intervent.diff$Taxa, 4, -1)
    species.base.intervent.diff$Taxa <- gsub("_", " ", species.base.intervent.diff$Taxa)
    
    rm(ppi.base.intervent.lefse.species, ppi.base.intervent.statistical.species, ppi.base.intervent.species.intersect)
  }
  
  { ### genus level ----
    
    ppi.base.intervent.statistical.genus <- paired.test(data = filtered.abun, group = "PPI", timepoint = c("Baseline", "Intervention"), level = "g", type = "gut")

    ppi.base.intervent.lefse.genus <- lefse.analysis(data = filtered.abun, group = "PPI", timepoint = c("Baseline", "Intervention"), level = "Genus", type = "gut")

    # select intersection of the wilcox and lefse result
    ppi.base.intervent.genus.intersect <- ppi.base.intervent.statistical.genus[ppi.base.intervent.statistical.genus$Taxa %in% ppi.base.intervent.lefse.genus$Taxa, ]
    ppi.base.intervent.genus.intersect$group <- "PPI"
    ppi.base.intervent.genus.intersect$level <- "Genus"

    genus.base.intervent.diff <- rbind(
      ppi.base.intervent.genus.intersect
    )
    genus.base.intervent.diff$Taxa <- str_sub(genus.base.intervent.diff$Taxa, 4, -1)
    genus.base.intervent.diff$Taxa <- gsub("_", " ", genus.base.intervent.diff$Taxa)
    
    rm(ppi.base.intervent.lefse.genus, ppi.base.intervent.statistical.genus, ppi.base.intervent.genus.intersect)
  }
  
  
}


# 3. Metacardis deconfound result -----------------------------------------
external.data <- as.data.frame(read_xlsx("00_variables/Bork_Nature_2021.xlsx", sheet = 2))
external.data <- external.data[str_detect(external.data$`Feature space`, "mOTU") & grepl("^PPI", external.data$Effector), ]
external.data <- external.data[!grepl("unclassified", external.data$`Feature display name`), ]
external.data <- external.data[!grepl("sp.", external.data$`Feature display name`), ]
external.data <- external.data[external.data$`Feature space` %in% c("mOTU, species","mOTU, genus"),]
external.data$`Feature space` <- gsub("mOTU, ", "", external.data$`Feature space`)
# external.data <- external.data[external.data$`Confounder status`=="SD",]
external.data <- external.data[, 3:ncol(external.data)]
colnames(external.data)[1:5] <- c("Disease", "Level", "Taxa", "FDR", "Effect.size")
external.data <- external.data[, c(1, 2, 3, 5)]
external.data <- external.data[external.data$Taxa != "Archaea", ]
colnames(external.data)[c(1, 4)] <- c("group", "cliff.d")
external.data <- na.omit(external.data)
external.data <- aggregate(external.data[,-c(1:3)],by=list(external.data$Level,external.data$Taxa),mean)
interest.species <- c(
  "Veillonella parvula", "Streptococcus sanguinis", "Streptococcus salivarius", "Streptococcus parasanguinis",
  "Streptococcus mutans", "Streptococcus australis", "Streptococcus anginosus", "Lactobacillus delbrueckii",
  "Haemophilus parainfluenzae", "Bifidobacterium longum", "Bifidobacterium dentium"
)
interest.genus <- c("Veillonella", "Streptococcus", "Rothia", "Lactobacillus", "Haemophilus", "Bacilli", "Actinomyces")
bork.list <- c(interest.species, interest.genus)
external.data <- external.data[external.data$Group.2 %in% bork.list, ]
external.data$Group <- "Forslund et al., 2021"
colnames(external.data)[1:3] <- c("Level","Taxa","Cliff.d")
external.data$Level <- Hmisc::capitalize(external.data$Level)
our.diff <- rbind(species.base.intervent.diff,genus.base.intervent.diff)
our.diff <- our.diff[,c(9,1,4,8)]
colnames(our.diff) <- colnames(external.data)
our.diff$Cliff.d <- our.diff$Cliff.d*(-1)

plot.data <- rbind(our.diff,external.data)
plot.data$color <- ifelse(plot.data$Cliff.d<0,"non-PPI use","PPI use")

plot.data$Group <- factor(plot.data$Group,levels=c("PPI","Forslund et al., 2021"))

plot.data.order <- plot.data %>% group_by(Taxa) %>% mutate(freq=n()) %>% arrange(freq)
# plot.data.order <- plot.data.order[,c(1,2,3,6)]
plot.data.order <- plot.data.order %>% arrange(freq,Group,Cliff.d)
plot.data.order <- plot.data.order[!duplicated(plot.data.order$Taxa),]
plot.data$Taxa <- factor(plot.data$Taxa,levels=plot.data.order$Taxa)

p1 <- plot.data %>% ggplot(aes(x=Cliff.d,y=Taxa,fill=color))+
  geom_col(width=0.3)+
  facet_grid(Level~Group,scales = "free_y",space = "free_y")+
  theme_bw()+
  scale_fill_manual(name="Group",values = c("#EF2326","#3952A4"),labels=c("non-PPI use","PPI use"))+
  theme(
    strip.background = element_rect(color = "black", fill = "#FFFFFF"),
    text = element_text(family = "serif"),
    axis.text.y = element_text(face = "italic")
  ) +
  labs(x="Effect size (Cliff's delta)")

p1
ggsave(p1,filename="02_Figure2/Figure2.pdf",width=8,height=8)

write_xlsx(plot.data,"SuppTable10.xlsx")







