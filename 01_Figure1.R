# 0. load packages --------------------------------------------------------
library(tidyverse) # Easily Install and Load the 'Tidyverse'
library(readxl) # Read Excel Files
library(writexl) # Export Data Frames to Excel 'xlsx' Format
library(RColorBrewer) # ColorBrewer Palettes
library(ggpubr) # 'ggplot2' Based Publication Ready Plots # 'ggplot2' Based Publication Ready Plots
library(patchwork) # The Composer of Plots
library(SIAMCAT) # Statistical Inference of Associations between Microbial Communities And host phenoTypes
library(ggsci) # Scientific Journal and Sci-Fi Themed Color Palettes for 'ggplot2'

# 0. lazy load ------------------------------------------------------------
load("00_variables/sample.list.RData")
load("00_variables/filtered.abun.RData")
metadata <- as.data.frame(read_xlsx("00_variables/metadata_newest.xlsx"))

# 0. prepare directory and load data --------------------------------------
# dir.create("00_supp_tables", showWarnings = F)
dir.create("01_Figure1", showWarnings = F)
dir.create("01_Figure1/01_input", showWarnings = F)

# 2. B: Oral and gut composition diff ----------------------------------------
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
        cohen <- effsize::cohen.d(subset$time1, subset$time2, paired = T)
        # cohen <- effsize::cliff.delta(subset$time1, subset$time2)
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
    h2ra.base.intervent.statistical.species <- paired.test(data = filtered.abun, group = "H2RA", timepoint = c("Baseline", "Intervention"), level = "s", type = "gut")

    ppi.base.intervent.lefse.species <- lefse.analysis(data = filtered.abun, group = "PPI", timepoint = c("Baseline", "Intervention"), level = "Species", type = "gut")
    h2ra.base.intervent.lefse.species <- lefse.analysis(data = filtered.abun, group = "H2RA", timepoint = c("Baseline", "Intervention"), level = "Species", type = "gut")

    # select intersection of the wilcox and lefse result
    ppi.base.intervent.species.intersect <- ppi.base.intervent.statistical.species[ppi.base.intervent.statistical.species$Taxa %in% ppi.base.intervent.lefse.species$Taxa, ]
    ppi.base.intervent.species.intersect$group <- "PPI"
    ppi.base.intervent.species.intersect$level <- "Species"
    h2ra.base.intervent.species.intersect <- h2ra.base.intervent.statistical.species[h2ra.base.intervent.statistical.species$Taxa %in% h2ra.base.intervent.lefse.species$Taxa, ]
    h2ra.base.intervent.species.intersect$group <- "H2RA"
    h2ra.base.intervent.species.intersect$level <- "Species"

    species.base.intervent.diff <- rbind(
      ppi.base.intervent.species.intersect, h2ra.base.intervent.species.intersect
    )
    species.base.intervent.diff$Taxa <- str_sub(species.base.intervent.diff$Taxa, 4, -1)
    species.base.intervent.diff$Taxa <- gsub("_", " ", species.base.intervent.diff$Taxa)

    rm(ppi.base.intervent.lefse.species, ppi.base.intervent.statistical.species, ppi.base.intervent.species.intersect)
    rm(h2ra.base.intervent.lefse.species, h2ra.base.intervent.statistical.species, h2ra.base.intervent.species.intersect)
  }

  { ### genus level ----

    ppi.base.intervent.statistical.genus <- paired.test(data = filtered.abun, group = "PPI", timepoint = c("Baseline", "Intervention"), level = "g", type = "gut")
    h2ra.base.intervent.statistical.genus <- paired.test(data = filtered.abun, group = "H2RA", timepoint = c("Baseline", "Intervention"), level = "g", type = "gut")

    ppi.base.intervent.lefse.genus <- lefse.analysis(data = filtered.abun, group = "PPI", timepoint = c("Baseline", "Intervention"), level = "Genus", type = "gut")
    h2ra.base.intervent.lefse.genus <- lefse.analysis(data = filtered.abun, group = "H2RA", timepoint = c("Baseline", "Intervention"), level = "Genus", type = "gut")

    # select intersection of the wilcox and lefse result
    ppi.base.intervent.genus.intersect <- ppi.base.intervent.statistical.genus[ppi.base.intervent.statistical.genus$Taxa %in% ppi.base.intervent.lefse.genus$Taxa, ]
    ppi.base.intervent.genus.intersect$group <- "PPI"
    ppi.base.intervent.genus.intersect$level <- "Genus"
    h2ra.base.intervent.genus.intersect <- h2ra.base.intervent.statistical.genus[h2ra.base.intervent.statistical.genus$Taxa %in% h2ra.base.intervent.lefse.genus$Taxa, ]
    h2ra.base.intervent.genus.intersect$group <- "H2RA"
    h2ra.base.intervent.genus.intersect$level <- "Genus"

    genus.base.intervent.diff <- rbind(
      ppi.base.intervent.genus.intersect, h2ra.base.intervent.genus.intersect
    )
    genus.base.intervent.diff$Taxa <- str_sub(genus.base.intervent.diff$Taxa, 4, -1)
    genus.base.intervent.diff$Taxa <- gsub("_", " ", genus.base.intervent.diff$Taxa)

    rm(ppi.base.intervent.lefse.genus, ppi.base.intervent.statistical.genus, ppi.base.intervent.genus.intersect)
    rm(h2ra.base.intervent.lefse.genus, h2ra.base.intervent.statistical.genus, h2ra.base.intervent.genus.intersect)
  }

  { ### export the content of supplementary table ----

    supp.table.1 <- rbind(species.base.intervent.diff, genus.base.intervent.diff)
    supp.table.1 <- supp.table.1[, c(1, 2, 3, 4, 8, 9)]
    colnames(supp.table.1) <- c("Taxon", "Statistical testing method", "P value", "Cohen's d", "Group", "Taxon rank")
    write_xlsx(supp.table.1, "SuppTable8.xlsx")
  }
}

{ ## 2.2 Plot ----

  { ### p1: heatmap of diff taxa ----

    plot.data.1 <- rbind(species.base.intervent.diff, genus.base.intervent.diff)
    plot.data.1 <- plot.data.1[, c(1, 4, 8, 9)]

    # taxa order by group and cohen.d
    x.axis.order <- plot.data.1 %>%
      group_by(Taxa) %>%
      mutate(num = n()) %>%
      arrange(desc(num), group, cohen.d) %>%
      select(Taxa) %>%
      unique()

    plot.data.1$Taxa <- factor(plot.data.1$Taxa, levels = x.axis.order$Taxa)
    plot.data.1$group <- factor(plot.data.1$group, levels = rev(c("PPI", "H2RA")))

    p1 <- plot.data.1 %>% ggplot(aes(x = Taxa, y = group, fill = cohen.d)) +
      geom_tile(aes(width = 0.6, height = 0.6)) +
      scale_fill_gradient2(name = "Cohen's d", low = "red", high = "#006837", mid = "white") +
      facet_grid(. ~ level, space = "free", scales = "free") +
      labs(x = "", y = "") +
      theme_bw() +
      theme(
        strip.background.x = element_rect(color = "black", fill = "white", linetype = "solid"),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle=90),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        text = element_text(family = "serif")
      )
    p1
  }

  { ### p2: oral bacteria in eHOMD ----

    ehomd <- as.data.frame(read_xlsx("00_variables/oral_bacteria_in_eHOMD.xlsx"))
    ehomd$eHOMD[ehomd$eHOMD == 0] <- "No"
    ehomd$eHOMD[ehomd$eHOMD == 1] <- "Yes"
    ehomd$Group <- "eHOMD"
    ehomd$Taxon <- factor(ehomd$Taxon, levels = x.axis.order$Taxa)

    p2 <- ehomd %>% ggplot(aes(x = Taxon, y = Group, fill = eHOMD)) +
      geom_tile(aes(width = 0.6, height = 0.6)) +
      scale_fill_manual(name = "Oral bacteria in eHOMD", values = c("white", "#CBB1D5")) +
      facet_grid(. ~ Rank, scales = "free", space = "free") +
      theme_bw() +
      labs(x = "", y = "") +
      theme(
        strip.background.x = element_rect(color = "black", fill = "white", linetype = "solid"),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        text = element_text(family = "serif")
      )
    p2
  }

  { ### p3: Taxa number statistics ----

    plot.data.2 <- plot.data.1
    plot.data.2$signal <- ifelse(plot.data.2$cohen.d > 0, "Baseline", "Intervention")
    plot.data.2 <- as.data.frame(table(plot.data.2[plot.data.2$level == "Species", c("group", "signal")]))
    plot.data.2$group <- as.character(plot.data.2$group)
    plot.data.2$signal <- as.character(plot.data.2$signal)
    plot.data.2 <- rbind(plot.data.2, c("eHOMD", "Yes", length(ehomd$eHOMD[ehomd$eHOMD == "Yes" & ehomd$Rank == "Species"])))
    plot.data.2$Freq <- as.numeric(plot.data.2$Freq)

    color <- c("#E67764", "#7FA493", "#CBB1D5")
    names(color) <- c("Intervention", "Baseline", "Yes")

    p3 <- plot.data.2 %>% ggplot(aes(x = group, y = Freq, fill = signal)) +
      geom_bar(stat = "identity", position = "stack", width = 0.3) +
      theme_bw() +
      labs(x = "", y = "Number of species") +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "right",
        text = element_text(family = "serif")
      ) +
      scale_y_continuous(position="right")+
      scale_fill_manual(name = "Enriched in", values = color) +
      coord_flip()
    p3
  }

  { ### p4: Prevalence change ----

    plot.data.3 <- filtered.abun[["gut"]]
    plot.data.3 <- separate(plot.data.3, "Taxa", into = c("k", "p", "c", "o", "f", "g", "s", "t"), sep = "\\|")
    plot.data.3.genus <- aggregate(plot.data.3[, -c(1:8)], by = list(plot.data.3[, "g"]), sum)
    plot.data.3.genus$Group.1 <- str_sub(plot.data.3.genus$Group.1, 4, -1)
    plot.data.3.genus$Group.1 <- gsub("_", " ", plot.data.3.genus$Group.1)
    plot.data.3.species <- aggregate(plot.data.3[, -c(1:8)], by = list(plot.data.3[, "s"]), sum)
    plot.data.3.species$Group.1 <- str_sub(plot.data.3.species$Group.1, 4, -1)
    plot.data.3.species$Group.1 <- gsub("_", " ", plot.data.3.species$Group.1)

    plot.data.3.species <- plot.data.3.species[plot.data.3.species$Group.1 %in% plot.data.1$Taxa, ]
    plot.data.3.genus <- plot.data.3.genus[plot.data.3.genus$Group.1 %in% plot.data.1$Taxa, ]
    plot.data.3.species <- gather(plot.data.3.species, sample, abundance, -Group.1)
    plot.data.3.genus <- gather(plot.data.3.genus, sample, abundance, -Group.1)

    plot.data.3.species <- merge(plot.data.3.species, sample.list[, c(1, 4, 5)], by = "sample")
    plot.data.3.genus <- merge(plot.data.3.genus, sample.list[, c(1, 4, 5)], by = "sample")
    plot.data.3.species <- plot.data.3.species[(plot.data.3.species$group %in% c("PPI", "H2RA")) & (plot.data.3.species$timepoint %in% c("Baseline", "Intervention")), ]
    plot.data.3.genus <- plot.data.3.genus[(plot.data.3.genus$group %in% c("PPI", "H2RA")) & (plot.data.3.genus$timepoint %in% c("Baseline", "Intervention")), ]
    plot.data.3.species$level <- "Species"
    plot.data.3.genus$level <- "Genus"
    plot.data.3 <- rbind(plot.data.3.species, plot.data.3.genus)
    colnames(plot.data.3)[2] <- "taxa"

    group.stats <- sample.list[sample.list$type == "gut", ] %>%
      group_by(group, timepoint) %>%
      summarise(num = n())
    plot.data.3 <- plot.data.3 %>%
      group_by(taxa, group, timepoint, level) %>%
      summarise(freq = sum(abundance != 0))

    plot.data.3 <- merge(plot.data.3, group.stats, by = c("group", "timepoint"))
    plot.data.3$prevalence <- plot.data.3$freq / plot.data.3$num * 100

    plot.data.3 <- plot.data.3[, -c(5, 6)]
    plot.data.3 <- spread(plot.data.3, timepoint, prevalence)
    plot.data.3$prev.change <- plot.data.3$Intervention - plot.data.3$Baseline

    plot.data.3$taxa <- factor(plot.data.3$taxa, levels = x.axis.order$Taxa)

    p4 <- plot.data.3 %>% ggplot(aes(x = taxa, y = prev.change, color = group)) +
      geom_point(size = 3.5) +
      scale_color_manual(name = "Group", labels = c("H2RA", "PPI"), values = c("#709AE1", "#FED439")) +
      theme_bw() +
      facet_grid(. ~ level, scales = "free", space = "free") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme(
        strip.background.x = element_rect(color = "black", fill = "white", linetype = "solid"),
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(face = "italic", angle = 90, hjust = 1, family = "serif"),
        legend.position = "top"
      ) +
      scale_y_continuous(labels = scales::percent_format(scale = 1)) +
      labs(y = "Prevalence change (Intervention - Baseline)", x = "")
    p4
  }

  p5 <- ((((p1 / p2) | p3) + plot_layout(widths = c(2, 0.2))) / ((p4 | p3) + plot_layout(widths = c(2, 0.2)))) + plot_layout(heights = c(1, 0.3))
  p5
  ggsave(p5, filename = "01_Figure1/Figure1B.pdf", width = 12, height = 8)
}


# 3. C: Abundance of oral species in the gut --------------------------------------------------------------------
oral.abun <- filtered.abun[["oral"]]
oral.abun <- separate(oral.abun, "Taxa", into = c("k", "p", "c", "o", "f", "g", "s", "t"), sep = "\\|")
oral.abun <- aggregate(oral.abun[, -c(1:8)], by = list(oral.abun[, "s"]), sum)
colnames(oral.abun)[1] <- "Taxa"
oral.abun$Prevalence <- apply(oral.abun[,-1],1,function(x){sum(x!=0)/length(x)*100})
oral.abun <- oral.abun[oral.abun$Prevalence>10,]

gut.abun <- filtered.abun[["gut"]]
gut.abun <- separate(gut.abun, "Taxa", into = c("k", "p", "c", "o", "f", "g", "s", "t"), sep = "\\|")
gut.abun <- aggregate(gut.abun[, -c(1:8)], by = list(gut.abun[, "s"]), sum)
colnames(gut.abun)[1] <- "Taxa"

ppi.gut.abun <- gut.abun[,c("Taxa",sample.list$sample[sample.list$group=="PPI"&sample.list$type=="gut"])]
ppi.gut.abun <- ppi.gut.abun[rowSums(ppi.gut.abun[,-1])!=0,]

h2ra.gut.abun <- gut.abun[,c("Taxa",sample.list$sample[sample.list$group=="H2RA"&sample.list$type=="gut"])]
h2ra.gut.abun <- h2ra.gut.abun[rowSums(h2ra.gut.abun[,-1])!=0,]

intersect.species.ppi <- intersect(oral.abun$Taxa,ppi.gut.abun$Taxa)
intersect.species.h2ra <- intersect(oral.abun$Taxa,h2ra.gut.abun$Taxa)

ppi.gut.abun <- ppi.gut.abun[ppi.gut.abun$Taxa %in% intersect.species.ppi,]
h2ra.gut.abun <- h2ra.gut.abun[h2ra.gut.abun$Taxa %in% intersect.species.h2ra,]

ppi.gut.abun.1 <- gather(ppi.gut.abun,Sample,Abundance,-Taxa)
h2ra.gut.abun.1 <- gather(h2ra.gut.abun,Sample,Abundance,-Taxa)

{### export the content of supplementary table ----
  
  table4 <- rbind(ppi.gut.abun.1,h2ra.gut.abun.1)
  table4 <- merge(table4,sample.list[,c(1,4,5)],by.x="Sample",by.y="sample")
  colnames(table4) <- c("Sample","Taxa","Abundance in the gut","Group","Timepoint")
  write_xlsx(table4,"SuppTable4.xlsx")
  
}

plot.data <- table4 %>% group_by(Sample) %>% summarise(Abundance=sum(`Abundance in the gut`))
plot.data <- merge(plot.data,sample.list[,c(1,3,4,5)],by.x="Sample",by.y="sample")

plot.data$group <- factor(plot.data$group,levels=c("PPI","H2RA"))

p1 <- plot.data %>% ggplot(aes(x=timepoint,y=Abundance,fill=timepoint))+
  geom_boxplot(width=0.3)+
  geom_signif(comparisons = list(c("Baseline","Intervention")),map_signif_level = T,step_increase = 0.2)+
  facet_grid(.~group)+
  theme_bw()+
  scale_fill_manual(name = "Group", labels = c("H2RA", "PPI"), values = c("#709AE1", "#FED439")) +
  theme(
    strip.background.x = element_rect(color = "black", fill = "white", linetype = "solid"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(family = "serif")  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(y = "Accumulated abundance of the oral bacteria in the gut", x = "")
p1

ggsave(p1,filename="01_Figure1/Figure1C.pdf",width=6,height=6)

# 统计每个人增加了多少丰度
plot.data.1 <- spread(plot.data[,-1],timepoint,Abundance)

wilcox.test(plot.data.1$Baseline[plot.data.1$group=="PPI"],plot.data.1$Intervention[plot.data.1$group=="PPI"],paired = T)
wilcox.test(plot.data.1$Baseline[plot.data.1$group=="H2RA"],plot.data.1$Intervention[plot.data.1$group=="H2RA"],paired = T)

plot.data.1$increase <- plot.data.1$Intervention-plot.data.1$Baseline
median(plot.data.1$increase[plot.data.1$group=="PPI"])
median(plot.data.1$increase[plot.data.1$group=="H2RA"&plot.data.1$id!="H15"])


# 4. Supp Fig (abandon) -------------------------------------------------------------
abun.data <- filtered.abun[["gut"]]
rownames(abun.data) <- abun.data$Taxa
abun.data <- abun.data[,-1]
shannon <- as.data.frame(vegan::diversity(t(abun.data), index = "shannon", base = exp(1)))
colnames(shannon) <- "shannon"
simpson.gini <- as.data.frame(vegan::diversity(t(abun.data), index = "simpson"))
colnames(simpson.gini) <- "simpson"

shannon <- merge(shannon, sample.list, by.x = "row.names", by.y = "sample")
simpson <- merge(simpson.gini, sample.list, by.x = "row.names", by.y = "sample")

shannon$timepoint <- factor(shannon$timepoint,levels=c("Baseline","Intervention"))
p1 <- shannon %>% subset(group=="PPI") %>% ggplot(aes(x = timepoint, y = shannon, fill = timepoint)) +
  geom_boxplot(width = 0.3) +
  scale_fill_lancet(name = "Timepoint") +
  theme_bw() +
  labs(x = "Timepoint", y = "Shannon") +
  geom_signif(comparisons = list(c("Baseline","Intervention")), map_signif_level = T, step_increase = 0.2,test = "wilcox.test")
p1

p2 <- shannon %>% subset(group=="H2RA") %>% ggplot(aes(x = timepoint, y = shannon, fill = timepoint)) +
  geom_boxplot(width = 0.3) +
  scale_fill_lancet(name = "Timepoint") +
  theme_bw() +
  labs(x = "Timepoint", y = "Shannon") +
  geom_signif(comparisons = list(c("Baseline","Intervention")), map_signif_level = T, step_increase = 0.2,test = "wilcox.test")
p2

p1 <- simpson %>% subset(group=="PPI") %>% ggplot(aes(x = timepoint, y = simpson, fill = timepoint)) +
  geom_boxplot(width = 0.3) +
  scale_fill_lancet(name = "Timepoint") +
  theme_bw() +
  labs(x = "Timepoint", y = "Simpson") +
  geom_signif(comparisons = list(c("Baseline","Intervention")), map_signif_level = T, step_increase = 0.2,test = "wilcox.test")
p1

p2 <- simpson %>% subset(group=="H2RA") %>% ggplot(aes(x = timepoint, y = simpson, fill = timepoint)) +
  geom_boxplot(width = 0.3) +
  scale_fill_lancet(name = "Timepoint") +
  theme_bw() +
  labs(x = "Timepoint", y = "Simpson") +
  geom_signif(comparisons = list(c("Baseline","Intervention")), map_signif_level = T, step_increase = 0.2,test = "wilcox.test")
p2

# 5. D: Modeling -------------------------------------------------------------

siamcat.model <- function(abundance, level, meta, case, method, feature.num = 50, feature.diff = c(), filename) {
  abundance <- separate(abundance, "Taxa", into = c("k", "p", "c", "o", "f", "g", "s", "t"), sep = "\\|")
  abundance <- aggregate(abundance[, -c(1:8)], by = list(abundance[, level]), sum)
  colnames(abundance)[1] <- "Taxa"

  rownames(abundance) <- abundance$Taxa
  abundance <- abundance[, -1]
  abundance <- abundance / 100

  meta <- merge(meta, metadata[, c(1, 3, 4, 7, 9, 15)], by.x = "id", by.y = "ID")
  rownames(meta) <- meta$sample
  meta <- meta[, -c(1, 2)]
  abundance <- abundance[, colnames(abundance) %in% rownames(meta)]

  # create label
  label <- create.label(meta = meta, label = "timepoint", case = case)
  # create object
  obj <- siamcat(feat = abundance, label = label, meta = meta)

  # filter(按照丰度进行过滤)
  obj <- filter.features(obj, filter.method = "abundance", cutoff = 1e-04)

  # normalize
  obj <- normalize.features(obj, norm.method = "log.std", norm.param = list(log.n0 = 1e-06, sd.min.q = 0.1), verbose = 2)

  # cross validation(交叉验证)
  obj <- create.data.split(obj, num.folds = 5, num.resample = 3)

  # model training
  obj <- train.model(obj, method = method)

  # model prediction
  obj <- make.predictions(obj)

  # model evaluation
  obj <- evaluate.predictions(obj)

  # 选择重要的特征重新建模
  # feature_weight <- feature_weights(obj)
  # feature_weight <- arrange(feature_weight,desc(feature_weight$mean.weight))
  # feature_select <- head(rownames(feature_weight),feature.num)
  #
  # feat.select.abun <- abundance[rownames(abundance) %in% feature_select,]
  if (length(feature.diff) != 0) {
    diff.select.abun <- abundance[rownames(abundance) %in% feature.diff, ]
  }

  # # 重新建模
  # select <- siamcat(feat=abundance,label=label,meta=meta)
  # select <- filter.features(select, filter.method = "abundance", cutoff = 0.001)
  # select <- normalize.features(select, norm.method = "log.std", norm.param = list(log.n0 = 1e-06, sd.min.q = 0.1), verbose = 2)
  # select <-  create.data.split(select, num.folds = 10, num.resample = 5)
  # select <- train.model(select, method = method,perform.fs=T,
  #                       param.fs = list(no_features = feature.num, method = "AUC", direction="absolute"))
  # select <- make.predictions(select)
  # select <-  evaluate.predictions(select)



  if (length(feature.diff) != 0) {
    diff <- siamcat(feat = diff.select.abun, label = label, meta = meta)
    diff <- filter.features(diff, filter.method = "abundance", cutoff = 0.001)
    diff <- normalize.features(diff, norm.method = "log.std", norm.param = list(log.n0 = 1e-06, sd.min.q = 0.1), verbose = 2)
    diff <- create.data.split(diff, num.folds = 5, num.resample = 3)
    diff <- train.model(diff, method = method)
    diff <- make.predictions(diff)
    diff <- evaluate.predictions(diff)

    model.evaluation.plot(
      "All features" = obj,
      # 'Select features'=select,
      "Diff features" = diff,
      fn.plot = filename,
      show.all = F
    )

    return(list(obj, select, feature_select, diff))
  } else {
    model.evaluation.plot(
      "All features" = obj,
      # 'Select features'=select,
      fn.plot = filename,
      show.all = F
    )

    # return(list(obj, select, feature_select))
    return(list(obj, select))
  }
}

{ ## 4.1 建模 ----

  { ### PPI ----

    ppi.species.siamcat <- siamcat.model(
      abundance = filtered.abun[["gut"]],
      level = "s",
      meta = sample.list[sample.list$type == "gut" & sample.list$group == "PPI" & sample.list$timepoint %in% c("Baseline", "Intervention"), c("sample", "timepoint", "id")],
      case = "Intervention",
      method = "randomForest",
      filename = "01_Figure1/PPI.species.siamcat.pdf"
    )
    save(ppi.species.siamcat, file = "00_variables/ppi.species.siamcat.RData")

    model.evaluation.plot(
      "All features" = ppi.species.siamcat[[1]],
      "All features" = ppi.species.siamcat[[1]],
      fn.plot = "01_Figure1/PPI.species.siamcat.pdf",
      colours = c("#FED439", "#FED439"),
      show.all = F
    )
  }

  { ### H2RA ----

    h2ra.species.siamcat <- siamcat.model(
      abundance = filtered.abun[["gut"]],
      level = "s",
      meta = sample.list[sample.list$type == "gut" & sample.list$group == "H2RA" & sample.list$timepoint %in% c("Baseline", "Intervention"), c("sample", "timepoint", "id")],
      case = "Intervention",
      method = "randomForest",
      filename = "01_Figure1/H2RA.species.siamcat.pdf"
    )

    save(h2ra.species.siamcat, file = "00_variables/h2ra.species.siamcat.RData")

    model.evaluation.plot(
      "All features" = h2ra.species.siamcat[[1]],
      "All features" = h2ra.species.siamcat[[1]],
      fn.plot = "01_Figure1/h2ra.species.siamcat.pdf",
      colours = c("#709AE1", "#709AE1"),
      show.all = F
    )
  }
}

{ ## 4.2 重要特征画图 ----

  { ### PPI 特征重要性画图 ----

    ppi.feat.weight <- feature_weights(ppi.species.siamcat[[1]]) %>% as.data.frame()
    ppi.feat.weight <- ppi.feat.weight[ppi.feat.weight$percentage >= 0.5, ]
    ppi.feat.weight <- ppi.feat.weight %>% arrange(median.rel.weight)
    ppi.feat.weight$Taxa <- rownames(ppi.feat.weight)
    # ppi.feat.weight$Taxa <- factor(ppi.feat.weight$Taxa,levels = ppi.feat.weight$Taxa)
    ppi.feat.weight$Taxa <- str_sub(ppi.feat.weight$Taxa, 4, -1)
    ppi.feat.weight$Taxa <- gsub("_", " ", ppi.feat.weight$Taxa)

    ppi.feat.weight$diff <- apply(ppi.feat.weight, 1, function(x) {
      if (x[8] %in% species.base.intervent.diff[species.base.intervent.diff$group == "PPI", 1]) {
        return("Yes")
      } else {
        return("No")
      }
    })
    ppi.feat.weight$color[ppi.feat.weight$diff == "Yes"] <- "#E67764"
    ppi.feat.weight$color[ppi.feat.weight$diff == "No"] <- "black"
    ppi.feat.weight$color[ppi.feat.weight$Taxa == "Turicibacter sanguinis"] <- "#7FA493"
    ppi.feat.weight <- ppi.feat.weight %>% arrange(desc(median.rel.weight))
    ppi.feat.weight <- head(ppi.feat.weight, 20)
    ppi.feat.weight$Taxa <- factor(ppi.feat.weight$Taxa, levels = rev(ppi.feat.weight$Taxa))

    # p1：特征重要性排序
    p1 <- ppi.feat.weight %>% ggplot(aes(x = median.rel.weight, y = Taxa, fill = color)) +
      geom_col(width = 0.5) +
      theme_bw() +
      labs(x = "Median relative feat. weight", y = "") +
      scale_fill_manual(name = "Differential abundant in", values = c("#7FA493", "#E67764", "black"), labels = c("Baseline", "Intervention", "None")) +
      theme(
        axis.text.y = element_text(face = "italic", family = "serif"),
        text = element_text(family = "serif"),
        legend.position = "top"
      ) +
      xlim(0, 0.15)
    p1

    # p2: Abundance 箱线图
    plot.data <- filtered.abun[["gut"]]
    plot.data <- separate(plot.data, "Taxa", into = c("k", "p", "c", "o", "f", "g", "s", "t"), sep = "\\|")
    plot.data.species <- aggregate(plot.data[, -c(1:8)], by = list(plot.data[, "s"]), sum)
    plot.data.species$Group.1 <- str_sub(plot.data.species$Group.1, 4, -1)
    plot.data.species$Group.1 <- gsub("_", " ", plot.data.species$Group.1)
    plot.data.species <- plot.data.species[plot.data.species$Group.1 %in% ppi.feat.weight$Taxa, ]
    plot.data.species <- gather(plot.data.species, sample, abundance, -Group.1)
    plot.data.species <- merge(plot.data.species, sample.list[, c(1, 4, 5)], by = "sample")
    plot.data.species <- plot.data.species[(plot.data.species$group %in% c("PPI")) & (plot.data.species$timepoint %in% c("Baseline", "Intervention")), ]
    colnames(plot.data.species)[2] <- "Taxa"
    plot.data.species$Taxa <- factor(plot.data.species$Taxa, levels = levels(ppi.feat.weight$Taxa))

    p2 <- plot.data.species %>% ggplot(aes(x = abundance, y = Taxa, fill = timepoint)) +
      geom_boxplot(width = 0.8, lwd = 0.3) +
      scale_x_log10(labels = scales::percent_format(scale = 1), limits = c(0.0001, 10)) +
      theme_bw() +
      scale_fill_manual(name = "Timepoint", values = c("#7FA493", "#E67764"), labels = c("Baseline", "During1")) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family = "serif"),
        legend.position = "top"
      ) +
      labs(x = "Abundance (Log10-scaled)", y = "")

    p2

    # p3: fold change
    plot.data.species <- plot.data.species %>%
      group_by(Taxa, group, timepoint) %>%
      summarise(mean.abun = mean(abundance))
    plot.data.species <- spread(plot.data.species, timepoint, mean.abun)
    plot.data.species$fc <- plot.data.species$Intervention / plot.data.species$Baseline
    plot.data.species <- merge(plot.data.species, ppi.feat.weight[, c(8, 10)], by = "Taxa")
    plot.data.species$Taxa <- factor(plot.data.species$Taxa, levels = levels(ppi.feat.weight$Taxa))

    p3 <- plot.data.species %>% ggplot(aes(x = fc, y = Taxa, fill = color)) +
      geom_col(width = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.1), limits = c(0.05, 4500)) +
      scale_fill_manual(name = "Timepoint", values = c("#7FA493", "#E67764", "black"), labels = c("Baseline", "Intervention", "None")) +
      theme_bw() +
      # geom_vline(xintercept = 1, linetype = "dashed") +
      theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top",
        text = element_text(family = "serif")
      ) +
      labs(x = "Fold change (Log10-scaled)", y = "")
    p3

    supp.table.2.ppi <- merge(ppi.feat.weight, plot.data.species, by = "Taxa")
    supp.table.2.ppi <- supp.table.2.ppi[, c(11, 1, 6, 14)]
    colnames(supp.table.2.ppi) <- c("Group", "Taxa", "Meidan relative weight", "Fold change")
  }

  { ### h2ra 特征重要性画图 ----

    h2ra.feat.weight <- feature_weights(h2ra.species.siamcat[[1]]) %>% as.data.frame()
    h2ra.feat.weight <- h2ra.feat.weight[h2ra.feat.weight$percentage >= 0.5, ]
    h2ra.feat.weight <- h2ra.feat.weight %>% arrange(median.rel.weight)
    h2ra.feat.weight$Taxa <- rownames(h2ra.feat.weight)
    # h2ra.feat.weight$Taxa <- factor(h2ra.feat.weight$Taxa,levels = h2ra.feat.weight$Taxa)
    h2ra.feat.weight$Taxa <- str_sub(h2ra.feat.weight$Taxa, 4, -1)
    h2ra.feat.weight$Taxa <- gsub("_", " ", h2ra.feat.weight$Taxa)

    h2ra.feat.weight$diff <- apply(h2ra.feat.weight, 1, function(x) {
      if (x[8] %in% species.base.intervent.diff[species.base.intervent.diff$group == "H2RA", 1]) {
        return("Yes")
      } else {
        return("No")
      }
    })
    h2ra.feat.weight$color[h2ra.feat.weight$diff == "Yes"] <- "#E67764"
    h2ra.feat.weight$color[h2ra.feat.weight$diff == "No"] <- "black"
    h2ra.feat.weight <- h2ra.feat.weight %>% arrange(desc(median.rel.weight))
    h2ra.feat.weight <- head(h2ra.feat.weight, 20)
    h2ra.feat.weight$Taxa <- factor(h2ra.feat.weight$Taxa, levels = rev(h2ra.feat.weight$Taxa))

    # p4：特征重要性排序
    p4 <- h2ra.feat.weight %>% ggplot(aes(x = median.rel.weight, y = Taxa, fill = color)) +
      geom_col(width = 0.5) +
      theme_bw() +
      labs(x = "Median relative feat. weight", y = "") +
      scale_fill_manual(name = "Differential abundant in", values = c("#E67764", "black"), labels = c("Intervention", "None")) +
      theme(
        axis.text.y = element_text(face = "italic", family = "serif"),
        text = element_text(family = "serif"),
        legend.position = "top"
      ) +
      xlim(0, 0.15)
    p4

    # p2: Abundance 箱线图
    plot.data <- filtered.abun[["gut"]]
    plot.data <- separate(plot.data, "Taxa", into = c("k", "p", "c", "o", "f", "g", "s", "t"), sep = "\\|")
    plot.data.species <- aggregate(plot.data[, -c(1:8)], by = list(plot.data[, "s"]), sum)
    plot.data.species$Group.1 <- str_sub(plot.data.species$Group.1, 4, -1)
    plot.data.species$Group.1 <- gsub("_", " ", plot.data.species$Group.1)
    plot.data.species <- plot.data.species[plot.data.species$Group.1 %in% h2ra.feat.weight$Taxa, ]
    plot.data.species <- gather(plot.data.species, sample, abundance, -Group.1)
    plot.data.species <- merge(plot.data.species, sample.list[, c(1, 4, 5)], by = "sample")
    plot.data.species <- plot.data.species[(plot.data.species$group %in% c("H2RA")) & (plot.data.species$timepoint %in% c("Baseline", "Intervention")), ]
    colnames(plot.data.species)[2] <- "Taxa"
    plot.data.species$Taxa <- factor(plot.data.species$Taxa, levels = levels(h2ra.feat.weight$Taxa))

    p5 <- plot.data.species %>% ggplot(aes(x = abundance, y = Taxa, fill = timepoint)) +
      geom_boxplot(width = 0.8, lwd = 0.3) +
      scale_x_log10(labels = scales::percent_format(scale = 1), limits = c(0.0001, 10)) +
      theme_bw() +
      scale_fill_manual(name = "Timepoint", values = c("#7FA493", "#E67764"), labels = c("Baseline", "During1")) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family = "serif"),
        legend.position = "top"
      ) +
      labs(x = "Abundance (Log10-scaled)", y = "")

    p5

    # p3: fold change
    plot.data.species <- plot.data.species %>%
      group_by(Taxa, group, timepoint) %>%
      summarise(mean.abun = mean(abundance))
    plot.data.species <- spread(plot.data.species, timepoint, mean.abun)
    plot.data.species$fc <- plot.data.species$Intervention / plot.data.species$Baseline
    plot.data.species <- merge(plot.data.species, h2ra.feat.weight[, c(8, 10)], by = "Taxa")
    plot.data.species$Taxa <- factor(plot.data.species$Taxa, levels = levels(h2ra.feat.weight$Taxa))

    p6 <- plot.data.species %>% ggplot(aes(x = fc, y = Taxa, fill = color)) +
      geom_col(width = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.1), limits = c(0.05, 4500)) +
      scale_fill_manual(name = "Timepoint", values = c("#E67764", "black"), labels = c("Intervention", "None")) +
      theme_bw() +
      # geom_vline(xintercept = 1, linetype = "dashed") +
      theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top",
        text = element_text(family = "serif")
      ) +
      labs(x = "Fold change (Log10-scaled)", y = "")
    p6

    supp.table.2.h2ra <- merge(h2ra.feat.weight, plot.data.species, by = "Taxa")
    supp.table.2.h2ra <- supp.table.2.h2ra[, c(11, 1, 6, 14)]
    colnames(supp.table.2.h2ra) <- c("Group", "Taxa", "Meidan relative weight", "Fold change")
  }
  p7 <- (p1 + p2 + p3) / (p4 + p5 + p6)
  p7
  ggsave(p7, filename = "01_Figure1/Figure1D.pdf", width = 12, height = 16)

  supp.table.2 <- rbind(supp.table.2.ppi, supp.table.2.h2ra)
  write_xlsx(supp.table.2, "00_supp_tables/SuppTable-Feature_weight.xlsx")
}

