# 0. load packages --------------------------------------------------------
library(tidyverse)
library(readxl)
library(writexl)
library(ggsci)
library(RColorBrewer)
library(ggpubr)
library(patchwork)

# 0. lazy load ------------------------------------------------------------
load("00_variables/filtered.abun.RData")
load("00_variables/sample.list.RData")
load("00_variables/same.indi.time.nGD.RData")
load("00_variables/raw.abun.RData")
load("00_variables/total.nGD.RData")
group.color <- c("#709AE1", "#FED439")
names(group.color) <- c("H2RA", "PPI")


timepoint.color <- c("#7FA493", "#E67764")
names(timepoint.color) <- c("Baseline", "Intervention")


# 1. prepare directory and load data --------------------------------------
dir.create("03_Figure3/", showWarnings = F)
# dir.create("03_Figure3/01_input", showWarnings = F)


# 2. A: oral & gut prevalence ordered by abundance ------------------------
calculate.prevalence <- function(abun, sample.list, group, timepoint) {
  gut.abun <- abun[["gut"]]
  oral.abun <- abun[["oral"]]
  gut.sample <- sample.list[sample.list$type == "gut" & (sample.list$group %in% group) & (sample.list$timepoint %in% timepoint), ]
  oral.sample <- sample.list[sample.list$type == "oral" & (sample.list$group %in% group) & (sample.list$timepoint %in% timepoint), ]

  sample <- merge(gut.sample, oral.sample, by = c("id", "timepoint"))
  sample <- sample[, c(1, 2, 5, 3, 7)]
  colnames(sample) <- c("id", "timepoint", "group", "gut", "oral")

  gut.abun <- gut.abun[, c("Taxa", sample$gut)]
  oral.abun <- oral.abun[, c("Taxa", sample$oral)]

  # prevalence 计算
  gut.abun$prev <- apply(gut.abun[, -1], 1, function(x) {
    sum(x != 0) / (length(x)) * 100
  })
  gut.abun$type <- "gut"
  gut.abun <- gut.abun[, c(1, ncol(gut.abun), (ncol(gut.abun) - 1), 2:(ncol(gut.abun) - 2))]
  oral.abun$prev <- apply(oral.abun[, -1], 1, function(x) {
    sum(x != 0) / (length(x)) * 100
  })
  oral.abun$type <- "oral"
  oral.abun <- oral.abun[, c(1, ncol(oral.abun), (ncol(oral.abun) - 1), 2:(ncol(oral.abun) - 2))]

  total.prev <- rbind(gut.abun[, 1:3], oral.abun[, 1:3])
  total.prev <- spread(total.prev, type, prev, fill = 0)

  total.prev$Type <- ""
  for (i in 1:nrow(total.prev)) {
    a <- total.prev[i, 2]
    b <- total.prev[i, 3]

    if (a >= 10 & b >= 10) {
      total.prev[i, 4] <- "Prevalent in both sites"
    } else if (a >= 10 & b < 10) {
      total.prev[i, 4] <- "Prevalent in gut"
    } else if (a < 10 & b >= 10) {
      total.prev[i, 4] <- "Prevalent in oral cavity"
    } else {
      total.prev[i, 4] <- "Not prevalent in either sites"
    }
  }

  total.prev <- total.prev[total.prev$Type != "Not prevalent in either sites", ]
  total.prev$Type <- factor(total.prev$Type, c("Prevalent in oral cavity", "Prevalent in both sites", "Prevalent in gut"))
  total.prev <- total.prev %>% arrange(Type)

  # 按照丰度平均值从大到小排序
  gut.abun$mean.abun <- rowMeans(gut.abun[, -c(1:3)])
  gut.mean.abun <- gut.abun[, c(1, 2, ncol(gut.abun))]
  oral.abun$mean.abun <- rowMeans(oral.abun[, -c(1:3)])
  oral.mean.abun <- oral.abun[, c(1, 2, ncol(oral.abun))]

  mean.abun.1 <- oral.mean.abun[oral.mean.abun$Taxa %in% as.character(total.prev$Taxa[total.prev$Type == "Prevalent in oral cavity"]), ]
  mean.abun.2 <- gut.mean.abun[gut.mean.abun$Taxa %in% as.character(total.prev$Taxa[total.prev$Type == "Prevalent in gut"]), ]
  mean.abun.3 <- gut.mean.abun[gut.mean.abun$Taxa %in% as.character(total.prev$Taxa[total.prev$Type == "Prevalent in both sites"]), ]
  mean.abun.1 <- mean.abun.1 %>% arrange(desc(mean.abun))
  mean.abun.2 <- mean.abun.2 %>% arrange(desc(mean.abun))
  mean.abun.3 <- mean.abun.3 %>% arrange(desc(mean.abun))

  baseline.order <- c(mean.abun.1$Taxa, mean.abun.3$Taxa, mean.abun.2$Taxa)
  total.prev$Taxa <- factor(total.prev$Taxa, levels = baseline.order)
  # save(total.prev,file="00_variables/oral.gut.prevalence.RData")

  total.prev.long <- gather(total.prev[, -4], site, prevalence, -Taxa)
  total.prev.long$prevalence <- as.numeric(total.prev.long$prevalence)
  total.prev.long$Taxa <- factor(total.prev.long$Taxa, levels = levels(total.prev$Taxa))
  total.prev.long$prevalence[total.prev.long$site == "gut"] <- total.prev.long$prevalence[total.prev.long$site == "gut"] * (-1)

  color <- c("#CCABDB", "#6E6F78", "#FAAC69")
  names(color) <- c("Prevalent in oral cavity", "Prevalent in both sites", "Prevalent in gut")

  label.x <- mean.abun.1[floor(nrow(mean.abun.1) / 2), 1]
  label.y <- mean.abun.3[floor(nrow(mean.abun.3) / 2), 1]
  label.z <- mean.abun.2[floor(nrow(mean.abun.2) / 2), 1]

  # p1: oral & gut prevalence
  p1 <- total.prev.long %>% ggplot(aes(x = Taxa, y = site, fill = prevalence)) +
    geom_tile() +
    geom_linerange(aes(xmin = mean.abun.1[1, 1], xmax = mean.abun.1[nrow(mean.abun.1), 1], y = 3), colour = color["Prevalent in oral cavity"], linewidth = 1.5, linetype = "solid") +
    geom_linerange(aes(xmin = mean.abun.3[1, 1], xmax = mean.abun.3[nrow(mean.abun.3), 1], y = 3), colour = color["Prevalent in both sites"], linewidth = 1.5, linetype = "solid") +
    geom_linerange(aes(xmin = mean.abun.2[1, 1], xmax = mean.abun.2[nrow(mean.abun.2), 1], y = 3), colour = color["Prevalent in gut"], linewidth = 1.5, linetype = "solid") +
    annotate("text", x = label.x, y = 3.3, label = paste0("Prevalent in oral cavity (", nrow(mean.abun.1), ")")) +
    annotate("text", x = label.y, y = 3.3, label = paste0("Prevalent in both sites (", nrow(mean.abun.3), ")")) +
    annotate("text", x = label.z, y = 3.3, label = paste0("Prevalent in gut (", nrow(mean.abun.2), ")")) +
    scale_fill_gradient2(low = color["Prevalent in gut"], mid = "white", high = color["Prevalent in oral cavity"], midpoint = 0, space = "Lab", limits = c(-100, 100)) +
    scale_y_discrete(labels = c("Gut prevalence", "Oral prevalence")) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      text = element_text(family = "serif")
    )
  p1

  # p2: scale bar
  oral.prev.scale <- ggplot(data.frame(x = seq(0, 1, by = 0.01)), aes(x = x, y = 1, fill = x)) +
    geom_tile() +
    scale_fill_continuous(low = "white", high = "#CCABDB") +
    # coord_flip() +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    scale_x_continuous(labels = scales::percent_format(scale = 100))
  oral.prev.scale

  gut.prev.scale <- ggplot(data.frame(x = seq(0, 1, by = 0.01)), aes(x = x, y = 1, fill = x)) +
    geom_tile() +
    scale_fill_continuous(low = "white", high = "#FAAC69") +
    # coord_flip() +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    scale_x_continuous(labels = scales::percent_format(scale = 100))
  gut.prev.scale

  return(list(plot = p1, oral.scale = oral.prev.scale, gut.scale = gut.prev.scale, order = baseline.order,prev.data=total.prev))
}

a.prevalence.plot <- calculate.prevalence(
  abun = filtered.abun,
  sample.list = sample.list,
  group = c("H2RA", "PPI"),
  timepoint = c("Baseline")
)
a.prevalence.plot[[1]]
a.prevalence.plot[[2]]
a.prevalence.plot[[3]]

ggsave(a.prevalence.plot[[2]], filename = "03_Figure3/Figure3A_oral.prev.scalebar.pdf", width = 1.5, height = 1)
ggsave(a.prevalence.plot[[3]], filename = "03_Figure3/Figure3A_gut.prev.scalebar.pdf", width = 1.5, height = 1)

supp.table <- a.prevalence.plot[[5]]
colnames(supp.table) <- c("Taxa","Prevalence in gut","Prevalence in oral cavity","Type")
write_xlsx(supp.table,"SuppTable11.xlsx")

# 3. A: oral-to-gut transmission (left) ------------------------------------------
# same.indi.time <- total.nGD[total.nGD$same_individual == "same_individual" & total.nGD$same_timepoint == "same_timepoint" & total.nGD$same_type == "different_type", ]
same.indi.time <- same.indi.time[same.indi.time$nGD <= 0.01, ]
# same.indi.time$group.1[same.indi.time$group.1=="H2B"] <- "H2RA"
# same.indi.time$group.2[same.indi.time$group.2=="H2B"] <- "H2RA"
# same.indi.time$timepoint.1[same.indi.time$timepoint.1=="First"] <- "Intervention"
# same.indi.time$timepoint.2[same.indi.time$timepoint.2=="First"] <- "Intervention"
# same.indi.time$timepoint.2[same.indi.time$timepoint.2=="Last"] <- "Follow_up"
# same.indi.time$timepoint.1[same.indi.time$timepoint.1=="Last"] <- "Follow_up"
same.indi.time <- same.indi.time[same.indi.time$timepoint.1 %in% c("Baseline","Intervention"),]

a.trans.data <- same.indi.time[, c(1, 4, 6, 7)]
a.trans.data <- merge(a.trans.data, filtered.abun[["oral.species"]][, c("Taxa", "t")], by.x = "SGB", by.y = "t")
a.trans.data$Taxa <- factor(a.trans.data$Taxa, levels = a.prevalence.plot[["order"]])

a.trans.data.1 <- a.trans.data[, c(3, 4, 5)] %>% unique()
a.trans.data.1$exist <- 1
a.trans.data.1 <- spread(a.trans.data.1, timepoint.1, exist, fill = 0)
colnames(a.trans.data.1)[1] <- "Group"

add.data <- data.frame()
for (i in unique(a.trans.data.1$Group)) {
  species <- a.trans.data.1$Taxa[a.trans.data.1$Group == i] %>% as.character()
  not.exist <- a.prevalence.plot[["order"]][!(a.prevalence.plot[["order"]] %in% species)]

  temp <- data.frame(
    Group = rep(i, length(not.exist)),
    Taxa = not.exist,
    Baseline = rep(0, length(not.exist)),
    Intervention = rep(0, length(not.exist))
  )
  add.data <- rbind(add.data, temp)
}

a.trans.data.1 <- rbind(a.trans.data.1, add.data)
a.trans.data.1 <- gather(a.trans.data.1, timepoint, exist, -c(Group, Taxa))
a.trans.data.1$exist <- as.character(a.trans.data.1$exist)
a.trans.data.1$timepoint <- factor(a.trans.data.1$timepoint, levels = rev(c("Baseline", "Intervention")))
a.trans.data.1$Taxa <- factor(a.trans.data.1$Taxa, levels = a.prevalence.plot[["order"]])
a.trans.data.1$Group <- factor(a.trans.data.1$Group, levels = c("PPI", "H2RA"))
a.trans.data.1$color.group <- paste0(a.trans.data.1$timepoint, "_", a.trans.data.1$exist)

color <- c("white", "#7FA493", "white", "#E67764")
names(color) <- c("Baseline_0", "Baseline_1", "Intervention_0", "Intervention_1")

a.transmission.plot <- a.trans.data.1 %>% ggplot(aes(x = Taxa, y = timepoint, fill = color.group)) +
  geom_tile() +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = color) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    strip.background = element_rect(color = "black", fill = "#FFFFFF"),
    strip.text.x = element_blank(),
    text = element_text(family = "serif")
  )
a.transmission.plot

a.left <- a.prevalence.plot[["plot"]] / a.transmission.plot + plot_layout(heights = c(0.5, 1))
a.left
ggsave(a.left, filename = "03_Figure3/Figure3A_left.pdf", width = 10, height = 3)

# supp table 11
supp.table11 <- merge(a.trans.data.1,a.prevalence.plot[["prev.data"]][,c(1,4)],by="Taxa")
supp.table11 <- supp.table11[supp.table11$exist==1,]
supp.table11 <- supp.table11[,-c(4,5)]
supp.table11 <- separate(supp.table11,col=Taxa,into=c('k','p','c','o','f','g','s','t'),sep="\\|")
supp.table11 <- supp.table11[,c(7,9,10,11)]
supp.table11$s <- gsub("s__","",supp.table11$s)
supp.table11$s <- gsub("_"," ",supp.table11$s)
colnames(supp.table11) <- c("Taxa","Group","Timepoint","Prevalent in body sites")
write_xlsx(supp.table11,"SuppTable11.xlsx")

# 4. A: oral-to-gut transmission (mid) ------------------------------------------
stats <- same.indi.time %>%
  group_by(subjectID.1, timepoint.1, group.1) %>%
  summarise(Freq = n())

stats <- spread(stats,timepoint.1,Freq)
stats[is.na(stats)] <- 0
stats <- gather(stats,timepoint.1,Freq,-c(subjectID.1,group.1))
stats$timepoint.1 <- factor(stats$timepoint.1, levels = c("Baseline", "Intervention"))
stats$timepoint.1 <- factor(stats$timepoint.1, levels = rev(c("Baseline", "Intervention")))
stats$group.1 <- factor(stats$group.1, levels = c("PPI", "H2RA"))


timepoint.color <- c("#7FA493", "#E67764")
names(timepoint.color) <- c("Baseline", "Intervention")

a.mid <- stats %>% ggplot(aes(y = timepoint.1, x = Freq, fill = timepoint.1)) +
  geom_boxplot(width = 0.3) +
  theme_bw() +
  facet_grid(group.1 ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(name = "Timepoint", values = timepoint.color) +
  labs(x = "Number of species with oral-to-gut transmission", y = "", title = "") +
  geom_signif(comparison = list(c("Baseline", "Intervention")), step_increase = 0.1, map_signif_level = T, tip_length = 0.02) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(color = "black", fill = "#FFFFFF"),
    strip.text.x = element_blank(),
    text = element_text(family = "serif")
  )

a.mid
ggsave(a.mid, filename = "03_Figure3/Figure3A_mid.pdf", width = 5, height = 4)

stats$timepoint.1 <- factor(stats$timepoint.1,levels=c("Baseline","Intervention"))
supp3.a <- stats %>% ggplot(aes(x = group.1, y = Freq, fill = group.1)) +
  geom_boxplot(width = 0.3) +
  theme_bw() +
  facet_grid(.~timepoint.1, scales = "free_x", space = "free_x") +
  scale_fill_manual(name = "Group", values = group.color) +
  labs(y = "Number of species with oral-to-gut transmission", x = "", title = "") +
  geom_signif(comparison = list(c("PPI", "H2RA")), step_increase = 0.1, map_signif_level = T, tip_length = 0.02) +
  theme(
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(color = "black", fill = "#FFFFFF"),
    strip.text.y = element_blank(),
    text = element_text(family = "serif")
  )
supp3.a
ggsave(supp3.a, filename = "03_Figure3/SuppFig3.pdf", width = 6, height = 6)
wilcox.test(stats.1$Baseline[stats.1$group.1=="PPI"],
            stats.1$Baseline[stats.1$group.1=="H2RA"])
wilcox.test(stats.1$Intervention[stats.1$group.1=="PPI"],
            stats.1$Intervention[stats.1$group.1=="H2RA"])


stats.1 <- spread(stats,timepoint.1,Freq)
stats.1.ppi <- stats.1[stats.1$group.1=="PPI",]
wilcox.test(stats.1.ppi$Baseline,
            stats.1.ppi$Intervention,paired = T)
wilcox.test(stats.1$Baseline[stats.1$group.1=="H2RA"],
            stats.1$Intervention[stats.1$group.1=="H2RA"],paired = T)


# 5. A: oral-to-gut transmission (right) ------------------------------------------
abun.species.with.trans <- same.indi.time[, c(1, 2, 5, 6, 7)]
abun.species.with.trans <- merge(abun.species.with.trans, filtered.abun[["oral.species"]][, c("Taxa", "t")], by.x = "SGB", by.y = "t")
colnames(abun.species.with.trans) <- c("SGB", "sample", "subject", "group", "timepoint", "Taxa")
gut.abun <- raw.abun[["gut"]]
gut.abun.long <- gather(gut.abun, sample, abundance, -Taxa)
abun.species.with.trans <- merge(abun.species.with.trans, gut.abun.long, by = c("sample", "Taxa"))

abun.species.with.trans <- aggregate(abun.species.with.trans[, -c(1:6)], by = list(abun.species.with.trans$subject, abun.species.with.trans$group, abun.species.with.trans$timepoint), sum)
colnames(abun.species.with.trans) <- c("subject", "group", "timepoint", "total.abundance")

abun.species.with.trans <- spread(abun.species.with.trans,timepoint,total.abundance)
abun.species.with.trans[is.na(abun.species.with.trans)] <- 0
abun.species.with.trans <- gather(abun.species.with.trans,timepoint,total.abundance,-c(subject,group))

abun.species.with.trans$timepoint <- factor(abun.species.with.trans$timepoint, levels = rev(c("Baseline", "Intervention")))
abun.species.with.trans$group <- factor(abun.species.with.trans$group, levels = c("PPI", "H2RA"))

timepoint.color <- c("#7FA493", "#E67764")
names(timepoint.color) <- c("Baseline", "Intervention")

a.right <- abun.species.with.trans %>% ggplot(aes(x = total.abundance, y = timepoint, fill = timepoint)) +
  geom_boxplot(width = 0.3) +
  facet_grid(group ~ .) +
  theme_bw() +
  scale_fill_manual(name = "Timepoint", values = timepoint.color) +
  geom_signif(comparison = list(c("Baseline", "Intervention")), step_increase = 0.1, map_signif_level = T, tip_length = 0.02) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(color = "black", fill = "#FFFFFF"),
    strip.text.x = element_blank(),
    text = element_text(family = "serif")
  ) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Total abundance of the species with oral-to-gut transmission in the gut", y = "", title = "")

a.right
ggsave(a.right, filename = "03_Figure3/Figure3A_right.pdf", width = 5, height = 4)

abun.species.with.trans$timepoint <- factor(abun.species.with.trans$timepoint,levels=c("Baseline","Intervention"))
supp3.a <- abun.species.with.trans %>% ggplot(aes(x = group, y = total.abundance, fill = group)) +
  geom_boxplot(width = 0.3) +
  theme_bw() +
  facet_grid(.~timepoint, scales = "free_x", space = "free_x") +
  scale_fill_manual(name = "Group", values = group.color) +
  labs(y = "Total abundance of the species with oral-to-gut transmission in the gut", x = "", title = "") +
  geom_signif(comparison = list(c("PPI", "H2RA")), step_increase = 0.1, map_signif_level = T, tip_length = 0.02) +
  theme(
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(color = "black", fill = "#FFFFFF"),
    strip.text.y = element_blank(),
    text = element_text(family = "serif")
  )+
  scale_y_continuous(labels = scales::percent_format(scale = 1))
supp3.a
ggsave(supp3.a, filename = "03_Figure3/SuppFig3.2.pdf", width = 6, height = 6)

abun.species.with.trans.1 <- abun.species.with.trans
abun.species.with.trans.1 <- spread(abun.species.with.trans.1,timepoint,total.abundance)
wilcox.test(abun.species.with.trans.1$Baseline[abun.species.with.trans.1$group=="PPI"],
            abun.species.with.trans.1$Intervention[abun.species.with.trans.1$group=="PPI"],paired = T)
wilcox.test(abun.species.with.trans.1$Baseline[abun.species.with.trans.1$group=="H2RA"],
            abun.species.with.trans.1$Intervention[abun.species.with.trans.1$group=="H2RA"],paired = T)

wilcox.test(abun.species.with.trans.1$Baseline[abun.species.with.trans.1$group=="PPI"],
            abun.species.with.trans.1$Baseline[abun.species.with.trans.1$group=="H2RA"])
wilcox.test(abun.species.with.trans.1$Intervention[abun.species.with.trans.1$group=="PPI"],
            abun.species.with.trans.1$Intervention[abun.species.with.trans.1$group=="H2RA"])

a.bottom.right <- a.mid + a.right
ggsave(a.bottom.right, filename = "03_Figure3/Figure3A_bottom.right.pdf", width = 10, height = 4)


# supptable
supp.table12 <- merge(stats,abun.species.with.trans,by.x="subjectID.1",by.y="subject")
supp.table12 <- supp.table12[supp.table12$timepoint.1==supp.table12$timepoint,]
supp.table12 <- supp.table12[,c(1,2,3,4,7)]
colnames(supp.table12) <- c("SubjectID","Group","Timepoint","Number of transmissible species","Total abundance of transmissible species")
write_xlsx(supp.table12,"SuppTable11.xlsx")

# 6. B: specific species oral-to-gut transmission ----------------------------
{ # 6.1 统计各组中各个细菌的transmission prevalence ----

  b.stats.1 <- same.indi.time %>%
    group_by(SGB, timepoint.1, group.1) %>%
    summarise(Freq = n())

  calculate.transmission.prev <- function(group, timepoint, background) {
    stats1.baseline <- b.stats.1 %>%
      subset((timepoint.1 == timepoint) & (group.1 == group)) %>%
      group_by(SGB) %>%
      summarise(Sum = sum(Freq))
    stats1.baseline <- merge(stats1.baseline, filtered.abun[["oral.species"]][, c("Taxa", "s", "t")], by.x = "SGB", by.y = "t")
    stats1.baseline$Ratio <- stats1.baseline$Sum / background * 100
    stats1.baseline <- stats1.baseline %>% arrange(Ratio)
    stats1.baseline$s <- str_sub(stats1.baseline$s, 4, -1)
    stats1.baseline$s <- gsub("_", " ", stats1.baseline$s)
    stats1.baseline$s <- factor(stats1.baseline$s, levels = unique(stats1.baseline$s))
    stats1.baseline$timepoint.1 <- timepoint
    stats1.baseline$group.1 <- group
    return(stats1.baseline)
  }

  h2ra.baseline <- calculate.transmission.prev(group = "H2RA", timepoint = "Baseline", background = 26)
  h2ra.intervention <- calculate.transmission.prev(group = "H2RA", timepoint = "Intervention", background = 26)

  ppi.baseline <- calculate.transmission.prev(group = "PPI", timepoint = "Baseline", background = 23)
  ppi.intervention <- calculate.transmission.prev(group = "PPI", timepoint = "Intervention", background = 23)
  
  # control.baseline <- calculate.transmission.prev(group = "Control", timepoint = "Baseline", background = 27)
  # control.intervention <- calculate.transmission.prev(group = "Control", timepoint = "Intervention", background = 27)
  

  b.trans.prev <- rbind(
    h2ra.baseline, h2ra.intervention,
    ppi.baseline,  ppi.intervention
    # control.baseline,control.intervention
  )

  b.trans.prev$SGB <- gsub("t__", "", b.trans.prev$SGB)
  b.trans.prev$s <- paste0(b.trans.prev$s, " [", b.trans.prev$SGB, "]")
  colnames(b.trans.prev) <- c("SGB", "Freq", "Taxa", "Species", "Ratio", "Timepoint", "Group")
}

{ # 6.2 计算Intervention和Baseline在PPI和H2RA两组中的差异 ----

  b.intervention.base <- b.trans.prev[, -c(1:3)]
  b.intervention.base <- spread(b.intervention.base, Timepoint, Ratio, fill = 0)
  b.intervention.base$minus <- b.intervention.base$Intervention - b.intervention.base$Baseline
  b.intervention.base <- b.intervention.base[, c(1, 2, 5)]
  b.intervention.base <- spread(b.intervention.base, Group, minus, fill = 0)

  # 补充b.trans.prev中存在但不存在b.intervention.base中的细菌
  not.exist <- setdiff(b.trans.prev$Species, b.intervention.base$Species)
  add.data <- data.frame(
    Species = not.exist,
    H2RA = rep(0, length(not.exist)),
    PPI = rep(0, length(not.exist))
  )
  b.intervention.base <- rbind(b.intervention.base, add.data)
  b.intervention.base$Signal <- ifelse(b.intervention.base$H2RA > b.intervention.base$PPI, "PPI < H2RA", ifelse(b.intervention.base$H2RA < b.intervention.base$PPI, "PPI > H2RA", "PPI = H2RA"))
  b.intervention.base.plot.y.order <- b.intervention.base %>% arrange(Signal, PPI, H2RA)
  b.intervention.base <- b.intervention.base[, -ncol(b.intervention.base)]
  # save(b.intervention.base.plot.y.order,file = "00_variables/b.intervention.base.plot.y.order.RData")
  b.intervention.base <- gather(b.intervention.base, Group, Minus, -Species)
  b.intervention.base$Facet <- "Intervention - Baseline"

  b.intervention.base$Species <- factor(b.intervention.base$Species, levels = b.intervention.base.plot.y.order$Species)
  b.intervention.base$Minus <- as.numeric(b.intervention.base$Minus)
  
  text.color <- unique(b.intervention.base[,1,drop=F])

  prev.data <- a.prevalence.plot[["prev.data"]]
  prev.data <- separate(prev.data,Taxa,into=c("k","p","c","o","f","g","s","t"),sep="\\|")
  prev.data <- prev.data[,c(7,8,11)]
  prev.data$s <- gsub("s__","",prev.data$s)
  prev.data$s <- gsub("_"," ",prev.data$s)
  prev.data$t <- gsub("t__","",prev.data$t)
  prev.data$Species <- paste0(prev.data$s," [",prev.data$t,"]")
  prev.data <- unique(prev.data[,c(4,3)])
  text.color <- merge(text.color,prev.data,by="Species")
  
  color <- data.frame(Type=c("Prevalent in oral cavity", "Prevalent in both sites", "Prevalent in gut"),Color=c("#CCABDB", "#6E6F78", "#FAAC69"))
  text.color <- merge(text.color,color,by="Type")
  rownames(text.color) <- text.color$Species
  text.color <- text.color[levels(b.intervention.base$Species),]
  # text.color <- text.color[,-1]

  b.intervention.base.plot <- b.intervention.base %>% ggplot(aes(x = Minus, y = Species)) +
    geom_line() +
    geom_point(aes(color = Group, size = Group)) +
    scale_size_manual(values = c(4.5, 4.5)) +
    theme_bw() +
    facet_wrap(. ~ Facet) +
    scale_y_discrete(position = "right") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "Transmission prevalence change", y = "", title = "") +
    scale_color_manual(name = "Group", values = group.color) +
    scale_x_continuous(labels = scales::percent_format(scale = 1)) +
    theme(
      plot.title = element_text(hjust = 0.5),
      # axis.text.y = element_text(colour = total.y.order$color,face="bold")
      # axis.text.y = element_blank(),
      strip.background = element_rect(color = "black", fill = "#FFFFFF"),
      # axis.ticks.y = element_blank(),
      axis.text.y = element_text(face = "italic",colour = text.color$Color),
      text = element_text(family = "serif"),
      # legend.position = "none"
    )
  b.intervention.base.plot
}

{ # 6.3 可视化各组中各个细菌的transmission prevalence ----

  b.trans.prev$Species <- factor(b.trans.prev$Species, levels = b.intervention.base.plot.y.order$Species)
  b.trans.prev$Timepoint <- factor(b.trans.prev$Timepoint, levels = c("Baseline", "Intervention"))
  b.trans.prev$Group <- factor(b.trans.prev$Group, levels = c("PPI", "H2RA"))
  save(b.trans.prev,file="00_variables/transmission.species.RData")
  
  b.trans.prev.plot <- b.trans.prev %>% ggplot(aes(x = Ratio, y = Species, color = Timepoint, shape = Timepoint)) +
    geom_point(size = 3.5) +
    theme_bw() +
    facet_wrap(. ~ Group) +
    scale_x_continuous(labels = scales::percent_format(scale = 1)) +
    scale_color_manual(values = timepoint.color, name = "Timepoint") +
    scale_shape_manual(name = "Timepoint", values = c(17, 15, 16, 18)) +
    labs(y = "", title = "", x = "Transmission prevalence") +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.background = element_rect(color = "black", fill = "#FFFFFF"),
      # strip.text.x = element_blank(),
      # axis.text.y = element_text(face="italic"),
      text = element_text(family = "serif"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      # legend.position = "none"
    )
  b.trans.prev.plot

  #supp table 13
  supp.table13 <- b.trans.prev[,-c(2,3)]
  colnames(supp.table13) <- c("SGB","Species","Transmission prevalence","Timepoint","Group")
  write_xlsx(supp.table13,"SuppTable13.xlsx")
  
}

b.plot <- ( b.trans.prev.plot | b.intervention.base.plot ) + plot_layout(widths = c(1.2, 0.3))
b.plot
ggsave(b.plot, filename = "03_Figure3/Figure3B.pdf", width = 16, height = 8)
ggsave(b.plot, filename = "03_Figure3/SuppFig4.pdf", width = 16, height = 8)


# 7. C: species characteristics & transmission prev change ----------------
{ # 7.1 口腔中的丰度与transmission prevalence 之间的相关性 ----
  
  oral.abun <- filtered.abun[["oral"]]
  oral.abun <- gather(oral.abun, sample, abundance, -Taxa)
  oral.abun <- merge(oral.abun, sample.list[, -2], by = "sample")
  oral.abun <- merge(oral.abun, b.trans.prev[, c(3, 4)], by = "Taxa")
  
  ppi.oral.abun.base.intervention <- oral.abun[oral.abun$group == "PPI" & oral.abun$timepoint %in% c("Baseline", "During1"), ]
  ppi.oral.abun.base.intervention <- ppi.oral.abun.base.intervention %>%
    group_by(Species, id) %>%
    summarise(mean = mean(abundance))
  
  ppi.oral.abun.base.intervention <- spread(ppi.oral.abun.base.intervention, id, mean) %>% as.data.frame()
  rownames(ppi.oral.abun.base.intervention) <- ppi.oral.abun.base.intervention$Species
  ppi.oral.abun.base.intervention <- ppi.oral.abun.base.intervention[, -1]
  ppi.oral.abun.base.intervention$mean.abun <- rowMeans(ppi.oral.abun.base.intervention)
  ppi.oral.abun.base.intervention <- ppi.oral.abun.base.intervention[, ncol(ppi.oral.abun.base.intervention), drop = F]
  
  ppi.base.intervention.prev <- b.intervention.base[b.intervention.base$Group == "PPI", ]
  rownames(ppi.base.intervention.prev) <- ppi.base.intervention.prev$Species
  ppi.base.intervention.prev <- ppi.base.intervention.prev[, -c(1, 2, 4), drop = F]
  ppi.oral.abun.base.intervention <- ppi.oral.abun.base.intervention[rownames(ppi.base.intervention.prev), , drop = F]
  
  ppi.abun.prev.cor <- cor.test(ppi.oral.abun.base.intervention$mean.abun, ppi.base.intervention.prev$Minus, method = "spearman")
}

{ # 7.2 口腔中的丰度与transmission prevalence 之间的相关性 ----
  
  oral.abun <- filtered.abun[["oral"]]
  oral.abun <- gather(oral.abun, sample, abundance, -Taxa)
  oral.abun <- merge(oral.abun, sample.list[, -2], by = "sample")
  oral.abun <- merge(oral.abun, b.trans.prev[, c(3, 4)], by = "Taxa")
  
  h2ra.oral.abun.base.intervention <- oral.abun[oral.abun$group == "H2RA" & oral.abun$timepoint %in% c("Baseline", "During1"), ]
  h2ra.oral.abun.base.intervention <- h2ra.oral.abun.base.intervention %>%
    group_by(Species, id) %>%
    summarise(mean = mean(abundance))
  
  h2ra.oral.abun.base.intervention <- spread(h2ra.oral.abun.base.intervention, id, mean) %>% as.data.frame()
  rownames(h2ra.oral.abun.base.intervention) <- h2ra.oral.abun.base.intervention$Species
  h2ra.oral.abun.base.intervention <- h2ra.oral.abun.base.intervention[, -1]
  h2ra.oral.abun.base.intervention$mean.abun <- rowMeans(h2ra.oral.abun.base.intervention)
  h2ra.oral.abun.base.intervention <- h2ra.oral.abun.base.intervention[, ncol(h2ra.oral.abun.base.intervention), drop = F]
  
  h2ra.base.intervention.prev <- b.intervention.base[b.intervention.base$Group == "H2RA", ]
  rownames(h2ra.base.intervention.prev) <- h2ra.base.intervention.prev$Species
  h2ra.base.intervention.prev <- h2ra.base.intervention.prev[, -c(1, 2, 4), drop = F]
  h2ra.oral.abun.base.intervention <- h2ra.oral.abun.base.intervention[rownames(h2ra.base.intervention.prev), , drop = F]
  
  h2ra.abun.prev.cor <- cor.test(h2ra.oral.abun.base.intervention$mean.abun, h2ra.base.intervention.prev$Minus, method = "spearman")
}

{ # 7.3 细菌的其他指标与transmission prevalence change的相关性(PPI) ----
  
  traitar0 <- as.data.frame(read.table("00_variables/predictions_majority-vote.txt", sep = "\t", header = T))
  traitar1 <- traitar0[grepl("^[0-9]", traitar0$X), ]
  traitar1$SGB <- apply(traitar1, 1, function(x) {
    paste0("SGB", str_split(x[1], "_")[[1]][1])
  })
  traitar1 <- traitar1[, c(1, ncol(traitar1), 2:(ncol(traitar1) - 1))]
  traitar2 <- traitar0[grepl("SGB", traitar0$X), ]
  traitar2$SGB <- traitar2$X
  traitar2 <- traitar2[, c(1, ncol(traitar2), 2:(ncol(traitar2) - 1))]
  
  #根据gtdbtk结果补充Streptococcus parasanguinis(F04-H14A.bin.23)
  traitar3 <- traitar0[traitar0$X=="F04_H14A.bin.23",]
  traitar3$SGB <- "SGB8071"
  traitar3 <- traitar3[, c(1, ncol(traitar3), 2:(ncol(traitar3) - 1))]
  
  traitar <- rbind(traitar1,traitar2,traitar3)
  colnames(traitar)[1] <- "Species"
  
  b.trans.prev.1 <- b.trans.prev
  b.trans.prev.1$sgb <- apply(b.trans.prev.1, 1, function(x) {
    str_split(x[1], "_")[[1]][1]
  })
  
  traitar <- merge(unique(b.trans.prev.1[, c(8, 4)]), traitar, by.x = "sgb", by.y = "SGB")
  traitar <- traitar[, -c(1, 3)]
  traitar <- traitar[!duplicated(traitar$Species.x), ]
  rownames(traitar) <- traitar$Species.x
  traitar <- traitar[, -1]
  
  base.intervention.prev <- merge(ppi.base.intervention.prev, h2ra.base.intervention.prev, by = "row.names")
  rownames(base.intervention.prev) <- base.intervention.prev$Row.names
  base.intervention.prev <- base.intervention.prev[, -1]
  colnames(base.intervention.prev) <- c("PPI.change", "H2RA.change")
  not.exist <- rownames(base.intervention.prev)[!(rownames(base.intervention.prev) %in% rownames(traitar))]
  base.intervention.prev <- base.intervention.prev[rownames(traitar), , drop = F]
  
  
  
  traitar <- traitar[rownames(base.intervention.prev), ]
  traitar <- merge(traitar, ppi.oral.abun.base.intervention, by = "row.names")
  rownames(traitar) <- traitar$Row.names
  traitar <- traitar[, -1]
  traitar <- traitar[rownames(base.intervention.prev), ]
  
  # total species
  corr.result <- data.frame()
  for (j in colnames(base.intervention.prev)) {
    for (i in colnames(traitar)) {
      corr <- cor.test(traitar[, i], base.intervention.prev[, j], method = "spearman")
      
      corr.result <- rbind(corr.result, c(j, i, corr[["p.value"]], corr[["estimate"]]))
    }
  }
  
  colnames(corr.result) <- c("Group", "Index", "P.value", "Correlation")
  # corr.result <- corr.result[corr.result$P.value<0.1,] %>% na.omit()
  corr.result.cor <- spread(corr.result[, -3], Group, Correlation)
  corr.result.cor$H2RA.change <- as.numeric(corr.result.cor$H2RA.change)
  corr.result.cor$PPI.change <- as.numeric(corr.result.cor$PPI.change)
  rownames(corr.result.cor) <- corr.result.cor$Index
  corr.result.cor <- corr.result.cor[, -1, drop = F]
  corr.result.cor <- corr.result.cor %>% arrange(PPI.change, H2RA.change)
  corr.result.cor <- corr.result.cor[, c(2, 1)]
  
  corr.result.p <- spread(corr.result[, -4], Group, P.value)
  corr.result.p$H2RA.change <- as.numeric(corr.result.p$H2RA.change)
  corr.result.p$PPI.change <- as.numeric(corr.result.p$PPI.change)
  rownames(corr.result.p) <- corr.result.p$Index
  corr.result.p <- corr.result.p[, -1, drop = F]
  corr.result.p <- corr.result.p %>% arrange(PPI.change, H2RA.change)
  corr.result.p <- corr.result.p[, c(2, 1)]
  corr.result.p[corr.result.p < 0.05 & corr.result.p > 0.01] <- "*"
  corr.result.p[corr.result.p < 0.01] <- "**"
  corr.result.p[is.na(corr.result.p)] <- ""
  corr.result.p[corr.result.p > 0.05] <- ""
  corr.result.p[corr.result.p < 0.1 & corr.result.p > 0.05] <- "."
  corr.result.p <- corr.result.p[rownames(corr.result.cor), ]
  
  pdf(file = "03_Figure3/Figure3C.pdf", width = 12, height = 6)
  pheatmap::pheatmap(t(corr.result.cor),
                     cluster_rows = F,
                     cluster_cols = F,
                     na_col = "white",
                     cellwidth = 10,
                     cellheight = 10,
                     display_numbers = t(corr.result.p)
  )
  dev.off()
}

{ # 7.4 细菌的各项指标与transmission prevalence的相关性 ----
  
  base.prev <- b.trans.prev[,c(4:7)]
  base.prev <- base.prev[base.prev$Timepoint=="Baseline",]
  base.prev <- spread(base.prev,Group,Ratio,fill=0)
  rownames(base.prev) <- base.prev$Species
  base.prev <- base.prev[,-c(1,2)]
  
  traitar4 <- traitar[rownames(base.prev),]
  
  # total species
  corr.result <- data.frame()
  for (j in colnames(base.prev)) {
    for (i in colnames(traitar4)) {
      corr <- cor.test(traitar4[, i], base.prev[, j], method = "spearman")
      
      corr.result <- rbind(corr.result, c(j, i, corr[["p.value"]], corr[["estimate"]]))
    }
  }
  
  colnames(corr.result) <- c("Group", "Index", "P.value", "Correlation")
  # corr.result <- corr.result[corr.result$P.value<0.1,] %>% na.omit()
  corr.result.cor <- spread(corr.result[, -3], Group, Correlation)
  corr.result.cor$H2RA <- as.numeric(corr.result.cor$H2RA)
  corr.result.cor$PPI <- as.numeric(corr.result.cor$PPI)
  rownames(corr.result.cor) <- corr.result.cor$Index
  corr.result.cor <- corr.result.cor[, -1, drop = F]
  corr.result.cor <- corr.result.cor %>% arrange(PPI, H2RA)
  corr.result.cor <- corr.result.cor[, c(2, 1)]
  
  corr.result.p <- spread(corr.result[, -4], Group, P.value)
  corr.result.p$H2RA <- as.numeric(corr.result.p$H2RA)
  corr.result.p$PPI <- as.numeric(corr.result.p$PPI)
  rownames(corr.result.p) <- corr.result.p$Index
  corr.result.p <- corr.result.p[, -1, drop = F]
  corr.result.p <- corr.result.p %>% arrange(PPI, H2RA)
  corr.result.p <- corr.result.p[, c(2, 1)]
  corr.result.p[corr.result.p < 0.05 & corr.result.p > 0.01] <- "*"
  corr.result.p[corr.result.p < 0.01] <- "**"
  corr.result.p[is.na(corr.result.p)] <- ""
  corr.result.p[corr.result.p > 0.05] <- ""
  corr.result.p[corr.result.p < 0.1 & corr.result.p > 0.05] <- "."
  corr.result.p <- corr.result.p[rownames(corr.result.cor), ]
  
  pdf(file = "03_Figure3/Figure3C.pdf", width = 12, height = 6)
  pheatmap::pheatmap(t(corr.result.cor),
                     cluster_rows = F,
                     cluster_cols = F,
                     na_col = "white",
                     cellwidth = 10,
                     cellheight = 10,
                     display_numbers = t(corr.result.p)
  )
  dev.off()
  
  
}

{ # 7.5 多元线性回归 ----
  
  data <- merge(base.intervention.prev,traitar,by="row.names")
  
  #ppi
  data1 <- data[,-c(1,3)]
  
  # 多元
  data1 <- data1[,c("PPI.change",rownames(corr.result.p)[corr.result.p$PPI.change<0.05][1:9])]
  ppi.lm <- lm(PPI.change~.,data=data1)
  summary(ppi.lm)
  
  # 单元
  
  
  
}

# 8. D: PTR ---------------------------------------------------------------
{ # 8.1 整理PTR数据 ----

  gut.ptr <- read.table("00_variables/gut_PTR.txt", sep = "\t", header = T)
  colnames(gut.ptr)[1] <- "Species"
  gut.ptr <- gather(gut.ptr, sample, ptr, -Species)
  gut.ptr$sample <- gsub("\\.", "-", gut.ptr$sample)
  gut.ptr.sgb <- gut.ptr[grepl("^[0-9]", gut.ptr$Species), ]
  gut.ptr.sgb$sgb <- apply(gut.ptr.sgb, 1, function(x) {
    str_split(x[1], "_")[[1]][1]
  })
  gut.ptr.sgb$sgb <- paste0("SGB", gut.ptr.sgb$sgb)
  sgb.list <- filtered.abun[["oral.species"]]
  sgb.list$SGB <- apply(sgb.list, 1, function(x) {
    str_split(x[9], "_")[[1]][3]
  })
  gut.ptr.sgb <- merge(gut.ptr.sgb, sgb.list[, c("SGB", "s")], by.x = "sgb", by.y = "SGB")
  gut.ptr.sgb <- merge(gut.ptr.sgb, sample.list[, -c(2)], by = "sample")
  # gut.ptr.sgb <- na.omit(gut.ptr.sgb)
  gut.ptr.sgb$ptr[is.na(gut.ptr.sgb$ptr)] <- 0

  gut.ptr.sgb$s <- str_sub(gut.ptr.sgb$s, 4, -1)
  gut.ptr.sgb$s <- gsub("_", " ", gut.ptr.sgb$s)
  gut.ptr.sgb$sgb <- gsub("t__", "", gut.ptr.sgb$sgb)
  gut.ptr.sgb$s <- paste0(gut.ptr.sgb$s, " [", gut.ptr.sgb$sgb, "]")
}

{ # 8.4 查看在Intervention 和 Baseline 中有显著差别的细菌 ----

  find.ptr.diff.species <- function(group.select) {
    ptr.diff <- gut.ptr.sgb %>% subset(group %in% group.select)

    ptr.diff.species <- data.frame()
    for (i in group.select) {
      data <- ptr.diff[ptr.diff$group == i, ]
      id <- sample.list$id[sample.list$group == i] %>% unique()

      for (j in unique(data$Species)) {
        temp <- data[data$Species == j, c("id", "timepoint", "ptr")]
        if (length(unique(temp$timepoint)) == 2) {
          temp <- spread(temp, timepoint, ptr, fill = 0)
          wilcox <- wilcox.test(as.numeric(temp[, 2]), as.numeric(temp[, 3]), paired = T)

          if (!is.nan(wilcox[["p.value"]]) & wilcox[["p.value"]] < 0.05) {
            ptr.diff.species <- rbind(ptr.diff.species, c(i, j, wilcox[["p.value"]]))
          }
        } else {
          exist.time <- unique(temp$timepoint)
          not.exist.time <- timepoint.select[!(timepoint.select %in% exist.time)]
          ptr.diff.species <- rbind(ptr.diff.species, c(i, j, paste0("only exist in ", exist.time)))
        }
      }
    }
    colnames(ptr.diff.species) <- c("Group", "Species", "P.value")

    return(ptr.diff.species)
  }

  ptr.diff.intervention.base <- find.ptr.diff.species(group.select = c("H2RA", "PPI"))

  ptr.diff.intervention.base <- merge(ptr.diff.intervention.base, unique(gut.ptr.sgb[, c("Species", "s")]), by = "Species")
  ptr.diff.intervention.base <- ptr.diff.intervention.base %>% arrange(Group)

  ptr.diff.intervention.base.ptr <- gut.ptr.sgb[gut.ptr.sgb$Species %in% ptr.diff.intervention.base$Species, ]
  ptr.diff.intervention.base.ptr$s <- factor(ptr.diff.intervention.base.ptr$s, levels = unique(ptr.diff.intervention.base$s))
  ptr.diff.intervention.base.ptr$group <- factor(ptr.diff.intervention.base.ptr$group, levels = c("PPI", "H2RA"))

  plot.data <- ptr.diff.intervention.base.ptr[, c(4:8)]

  signif.result <- data.frame()
  for (i in unique(plot.data$group)) {
    for (j in unique(plot.data$s)) {
      temp <- plot.data[plot.data$group == i & plot.data$s == j, ]
      temp <- spread(temp, timepoint, ptr, fill = 0)
      wilcox <- wilcox.test(temp$Baseline, temp$Intervention, paired = T)
      signif.result <- rbind(signif.result, c(i, j, wilcox[["p.value"]]))
    }
  }
  colnames(signif.result) <- c("Group", "Species", "P.value")
  signif.result$sig[signif.result$P.value > 0.01 & signif.result$P.value < 0.05] <- "*"
  signif.result$sig[signif.result$P.value > 0.001 & signif.result$P.value < 0.01] <- "**"


  ptr.diff.intervention.base.plot <- plot.data %>%
    ggplot(aes(x = ptr, y = s, fill = timepoint)) +
    geom_boxplot(width = 0.6, lwd = 0.2, outlier.size = 1) +
    theme_bw() +
    facet_grid(. ~ group) +
    scale_x_log10() +
    scale_y_discrete(position = "right") +
    scale_fill_manual(name = "Timepoint", values = timepoint.color) +
    theme(
      strip.background = element_rect(color = "black", fill = "#FFFFFF"),
      axis.text.y = element_text(face = "italic"),
      text = element_text(family = "serif"),
      # legend.position = "none",
    ) +
    labs(y = "", x = "Peak-to-Trough Ratio (PTR, Log-scaled)") +
    geom_signif(comparisons = list(c("Baseline", "Intervention")))

  ptr.diff.intervention.base.plot
  # ggsave(ptr.diff.intervention.base.plot,filename="03_Figure3/Figure3C.pdf",width=7,height = 6)

  b.plot <- (b.trans.prev.plot | b.intervention.base.plot | ptr.diff.intervention.base.plot) + plot_layout(widths = c(2, 1, 2))
  b.plot
  ggsave(b.plot, filename = "03_Figure3/Figure3B.pdf", width = 18, height = 7)
}


