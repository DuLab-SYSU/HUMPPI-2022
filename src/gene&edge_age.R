###############################################
### New genes drive the evolution of human PPIN
library(ggsignif)
library(RColorBrewer)
library(tidyverse)
mytheme <- theme(panel.grid.major = element_line(colour=brewer.pal(9,"Pastel1")[9],
                                                 linetype = "longdash"),
                 panel.background = element_rect(fill='transparent', color="#000000"),
                 panel.border=element_rect(fill='transparent', color='black'),
                 axis.text = element_text(size = 16),
                 axis.title = element_text(size = 20),
                 legend.text = element_text(size = 13),
                 legend.title = element_text(size = 16),
                 legend.background = element_blank())

### Node and edge number on each age branch
node_info <- read.table('../intermediate/nodes_info.txt', header = TRUE, sep = '\t')
edge_info <- read.table('../intermediate/edges_info.txt', header = TRUE, sep = '\t')

library(plyr)
node_age_num <- ddply(node_info[node_info$Branch!='-' & node_info$Degree!='-',], 
                            .(Branch), nrow)
edge_age_num <- ddply(edge_info[edge_info$Edge_age!='-',], .(Edge_age), nrow)
edge_age_num$Edge_age <- as.numeric(as.character(edge_age_num$Edge_age))

p_node_age <- ggplot(node_age_num, aes(factor(Branch, levels = rev(0:13)), V1)) +
  geom_line(size = 1.6, group = 1, colour = "#E69F00") +
  scale_y_log10() +
  scale_x_discrete(expand = c(0.02,0.02)) +
  labs(x="Age Branch",y="Node Number") +
  geom_text(label=node_age_num$V1, vjust=0.5, hjust=0.1, color = brewer.pal(8,"Dark2")[8]) +
  coord_flip() +
  mytheme +
  theme(legend.title = element_blank(),legend.position = c(0.9, 0.85),
        axis.text.x = element_text(angle = 60, vjust=1, hjust=1))

p_edge_age <- ggplot(edge_age_num, aes(factor(Edge_age, levels = rev(0:13)), V1)) +
  geom_line(size = 1.6, group = 1, colour = "#56B4E9") +
  scale_y_log10() +
  scale_x_discrete(expand = c(0.02,0.02)) +
  labs(x="Age Branch",y="Edge Number") +
  geom_text(label=edge_age_num$V1, vjust=0.5, hjust=0.1, color = brewer.pal(8,"Dark2")[8]) +
  coord_flip() +
  mytheme +
  theme(legend.title = element_blank(),legend.position = c(0.9, 0.85),
        axis.text.x = element_text(angle = 60, vjust=1, hjust=1))

# Gene number for Different category on each age branch
node_age_num2 <- ddply(node_info[node_info$Branch!='-' & node_info$Degree!='-',], 
                      .(Branch, Gene_class), nrow)

p_node_age2 <- ggplot(node_age_num2, 
                      aes(factor(Branch, levels = rev(0:13)), V1, group = Gene_class,
                          color = factor(Gene_class, levels = c('IRG&VTG', 'IRG', 'VTG', 'Other')))) +
  scale_y_log10() +
  geom_line(size = 1.6) +
  scale_color_manual(values = c('#FF6666', '#00CCCC', '#FFCC00', '#999999')) +
  scale_x_discrete(expand = c(0.02,0.02)) +
  labs(x="Age Branch",y="Node Number",color = "") +
  geom_text(label=node_age_num2$V1, vjust=0.5, hjust=0.1, aes(color = Gene_class)) +
  coord_flip() +
  mytheme +
  theme(legend.title = element_blank(),legend.position = c(0.8, 0.15),
        axis.text.x = element_text(angle = 60, vjust=1, hjust=1))

p_s6bcd<- ggpubr::ggarrange(p_node_age, p_node_age2, p_edge_age, ncol=3, 
                         labels = c("B", "C", "D"))
ggsave(p_s6bcd,
       filename = "../result/Fig_S6BCD.pdf",
       width = 7.4,
       height = 8.2)


### new genes integrated into PPIN
gNetwork_age <- node_info[node_info$Branch!='-' & node_info$Degree!='-',]
gNetwork_age$Degree <- as.numeric(as.character(gNetwork_age$Degree))
gNetwork_age$Eigenvector_centrality <- as.numeric(as.character(gNetwork_age$Eigenvector_centrality))
gNetwork_age$Clustering_coefficient <- as.numeric(as.character(gNetwork_age$Clustering_coefficient))
gNetwork_age$Assortativity <- as.numeric(as.character(gNetwork_age$Assortativity))

branch_time <- data.frame(Branch = 0:13,
                          Node_time = c(454.6, 361.2, 324.5, 220.2, 176.1, 104.7, 97.4,
                                        91, 44.2, 29.6, 18.8, 15.1, 6.1, 0))
divergence_time <- c(500)
for (i in 1:nrow(branch_time)-1) {
  divergence_time <- c(divergence_time, 
                       (branch_time[i,2] + branch_time[i+1,2])/2)
}

branch_time$Divergence_time <- divergence_time
# write.table(branch_time,
#             '~/Documents/Computational_Biology/Data/Gene_age/gene_list/branch_time.txt',
#             quote = FALSE, sep = "\t", row.names = FALSE)

gNetwork_age <- merge(gNetwork_age, branch_time, by.x = "Branch", by.y = "Branch")
gNetwork_age$Link_rate <- gNetwork_age$Degree/gNetwork_age$Divergence_time

gNetwork_age2 <- node_info[node_info$Branch!='-' & node_info$Degree!='-',]
gNetwork_age2 <- gNetwork_age2[gNetwork_age2$Specificity!='-',]
gNetwork_age2$Specificity <- as.numeric(as.character(gNetwork_age2$Specificity))

gNetwork_age2 <- merge(gNetwork_age2, branch_time, by.x = "Branch", by.y = "Branch")

## Comparsion of degree, centrality and expression specificity across 4 categories of genes
my_comparsion <- list()
for (i in 1:6) {
  my_comparsion[[i]] <- combn(unique(as.character(gNetwork_age$Gene_class)), 2)[,i]
}

p1 <- ggplot(gNetwork_age, aes(factor(Gene_class, levels = c('IRG&VTG', 'IRG', 'VTG', 'Other')), 
                               Degree)) +
  geom_boxplot(aes(color = Gene_class)) +
  scale_color_manual(limits = c('IRG&VTG', 'IRG', 'VTG', 'Other'),
                     values = c('#FF6666', '#00CCCC', '#FFCC00', '#999999')) +
  geom_signif(comparisons = my_comparsion,
              map_signif_level = TRUE,
              textsize = 5, test = wilcox.test, step_increase = .08) +
  labs(x="Gene Categories",
       y="Degree in the Human PPIN") +
  mytheme +
  theme(plot.margin = unit(x=c(6,6,6,6),units = 'mm'),
        legend.title = element_blank(),legend.text = element_text(size = 13),
        legend.position = "top")

p2 <- ggplot(gNetwork_age, aes(factor(Gene_class, levels = c('IRG&VTG', 'IRG', 'VTG', 'Other')),
                               Eigenvector_centrality)) +
  geom_boxplot(aes(color = Gene_class)) +
  scale_color_manual(limits = c('IRG&VTG', 'IRG', 'VTG', 'Other'),
                     values = c('#FF6666', '#00CCCC', '#FFCC00', '#999999')) +
  geom_signif(comparisons = my_comparsion,
              map_signif_level = TRUE,
              textsize = 5, test = wilcox.test, step_increase = .08) +
  labs(x="Gene Categories",
       y="Eigenvector Centrality in the Human PPIN") +
  mytheme +
  theme(plot.margin = unit(x=c(6,6,6,6),units = 'mm'),
        legend.title = element_blank(),legend.text = element_text(size = 13),
        legend.position = "top")

p3 <- ggplot(gNetwork_age2, aes(factor(Gene_class, levels = c('IRG&VTG', 'IRG', 'VTG', 'Other')), 
                                Specificity)) +
  geom_boxplot(aes(color = Gene_class)) +
  scale_color_manual(limits = c('IRG&VTG', 'IRG', 'VTG', 'Other'),
                     values = c('#FF6666', '#00CCCC', '#FFCC00', '#999999')) +
  geom_signif(comparisons = my_comparsion,
              map_signif_level = TRUE,
              textsize = 5, test = wilcox.test, step_increase = .08) +
  labs(x="Gene Categories",
       y="Specificity in the Human PPIN") +
  mytheme +
  theme(plot.margin = unit(x=c(6,6,6,6),units = 'mm'),
        legend.title = element_blank(),legend.text = element_text(size = 13),
        legend.position = "top")

p_s7 <- ggpubr::ggarrange(p1, p2, p3, 
                  nrow = 1, labels = c("A", "B", "C"))
ggsave(p_s7,
       filename = "../result/Fig_S7.pdf",
       width = 15.5,
       height = 5.8)

# Average degree for different gene category on evolutionay branch
degree_age <- Rmisc::summarySE(gNetwork_age, measurevar="Degree", 
                               groupvars = c("Gene_class", "Divergence_time"))

pd <- position_dodge(0.1)
library(ggpmisc)

p_degree <- ggplot(degree_age, aes(Divergence_time, Degree, 
                       color = factor(Gene_class, levels = c('IRG&VTG', 'IRG', 'VTG', 'Other')))) + 
  geom_point(position = pd, size=5.5, alpha = 0.8) +
  scale_color_manual(values = c('#FF6666', '#00CCCC', '#FFCC00', '#999999')) +
  geom_errorbar(aes(ymin=Degree-se,ymax=Degree+se),width=.4, position=pd) +
  scale_x_reverse() +
  ylim(0, 50) +
  labs(x="Divergence time (Million Year Ago)",
       y="Average Degree in the Human PPIN") +
  mytheme +
  theme(plot.margin = unit(x=c(6,6,6,6),units = 'mm'),
        legend.title = element_blank(),legend.text = element_text(size = 13),
        legend.position = "top") +
  stat_smooth(formula = y ~ x, method = "lm", se=FALSE) +
  stat_fit_glance(aes(label = paste("R-squared = ", signif(..r.squared.., digits = 3),
                                    ", p = ", signif(..p.value.., digits = 4), sep = "")),
                  label.x = 0.1, size = 4)

# Average eigenvector_centrality for different gene category on evolutionay branch
eigenvector_age <- Rmisc::summarySE(gNetwork_age, measurevar="Eigenvector_centrality", 
                                groupvars = c("Gene_class", "Divergence_time"))

p_eigenvector <- ggplot(eigenvector_age, aes(Divergence_time, Eigenvector_centrality, 
                                     color = factor(Gene_class, levels = c('IRG&VTG', 'IRG', 'VTG', 'Other')))) + 
  geom_point(position = pd, size=5.5, alpha = 0.8) +
  scale_color_manual(values = c('#FF6666', '#00CCCC', '#FFCC00', '#999999')) +
  geom_errorbar(aes(ymin=Eigenvector_centrality-se,ymax=Eigenvector_centrality+se),width=.4, position=pd) +
  scale_x_reverse() +
  ylim(0, 0.013) +
  labs(x="Divergence time (Million Year Ago)",
       y="Average Eigenvector Centrality in the Human PPIN") +
  mytheme +
  theme(plot.margin = unit(x=c(6,6,6,6),units = 'mm'),
        legend.title = element_blank(),legend.text = element_text(size = 13),
        legend.position = "top") +
  stat_smooth(formula = y ~ x, method = "lm", se=FALSE) +
  stat_fit_glance(aes(label = paste("R-squared = ", signif(..r.squared.., digits = 3),
                                    ", p = ", signif(..p.value.., digits = 4), sep = "")),
                  label.x = 0.1, size = 4)

### (Tissue-specific expression) gene pleiotropy
specificity_age <- Rmisc::summarySE(gNetwork_age2, measurevar="Specificity", 
                                     groupvars = c("Gene_class", "Divergence_time"))

p_specificity <- ggplot(specificity_age, aes(Divergence_time, Specificity, 
                                               color = factor(Gene_class, levels = c('IRG&VTG', 'IRG', 'VTG', 'Other')))) + 
  geom_point(position = pd, size=5.5, alpha = 0.8) +
  scale_color_manual(values = c('#FF6666', '#00CCCC', '#FFCC00', '#999999')) +
  geom_errorbar(aes(ymin=Specificity-se,ymax=Specificity+se),width=.4, position=pd) +
  scale_x_reverse() +
  ylim(0.15, 1) +
  labs(x="Divergence time (Million Year Ago)",
       y="Average Specificity in the Human PPIN") +
  mytheme +
  theme(plot.margin = unit(x=c(6,6,6,6),units = 'mm'),
        legend.title = element_blank(),legend.text = element_text(size = 13),
        legend.position = "top") +
  stat_smooth(formula = y ~ x, method = "lm", se=FALSE) +
  stat_fit_glance(aes(label = paste("R-squared = ", signif(..r.squared.., digits = 3),
                                    ", p = ", signif(..p.value.., digits = 4), sep = "")),
                  label.x = 0.1, size = 4)
  
### Communities evolution (or edges evolution)
## Edge assigns to intra- or inter-communities 
edge_info <- edge_info[edge_info$Category!='-',]

# All genes
edge_info2 <- edge_info[edge_info$Category!='-' & edge_info$Edge_age!='-',]

require(plyr)
edge_num <- ddply(edge_info2, .(Category, Edge_age), nrow)
edge_num <- dplyr::rename(edge_num, Number=V1)
edge_num <- merge(edge_num, branch_time, by.x = 'Edge_age', by.y = 'Branch')

Frequency <- c()
for (i in 1:nrow(edge_num)) {
  age <- edge_num[i,1]
  freq <- edge_num[i,3]/sum(edge_num[edge_num$Edge_age==age,3])
  freq <- round(freq, 4)
  Frequency <- c(Frequency, freq)
}

edge_num$Frequency <- Frequency

P_intra_freq <- ggplot(subset(edge_num, Category == 'Intra-'), 
                       aes(Divergence_time, Frequency*100)) +
  geom_point(size=5.5, color = '#00A087B2') +
  scale_x_reverse() +
  labs(x="Divergence time (Million Year Ago)",
       y="Proportion of links located in\nintra-communities (%)") +
  mytheme +
  theme(plot.margin = unit(x=c(6,6,6,6),units = 'mm')) +
  stat_smooth(formula = y ~ x, method = "lm", se=FALSE, color = '#E64B35B2') +
  stat_fit_glance(aes(label = paste("R-squared = ", signif(..r.squared.., digits = 3),
                                    ", p = ", signif(..p.value.., digits = 4), sep = "")),
                  label.x = 0.1, size = 4)

# Different gene category (2021-09-22)
edge_num2 <- ddply(edge_info2, .(Category, Edge_age, Edge_class), nrow)
Frequency <- c()
for (i in 1:nrow(edge_num2)) {
  age <- edge_num2[i,2]
  Class <- edge_num2[i,3]
  freq <- edge_num2[i,4]/
    sum(edge_num2[edge_num2$Edge_age==age & edge_num2$Edge_class==Class,4])
  freq <- round(freq, 4)
  Frequency <- c(Frequency, freq)
}

edge_num2$Frequency <- Frequency
edge_num2 <- merge(edge_num2, branch_time, by.x = 'Edge_age', by.y = 'Branch')

P_intra_freq2 <- ggplot(subset(edge_num2, Category == 'Intra-'), 
       aes(Divergence_time, Frequency*100, 
           color = factor(Edge_class, levels = c('IRE&VTE', 'IRE', 'VTE', 'Other')))) +
  geom_point(position = pd, size=5.5, alpha = 0.8) +
  scale_color_manual(values = c('#FF6666', '#00CCCC', '#FFCC00', '#999999')) +
  scale_x_reverse() +
  labs(x="Divergence time (Million Year Ago)",
       y="Proportion of links located in\nintra-communities (%)") +
  mytheme +
  theme(plot.margin = unit(x=c(6,6,6,6),units = 'mm'),
        legend.title = element_blank(),legend.text = element_text(size = 13),
        legend.position = "top") +
  stat_smooth(formula = y ~ x, method = "lm", se=FALSE) +
  stat_fit_glance(aes(label = paste("R-squared = ", signif(..r.squared.., digits = 3),
                                    ", p = ", signif(..p.value.., digits = 4), sep = "")),
                  label.x = 0.1, size = 4)

p_4abcdf <- ggpubr::ggarrange(p_degree, p_eigenvector, p_specificity,
                  P_intra_freq, P_intra_freq2, 
                  ncol=3, nrow = 2, labels = c("A", "B", "C", "D", "F"))
ggsave(p_4abcdf,
       filename = "../result/Fig_4ABCDF.pdf",
       width = 12.8,
       height = 8.2)

# IRE&VTE vs all edges (hypergeometric test)
edge_num3 <- edge_info2 %>%
  group_by(Edge_age) %>%
  dplyr::summarise(Number = n()) %>%
  dplyr::rename(Total_num = Number)

edge_num3 <- edge_info2 %>%
  group_by(Category, Edge_age) %>%
  dplyr::summarise(Number = n()) %>%
  filter(Category == "Intra-") %>%
  dplyr::rename(Total_Intra_num = Number) %>%
  inner_join(edge_num3, by = "Edge_age")

IRE_VTE_num <- edge_info2 %>%
  filter(Edge_class == "IRE&VTE") %>% 
  group_by(Edge_age) %>%
  dplyr::summarise(Number = n()) %>%
  dplyr::rename(IRE_VTE_num = Number)

IRE_VTE_num <- edge_info2 %>%
  filter(Edge_class == "IRE&VTE") %>% 
  group_by(Category, Edge_age) %>%
  dplyr::summarise(Number = n()) %>%
  filter(Category == "Intra-") %>%
  dplyr::rename(IRE_VTE_Intra_num = Number) %>%
  inner_join(IRE_VTE_num, by = "Edge_age")

edge_num3 <- edge_num3 %>%
  inner_join(IRE_VTE_num, by = c("Edge_age", "Category"))

edge_num3$P_value <- phyper(edge_num3$IRE_VTE_Intra_num-1, 
                            edge_num3$Total_Intra_num,
                            edge_num3$Total_num, 
                            edge_num3$IRE_VTE_num, 
                            lower.tail = TRUE)

signif_fun <- function(x){
  s <- c()
  for (i in x) {
    ifelse(i < 0.001, sn <- "***", 
           ifelse(i < 0.01, sn <- "**", 
                  ifelse(i < 0.05, sn <- "*", 
                         sn <- "ns")))
    s <- c(s, sn)
  }
  s
}

edge_num3 <- edge_num3 %>%
  mutate(Signif = signif_fun(P_value))


