################################################
### Topological and functional analysis of VTGs

### Specific and pan targets
gene_virus_num <- read.table("../intermediate/viralNum_target_a_gene.txt",
                             header = TRUE, sep = "\t", quote = "")
node_info <- read.table('../intermediate/nodes_info.txt', 
                        header = TRUE, sep = '\t')

vtg_info <- merge(gene_virus_num, node_info, by = 'Gene')

library(ggpubr)
p_degree <- ggplot(subset(vtg_info, Degree != '-'), 
                   aes(Virally_targeted, as.numeric(as.character(Degree)))) +
  geom_boxplot()+
  scale_x_discrete(limits = c('Specific', 'Pan')) +
  stat_compare_means(comparisons = list(c('Specific', 'Pan'))) +
  labs(x = '', y = 'Degree') +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust=1, hjust=1))

p_eigenvector <- ggplot(subset(vtg_info, Degree != '-'), 
                        aes(Virally_targeted, as.numeric(as.character(Eigenvector_centrality)))) +
  geom_boxplot()+
  scale_x_discrete(limits = c('Specific', 'Pan')) +
  stat_compare_means(comparisons = list(c('Specific', 'Pan'))) +
  labs(x = '', y = 'Eigenvector centrality') +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust=1, hjust=1))

p_3b <- ggpubr::ggarrange(p_degree, p_eigenvector, nrow=1, labels = c("B", ""))
ggsave(p_3b,
       filename = "../result/Fig_3B.pdf",
       width = 5.4,
       height = 5.1)

# Gene ontology over-representation analysis
library(clusterProfiler)
bp <- read.table('../data/MSigDB7.4/c5.general_bp.v7.4.symbols_GO.ID.txt',
                 header = FALSE, sep = '\t', quote = "")
colnames(bp) <- c("GO_ID", 'Description', 'Gene')
bp2 <- bp[bp$Gene %in% node_info$Gene,]

jaccard <- function(x,y){
  jaccard_index <- length(intersect(x, y))/length(union(x, y))
  jaccard_index
}

pan_genes <- vtg_info[vtg_info$Virally_targeted == 'Pan',1]
pan_go <- enricher(pan_genes,
                   pvalueCutoff = 0.01,
                   pAdjustMethod = "BH",
                   universe=as.character(node_info$Gene),
                   qvalueCutoff = 0.01,
                   minGSSize = 1,
                   maxGSSize = 50000,
                   TERM2GENE= bp2[,2:3])

pan_go2 <- pan_go@result
pan_go2 <- pan_go2[pan_go2$pvalue < 0.01 & pan_go2$p.adjust < 0.01,c(2:6,9,8)]

if(nrow(pan_go2)>0){
  reduRow <- c()
  for(n in 1:nrow(pan_go2)){
    for(m in 1:nrow(pan_go2)){
      if(n != m){
        list1 = strsplit(pan_go2$geneID[n], split = '/', fixed = T)[[1]]
        list2 = strsplit(pan_go2$geneID[m], split = '/', fixed = T)[[1]]
        ji = jaccard(list1, list2)
        if(ji >= 0.9){
          if(pan_go2$p.adjust[n] > pan_go2$p.adjust[m]){
            reduRow <- c(reduRow, n)
          } else {
            reduRow <- c(reduRow, m)
          }
        }
      }
    }
  }
  
  if(length(reduRow)>0){
    pan_go3 <- pan_go2[-unique(reduRow),]
  }
}

pan_go3 <- pan_go3 %>% 
  mutate(GeneRatio = paste0("(", GeneRatio, ")"), BgRatio = paste0("(", BgRatio, ")"))
write.table(pan_go3, "../result/panVTG_GO_bp.txt",
            quote = FALSE, sep = '\t', row.names = FALSE)

pan_go4 <- pan_go3[1:20,]

p_go <- ggplot(pan_go4, aes('', Description, size = -log10(p.adjust))) +
  geom_point(color = '#FC4E07') + 
  scale_y_discrete(limits = pan_go4[,1][rev(order(pan_go4$p.adjust[1:20]))]) +
  scale_size_continuous(range = c(2,8)) +
  labs(x = '', y = '') +
  mytheme +
  theme(axis.ticks.x = element_blank())

ggsave(p_go,
       filename = "../result/Fig_3C.pdf",
       width = 9.6,
       height = 7)

# overlap between pan VTGs and essential genes
essen_genes <- read.table("../data/essential_genes.txt", 
                          header = FALSE, quote = "", sep = "\t")
essen_genes2 <- intersect(essen_genes$V1, node_info$Gene)
library(Vennerable)

pdf(file = "../result/Fig_S3C.pdf")
pe_plot <- Venn(list("Pan VTGs" = pan_genes,
                     "Essential genes" = essen_genes2))
plot(pe_plot, doWeight = T)
dev.off()

phyper(length(intersect(essen_genes2, pan_genes))-1, 
       length(essen_genes2), 
       17402, 
       length(pan_genes),
       lower.tail = F)

### Virally-targeted genes
vtgs <- read.table("../result/virally-targeted_genes29.txt", 
                   header = TRUE, sep = '\t', quote = "")

node_centrality <- node_info %>%
  filter(Degree != '-') %>%
  dplyr::select(Gene, Degree, Eigenvector_centrality) %>%
  mutate(Degree = as.numeric(Degree),
         Eigenvector_centrality = as.numeric(Eigenvector_centrality))

vtg_net <- merge(vtgs, node_centrality, by = 'Gene')

vtg_num_net <- data.frame(table(vtg_net$Abbre))
colnames(vtg_num_net) <- c("Abbre", "Number")

random_centrality <- data.frame(Number = numeric(), Abbre = character(),
                                Degree = numeric(), Eigenvector_centrality = numeric())
for (v in vtg_num_net$Abbre) {
  random_node <- data.frame(Number = 1:10000, Abbre = v)
  for (attr in colnames(node_centrality)[2:3]) {
    var <- assign(attr, c())
    
    for (i in 1:10000) {
      set.seed(i)
      node_samp <- node_centrality[sample(1:nrow(node_centrality),
                                          size = vtg_num_net[vtg_num_net$Abbre == v, 2], replace = FALSE),]
      var <- c(var, mean(as.numeric(as.character(node_samp[attr][,1]))))
    }
    
    random_node[attr] <- var
  }
  
  random_centrality <- rbind(random_centrality, random_node)
}

random_centrality <- vtg_net %>%
  dplyr::select(Abbre, Baltimore) %>%
  distinct() %>%
  right_join(random_centrality, by = "Abbre")

random_centrality2 <- random_centrality %>%
  group_by(Abbre) %>%
  dplyr::summarise(Degree_mean = mean(Degree),
                   Degree_sd = sd(Degree),
                   Eigenvector_centrality_mean = mean(Eigenvector_centrality),
                   Eigenvector_centrality_sd = sd(Eigenvector_centrality))

vtg_centrality <- vtg_net %>%
  group_by(Abbre, Baltimore) %>%
  dplyr::summarise(Degree = mean(Degree),
                   Eigenvector_centrality = mean(Eigenvector_centrality)) %>%
  inner_join(random_centrality2, by = "Abbre")

vtg_centrality <- vtg_centrality %>%
  mutate(Degree_Z = (Degree - Degree_mean)/Degree_sd,
         Degree_P = 2*pnorm(-abs(Degree_Z)),
         Eigenvector_centrality_Z = (Eigenvector_centrality - 
                                       Eigenvector_centrality_mean)/Eigenvector_centrality_sd,
         Eigenvector_centrality_P = 2*pnorm(-abs(Eigenvector_centrality_Z)))

signif_fun <- function(x){
  ifelse(x < 0.001, sn <- "***", 
         ifelse(x < 0.01, sn <- "**", 
                ifelse(x < 0.05, sn <- "*", 
                       sn <- "ns")))
  sn
}

vtg_centrality <- vtg_centrality %>%
  mutate(Degree_Signif = signif_fun(Degree_P),
         Eigenvector_centrality_Signif = signif_fun(Eigenvector_centrality_P))

p_v_degree <- ggplot(random_centrality, aes(Abbre, Degree, fill = Baltimore)) +
  geom_boxplot(alpha = 0.75, outlier.shape = 1) +
  scale_x_discrete(limits = vtg_centrality$Abbre[rev(order(vtg_centrality$Degree))]) +
  scale_fill_manual(values = c(brewer.pal(9, 'Set1')[c(1,2,3,7)])) +
  geom_point(data = vtg_centrality, aes(Abbre, Degree), size = 2, 
             color = "red") +
  geom_text(data = vtg_centrality, aes(x = Abbre, y = Degree + 0.5, label = Degree_Signif)) + 
  labs(x="",y="Degree") +
  mytheme +
  theme(legend.title = element_blank(),legend.text = element_text(size = 13),
        legend.position = "top", legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))

p_v_eigenvector <- ggplot(random_centrality, aes(Abbre, 1000*Eigenvector_centrality, fill = Baltimore)) +
  geom_boxplot(alpha = 0.75, outlier.shape = 1) +
  scale_x_discrete(limits = vtg_centrality$Abbre[rev(order(vtg_centrality$Degree))]) +
  scale_fill_manual(values = c(brewer.pal(9, 'Set1')[c(1,2,3,7)])) +
  geom_point(data = vtg_centrality, aes(Abbre, 1000*Eigenvector_centrality), size = 2, 
             color = "red") +
  geom_text(data = vtg_centrality, aes(x = Abbre, y = 1000*Eigenvector_centrality + 0.2, 
                                       label = Eigenvector_centrality_Signif)) + 
  labs(x="",y="Eigenvector centrality (10-3)") +
  mytheme +
  theme(legend.title = element_blank(),legend.text = element_text(size = 13),
        legend.position = "top", legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))

p_3a <- ggpubr::ggarrange(p_v_degree, p_v_eigenvector, nrow=2, labels = c("A", ""))
ggsave(p_3a,
       filename = "../result/Fig_3A.pdf",
       width = 7.2,
       height = 8.5)


### Virus-specific targeted genes and Network topological properties
vtgs2 <- vtgs %>%
  dplyr::select(Gene, Abbre) %>%
  mutate(Number = 1) %>%
  pivot_wider(names_from = Abbre, values_from = Number)
vtgs2[is.na(vtgs2)] <- 0
vtgs2$Number <- apply(vtgs2[,2:ncol(vtgs2)], 1,
                      function(x) sum(x == 1))

vtgs2 <- vtgs2 %>%
  filter(Number == 1) %>%
  dplyr::select(-Number)

virus <- colnames(vtgs2)[-1]
vtgs2$Label <- apply(vtgs2[,2:ncol(vtgs2)], 1,
                     function(x) virus[which(x == 1)])

vtg_num <- vtgs2 %>%
  group_by(Label) %>%
  dplyr::summarise(Number = n())

vtg_num <- vtgs %>%
  dplyr::select(Baltimore, Family, Abbre) %>%
  distinct() %>%
  right_join(vtg_num, by = c("Abbre" = "Label"))

vtg_net <- node_info %>%
  dplyr::select(Gene, Degree, Eigenvector_centrality) %>%
  filter(Degree != "-") %>%
  inner_join(vtgs2[,c(1, ncol(vtgs2))], by = "Gene") %>%
  left_join(vtgs[,3:5], by = c("Label" = "Abbre")) %>%
  distinct() %>%
  dplyr::rename(Abbre = Label) %>% 
  mutate(Degree = as.numeric(Degree),
         Eigenvector_centrality = as.numeric(Eigenvector_centrality))

vtg_num_net2 <- data.frame(table(vtg_net$Abbre))
colnames(vtg_num_net2) <- c("Abbre", "Number")
vtg_num_net2 <- vtg_net %>%
  dplyr::select(Abbre, Baltimore) %>%
  distinct() %>%
  inner_join(vtg_num_net2, by = "Abbre")

p_vtg_num <- ggplot(vtg_num_net2, aes(Abbre, Number, fill = Baltimore)) +
  geom_bar(stat = 'identity', width = .8, alpha = .75) +
  geom_text(aes(x = Abbre, y = Number+10, label = Number, angle = 270)) +
  scale_x_discrete(limits = vtg_num_net2[,1][order(vtg_num_net2$Number)]) +
  scale_y_continuous(expand = c(0.01, 0.03)) +
  scale_fill_manual(values = c(brewer.pal(9, 'Set1')[c(1,2,3,7)])) +
  labs(x = "", y = "Gene Number") +
  coord_flip() +
  mytheme +
  theme(legend.position = c(0.8, 0.8))

vtg_random_centrality <- data.frame(Number = numeric(), Abbre = character(),
                                    Degree = numeric(), Eigenvector_centrality = numeric())
for (v in vtg_num_net2$Abbre) {
  random_node <- data.frame(Number = 1:10000, Abbre = v)
  for (attr in colnames(node_centrality)[2:3]) {
    var <- assign(attr, c())
    
    for (i in 1:10000) {
      set.seed(i)
      node_samp <- node_centrality[sample(1:nrow(node_centrality),
                                          size = vtg_num_net2[vtg_num_net2$Abbre == v, 3], replace = FALSE),]
      var <- c(var, mean(as.numeric(as.character(node_samp[attr][,1]))))
    }
    
    random_node[attr] <- var
  }
  
  vtg_random_centrality <- rbind(vtg_random_centrality, random_node)
}

vtg_random_centrality <- vtg_net %>%
  dplyr::select(Abbre, Baltimore) %>%
  distinct() %>%
  right_join(vtg_random_centrality, by = "Abbre")

vtg_random_centrality2 <- vtg_random_centrality %>%
  group_by(Abbre) %>%
  dplyr::summarise(Degree_mean = mean(Degree),
                   Degree_sd = sd(Degree),
                   Eigenvector_centrality_mean = mean(Eigenvector_centrality),
                   Eigenvector_centrality_sd = sd(Eigenvector_centrality))

vtg_centrality2 <- vtg_net %>%
  group_by(Abbre, Baltimore) %>%
  dplyr::summarise(Degree = mean(Degree),
                   Eigenvector_centrality = mean(Eigenvector_centrality)) %>%
  inner_join(vtg_random_centrality2, by = "Abbre")

vtg_centrality2 <- vtg_centrality2 %>%
  mutate(Degree_Z = (Degree - Degree_mean)/Degree_sd,
         Degree_P = 2*pnorm(-abs(Degree_Z)),
         Eigenvector_centrality_Z = (Eigenvector_centrality - 
                                       Eigenvector_centrality_mean)/Eigenvector_centrality_sd,
         Eigenvector_centrality_P = 2*pnorm(-abs(Eigenvector_centrality_Z)))

vtg_centrality2 <- vtg_centrality2 %>%
  mutate(Degree_Signif = signif_fun(Degree_P),
         Eigenvector_centrality_Signif = signif_fun(Eigenvector_centrality_P))

p_vtg_d <- ggplot(vtg_random_centrality, aes(Abbre, Degree, fill = Baltimore)) +
  geom_boxplot(alpha = 0.75, outlier.shape = 1) +
  scale_x_discrete(limits = vtg_centrality2$Abbre[order(vtg_centrality2$Degree)]) +
  scale_fill_manual(values = c(brewer.pal(9, 'Set1')[c(1,2,3,7)])) +
  geom_point(data = vtg_centrality2, aes(Abbre, Degree), size = 2, 
             color = "red") +
  geom_text(data = vtg_centrality2, aes(x = Abbre, y = Degree + 0.5, label = Degree_Signif)) + 
  labs(x="",y="Degree") +
  coord_flip() +
  mytheme +
  theme(legend.position = c(0.8, 0.8))

p_vtg_e <- ggplot(vtg_random_centrality, aes(Abbre, 1000*Eigenvector_centrality, fill = Baltimore)) +
  geom_boxplot(alpha = 0.75, outlier.shape = 1) +
  scale_x_discrete(limits = vtg_centrality2$Abbre[order(vtg_centrality2$Degree)]) +
  scale_fill_manual(values = c(brewer.pal(9, 'Set1')[c(1,2,3,7)])) +
  geom_point(data = vtg_centrality2, aes(Abbre, 1000*Eigenvector_centrality), size = 2, 
             color = "red") +
  geom_text(data = vtg_centrality2, aes(x = Abbre, y = 1000*Eigenvector_centrality + 0.2, 
                                        label = Eigenvector_centrality_Signif)) + 
  labs(x="",y="Eigenvector centrality (10-3)") +
  coord_flip() +
  mytheme +
  theme(legend.position = c(0.8, 0.8))

p_specific <- ggpubr::ggarrange(p_vtg_num, p_vtg_d, p_vtg_e, nrow=1,
                  labels = c("A", "B", ""))
ggsave(p_specific,
       filename = "../intermediate/viral_specific_targets.pdf",
       width = 15.7,
       height = 8.2)

## Functional enrichment (2021-12-27)
# GO
bp <- read.table('../data/MSigDB7.4/c5.bp.v7.4.symbols_GO.ID.txt',
                 header = FALSE, sep = '\t', quote = "")
colnames(bp) <- c("GO_ID", 'Description', 'Gene')
bp2 <- bp[bp$Gene %in% node_info$Gene,]

vtg_go <- data.frame(Abbr = character(), Description = character(), GeneRatio = character(), 
                     BgRatio = character(), pvalue = character(), p.adjust = character(), 
                     qvalue = character(), geneID = character())
vtgs2 <- as.data.frame(vtgs2)

for (i in unique(vtgs2$Label)){
  go <- enricher(vtgs2[vtgs2$Label == i,1],
                 pvalueCutoff = 0.01,
                 pAdjustMethod = "BH",
                 universe=as.character(node_info$Gene),
                 qvalueCutoff = 0.01,
                 minGSSize = 1,
                 maxGSSize = 50000,
                 TERM2GENE= bp2[,2:3])
  
  go2 <- go@result
  go2 <- go2[go2$pvalue < 0.01 & go2$p.adjust < 0.01,-1]
  go2 <- go2[,1:7]
  go2$Abbr <- rep(i, nrow(go2))
  go2 <- go2[,c(8,1:7)]
  
  if(nrow(go2)>0){
    reduRow <- c()
    for(n in 1:nrow(go2)){
      for(m in 1:nrow(go2)){
        if(n != m){
          list1 = strsplit(go2$geneID[n], split = '/', fixed = T)[[1]]
          list2 = strsplit(go2$geneID[m], split = '/', fixed = T)[[1]]
          ji = jaccard(list1, list2)
          if(ji >= 0.99){
            if(go2$p.adjust[n] > go2$p.adjust[m]){
              reduRow <- c(reduRow, n)
            } else {
              reduRow <- c(reduRow, m)
            }
          }
        }
      }
    }
    
    if(length(reduRow)>0){
      go3 <- go2[-unique(reduRow),]
      vtg_go <- rbind(vtg_go, go3)
    }else{
      vtg_go <- rbind(vtg_go, go2)
    }
  }
  print(i)
}

vtg_go <- vtg_go %>% 
  mutate(GeneRatio = paste0("(", GeneRatio, ")"), BgRatio = paste0("(", BgRatio, ")"))
write.table(vtg_go[,c(1:4,6)], "../result/viral_specificTarget_GO_bp",
            sep = "\t", 
            row.names = FALSE, quote = FALSE)

p_go_fun <- function(v) {
  viral_go <- vtg_go[vtg_go$Cluster == v,]
  p <- ggplot(viral_go, 
              aes('', Description, size = -log10(p.adjust))) +
    geom_point(color = '#FC4E07', aes(size = )) +
    scale_y_discrete(limits = viral_go[,2][rev(order(viral_go$p.adjust))]) +
    labs(x = '', y = '', title = v) +
    mytheme +
    theme(axis.ticks.x = element_blank())
  p
}