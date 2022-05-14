### Node (gene) information

# Protein-coding genes 
genes <- read.table('../data/protein_coding_genes.txt', header = FALSE)
colnames(genes) <- c('Gene')

# Immune-related genes (IRGs)
irgs <- read.table('../data/IRGs/all_IRGs.txt', header = FALSE)
colnames(irgs) <- c('Gene')

# Virally-targeted genes (VTGs)
vtgs <- read.table('../intermediate/viralNum_target_a_gene.txt',
                   header = TRUE, sep = '\t', quote = "")

# Node topological properties
centrality <- read.table('../intermediate/Gene_centrality.txt',
                         header = TRUE, sep = '\t')

# Cluster of each gene in
clusters <- read.table('../result/Module_members.txt',
                       header = TRUE, sep = '\t')
clusters <- clusters[1:9004,]


# Gene age
age <- read.table('../data/human_name_age.txt',
                  header = TRUE, sep = '\t')

# Combine all the gene information 
genes$Gene_class <- apply(genes[,1,drop = FALSE], 1,
                               function(x) ifelse((x %in% irgs$Gene) && (x %in% vtgs$Gene),
                                                  'IRG&VTG', ifelse(x %in% irgs$Gene,
                                                                    'IRG', ifelse(x %in% vtgs$Gene,
                                                                                  'VTG', 'Other'))))

genes <- merge(genes, age, by.x = 'Gene', by.y = 'Gene', all.x = TRUE)
genes <- merge(genes, centrality, by.x = 'Gene', by.y = 'Gene', all.x = TRUE)
genes <- merge(genes, clusters, by.x = 'Gene', by.y = 'Gene', all.x = TRUE)
genes[is.na(genes)] <- "-"

tissue_exp <- read.table("../data/GTEx_Analysis_gene_tissue-specific_expression.gct",
                         header = TRUE, sep = '\t')
tissue_exp <- tissue_exp[,c(1,44:46)]

node_info <- merge(genes, tissue_exp, by.x = "Gene", by.y = "Description", all.x = TRUE)
node_info[is.na(node_info)] <- "-"
write.table(node_info,
            '../intermediate/nodes_info.txt',
            quote = FALSE, sep = "\t", row.names = FALSE)

### The landscape of human protein interactome
## The number of nodes and interactions in different HT experiments
Datasets <- c('Lit-2020', 'Luck et al. 2020', 'BioPlex 2.0', 'BioPlex 1.0', 
              'Hein et al. 2015', 'Rolland et al. 2014', 'Drew et al. 2017')
Methods <- c('Literature-curation', 'Y2H', 'AP-MS', 'AP-MS', 'AP-MS', 'Y2H', 'Literature-curation')
Interactions <- c(78261, 52569, 56553, 23744, 28504, 13944,  56735)
Genes <- c(11202, 8275, 10961, 7668, 5462, 4303,  7777)

interaction_num <- data.frame(Datasets = Datasets, Methods = Methods,
                              Interactions = Interactions, Genes = Genes)

p1 <- ggplot(interaction_num, aes(Genes, Interactions, 
                                  color = factor(Methods, levels = c('AP-MS', 'Y2H', 'Literature-curation')))) +
  geom_point(aes(size = Interactions)) +
  scale_size_continuous(range = c(4,30)) +
  labs(x = 'Gene count', y = 'Interaction count', color = "Methods") +
  xlim(3800,12000) +
  ylim(10000, 90000) +
  geom_text(aes(label = Datasets), hjust=0.5, vjust=-4) +
  #coord_fixed() +
  mytheme +
  theme(legend.position = c(0.15, 0.7),
        plot.margin = unit(x=c(6,6,6,6),units = 'mm'))

## ## Power-law degree distributions (2020-11-04)
node_info <- read.table('../intermediate/nodes_info.txt', 
                        sep = '\t', header = TRUE)
Degree <- subset(node_info, Degree != '-')

link_num <- data.frame(table(log10(as.numeric(as.character(Degree$Degree)))))

link_num$fraction <- log10(link_num$Freq/sum(link_num$Freq))
link_num$Var1 <- as.numeric(as.character(link_num$Var1))

p2 <- ggplot(link_num,aes(Var1, fraction)) + 
  geom_point(color='#999999', size = 3) +
  scale_y_continuous(limits = c(-4.6,0)) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  labs(x="Number of links (log10)",y="Fraction of genes (log10)") +
  mytheme +
  theme(plot.margin = unit(x=c(6,6,6,6),units = 'mm')) +
  geom_smooth(method = "lm", fill = NA, color='red', formula = y ~ x) +
  geom_text(x=2.5, y=-0.5, label = expression(italic("P(k)~k")^ italic("-1.7")), size = 6)

p_1bc <- ggpubr::ggarrange(p1, p2, nrow=1, labels = c("B", "C"))
ggsave(p_1bc,
       filename = "../result/Fig_1BC.pdf",
       width = 16.5,
       height = 8.2)

### Network properties randomly selected gene subsets
node_info <- read.table('../intermediate/nodes_info.txt', 
                        sep = '\t', header = TRUE)
node_centrality <- subset(node_info, Degree != '-')

IRG <- subset(node_centrality, Gene_class == 'IRG&VTG' | Gene_class == 'IRG')

random_node <- data.frame(Number = 1:10000)
for (attr in colnames(node_centrality)[4:9]) {
  var <- assign(attr, c())
  
  for (i in 1:10000) {
    set.seed(i)
    node_samp <- node_centrality[sample(1:nrow(node_centrality), 
                                        size = nrow(IRG), replace = FALSE),]
    var <- c(var, mean(as.numeric(as.character(node_samp[attr][,1]))))
  }
  
  random_node[attr] <- var
}

# statistics
annotation <- data.frame(Vars = c('Observed', 'Z-Score', 'P-Value', 'Mean', 'Sd'))
for (attr in colnames(node_centrality)[4:9]) {
  Obs <- mean(as.numeric(as.character(IRG[attr][,1])))
  Mean <- mean(random_node[attr][,1])
  Sd <- sd(random_node[attr][,1])
  Zsc <- (Obs - Mean)/Sd
  Pval <- 2*pnorm(-abs(Zsc))
  var <- assign(attr, c(Obs, Zsc, Pval, Mean, Sd))
  
  annotation[attr] <- var
}


# plots
p_k <- ggplot(random_node, aes(Degree)) +
  geom_histogram(bins = 80, aes(y=..density..), fill = brewer.pal(8,"Set2")[8],
                 color = brewer.pal(9,"Set1")[9]) +
  stat_function(aes(Degree),fun=dnorm, colour="steelblue", 
                xlim = c(min(random_node$Degree), max(random_node$Degree)),
                args = list(mean = mean(random_node$Degree), sd = sd(random_node$Degree))) +
  geom_point(aes(x = mean(as.numeric(as.character(IRG$Degree))), y=0), 
             color = 'red') +
  geom_segment(aes(x = mean(as.numeric(as.character(IRG$Degree))), y=0.25,
                   xend = mean(as.numeric(as.character(IRG$Degree))), yend=0.02),
               arrow = arrow(length = unit(0.3, "cm")), color = 'red') +
  scale_y_continuous(expand = c(0.005,0.005)) +
  annotate("text", x = 22, y = 0.3, color = "red",
           label = "Observed: 24.115\nZ-Score: 26.016\nP-Value < 2.2e-16") +
  annotate("text", x = max(random_node$Degree)+1, y = 0.9, color = "steelblue",
           label = "Mean: 13.962\nSd: 0.390") +
  labs(x = 'Mean Degree', y = 'Density') +
  mytheme

p_C <- ggplot(random_node, aes(Clustering_coefficient)) +
  geom_histogram(bins = 30, aes(y=..density..), fill = brewer.pal(8,"Set2")[8],
                 color = brewer.pal(9,"Set1")[9]) +
  stat_function(aes(Clustering_coefficient),fun=dnorm, colour="steelblue", 
                xlim = c(min(random_node$Clustering_coefficient), 
                         max(random_node$Clustering_coefficient)),
                args = list(mean = mean(random_node$Clustering_coefficient), 
                            sd = sd(random_node$Clustering_coefficient))) +
  geom_point(aes(x = mean(as.numeric(as.character(IRG$Clustering_coefficient))), y=0), 
             color = 'red') +
  geom_segment(aes(x = mean(as.numeric(as.character(IRG$Clustering_coefficient))), y=25,
                   xend = mean(as.numeric(as.character(IRG$Clustering_coefficient))), yend=2),
               arrow = arrow(length = unit(0.3, "cm")), color = 'red') +
  scale_y_continuous(expand = c(0.005,0.005)) +
  annotate("text", x = 0.151, y = 32, color = "red",
           label = "Observed: 0.153\nZ-Score: -2.344\nP-Value < 2.2e-16") +
  annotate("text", x = 0.17, y = 110, color = "steelblue",
           label = "Mean: 0.160\nSd: 0.003") +
  labs(x = 'Mean Clustering Coefficient', y = 'Density') +
  mytheme

p_NC <- ggplot(random_node, aes(Assortativity)) +
  geom_histogram(bins = 30, aes(y=..density..), fill = brewer.pal(8,"Set2")[8],
                 color = brewer.pal(9,"Set1")[9]) +
  stat_function(aes(Assortativity),fun=dnorm, colour="steelblue", 
                xlim = c(min(random_node$Assortativity), 
                         max(random_node$Assortativity)),
                args = list(mean = mean(random_node$Assortativity), 
                            sd = sd(random_node$Assortativity))) +
  geom_point(aes(x = mean(as.numeric(as.character(IRG$Assortativity))), y=0), 
             color = 'red') +
  geom_segment(aes(x = mean(as.numeric(as.character(IRG$Assortativity))), y=0.1,
                   xend = mean(as.numeric(as.character(IRG$Assortativity))), yend=0.005),
               arrow = arrow(length = unit(0.3, "cm")), color = 'red') +
  scale_y_continuous(expand = c(0.005,0.005)) +
  annotate("text", x = 75.5, y = 0.12, color = "red",
           label = "Observed: 74.787\nZ-Score: -6.624\nP-Value < 2.2e-16") +
  annotate("text", x = 82.5, y = 0.4, color = "steelblue",
           label = "Mean: 80.481\nSd: 0.860") +
  labs(x = 'Mean Assortativity', y = 'Density') +
  mytheme

p_BC <- ggplot(random_node, aes(Betweenness_centrality)) +
  geom_histogram(bins = 50, aes(y=..density..), fill = brewer.pal(8,"Set2")[8],
                 color = brewer.pal(9,"Set1")[9]) +
  stat_function(aes(Betweenness_centrality),fun=dnorm, colour="steelblue", 
                xlim = c(min(random_node$Betweenness_centrality), 
                         max(random_node$Betweenness_centrality)),
                args = list(mean = mean(random_node$Betweenness_centrality), 
                            sd = sd(random_node$Betweenness_centrality))) +
  geom_point(aes(x = mean(as.numeric(as.character(IRG$Betweenness_centrality))), y=0), 
             color = 'red') +
  geom_segment(aes(x = mean(as.numeric(as.character(IRG$Betweenness_centrality))), y=4000,
                   xend = mean(as.numeric(as.character(IRG$Betweenness_centrality))), yend=300),
               arrow = arrow(length = unit(0.3, "cm")), color = 'red') +
  scale_y_continuous(expand = c(0.005,0.005)) +
  annotate("text", x = 0.00045, y = 5200, color = "red",
           label = "Observed: 4.987e-04\nZ-Score: 17.427\nP-Value < 2.2e-16") +
  annotate("text", x = 0.0003, y = 22500, color = "steelblue",
           label = "Mean: 2.289e-04\nSd: 1.549e-05") +
  labs(x = 'Mean Betweenness Centrality', y = 'Density') +
  mytheme

p_x <- ggplot(random_node, aes(Eigenvector_centrality)) +
  geom_histogram(bins = 60, aes(y=..density..), fill = brewer.pal(8,"Set2")[8],
                 color = brewer.pal(9,"Set1")[9]) +
  stat_function(aes(Eigenvector_centrality),fun=dnorm, colour="steelblue", 
                xlim = c(min(random_node$Eigenvector_centrality), 
                         max(random_node$Eigenvector_centrality)),
                args = list(mean = mean(random_node$Eigenvector_centrality), 
                            sd = sd(random_node$Eigenvector_centrality))) +
  geom_point(aes(x = mean(as.numeric(as.character(IRG$Eigenvector_centrality))), y=0), 
             color = 'red') +
  geom_segment(aes(x = mean(as.numeric(as.character(IRG$Eigenvector_centrality))), y=800,
                   xend = mean(as.numeric(as.character(IRG$Eigenvector_centrality))), yend=50),
               arrow = arrow(length = unit(0.3, "cm")), color = 'red') +
  scale_y_continuous(expand = c(0.005,0.005)) +
  annotate("text", x = 0.0065, y = 1000, color = "red",
           label = "Observed: 7.118e-03\nZ-Score: 26.677\nP-Value < 2.2e-16") +
  annotate("text", x = 0.005, y = 3000, color = "steelblue",
           label = "Mean: 4.066e-03\nSd: 1.144e-04") +
  labs(x = 'Mean Eigenvector Centrality', y = 'Density') +
  mytheme

p_SP <- ggplot(random_node, aes(Closeness_centrality)) +
  geom_histogram(bins = 50, aes(y=..density..), fill = brewer.pal(8,"Set2")[8],
                 color = brewer.pal(9,"Set1")[9]) +
  stat_function(aes(Closeness_centrality),fun=dnorm, colour="steelblue", 
                xlim = c(min(random_node$Closeness_centrality), 
                         max(random_node$Closeness_centrality)),
                args = list(mean = mean(random_node$Closeness_centrality), 
                            sd = sd(random_node$Closeness_centrality))) +
  geom_point(aes(x = mean(as.numeric(as.character(IRG$Closeness_centrality))), y=0), 
             color = 'red') +
  geom_segment(aes(x = mean(as.numeric(as.character(IRG$Closeness_centrality))), y=150,
                   xend = mean(as.numeric(as.character(IRG$Closeness_centrality))), yend=10),
               arrow = arrow(length = unit(0.3, "cm")), color = 'red') +
  scale_y_continuous(expand = c(0.005,0.005)) +
  annotate("text", x = 0.285, y = 200, color = "red",
           label = "Observed: 0.287\nZ-Score: 20.486\nP-Value < 2.2e-16") +
  annotate("text", x = 0.277, y = 600, color = "steelblue",
           label = "Mean: 0.274\nSd: 5.981e-04") +
  labs(x = 'Mean Closeness Centrality', y = 'Density') +
  mytheme

# The neighborhood enrichment scores
safe_node_property <- read.table("../intermediate/SAFE/human_gobp-node_properties_annotation-highest.txt",
                                 header = TRUE, sep = "\t", quote = "")
domain_num <- safe_node_property %>%
  dplyr::select(Node.label, Domain..predominant., Neighborhood.score..max.1..min.0...predominant., 
                Total.number.of.enriched.domains) %>%
  dplyr::rename(Gene = Node.label,
                Predominant_domain = Domain..predominant.,
                Neighborhood_score = Neighborhood.score..max.1..min.0...predominant.,
                Domain_num = Total.number.of.enriched.domains)

safe_neighbor_score <- read.table("../intermediate/SAFE/human_gobp-neighborhood_scores_annotation-highest.txt",
                                  header = TRUE, sep = "\t", quote = "", skip = 3)
safe_neighbor_score$Function_num <- apply(safe_neighbor_score[,-1], 1, 
                                          function(x) sum(x >= 0.323))

safe_domain_fun <- dplyr::select(safe_neighbor_score, name, Function_num) %>%
  dplyr::rename(Gene = name) %>%
  inner_join(domain_num[,c(1,4)], by = "Gene")

IRG2 <- safe_domain_fun %>%
  left_join(node_centrality[,1:2], by = "Gene") %>%
  filter(Gene_class == 'IRG&VTG' | Gene_class == 'IRG')

random_node2 <- data.frame(Number = 1:10000)
for (attr in colnames(safe_domain_fun)[2:3]) {
  var <- assign(attr, c())
  
  for (i in 1:10000) {
    set.seed(i)
    node_samp <- safe_domain_fun[sample(1:nrow(safe_domain_fun), 
                                        size = nrow(IRG2), replace = FALSE),]
    var <- c(var, mean(as.numeric(as.character(node_samp[attr][,1]))))
  }
  
  random_node2[attr] <- var
}

annotation2 <- data.frame(Vars = c('Observed', 'Z-Score', 'P-Value', 'Mean', 'Sd'))
for (attr in colnames(IRG2)[2:3]) {
  Obs <- mean(as.numeric(as.character(IRG2[attr][,1])))
  Mean <- mean(random_node2[attr][,1])
  Sd <- sd(random_node2[attr][,1])
  Zsc <- (Obs - Mean)/Sd
  #Alt <- ifelse(Zsc > 0, 'greater', 'less')
  #print(Alt)
  #Test <- t.test(random_node[attr][,1], mu = Obs)
  #Pval <- Test$p.value
  Pval <- 2*pnorm(-abs(Zsc))
  var <- assign(attr, c(Obs, Zsc, Pval, Mean, Sd))
  
  annotation2[attr] <- var
}

p_fun <- ggplot(random_node2, aes(Function_num)) +
  geom_histogram(bins = 80, aes(y=..density..), fill = brewer.pal(8,"Set2")[8],
                 color = brewer.pal(9,"Set1")[9]) +
  stat_function(aes(Function_num),fun=dnorm, colour="steelblue", 
                xlim = c(min(random_node2$Function_num), max(random_node2$Function_num)),
                args = list(mean = mean(random_node2$Function_num), 
                            sd = sd(random_node2$Function_num))) +
  geom_point(aes(x = mean(as.numeric(as.character(IRG2$Function_num))), y=0), 
             color = 'red') +
  geom_segment(aes(x = mean(as.numeric(as.character(IRG2$Function_num))), y=0.2,
                   xend = mean(as.numeric(as.character(IRG2$Function_num))), yend=0.02),
               arrow = arrow(length = unit(0.3, "cm")), color = 'red') +
  scale_y_continuous(expand = c(0.005,0.005)) +
  annotate("text", x = 28, y = 0.3, color = "red",
           label = "Observed: 30.861\nZ-Score: 24.045\nP-Value = 9.459e-128") +
  annotate("text", x = max(random_node2$Function_num)+0.01, y = 0.5, color = "steelblue",
           label = "Mean: 14.957\nSd: 0.661") +
  labs(x = 'Function Number', y = 'Density') +
  mytheme

p_2bcdef <- ggpubr::ggarrange(p_k, p_x, p_C, p_NC, p_fun,
                  nrow=2, ncol = 3, labels = c("B", "C", "D", "E", "F"))
ggsave(p_2bcdef,
       filename = "../result/Fig_2BCDEF.pdf",
       width = 14.2,
       height = 8.2)
