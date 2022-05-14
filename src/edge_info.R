### Edge information (2021-09-21)
# HUMPPI-2022
humppi2022 <- read.table('../intermediate/integrated_ppi.txt',
                      sep = '\t', header = FALSE)
colnames(humppi2022) <- c("Gene_A", "Gene_B")

node_info <- read.table('../intermediate/nodes_info.txt', 
                        header = TRUE, sep = '\t')
gene_class <- node_info[,1:2]

# Cluster of each gene in
clusters <- read.table('../result/Module_members.txt',
                       header = TRUE, sep = '\t')
clusters <- clusters[1:9004,]

# Gene age
age <- read.table('../data/human_name_age.txt',
                  header = TRUE, sep = '\t')

# Combine all the gene information 
humppi2022 <- merge(humppi2022, clusters, by.x = 'Gene_A', by.y = 'Gene', all.x = TRUE)
colnames(humppi2022) <- c('Gene_A', 'Gene_B', 'Module_A')
humppi2022 <- merge(humppi2022, clusters, by.x = 'Gene_B', by.y = 'Gene', all.x = TRUE)
humppi2022 <- dplyr::rename(humppi2022, Module_B=Module)
humppi2022 <- merge(humppi2022, age, by.x = 'Gene_A', by.y = 'Gene', all.x = TRUE)
humppi2022 <- dplyr::rename(humppi2022, Branch_A=Branch)
humppi2022 <- merge(humppi2022, age, by.x = 'Gene_B', by.y = 'Gene', all.x = TRUE)
humppi2022 <- dplyr::rename(humppi2022, Branch_B=Branch)
humppi2022 <- merge(humppi2022, gene_class, by.x = 'Gene_A', by.y = 'Gene', all.x = TRUE)
humppi2022 <- dplyr::rename(humppi2022, Class_A=Gene_class)
humppi2022 <- merge(humppi2022, gene_class, by.x = 'Gene_B', by.y = 'Gene', all.x = TRUE)
humppi2022 <- dplyr::rename(humppi2022, Class_B=Gene_class)
humppi2022 <- humppi2022[,c(2,1,3:8)]
humppi2022[is.na(humppi2022)] <- "-"

# Edge age and category
myFun <- function(humppi2022, c1, c2) {
  ifelse(humppi2022[c1]=='-' | humppi2022[c2]=='-', '-', 
         ifelse(humppi2022[c1]==humppi2022[c2], 'Intra-', 'Inter-'))
}

humppi2022$Category <- apply(humppi2022, 1, myFun, c1 = 'Module_A', c2 = 'Module_B')

myFun2 <- function(humppi2022, c1, c2) {
  ifelse(humppi2022[c1]=='-' | humppi2022[c2]=='-', '-', 
         max(as.numeric(as.character(humppi2022[c1])), as.numeric(as.character(humppi2022[c2]))))
}

humppi2022$Edge_age <- apply(humppi2022, 1, myFun2, c1 = 'Branch_A', c2 = 'Branch_B')

Edge_class <- c()
for (row in 1:nrow(humppi2022)) {
  if ((humppi2022[row,7] == "IRG&VTG" | humppi2022[row,8] == "IRG&VTG") |
      (humppi2022[row,7] == "IRG" & humppi2022[row,8] == "VTG") |
      (humppi2022[row,7] == "VTG" & humppi2022[row,8] == "IRG")) {
    Edge_class <- c(Edge_class, "IRE&VTE")
  }
  else if (humppi2022[row,7] == "IRG" | humppi2022[row,8] == "IRG") {
    Edge_class <- c(Edge_class, "IRE")
  }
  else if (humppi2022[row,7] == "VTG" | humppi2022[row,8] == "VTG") {
    Edge_class <- c(Edge_class, "VTE")
  } 
  else {
    Edge_class <- c(Edge_class, "Other")
  }
}

humppi2022$Edge_class <- Edge_class

write.table(humppi2022,
            '../intermediate/edges_info.txt',
            quote = FALSE, sep = "\t", row.names = FALSE)
