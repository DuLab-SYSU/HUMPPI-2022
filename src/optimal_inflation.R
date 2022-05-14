library(ggplot2)
library(RColorBrewer)
mytheme <- theme(panel.grid.major = element_line(colour=brewer.pal(12,"Set3")[9],
                                                 linetype = "longdash"),
                 panel.background = element_rect(fill='transparent', color="#000000"),
                 panel.border=element_rect(fill='transparent', color='black'),
                 axis.text = element_text(size = 16),
                 axis.title = element_text(size = 20),
                 legend.text = element_text(size = 13),
                 legend.background = element_blank(),
                 plot.margin = unit(x=c(6,6,6,6),units = 'mm'))


## Functional similarity
funSim <- read.table("../intermediate/module_funSim.txt", 
                     header = FALSE, sep = '\t')
colnames(funSim) <- c("Inflation", "Module", "Fun_Sim")

funSim2 <- Rmisc::summarySE(funSim, measurevar="Fun_Sim",
                            groupvars = c("Inflation"))

p_sim <- ggplot(funSim2, aes(Inflation, Fun_Sim)) +
  geom_point(size = .5, color = "steelblue") +
  geom_line(group = 1, color = "steelblue") +
  geom_errorbar(aes(ymin=Fun_Sim-se,ymax=Fun_Sim+se),
                width=.01, size = .05, color = "blue") +
  geom_vline(xintercept = funSim2[funSim2$Fun_Sim == max(funSim2$Fun_Sim),1], 
             linetype = 5, color = "red") +
  annotate("text", x = 2.7, y = 0.024, color = "red",
           label = "Inflation = 2.28") +
  scale_x_continuous(expand = c(0.01,0.01), breaks = seq(1.5, 5, 0.5)) +
  scale_y_continuous(breaks = seq(0.02, 0.04, 0.002)) +
  labs(x = "Inflation parameter", y = "Mean similarity value") +
  mytheme

## Module size distribution
modules <- read.table("../intermediate/batch_MCL_out/out.HUMPPI2022.I228", 
                       header = FALSE, sep = " ")
modules$Size <- apply(modules[,1,drop = FALSE], 1,
                      function(x) length(stringr::str_split(x, "\t")[[1]]))
modules <- dplyr::rename(modules, Member = V1)
p_size <- ggplot(modules, aes(Size)) +
  geom_histogram(bins = 75, fill = "steelblue", alpha = .85) +
  scale_x_continuous(breaks = seq(0, 75, 15), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = "Modules size", y = "modules number") +
  mytheme

library(cowplot)
ggdraw() +
  draw_plot(p_sim, x=0, y=0, width = .55, height = 1) +
  draw_plot(p_size, x=0.55, y=0, width = .45, height = 1) +
  draw_plot_label(label = c("A", "B"),
                  x = c(0, 0.55), y = c(1, 1), size = 20)

modules$Size_region <- apply(modules[,2,drop = FALSE], 1,
                              function(x) ifelse(x == 1, 1, ifelse(
                                x == 2, 2, ">=3")))

module_size <- data.frame(t(table(modules$Size_region)))
module_size$label_value <- paste0(module_size$Freq, '(', 
                                   round(module_size$Freq/sum(module_size$Freq) * 100, 1), '%)')
label_value <- paste0(module_size$Freq, ' (', 
                      round(module_size$Freq/sum(module_size$Freq) * 100, 1), '%)')
p_pie <- ggplot(module_size, aes("", Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = '', y = '', title = '', fill = "Module size") +
  coord_polar(theta = 'y') +
  geom_text(y = module_size$Freq/2 + 
              c(0, cumsum(module_size$Freq)[-length(module_size$Freq)]),
            x = sum(module_size$Freq)/2300, label = label_value,
            color = "white") +
  scale_fill_manual(limits = c(1, 2, ">=3"),
                    values = c('#45BF86', '#0378A6', '#D97652')) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  theme_void()


library(cowplot)
p_s2 <- ggdraw() +
  draw_plot(p_sim, x=0, y=0, width = .5, height = 1) +
  draw_plot(p_size, x=0.5, y=0, width = .5, height = 1) +
  draw_plot(p_pie, x=0.65, y=0.45, width = .3, height = .6) +
  draw_plot_label(label = c("A", "B"),
                  x = c(0, 0.5), y = c(1, 1), size = 20)

ggsave(p_s2,
       filename = "../result/Fig_S2.pdf",
       width = 15,
       height = 6)
