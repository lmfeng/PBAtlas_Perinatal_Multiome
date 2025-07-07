samples <- read.table("/PBA/original_data/1129_sample_rnaandprotein_175.txt", row.names = 1, header=TRUE, sep = '\t')
samples$Days[samples$Age == "E76"] <- 76
samples$Days[samples$Age == "E85"] <- 85
samples$Days[samples$Age == "E94"] <- 94
samples$Days[samples$Age == "E104"] <- 104
samples$Days[samples$Age == "E109"] <- 109
samples$Days[samples$Age == "P0"] <- 115
samples$Days[samples$Age == "P3"] <- 118
samples$Days[samples$Age == "P30"] <- 145
library(ggplot2)
library(dplyr)
samples <- samples %>%  mutate(Ratio = Brain_weight / Samples_weight)
age_levels <- c("E76", "E85", "E94", "E104", "E109", "P0", "P3", "P30")
samples$Age_numeric <- as.numeric(factor(samples$Age, levels = age_levels, ordered = TRUE))
table(samples$Age_numeric)
str(samples$Age_numeric)

pdf("/PBA/MF1_B.curve.pdf",height = 6,width = 6)
p <- ggplot(data = samples, aes(x = Days, y = Brain_weight, color = Age)) +
  geom_point(aes(shape = Sex), size = 4) +
  scale_color_manual(values = c('#f3993a','#ffee6f','#add5a2','#7a7b78','#9933cc','#0066ff','#33cccc','#ff66cc'),limits = age_levels) +
  stat_smooth(method = 'loess',color = "#83A0be",size = 1,fill = "#83A0be",span = 0.65) +
  scale_x_continuous(breaks = samples$Days, labels = samples$Age) +
  theme_minimal() +
  labs(x = "Age (days)",y ="Brain weight(g)",title = "Ratio vs Age",color = "Age",shape = "Sex") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.6), 
        axis.ticks = element_line(colour = "black", size = 0.6), 
        axis.text = element_text(colour = "black", size = 10)) + 
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2,override.aes = list(size = 3)))
print(p)


p <- ggplot(data = samples, aes(x = Days, y = Ratio, color = Age)) +
  geom_point(aes(shape = Sex), size = 4) +
  scale_color_manual(values = c('#f3993a','#ffee6f','#add5a2','#7a7b78','#9933cc','#0066ff','#33cccc','#ff66cc'), 
                     limits = age_levels) +
  stat_smooth(method = 'loess', color = "#83A0be", size = 0.6, formula = 'y ~ x', fill = "#83A0be",span = 0.65) + 
  scale_x_continuous(breaks = samples$Days, labels = samples$Age) +
  theme_minimal() +
  labs(x = "Age (days)", y = "Brain weight/Body weight", 
       title = "", 
       color = "Age", 
       shape = "Sex") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", size = 0.6),
        axis.ticks = element_line(colour = "black", size = 0.6),
        axis.text = element_text(colour = "black", size = 10)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2,override.aes = list(size = 3)))
print(p)

p <- ggplot(data = samples, aes(x = Days, y = Samples_weight, color = Age)) +
  geom_point(aes(shape = Sex), size = 4) +
  scale_color_manual(values = c('#f3993a','#ffee6f','#add5a2','#7a7b78','#9933cc','#0066ff','#33cccc','#ff66cc'),
                     limits = age_levels) +
  stat_smooth(method = 'loess', color = "#83A0be", size = 1, fill = "#83A0be", span = 0.65) +
  scale_x_continuous(breaks = samples$Days, labels = samples$Age) +
  theme_minimal() +
  labs(x = "Age (days)", y = "Body weight(kg)",
       title = "Brain weight", 
       color = "Age", 
       shape = "Sex") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", size = 0.6),
        axis.ticks = element_line(colour = "black", size = 0.6),
        axis.text = element_text(colour = "black", size = 10)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, override.aes = list(size = 3)))

print(p)
dev.off()

