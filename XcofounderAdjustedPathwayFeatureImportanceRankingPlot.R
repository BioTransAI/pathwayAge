library(ggplot2)
library(cowplot)
library(Metrics)
library(dplyr)
library(ggsignif)
library(ggtext)



data <- read.csv("./discovery3KGOSubRankWithGeneCount.csv")
custom_order <- data$Description
textSize <- 12
## Create a scatter plot

p <- ggplot(data, aes(x= RhoAbs, y= Description, size = GeneCount, color = Tag)) +
  geom_point(alpha=0.9) +  # Set point color based on y values
  scale_color_manual(values = c("developmental process" = "#7FB3D5", "metabolic process" = "#BC8B56", "response to stimulus"= "#F6E758", "cellular process" ="#6D936E",
                                    "biological regulation"= "#949b9b", "localization"= "#8c7e78", "multicellular organismal process" = "#3c4444")) +  # Customize y-axis breaks
  ylab("Biological Pathways")+
  scale_y_discrete(limits = custom_order)+
  theme(
    axis.text.y = element_markdown(color = data$color,size = textSize, face="bold"),
    axis.text.x = element_text(size = textSize, face="bold"),
    legend.text=element_text(size=16),
    legend.title = element_blank(),
    axis.title=element_text(size=16,face="bold")  
    )  # Match y-axis text color to point color
getwd()
# print(plot_grid (p, labels = c("")), label_size = 25, ncol=1)
ggsave("./control3000GOSubRankWithGeneCount.png", width = 35, height = 20, units = "cm")
