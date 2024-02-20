#public data analysis

library(readxl)
library(corrplot)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(pastecs)
library(ggpubr)
library(car)
library(dplyr)
library(lattice)
library(grid)



raw_data <- read_excel("public data/mortality_risk_clinical/Mortality_incidence_sociodemographic_and_clinical_data_in_COVID19_patients.xlsx")

data_no_ill <- subset(raw_data, !apply(raw_data[, c(11, 12, 13, 14, 15, 16, 19, 20, 21, 22, 23, 24, 25, 26)] == 1, 1, any))

data_no_ill <- subset(raw_data, !apply(raw_data[, c(11, 12, 13, 14, 15, 16, 19, 20, 21, 22, 23, 24, 25, 26)] == 1, 1, any))

data_no_dia <- subset(data_no_ill, !apply(data_no_ill[, c(17, 18)] == 1, 1, any))

sum(data_no_dia$`DM Complicated`)

data_with_glu <- subset(data_no_dia, Glucose != 0)

ggplot(data_with_glu, aes(x = factor(Death), y = Glucose)) +
  geom_boxplot() +
  labs(title = "Boxplot of Continuous Variable by Category", x = "Category", y = "Continuous Variable")


t.test(subset(data_with_glu, Death = 0)$Glucose, subset(data_with_glu, Death = 1)$Glucose)
