
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
library(knitr)
library(data.table)
library(psych)
library(EnvStats)
library(glmnet)

library(bestglm)
library(caret)
library(pROC)
library(glmnet)


##############################
#Data preparation

stats_new <- read_excel("stats_new.xlsx")

analize = stats_new[c('Lytis', 'Amzius', 'Atvykimo glu', 'Atv RITS glu',
                        'Dexametazonas', 'RITS valandos', 'Lovadieniai', 'Exitus RITS', 
                        'SAPS', 'SOFA', 'APACHEII', 'ISARIC4C','Poliligotumas',
                        'CD','LILGFG60', 'PaO2FiO2', 'DPV10', 'PIT', 'Vazopresoriai10')]

colnames(analize) <- c('Lytis', 'Amzius', 'Atvykimo_glu', 'Atv_RITS_glu',
                       'Dexametazonas', 'RITS_valandos', 'Lovadieniai', 'Exitus_RITS', 
                       'SAPS', 'SOFA', 'APACHEII', 'ISARIC4C','Poliligotumas',
                       'CD','LILGFG60', 'PaO2FiO2', 'DPV10', 'PIT', 'Vazopresoriai10')

analize$PaO2FiO2 <- as.numeric(analize$PaO2FiO2)

analize_nona <- na.omit(analize)

##############################
#Descriptive Stats

options(digits=3)
summary_df = stat.desc(analize_nona, basic=TRUE, desc=TRUE)
summary_df = round(summary_df, 2)
summary_reduced = summary_df[c(4,5,6,8,9,10,13),]
t_summary_reduced = transpose(summary_reduced)
colnames(t_summary_reduced) = rownames(summary_reduced)
rownames(t_summary_reduced) = colnames(summary_reduced)

kable(t_summary_reduced, caption = 'Summary Statistics')


pdf('summary_stats.pdf')
grid.table(t_summary_reduced)
dev.off()

scaled.dat <- scale(analize_nona)

boxplot(scaled.dat[,c('Amzius', 'Atvykimo_glu','Atv_RITS_glu','RITS_valandos','Lovadieniai',
                        'SAPS', 'SOFA', 'APACHEII', 'ISARIC4C', 'PaO2FiO2')])

boxplot(analize_nona[,c('Amzius', 'Atvykimo_glu','Atv_RITS_glu','RITS_valandos','Lovadieniai',
                     'SAPS', 'SOFA', 'APACHEII', 'ISARIC4C', 'PaO2FiO2')])
#Many outliers

##############################
#Check normal distribution

analize_nona_factors <- analize_nona
fac_colnames = c('Lytis', 'Dexametazonas', 'Exitus_RITS', 'Poliligotumas', 'CD', 'LILGFG60', 'DPV10', 'PIT','Vazopresoriai10')

analize_nona_factors[,fac_colnames] <- lapply(analize_nona_factors[,fac_colnames], factor)

nums <- unlist(lapply(analize_nona_factors, is.numeric), use.names = FALSE)

shap = lapply(analize_nona_factors[nums], shapiro.test)


##############################

#Check the outliers

pairs(analize_nona_factors[nums], cex.labels = 2)

outliers <- rosnerTest(analize_nona_factors$Atvykimo_glu, k=7)
outliers$all.stats

##############################
#Scatter Plots by CD

cd_colours = c('0' = 'blue', '1' = 'red')
plot(analize_nona_factors[nums], col = cd_colours[as.character(analize_nona$CD)],
     cex.labels=2)


#Netiklsinga
cor_matrix <- cor(analize_nona, method = 'spearman')

corrplot(cor_matrix, method='color')

#Tikslinga
cor_matrix <- cor(analize_nona_factors[nums], method = 'spearman')

corrplot(cor_matrix, method='color')

##############################
#Scatter Plots

pairs.panels(analize_nona_factors[nums], smooth=FALSE, method='spearman', ellipses=FALSE,
             lm=TRUE,cor=TRUE,jiggle=FALSE, pch=20,cex.cor=3, cex=2, stars=TRUE)

pairs(analize_nona_factors[nums], )


##############################
#Scatter Plots with correlation p-value

reshaped_df <- reshape2::melt(analize_nona, id.vars = "Atvykimo_glu")

plot <- ggplot(reshaped_df, aes(x = value, y = Atvykimo_glu), color = factor(CD)) +
  geom_point(alpha = 0.5) +
  ggtitle("Atvykimo gliukozė ~ Visi") +
  facet_wrap(~variable, scales = "free") +
  scale_color_manual(values = custom_colors) +
  stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "top") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

print(plot)

###################################
#Scatter Plots with correlation p-value and separated by CD

reshaped_df <- reshape2::melt(analize_nona, id.vars = c("Atvykimo_glu", "CD"))

custom_colors <- c("blue", "red")

plot <- ggplot(reshaped_df, aes(x = Atvykimo_glu, y = value, color = factor(CD))) +
  geom_point(alpha = 0.5) +
  ggtitle("Atvykimo_glu ~ Visi") +
  facet_wrap(~variable, scales = "free") +
  scale_color_manual(values = custom_colors) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "top") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

print(plot)

##############################
#Scatter Plot only cd=0 with correlation p-value

reshaped_df <- reshape2::melt(analize_nona[which(analize_nona$CD == 0),], id.vars = "Atvykimo_glu")

plot <- ggplot(reshaped_df, aes(x = Atvykimo_glu, y = value), col = 'blue') +
  geom_point(alpha = 0.5, color='blue') +
  ggtitle("Atvykimo_glu vs All") +
  facet_wrap(~variable, scales = "free") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "top") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

print(plot)

##############################
#Scatter Plot only cd=1 with correlation p-value

reshaped_df <- reshape2::melt(analize_nona[which(analize_nona$CD == 1),], id.vars = "Atvykimo_glu")

plot <- ggplot(reshaped_df, aes(x = Atvykimo_glu, y = value), col = 'blue') +
  geom_point(alpha = 0.5, color='red') +
  ggtitle("Atvykimo_glu vs All") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~variable, scales = "free") +
  stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "top") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

print(plot)

###################################
#Scatter Plots with correlation p-value and separated by Exitus where CD =0

reshaped_df <- reshape2::melt(analize_nona[which(analize_nona$CD == 0),], id.vars = c("Atvykimo_glu", "Exitus_RITS"))

custom_colors <- c("blue", "red")

plot <- ggplot(reshaped_df, aes(x = Atvykimo_glu, y = value, color = factor(Exitus_RITS))) +
  geom_point(alpha = 0.5) +
  ggtitle("Atvykimo_glu vs All") +
  facet_wrap(~variable, scales = "free") +
  scale_color_manual(values = custom_colors) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "top") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

print(plot)

###################################
#Cor.test
analize_nona_nocd = analize_nona[which(analize_nona$CD==0),][-'CD']
analize_nona_nocd = subset(analize_nona_nocd, select = -CD)
cor.test(analize_nona_nocd$Atvykimo_glu)

##############################
#Anova Test

analize_nona_nocd_factors = analize_nona_factors[which(analize_nona_factors$CD == 0),]
analize_nona_nocd_factors = subset(analize_nona_nocd_factors, select = -CD)


leveneTest(analize_nona_nocd_factors$Atvykimo_glu, analize_nona_nocd_factors$Lytis)
leveneTest(analize_nona_nocd_factors$Atvykimo_glu, analize_nona_nocd_factors$Dexametazonas)
leveneTest(analize_nona_nocd_factors$Atvykimo_glu, analize_nona_nocd_factors$Exitus_RITS)
leveneTest(analize_nona_nocd_factors$Atvykimo_glu, analize_nona_nocd_factors$Poliligotumas)
leveneTest(analize_nona_nocd_factors$Atvykimo_glu, analize_nona_nocd_factors$LILGFG60)
leveneTest(analize_nona_nocd_factors$Atvykimo_glu, analize_nona_nocd_factors$DPV10)
leveneTest(analize_nona_nocd_factors$Atvykimo_glu, analize_nona_nocd_factors$PIT)
leveneTest(analize_nona_nocd_factors$Atvykimo_glu, analize_nona_nocd_factors$Vazopresoriai10)

mod = aov(Atvykimo_glu ~ x)

###################################
#Exploring sugar levels when accounting for CD

bwplot(Atvykimo_glu ~ CD|Exitus_RITS, data=analize_nona_factors,aspect = 2,
         ylab="Atvykimo_glu", xlab="Exitus",
         main = "Boxplots of glu (left) for CD presence (top) by death (bottom)")

bwplot(Atvykimo_glu ~ Exitus_RITS, data=analize_nona_nocd_factors,aspect = 2,
       ylab="Atvykimo_glu", xlab="Exitus",
       main = "Boxplots of glu (left) for CD presence (top) by death (bottom)")

#ANOVA
ggboxplot(analize_nona_factors, x = "CD", y = "Atvykimo_glu", color = "Exitus_RITS",
          palette = c("#00AFBB", "#E7B800"))

aov_glu = aov(Atvykimo_glu ~ Exitus_RITS + CD, data=analize_nona_factors)
summary(aov_glu)
TukeyHSD(aov_glu, conf.level = 0.90)

#######################################
#t-tests

p_exitus = t.test(Atvykimo_glu ~ Exitus_RITS, data = analize_nona_nocd_factors, var.equal=FALSE)
p_dex = t.test(Atvykimo_glu ~ Dexametazonas, data = analize_nona_nocd_factors, var.equal=TRUE)
p_dpv = t.test(Atvykimo_glu ~ DPV10, data = analize_nona_nocd_factors, var.equal=FALSE)
p_vaz = t.test(Atvykimo_glu ~ Vazopresoriai10, data = analize_nona_nocd_factors, var.equal=FALSE)
p_pit = t.test(Atvykimo_glu ~ PIT, data = analize_nona_nocd_factors, var.equal=FALSE)

p_lova = cor.test(analize_nona_nocd_factors$Atvykimo_glu, analize_nona_nocd_factors$Lovadieniai, method = 'pearson')
p_isa = cor.test(analize_nona_nocd_factors$Atvykimo_glu, analize_nona_nocd_factors$ISARIC4C, method = 'pearson')
p_pa = cor.test(analize_nona_nocd_factors$Atvykimo_glu, analize_nona_nocd_factors$PaO2FiO2, method = 'pearson')

p_adjusted = p.adjust(c(p_exitus$p.value, p_dex$p.value, p_dpv$p.value, p_vaz$p.value, p_pit$p.value, p_lova$p.value, p_isa$p.value, p_pa$p.value))

p_adjusted
#########################################
#Tukey
TukeyHSD(aov, conf.level = 0.90)

#######################################
#other linear models

samp <- sample(1:nrow(analize_nona_nocd_factors), 0.8*164)
train = analize_nona_nocd_factors[samp,]
test = analize_nona_nocd_factors[-samp,]
x_train <- data.matrix(train[-8])
x_test <- data.matrix(test[-8])
cv_lasso <- cv.glmnet(x_train, y = train$Exitus_RITS, alpha = 1, family='binomial')
best_lambda <- cv_lasso$lambda.min
plot(cv_lasso)
lasso_best <- glmnet(x_train, y = as.factor(train$Exitus_RITS), alpha = 1, family = 'binomial', lambda = best_lambda)
coef(lasso_best)

lasso_pred <- predict(lasso_best, newx = x_test, type = "response")
optimal_lasso <- round(lasso_pred)
confusionMatrix(as.factor(test$Exitus_RITS), as.factor(optimal_lasso[,1]))
                
                
log_exit = glm(Exitus_RITS ~ . - Atv_RITS_glu, data = analize_nona, family=binomial(link='logit'))

summary(log_exit)

best.logit_AIC <- bestglm(train, IC = "AIC", family = binomial, method = "exhaustive")
##############################
#Looking for Multicollinearity

vif(linearm)

##############################
#Scatter plots for binomial data
binomial_scatter <- function(data, cont, bin) {
  scatter_plot <- ggplot(data, aes(x = {{cont}}, y = {{bin}})) +
    geom_point() +
    stat_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, color = "blue") +
    labs(title = "Scatter Plot with Logistic Regression Line", x = "Continuous Variable", y = "Binomial Variable")
  print(scatter_plot)
}

binomial_scatter(analize_nona, Atvykimo_glu, Poliligotumas)



##############################
#boxplots

boxplot_plot <- ggplot(analize_nona, aes(x = factor(Poliligotumas), y = Atvykimo_glu)) +
  geom_boxplot() +
  labs(title = "Boxplot of Continuous Variable by Category", x = "Category", y = "Continuous Variable")

# Display the plot
print(boxplot_plot)



#



