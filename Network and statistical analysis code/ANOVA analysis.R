#############################################
########reed and accommodate the data########
#############################################
library(WRS2)
library(ggpubr)
library(ggplot2)
library(afex)
library(dplyr)
library(tidyverse)
library(rstatix)
#Data for all assays
#run the specific data each time to run the analysis for the selected behavior 
AnovBM <- read.csv("Bombeo - Anova.csv")#pharyngeal pumping behavior 
AnovBM <- read.csv("Coleteo - Anova.csv") #swiming behavior
AnovBM <- ByC



AnovBM <- AnovBM %>%
  gather(key = "time", value = "score", Dia.1, Dia.3, Dia.5, Dia.7, Dia.9) %>%
  convert_as_factor(Gusano, time)
AnovBM$Cepa <- as.factor(AnovBM$Cepa)
head(AnovBM)


##################################
#########Model ASSUMPTIONS########
##################################
AnovBM%>%
  group_by(Cepa, time) %>%
  get_summary_stats(score, type = "mean_sd")

##visualization
bxp<- ggboxplot(AnovBM, "time", "score", color = "Cepa", palette= "jco", notch = TRUE,  add = "jitter")

bxp

# outliers 
AnovBM%>%group_by(time, Cepa) %>%
  identify_outliers(score)

##residuals normal distribution
AnovBM%>%
  group_by(time, Cepa) %>%
  shapiro_test(score)

ggqqplot(AnovBM, "score", ggtheme = theme_bw()) +facet_grid (time ~ Cepa)

#HOMogeneity of variace
AnovBM%>%  group_by(time)%>% levene_test(score~ Cepa)

### COVARIANCES
box_m(AnovBM[, "score", drop = FALSE],AnovBM$Cepa)

##### Three-way mixed ANOVA: 2 between- and 1 within-subjects factors
#two way 
AnovBM %>% sample_n_by(Cepa, size = 1)

res.aov <- anova_test(
  data = AnovBM, dv = score, wid = Gusano,
  within = time, between = c(Cepa)
)
res.aov
get_anova_table(res.aov)


########Post-hoc test#######
######Effect of all strains in time#######
one.way <- AnovBM %>%
  group_by(time) %>%
  anova_test(dv = score, wid = Gusano, between = Cepa) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
#######Comparison between strains######
pwc <- AnovBM %>%
  group_by(time) %>%
  pairwise_t_test(score ~ Cepa, p.adjust.method = "bonferroni")
print(pwc, n = 40)

#### Time effect on each strain#####
one.way2 <- AnovBM %>%
  group_by(Cepa) %>%
  anova_test(dv = score, wid = Gusano, within = time) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

#####Diference between the day for each strain######
pwc2 <- AnovBM %>%
  group_by(Cepa) %>%
  pairwise_t_test(
    score ~ time, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) %>%
  select(-df, -statistic, -p) # Remove details
print(pwc2, n = 40)



####P values in boxplots####
pwc <- pwc %>% add_xy_position(x = "time")
bxp + 
  stat_pvalue_manual(pwc) +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))

###P VALUES in boxplots, comparision between days in the same strains###
pwc2 <- pwc2 %>% add_xy_position(x = "time")
bxp + 
  stat_pvalue_manual(pwc2) +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc2))

###############################################################
####################Non parametric Anova#######################
###############################################################
####Mixed Anova####
bwtrim(score ~ Cepa*time, id = Gusano, data = AnovBM)
#####trimming 0, analysis######
bwtrim(score ~ Cepa*time, id = Gusano, data = AnovBM, tr = 0)
###M huber stimator#####
a<-sppba(score ~ Cepa*time, id = Gusano, data = AnovBM)
b<-sppbb(score ~ Cepa*time, id = Gusano, data = AnovBM)
bi<-sppbi(score ~ Cepa*time, id = Gusano, data = AnovBM)
a
b
bi

