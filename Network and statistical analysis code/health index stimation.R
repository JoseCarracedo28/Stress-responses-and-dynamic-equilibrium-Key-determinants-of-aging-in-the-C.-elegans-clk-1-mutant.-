library(WRS2)
library(ggpubr)
library(ggplot2)
library(afex)
library(dplyr)
library(tidyverse)
library(rstatix)
library(scales)




############load the data#############
BombeoIS <- read.csv("Bombeo  - Anova.csv") #Data for pharingeal pumpim behavior
ColeteoIS <- read.csv("Coleteo - Anova.csv") #Data for swiming behavior
######################
###Pumping data####
######################
BombeoIS <- BombeoIS %>%
  gather(key = "time", value = "score", Dia.1, Dia.3, Dia.5, Dia.7, Dia.9) %>%
  convert_as_factor(Gusano, time)
BombeoIS$Cepa <- as.factor(BombeoIS$Cepa)
head(BombeoIS) #processing de data, selecting the time in days

##########Change data to a relative scale 
BombeoIS$score_normal <- rescale(BombeoIS$score)
BombeoIS


BombN2 <- subset(BombeoIS, Cepa == 'N2') #subset the pumping data by strain
write.csv(BombN2, "BombeoRS_N2.csv ", row.names=FALSE)
Bombclk1 <- subset(BombeoIS, Cepa == 'clk1')
write.csv(Bombclk1, "BombeoRS_CLK1.csv ", row.names=FALSE)
Bombaak2 <- subset(BombeoIS, Cepa == 'aak2')
write.csv(Bombaak2, "BombeoRS_AAK2.csv ", row.names=FALSE)
Bombdoble <- subset(BombeoIS, Cepa == 'doble')
write.csv(Bombdoble, "BombeoRS_DOBLE.csv ", row.names=FALSE)

#######################
###Swiming data########
#######################
ColeteoIS <- ColeteoIS %>%
  gather(key = "time", value = "score", Dia.1, Dia.3, Dia.5, Dia.7, Dia.9) %>%
  convert_as_factor(Gusano, time)
ColeteoIS$Cepa <- as.factor(ColeteoIS$Cepa)
head(ColeteoIS) #processing de data, selecting the time in days

##########Change data to a relative scale 
ColeteoIS$score_normal <- rescale(ColeteoIS$score)
ColeteoIS


colN2 <- subset(ColeteoIS, Cepa == 'N2')#subset the swiming data by strain
write.csv(colN2, "ColeteoRS_N2.csv ", row.names=FALSE)
colclk1 <- subset(ColeteoIS, Cepa == 'clk1')
write.csv(colclk1, "ColeteoRS_CLK1.csv ", row.names=FALSE)
colaak2 <- subset(ColeteoIS, Cepa == 'aak2')
write.csv(colaak2, "ColeteoRS_AAK2.csv ", row.names=FALSE)
coldoble <- subset(ColeteoIS, Cepa == 'doble')
write.csv(coldoble, "ColeteoRS_DOBLE.csv ", row.names=FALSE)




####################################
#########combine data###########
####################################
##Pumping & swiming##
suma <- BombeoIS$score_normal + ColeteoIS$score_normal
ByC <- data.frame(BombeoIS$Gusano, BombeoIS$Cepa, BombeoIS$time, suma)
colnames(ByC)[1]  <- "Gusano" #worm
colnames(ByC)[2]  <- "Cepa" #strain
colnames(ByC)[3]  <- "time" #time
colnames(ByC)[4]  <- "score" #score (value)
ByC

ByCN2 <- subset(ByC, Cepa == 'N2')
write.csv(ByCN2, "ByCRS_N2.csv ", row.names=FALSE)
ByCclk1 <- subset(ByC, Cepa == 'clk1')
write.csv(ByCclk1, "ByCRS_CLK1.csv ", row.names=FALSE)
ByCaak2 <- subset(ByC, Cepa == 'aak2')
write.csv(ByCaak2, "ByCRS_AAK2.csv ", row.names=FALSE)
ByCdoble <- subset(ByC, Cepa == 'doble')
write.csv(ByCdoble, "ByCRS_DOBLE.csv ", row.names=FALSE)

##################################################
###############Pumpimg analysis#################
##################################################

#########################
####pumping raw data####
#N2#
integBN2 <- function(x) {107.5 + sin(x) - 4.312*x} #Set the function that fits the data
integrate(integBN2, lower = 1, upper = 9) #get the definite integral to define the area
#CLK-1#
integBclk1 <- function(x) {95.59 - 10.11*x - 33.9*cos(0.499 + x)}
integrate(integBclk1, lower = 1, upper = 9)
#AAK-2#
integBAAK2 <- function(x) {100.9 - 3.014*x*cos(3.551*x)}
integrate(integBAAK2, lower = 1, upper = 9)
#DOuble mutant#
integBdoble <- function(x) {41.81 - 20.9*cos(1.985 - 1237*x)}
integrate(integBdoble, lower = 1, upper = 9,  subdivisions=2000)

Indice_completo_bombeo <- c(688.9714, 396.6465, 802.4906, 334.4834)

names(Indice_completo_bombeo) <- c("N2", "clk-1", "aak-2", "clk-1;aak-2")
Indice_completo_bombeo

Indice_completo_bombeo_nor <- rescale(Indice_completo_bombeo, to = c(0, 100))
Indice_completo_bombeo_nor

#################################
####Pumpimg in relative scale####
#N2#
nor_integBN2 <- function(x) {0.7226 + 0.0002856*x^2 - 0.03398*x}
integrate(nor_integBN2, lower = 1, upper = 9)
#CLK-1#
nor_integBclk1 <- function(x) {(0.9336 - cos(0.1651 + x))/x}
integrate(nor_integBclk1, lower = 1, upper = 9)
#AAK-2#
nor_integBAAK2 <- function(x) {0.6006 - 0.0998*cos(2.16 + 6.771*x)}
integrate(nor_integBAAK2, lower = 1, upper = 9)
#double mutant#
nor_integBdoble <- function(x) {0.2473 + 0.148*sin(0.4097 + 0.4767*x)}
integrate(nor_integBdoble, lower = 1, upper = 9)

Indice_normal_bombeo <- c(4.490906, 2.447711, 4.807895, 2.1785)

names(Indice_normal_bombeo) <- c("N2", "clk-1", "aak-2", "clk-1;aak-2")
Indice_normal_bombeo



##################################################
##############Swiming Analysis#################
##################################################
#######################
####Swiming raw data####
#N2#
integCN2 <- function(x) {111.8 + 10.19*sin(0.7253 - x)}
integrate(integCN2, lower = 1, upper = 9)
#CLK-1#
integCclk1 <- function(x) {91.16 - x*sin(4.576*x) - 4.132*x}
integrate(integCclk1, lower = 1, upper = 9)
#AAK-2#
integCAAK2 <- function(x) {91.7-27.7*sin(4.81 + 0.277*x)}
integrate(integCAAK2, lower = 1, upper = 9)
#DOuble muntan#
integCdoble <- function(x) {37.58 - 19.65*cos(71.9*x)}
integrate(integCdoble, lower = 1, upper = 9,  subdivisions=2000)

Indice_completo_coleteo <- c(880.4303, 562.1465, 749.3613, 300.754)

names(Indice_completo_coleteo) <- c("N2", "clk-1", "aak-2", "clk-1;aak-2")
Indice_completo_coleteo

Indice_completo_coleteo_nor <- rescale(Indice_completo_coleteo, to = c(0, 100))
Indice_completo_coleteo_nor

#######################################
####Swiming data in a relative scale###
#N2#
nor_integCN2 <- function(x) {0.8081 + 0.07782*cos(8418 + x)}
integrate(nor_integCN2, lower = 1, upper = 9)
#CLK-1#
nor_integCclk1 <- function(x) {0.6554 + cos(1.567*x) - 0.0337*x}
integrate(nor_integCclk1, lower = 1, upper = 9)
#AAK-2#
nor_integCAAK2 <- function(x) {0.9028 - 0.04774*x}
integrate(nor_integCAAK2, lower = 1, upper = 9)
#DOuble mutant#
nor_integCdoble <- function(x) {9.561/(23.09 + x^2)}
integrate(nor_integCdoble, lower = 1, upper = 9)

Indice_normal_coleteo <- c(6.57374, 3.894832, 5.1328, 1.74142)

names(Indice_normal_coleteo) <- c("N2", "clk-1", "aak-2", "clk-1;aak-2")
Indice_normal_coleteo

###########################################################
##############Pumping and swiming analysis#################
###########################################################
#N2#
integCYBN2 <- function(x) {1.555 - 0.03881*x - 0.0545*sin(x)}
integrate(integCYBN2, lower = 1, upper = 9)
#CLK-1#
integCYBclk1 <- function(x) {1.276 - 0.09899*x - 0.2037*cos(0.5456 + x)}
integrate(integCYBclk1, lower = 1, upper = 9)
#AAK-2#
integCYBAAK2 <- function(x) {1.25 + 0.3095*cos(-31.13*x)}
integrate(integCYBAAK2, lower = 1, upper = 9)
#DOuble mutant#
integCYBdoble <- function(x) {0.4567 + 0.2968*cos(0.5064 - 0.3841*x)}
integrate(integCYBdoble, lower = 1, upper = 9,  subdivisions=2000)

Indice_completo_CYB <- c(10.8085, 6.476587, 9.997454, 3.894631)

names(Indice_completo_CYB) <- c("N2", "clk-1", "aak-2", "clk-1;aak-2")
Indice_completo_CYB


