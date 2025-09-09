#Boolean network configuration 2 analysis
#Load the required library
library(BoolNet)
library(BoolNetPerturb)
#Set the Boolean rules for the second configuration 
redclk12 <- tempfile(pattern = "testNet")
sink(redclk12)
cat("targets, factors\n")
cat("clk1, 0\n")
cat("ETC, clk1\n")
cat("pink1, !ETC\n")
cat("skn1, ros | !ETC\n")
cat("taf4, !clk1\n")
cat("met, !hif1\n")
cat("hlh11, !atfs1 & !mtor\n")
cat("unc51, ampk & !mtor\n")
cat("hlh30, !mtor\n")
cat("creb, taf4 & !crtc1\n")
cat("pmp, !clk1 & (ETC2 | !hsp60)\n")
cat("ATG, unc51 & hlh30\n")
cat("lipl4, hlh30\n")
cat("MTG, pink1 & skn1 & unc51\n")
cat("fzo1, creb\n")
cat("clpp1, pmp\n") 
cat("LDs, ATG\n")
cat("ETC2, fzo1 & MTG\n")
cat("atfs1, clpp1\n") 
cat("ros, !ETC & (!ampk | hif1)\n")
cat("hsp60, atfs1\n")
cat("hif1, ros & !ampk\n")
cat("atgl1, !hlh11\n")  
cat("lip, lipl4 & LDs & atgl1\n")
cat("betaox, lip\n") 
cat("ATP, (ETC2 & betaox) | ETC\n")
cat("crtc1, !ampk\n")
cat("mtor, !ampk & !unc51\n")
cat("ampk,  (ros & !hif1)\n") 
sink()

Nred2 <- loadNetwork(redclk12)
print(Nred2)
#get attractors#
atracctoresN2 <- getAttractors(Nred2)
plotAttractors(atracctoresN2)

#Test for the wildtype conditions (clk1nor2) and duoble mutant(doblemut2)
clk1nor2 <- fixGenes(Nred2, "clk1", 1)
doblemut2 <- fixGenes(Nred2, c("ampk","clk1"), c(0,0))

atracctoresclk112 <- getAttractors(clk1nor2)
plotAttractors(atracctoresclk112)

atracctoresdoble2 <- getAttractors(doblemut2)
plotAttractors(atracctoresdoble2)

##################################################
#############equations in desolve###############
##################################################
library(deSolve)
net.ode2 <- booleanToODE(Nred2, keep.input = TRUE)
out2 <- ode(func = net.ode2$func, 
           parms = net.ode2$parameters, 
           y = net.ode2$state, 
           times = seq(0, 20, 0.1)) #set of intersting states to start the simulation 
estados_de_interes <- list(c(0,1,0,0,0,1,1,0,1,1,1,1,0,1,1,0,0,0,0,0,1,1,0,1,1,0,0,0,0), 
                           c(0,1,1,0,0,1,1,0,1,1,1,1,0,1,1,0,0,0,0,0,1,1,0,1,1,0,0,0,0), 
                           c(0,1,0,1,0,1,1,0,1,1,1,1,0,1,1,0,0,0,0,0,1,1,0,1,1,0,0,0,0),
                           c(0,1,1,0,1,1,1,0,1,1,1,1,0,1,1,0,0,0,0,0,1,1,0,1,1,0,0,0,0),
                           c(0,1,0,0,0,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0),
                           c(0,1,0,1,0,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0),
                           c(0,1,1,1,0,1,0,0,1,1,1,1,0,0,1,1,1,0,0,0,1,1,1,1,1,0,0,0,0),
                           c(0,1,1,1,1,1,0,0,1,1,1,1,0,0,1,1,1,0,0,0,1,1,1,1,1,0,0,0,0),
                           c(0,1,0,1,0,1,1,0,1,1,1,1,0,0,1,1,1,0,0,0,1,1,1,1,1,0,0,0,0),
                           c(0,1,0,1,1,1,1,0,1,1,1,1,0,0,1,1,1,0,0,0,1,1,1,1,1,0,0,0,0),
                           c(0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1))

estados_de_interes
for(i in 1:11){print(estados_de_interes[[i]])}

names(estados_de_interes) <- c("Atractor 1", "Atractor 2", " Atractor 3", "Atractor 4", "Atractor 5", 
                               "Atractor 6", "Atractor 7", "Atractor 8", "Atractor 9", "Atractor 10", "Estado inicial")
estados_de_interes
library(ggplot2)
library(reshape2)
###########################################################################
##################Continuos model for the clk-1 mutant####################
##########################################################################
for(i in 1:11){state <- validateState(estados_de_interes[[i]], Nred2$genes)
net.ode2 <- booleanToODE(Nred2, keep.input = TRUE)
out2 <- ode(func = net.ode2$func, 
           parms = net.ode2$parameters, 
           y = state, 
           times = seq(0, 40, 0.01))
outggp2 <- melt(as.data.frame(as.matrix(out2)), id='time')
pdf(paste("/home/Estados config 2", names(estados_de_interes)[i] ,"pdf"))
print(ggplot(outggp2, aes(time, value, col=variable)) + 
        geom_line() +
        ggtitle(paste(names(estados_de_interes)[i])) +
        theme_bw())
dev.off()
}

state <- validateState(c(0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1), Nred2$genes)
net.ode2 <- booleanToODE(Nred2, keep.input = TRUE)
out2 <- ode(func = net.ode2$func, 
            parms = net.ode2$parameters, 
            y = state, 
            times = seq(0, 40, 0.01))
outggp2 <- melt(as.data.frame(as.matrix(out2)), id='time')
print(ggplot(outggp2, aes(time, value, col=variable)) + 
        geom_line() +
        ggtitle(paste(names(estados_de_interes)[i])) +
        theme_bw())
tail(out2)
