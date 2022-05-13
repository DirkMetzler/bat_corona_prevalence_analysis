## Load data libraries and data

library(glmmTMB)
library(DHARMa)
library(ape)
library(Rphylip)

data <- read.csv("BatCoV.csv",
                 stringsAsFactors = FALSE)
data$Species <-  paste(data$Bat_genus,data$Species_ep)
data$UniqueID <- as.factor(data$UniqueID)
data$IDstudy <- as.factor(data$IDstudy)
data$Clim <- sapply(strsplit(data$Climate, "_"),"[",1)
data$Clim <- as.factor(data$Clim)

### further below in this script we will read more data:
## from China_points.csv and
## 2015_ShiRabosky_tree.nex.gz

set.seed(123)

## fit a first model

zimod <- glmmTMB(cbind(N_pos,N_sampled-N_pos) ~ HF10  +  Richness + abs(LatDD) +
                     (1|Clim) + (1| Bat_genus / Species) + (1 | IDstudy),
                 ziformula= ~ HF10 +  Richness + abs(LatDD) +  (1|Clim) +
                     (1| Bat_genus / Species) + (1 | IDstudy),
                 data, family="binomial")
summary(zimod)

## function to simulate conditioned on randon effects; needed
## for calculation of quatile residuals with DHARMa functions
simulate.myglmmTMB <- function(obj,nsim=1) {
    ## print("Running simulation conditioned on random effects")
    dfr <- data.frame(row.names=row.names(obj$frame))
    zp <- predict(obj, re.form=NULL, type="zprob") ## zero inflation probability
    cp <- predict(obj, re.form=NULL, type="conditional") ## p parameter if p>0
    cn <- apply(obj$frame[[1]],1,sum)
    for(i in 1:nsim) {
        k <- (1-rbinom(length(cn),size=1,p=zp))*rbinom(length(cn),size=cn,p=cp)
        dfr[paste("sim_",i,sep="")] <- cbind(k,cn-k)
    }
    dfr
}

## We check the model with the simulation-based quantile residuals as calculated
## with the DHARMa package.  For this, we need a simulation function that
## simulates outcomes of the fitted model conditioned on the fitted values for
## the random effects.  (Detail regarding the R code: This function will only be
## used when we assign the class myglmmTMB as primary class to a
## fitted-model object.)

mzimod <- zimod  ## we copy the model and assing class myglmmTMB, sucht that 
                 ## zimod will call the origninal glmmTMB simulator
class(mzimod) <- c("myglmmTMB","glmmTMB")
simres <- simulateResiduals(fittedModel = mzimod)
plot(simres)

## As proposed by the warning message we use a bootstrap-based check for
## outliers:

testOutliers(zimod, type="bootstrap")

## We can also check wether the residuals still contain information on other
## variables that we don't have in our models:

kruskal.test(simres$scaledResiduals~data$Continent)
kruskal.test(simres$scaledResiduals~data$Country)
kruskal.test(simres$scaledResiduals~data$Feeding_gu)
kruskal.test(simres$scaledResiduals~data$sampletime)
 
## The p-values in the summary of the model in which HF10 have impact on the
## presence of the virus and of the infection rate, indicate that this impact
## is only significant for the infection rate. These p-values are based on
## approximations and we check them by simulation. For three different model
## comparisons we simulate data for the respective null model (including the
## estimated random effects), analyse the simulated data both with the null
## model and with the alternative model. We do this 1000 times and check in
## how many of these simulations the log(likelihood ratio) was as large as the
## log(likelihood ratio) for the real data.

## model with no human impact
zi0mod <- glmmTMB(cbind(N_pos,N_sampled-N_pos) ~ Richness +
                      (1|Clim) + abs(LatDD) +
                      (1| Bat_genus / Species) + (1 | IDstudy),
                  ziformula= ~ Richness +  (1|Clim) + abs(LatDD) +
                      (1| Bat_genus / Species) + (1 | IDstudy),
                  data, family="binomial")

## model with human impact on conditional model (prevalence) but not on zero
## inflation (presence)
zi0modB <- glmmTMB(cbind(N_pos,N_sampled-N_pos) ~ HF10 +  Richness  + abs(LatDD)
                   +  (1|Clim) + (1| Bat_genus / Species) + (1 | IDstudy),
                 ziformula= ~ Richness  + abs(LatDD) +  (1|Clim) +
                     (1| Bat_genus / Species) + (1 | IDstudy),
                 data, family="binomial")


## run the simulations or, if already done, load the simulation results from a
## data file

if(file.exists("lrtsimclim.RData")) {
    load("lrtsimclim.RData")
} else {
    class(zi0mod) <- c("myglmmTMB","glmmTMB")
    class(zi0modB) <- c("myglmmTMB","glmmTMB")
    llr <- numeric()
    llr0B <- numeric()
    for(i in 1:1000) {
        llr[i] <- NA
        while(is.na(llr[i])) {
            sim0dat <- simulate(zi0mod)
            s0 <- try(logLik(glmmTMB(as.matrix(sim0dat) ~ Richness   + abs(LatDD)+  (1|Clim) +
                                         (1| Bat_genus / Species) + (1 | IDstudy),
                                     ziformula= ~ Richness  + abs(LatDD) +  (1|Clim) +
                                         (1| Bat_genus / Species) + (1 | IDstudy),
                                     data, family="binomial")))
            s0B <- try(logLik(glmmTMB(as.matrix(sim0dat) ~ HF10   + abs(LatDD)+  Richness  
                                      +  (1|Clim) + (1| Bat_genus / Species) + (1 | IDstudy),
                                     ziformula= ~ Richness + abs(LatDD) +  (1|Clim) +
                                         (1| Bat_genus / Species) + (1 | IDstudy),
                                     data, family="binomial")))
            s1 <- try(logLik(glmmTMB(as.matrix(sim0dat) ~ HF10   + abs(LatDD)+  Richness  
                                     +  (1|Clim) + (1| Bat_genus / Species) + (1 | IDstudy),
                                     ziformula= ~ HF10 +  Richness  + abs(LatDD) +  (1|Clim) +
                                         (1| Bat_genus / Species) + (1 | IDstudy),
                                     data, family="binomial")))
            if(is.numeric(s1) & is.numeric(s0) & is.numeric(s0B)) {
                llr[i] <- s1-s0
                llr0B[i] <- s0B-s0
            } 
            cat("llr",i,llr[i],"\n",file="llr.txt",append=TRUE)
            cat("llr0B",i,llr0B[i],"\n",file="llr.txt",append=TRUE)
        }
    }
    llrB <- numeric()
    for(i in 1:1000) {
        llrB[i] <- NA
        while(is.na(llrB[i])) {
            sim0dat <- simulate(zi0modB)
            s0B <- try(logLik(glmmTMB(as.matrix(sim0dat) ~ HF10   + abs(LatDD)+  Richness 
                                      + (1|Clim) + (1| Bat_genus / Species) + (1 | IDstudy),
                                     ziformula= ~ Richness  + abs(LatDD) +  (1|Clim) +
                                         (1| Bat_genus / Species) + (1 | IDstudy),
                                     data, family="binomial")))
            s1 <- try(logLik(glmmTMB(as.matrix(sim0dat) ~ HF10   + abs(LatDD) +  Richness 
                                     + (1|Clim) + (1| Bat_genus / Species) + (1 | IDstudy),
                                     ziformula= ~ HF10 +  Richness  + abs(LatDD) +  (1|Clim) +
                                         (1| Bat_genus / Species) + (1 | IDstudy),
                                     data, family="binomial")))
            if(is.numeric(s0B) & is.numeric(s1)) {
                llrB[i] <- s1-s0B
            }
            cat("llrB",i,llrB[i],"\n",file="llr.txt",append=TRUE)
        }
    }
    save(llr,llrB,llr0B,file="lrtsimclim.RData")
}

hist(llr,xlim=c(0,13.1),60,
     main="Null model: HF10 no impact,
alternative model: HF10 has impact on virus presence and infection rate")
abline(v=logLik(zimod)-logLik(zi0mod), col="red", lwd=2)

hist(llr0B,xlim=c(0,13.1),60,main="Null model: HF10 no impact,
 alternative model: HF10 has impact on infection rate but not on presence")
 abline(v=logLik(zi0modB)-logLik(zi0mod), col="red", lwd=2)

hist(llrB,xlim=c(0,13.1),60,
     main="Null model: HF10 has impact on infection rate but not on presence,
alternative model: HF10 has impact on virus presence and infection rate")
abline(v=logLik(zimod)-logLik(zi0modB), col="red", lwd=2)

## We calculate p-values for these comparisons:

## Null model: HF50 no impact,
## alternative model: HF50 has impact on virus presence and infection rate.
## number of simulations (out of 1000) that led to a larger likelihood ratio:
sum(logLik(zimod)-logLik(zi0mod) <= llr)
## p-value:
(sum(logLik(zimod)-logLik(zi0mod) <= llr) + 1) / (1 + length(llr))

## Null model: HF50 no impact,
## alternative model: HF50 has impact on infection rate but not on presence.
## number of simulations (out of 1000) that led to a larger likelihood ratio:
sum(logLik(zi0modB)-logLik(zi0mod) <= llr0B)
## p-value:
(sum(logLik(zi0modB)-logLik(zi0mod) <= llr0B) + 1) / (1 + length(llr0B)) 

## Null model: HF50 has impact on infection rate but not on presence,
## alternative model: HF50 has impact on virus presence and infection rate.
## number of simulations (out of 1000) that led to a larger likelihood ratio:
sum(logLik(zimod)-logLik(zi0modB) <= llrB)
## p-value:
(sum(logLik(zimod)-logLik(zi0modB) <= llrB) + 1) / (1 + length(llrB)) 

## Comparison of AIC values (smaller is better):

AIC(zi0mod, zi0modB, zimod)

## The winner of the AIC comparison is the model zi0modB in which humas have
## an impact on the infection rate of infected populations but not on the
## probability that a population is infected.


### Testing variance components for species and genera
## We apply likelihood-ratio tests with chi-square approximation to test whether the
## variance components of species and genera lead to a significantly better fit.
## (Note that more accurate tests exist but are computationally more demanding; the
## chi-square approximation is not very accurate but rather on the conservative site.)

zi0modBwoS <- glmmTMB(cbind(N_pos,N_sampled-N_pos) ~ HF10 +  Richness  + abs(LatDD)
                   +  (1|Clim) + (1| Bat_genus) + (1 | IDstudy),
                 ziformula= ~ Richness  + abs(LatDD) +  (1|Clim) +
                     (1| Bat_genus) + (1 | IDstudy),
                 data, family="binomial")
zi0modBwoSG <- glmmTMB(cbind(N_pos,N_sampled-N_pos) ~ HF10 +  Richness  + abs(LatDD)
                   +  (1|Clim) + (1 | IDstudy),
                 ziformula= ~ Richness  + abs(LatDD) +  (1|Clim) +
                     + (1 | IDstudy),
                 data, family="binomial")
anova(zi0modBwoSG, zi0modBwoS, zi0modB)

zi0modBwoG <- glmmTMB(cbind(N_pos,N_sampled-N_pos) ~ HF10 +  Richness  + abs(LatDD)
                   +  (1|Clim) + (1| Species) + (1 | IDstudy),
                 ziformula= ~ Richness  + abs(LatDD) +  (1|Clim) +
                     (1| Species) + (1 | IDstudy),
                 data, family="binomial")
anova(zi0modBwoG, zi0modB)

## Accounting for between-species variation improves the model fit significantly,
## also when already accounting for between-genera variation.
## Also accounting for variation between genera does only lead to a significant
## improvement when not accounting for species variation.

#### Estimated effects according to the preferred model zi0modB

summary(zi0modB)
fixef(zi0modB)
ranef(zi0modB)$zi$Clim
ranef(zi0modB, condVar=TRUE)$cond$Clim

## Species with strong impact on infection probability that the population is
## uninfected; note that here negative values mean that the probability is
## larger that there are infected individuals in the population and larger
## values mean that the probability is higher that the population is not
## infected.

sz <- ranef(zi0modB)$zi$Species
sel <- abs(sz[["(Intercept)"]]) > 0.1
data.frame(species=rownames(sz)[sel],zi.ranef=sz[sel,])

## Species with strong impact on conditional model; not that the larger the values
## the higher the fraction of infected individual if there is an infection in the
## population. (And negative values mean: if there is an infection it tends to
## affect only few individuals.)

sc <- ranef(zi0modB)$cond$Species
sel <- abs(sc[["(Intercept)"]])>1
data.frame(species=rownames(sc)[sel],cond.ranef=sc[sel,])

## Wald confidence intervals for the fitted model parameters (fixed effects and
## standard deviations of random effects).

confint(zi0modB)

### Overall predictions of the model zi0modB for an average situation
## For the random factors climate, bat species, bat genus and study ID we take
## an average over the factors.
## For mammalian species richness and the distance to the equator (that is, the
## absolute value of latitude) we take a rounded means of the values in our dataset:

round(mean(data$Richness))
round(mean(abs(data$LatDD)))

nd <- data.frame(HF10=0:100/100,Richness=rep(round(mean(data$Richness)),101),
                 Clim=rep(NA,101),LatDD=rep(round(mean(abs(data$LatDD))),101),
                 Bat_genus=rep(NA,101),Species=rep(NA,101),IDstudy=rep(NA,101)
                 )
pr <- predict(zi0modB, t="response", newdata=nd, re.form=NULL,se.fit=TRUE)
## pr <- predict(zimod, t="response", newdata=nd, re.form=NULL,se.fit=TRUE)
plot(0:100/100,rep(1,101),t="l",ylim=c(0,max(pr$fit+pr$se)),
     ylab="Predicted proportion of infected individuals",
     main="Pred. infection rate +/- std. error",
     xlab="Human impact measure HF10")
## lines(0:100/100,pr$fit+pr$se,lty=2)
## lines(0:100/100,pr$fit-pr$se,lty=2)
pr <- predict(zi0modB, t="response", newdata=nd, re.form=NULL,se.fit=TRUE)
lines(0:100/100,pr$fit,col="red")
lines(0:100/100,pr$fit+pr$se,lty=2,col="red")
lines(0:100/100,pr$fit-pr$se,lty=2,col="red")
abline(h=0)

pr <- predict(zi0modB, t="zprob", newdata=nd, re.form=NULL,se.fit=TRUE)
plot(0:100/100,rep(1,101),t="l",ylim=c(0,max(1-pr$fit+pr$se)),
     ylab="Probability of presence of infection",
     main="Pred. infection presence probability +/- std. error",
     xlab="Human impact measure HF10")
lines(0:100/100,1-pr$fit,col="red")
lines(0:100/100,1-pr$fit+pr$se,lty=2,col="red")
lines(0:100/100,1-pr$fit-pr$se,lty=2,col="red")
abline(h=0)

pr <- predict(zi0modB, t="conditional", newdata=nd, re.form=NULL,se.fit=TRUE)
plot(0:100/100,rep(1,101),t="l",ylim=c(0,max(pr$fit+pr$se)),
     ylab="Prevalence in case of present infection",
     main="Pred. prevalence +/- std. error",
     xlab="Human impact measure HF10")
lines(0:100/100,pr$fit,col="red")
lines(0:100/100,pr$fit+pr$se,lty=2,col="red")
lines(0:100/100,pr$fit-pr$se,lty=2,col="red")
abline(h=0)

## Similar for each climate zone and mammals richness and the average distance to the equator in this climate zone:

ldd <- round(tapply(abs(data$LatDD), data$Clim, mean))
rich <- round(tapply(data$Richness, data$Clim, mean))

plot(0:100/100,rep(1,101),t="l",ylim=c(0,0.13),
     ylab="Predicted proportion of infected individuals",
     main="Pred. infection rate",
     xlab="Human impact measure HF10")
abline(h=0)
for(i in 1:length(levels(data$Clim))) { 
    nd <- data.frame(HF10=0:100/100,Richness=rep(rich[i],101),
                 Clim=rep(levels(data$Clim)[i],101),LatDD=rep(ldd[i],101),
                 Bat_genus=rep(NA,101),Species=rep(NA,101),IDstudy=rep(NA,101)
                 )
    pr <- predict(zi0modB, t="response", newdata=nd, re.form=NULL)
    lines(0:100/100, pr, col=i)
}
legend("bottomright", legend=levels(data$Clim), lty=1, col=1:4, bg="white")

plot(0:100/100,rep(1,101),t="l",ylim=c(0,0.9),
     ylab="Probability of presence of infection",
     main="Pred. infection presence probability",
     xlab="Human impact measure HF10")
abline(h=0)
for(i in 1:length(levels(data$Clim))) { 
    nd <- data.frame(HF10=0:100/100,Richness=rep(rich[i],101),
                 Clim=rep(levels(data$Clim)[i],101),LatDD=rep(ldd[i],101),
                 Bat_genus=rep(NA,101),Species=rep(NA,101),IDstudy=rep(NA,101)
                 )
    pr <- predict(zi0modB, t="zprob", newdata=nd, re.form=NULL)
    lines(0:100/100, 1-pr, col=i)
}
legend("bottomright", legend=levels(data$Clim), lty=1, col=1:4, bg="white")

plot(0:100/100,rep(1,101),t="l",ylim=c(0,0.18),
     ylab="Prevalence",
     main="Pred. prevalence in case of presence of the virus",
     xlab="Human impact measure HF10")
abline(h=0)
for(i in 1:length(levels(data$Clim))) { 
    nd <- data.frame(HF10=0:100/100,Richness=rep(rich[i],101),
                 Clim=rep(levels(data$Clim)[i],101),LatDD=rep(ldd[i],101),
                 Bat_genus=rep(NA,101),Species=rep(NA,101),IDstudy=rep(NA,101)
                 )
    pr <- predict(zi0modB, t="conditional", newdata=nd, re.form=NULL)
    lines(0:100/100, pr, col=i)
}
legend("bottomright", legend=levels(data$Clim), lty=1, col=1:4, bg="white")


## The next plots show predicted infection probabilities for combinations of climate
## zone and mammals richness, and the average distance to the equator in this climate zone.
## From top to bottom the lines refer to (mammalian) species richness values of 10, 50, 100,
## and 150.
## (For cool climate the latter two are shown in grey as these values are not realistic
## for cool climate.)
## Note that the predictions refer to an average species of our dataset (that is, species
## effect is not only averaged over the species in the repective climate zone).

ldd <- round(tapply(abs(data$LatDD), data$Clim, mean))
rich <- c(5, 30, 75, 125, 175)

(colo <- c(rgb(.745,.91,1.), rgb(.451,.698,1.), rgb(.0,.443,.996), 
           rgb(.0,.302,.655), rgb(.0,.145,.4511)))

plot(0:100/100,rep(0,101),t="l",ylim=c(0,13),
     ylab="Prevalence [%]",
     main="Cool climate",
     xlab="Human impact measure HF10")
for(i in 1:5) { 
    nd <- data.frame(HF10=0:100/100,Richness=rep(rich[i],101),
                 Clim="Cool",
                 LatDD=rep(ldd["Cool"],101),
                 Bat_genus=rep(NA,101),Species=rep(NA,101),IDstudy=rep(NA,101)
                 )
    pr <- predict(zi0modB, t="response", newdata=nd, re.form=NULL)
    lines(0:100/100, 100*pr, col=colo[i])
}
legend("topleft", pch=15, legend=rich, col=colo, title="Richness",
       bg="white")

## dev.print(device=svg, file="SVG-Figs/zi0Bavgclripredict-1.svg")

(colo <- c(rgb(0.91,0.745,0.996), rgb(0.875,0.447,1), rgb(.667,0,0.898), 
           rgb(.518,0,.663), rgb(.302,.0,.455)))

plot(0:100/100,rep(0,101),t="l",ylim=c(0,13),
     ylab="Prevalence [%]",
     main="Warm climate",
     xlab="Human impact measure HF10")
abline(h=0)
for(i in 1:5) { 
    nd <- data.frame(HF10=0:100/100,Richness=rep(rich[i],101),
                 Clim="Warm",
                 LatDD=rep(ldd["Warm"],101),
                 Bat_genus=rep(NA,101),Species=rep(NA,101),IDstudy=rep(NA,101)
                 )
    pr <- predict(zi0modB, t="response", newdata=nd, re.form=NULL)
    lines(0:100/100, 100*pr, col=colo[i])
}
legend("bottomright", pch=15, legend=rich, col=colo, title="Richness",
       bg="white")


## dev.print(device=svg, file="SVG-Figs/zi0Bavgclripredict-2.svg")

(colo <- c(rgb(1.,.655,.498), rgb(1.,.333,.0), rgb(.898,.298,0.), 
           rgb(.655,.22,.0), rgb(.455,.149,.0)))

plot(0:100/100,rep(0,101),t="l",ylim=c(0,13),
     ylab="Prevalence [%]",
     main="Sub-tropical climate",
     xlab="Human impact measure HF10")
abline(h=0)
for(i in 1:5) { 
    nd <- data.frame(HF10=0:100/100,Richness=rep(rich[i],101),
                 Clim="Sub",
                 LatDD=rep(ldd["Sub"],101),
                 Bat_genus=rep(NA,101),Species=rep(NA,101),IDstudy=rep(NA,101)
                 )
    pr <- predict(zi0modB, t="response", newdata=nd, re.form=NULL)
    lines(0:100/100, 100*pr, col=colo[i])
}
legend("bottomright", pch=15, legend=rich, col=colo, title="Richness",
       bg="white")

## dev.print(device=svg, file="SVG-Figs/zi0Bavgclripredict-3.svg")

(colo <- c(rgb(.906,.902,.0), rgb(1.,.667,.004), rgb(.902,.596,.0), 
           rgb(.659,.439,.0), rgb(.451,.298,.0)))


plot(0:100/100,rep(0,101),t="l",ylim=c(0,13),
     ylab="Prevalence [%]",
     main="Tropical climate",
     xlab="Human impact measure HF10")
abline(h=0)
for(i in 1:5) { 
    nd <- data.frame(HF10=0:100/100,Richness=rep(rich[i],101),
                 Clim="Tropical",
                 LatDD=rep(ldd["Tropical"],101),
                 Bat_genus=rep(NA,101),Species=rep(NA,101),IDstudy=rep(NA,101)
                 )
    pr <- predict(zi0modB, t="response", newdata=nd, re.form=NULL)
    lines(0:100/100, 100*pr, col=colo[i])
}
legend("bottomright", pch=15, legend=rich, col=colo, title="Richness",
       bg="white")

## dev.print(device=svg, file="SVG-Figs/zi0Bavgclripredict-4.svg")

## The next plots show predicted prevalences (conditioned on presence of the virus) 
## for combinations of climate
## zone and mammals richness, and the average distance to the equator in this climate zone.
## From top to bottom the lines refer to (mammalian) species richness values of 10, 50, 100,
## and 150.
## (For cool climate the latter two are shown in grey as these values are not realistic
## for cool climate.)
## Note that the predictions refer to an average species of our dataset (that is, species
## effect is not only averaged over the species in the repective climate zone).
## And the plot are in percent now.

ldd <- round(tapply(abs(data$LatDD), data$Clim, mean))
rich <- c(10, 50, 100, 150)
(colo <- c(rgb(0.4,0.4,1), rgb(0,0,0.8), rgb(.5,.5,.5), rgb(.5,.5,.5)))

plot(0:100/100, rep(1,101),t="l", ylim=c(0,18),
     ylab="Prevalence [%]",
     main="Cool climate",
     xlab="Human impact measure HF10")
abline(h=0)
for(i in 1:4) { 
    nd <- data.frame(HF10=0:100/100,Richness=rep(rich[i],101),
                 Clim="Cool",
                 LatDD=rep(ldd["Cool"],101),
                 Bat_genus=rep(NA,101),Species=rep(NA,101),IDstudy=rep(NA,101)
                 )
    pr <- predict(zi0modB, t="conditional", newdata=nd, re.form=NULL)
    lines(0:100/100, 100*pr, col=colo[i])
}
legend("topleft", pch=15, legend=rich, col=colo, title="Richness",
       bg="white")

## dev.print(device=svg, file="SVG-Figs/zi0Bavgclripredictpr-1.svg")

(colo <- c(rgb(0.8,0.4,1), rgb(0.6,0.2,0.9), rgb(.4,0,0.8), rgb(.3,0,.7)))

plot(0:100/100,rep(1,101),t="l",ylim=c(0,18),
     ylab="Prevalence [%]",
     main="Warm climate",
     xlab="Human impact measure HF10")
abline(h=0)
for(i in 1:4) { 
    nd <- data.frame(HF10=0:100/100,Richness=rep(rich[i],101),
                 Clim="Warm",
                 LatDD=rep(ldd["Warm"],101),
                 Bat_genus=rep(NA,101),Species=rep(NA,101),IDstudy=rep(NA,101)
                 )
    pr <- predict(zi0modB, t="conditional", newdata=nd, re.form=NULL)
    lines(0:100/100, 100*pr, col=colo[i])
}
legend("bottomright", pch=15, legend=rich, col=colo, title="Richness",
       bg="white")

## dev.print(device=svg, file="SVG-Figs/zi0Bavgclripredictpr-2.svg")


(colo <- c(rgb(1,0,0.7), rgb(1,0,0.5), rgb(.9,0,0.4), rgb(.7,0,.2)))

plot(0:100/100,rep(1,101),t="l",ylim=c(0, 18),
     ylab="Prevalence [%]",
     main="Sub-tropical climate",
     xlab="Human impact measure HF10")
abline(h=0)
for(i in 1:4) { 
    nd <- data.frame(HF10=0:100/100,Richness=rep(rich[i],101),
                 Clim="Sub",
                 LatDD=rep(ldd["Sub"],101),
                 Bat_genus=rep(NA,101),Species=rep(NA,101),IDstudy=rep(NA,101)
                 )
    pr <- predict(zi0modB, t="conditional", newdata=nd, re.form=NULL)
    lines(0:100/100, 100*pr, col=colo[i])
}
legend("bottomright", pch=15, legend=rich, col=colo, title="Richness",
       bg="white")

## dev.print(device=svg, file="SVG-Figs/zi0Bavgclripredictpr-3.svg")

(colo <- c(rgb(1,0.7,0.4), rgb(1,0.5,0.2), rgb(.8,0.4,0), rgb(.6,.2, 0)))

plot(0:100/100,rep(1,101),t="l",ylim=c(0, 18),
     ylab="Prevalence [%]",
     main="Tropical climate",
     xlab="Human impact measure HF10")
abline(h=0)
for(i in 1:4) { 
    nd <- data.frame(HF10=0:100/100,Richness=rep(rich[i],101),
                 Clim="Tropical",
                 LatDD=rep(ldd["Tropical"],101),
                 Bat_genus=rep(NA,101),Species=rep(NA,101),IDstudy=rep(NA,101)
                 )
    pr <- predict(zi0modB, t="conditional", newdata=nd, re.form=NULL)
    lines(0:100/100, 100*pr, col=colo[i])
}
legend("topleft", pch=15, legend=rich, col=colo, title="Richness",
       bg="white")

## dev.print(device=svg, file="SVG-Figs/zi0Bavgclripredictpr-4.svg")

## The following plots show for the species that are known to live in the respective (wider)
## area and are in the scope of our study how the probability of an individual to be infected
## depends on HF10 (according to our model prediction, of course).

r <- read.csv("China_points.csv",dec=",")
ep <- r[1:7] 
for(i in 8:46) {
    d <- r[r[i]!="",c(1:6,i)]
    names(d)[7] <- "Bat_species"
    ep <- rbind(ep,d)
}
ep$Clim <- sapply(sapply(as.character(ep$Climate), strsplit, " "),"[",1)
ep$Bat_genus <- sapply(sapply(as.character(ep$Bat_species), strsplit, "_"),"[",1)
ep$Species_ep <- sapply(sapply(as.character(ep$Bat_species), strsplit, "_"),"[",2)
ep$Species <-  paste(ep$Bat_genus,ep$Species_ep)
ep$place <- paste(ep$Country, ep$locality)

par(mfcol=c(1,2))
for(pl in unique(ep$place)) {
    plot(c(), main=pl, xlim=c(0,1), ylim=c(0,35),
         ylab="Proportion of infected individuals [%]",
         xlab="Human impact measure HF10")
    abline(h=0)
    sppresds <- numeric()
    for(sp in unique(ep$Species[ep$place==pl])) {
        i <- which(ep$place==pl & ep$Species==sp)[1]
        nd <- data.frame(HF10=0:100/100,Richness=rep(ep$Mammals[i],101),
                         Clim=rep(ep$Clim[i],101), LatDD=rep(ep$lat[i],101),
                         Bat_genus=rep(ep$Bat_genus[i],101), Species=rep(ep$Species[i],101), IDstudy=rep(NA,101)
                         )
        pr <- predict(zi0modB, t="response", newdata=nd, re.form=NULL,se.fit=TRUE)
        lines(0:100/100, 100*pr$fit, col="brown")
        sppresds[sp] <- pr$fit[50]
    }
    if( !(sp %in% data$Species) ) cat("Species ", sp," not found.\n")
    cat("\nLines in plot for ",pl, "refer to (from top to bottom):\n")
    cat(names(sppresds[order(-sppresds)]), sep=", ", fill=70)
    cat("\n")
}


### Predictions of zi0modB as a function of the parameter values

ranef(zi0modB)$zi$Clim
ranef(zi0modB)$cond$Clim


r <- 90
l <- 20
h <- 0.5
f <- function(x) 1/(1+exp(-x))
c1 <- ranef(zi0modB)$zi$Clim["Cool",] 
c2 <- ranef(zi0modB)$cond$Clim["Cool",]
fixef(zi0modB)
(fz <- fixef(zi0modB)$zi)
(fc <- fixef(zi0modB)$cond)
(1-f(fz[["(Intercept)"]]+fz[["Richness"]]*r+fz[["abs(LatDD)"]]*abs(l)+c1))*
    (f(fc[["(Intercept)"]] + fc[["HF10"]]*h+fc[["Richness"]]*r+fc[["abs(LatDD)"]]*abs(l)+c2))

## Let's check that the R predict command gives the same value for this model
predict(zi0modB, t="response", 
        newdata=data.frame(HF10=h,Richness=90,LatDD=20,Clim="Cool",
                           Bat_genus=NA,Species=NA,IDstudy=NA), re.form=NULL)


## For comparison also with rounded values as given above:
c1 <- 0.2108
c2 <- -0.4972
(1-f(-0.659-0.00000553*r -0.0125 * abs(l) + c1))*(f(-2.508 + 0.576*h -0.001198*r -0.00098*abs(l) + c2))


### Plots for HF10-effects on infection probability for various values of Richness and different Climate zones

## The next plots show predictions of model zi0modB for various combinations of
## parameter values.
## The predicted value is the probability that the infection is present in a population times
## the average prevalence in the case of an infection.
## Note that these predictions are perhaps a bit more stable than those of the next two sections as increasing the
## conditional prevalence can be partly compensated by an decrease of the probability of the presence of an
## infection in the population and vice versa.

## First, we vary HF10 and species Richness and keep the latitude fixed at 20 (or -20):

par(mfcol=c(1,4))
leg <- TRUE
for(cl in levels(data$Clim)) {
    pr <- predict(zi0modB, t="response",
                  newdata=data.frame(HF10=0:100/100,
                                     Clim=cl,
                                     Richness=90,
                                     LatDD=20,
                                     Bat_genus=NA, Species=NA, IDstudy=NA),
                  re.form=NULL)
    plot(0:100/100,pr,ylim=c(0,0.15),t="l",xlab="HF10",
         ylab="Predicted infection probability for average bat",
         main=paste(cl))
    abline(h=0)
    ri  <-  c(20,60,140,180)
    for(i in 1:4) {
        pr <- predict(zi0modB, t="response",
                      newdata=data.frame(HF10=0:100/100,
                                         Clim=cl,
                                         Richness=ri[i],
                                         LatDD=20,
                                         Bat_genus=NA, Species=NA, IDstudy=NA),
                      re.form=NULL)
        lines(0:100/100,pr,col=i+1)
    }
    if(leg) 
        legend("topleft",lty=1,col=c(2,3,1,4,5),
               legend=c(paste("Richness=", c(20,60,90,140,180))))
    leg <- FALSE
}

## Second, we vary HF10 and species latitude and keep species richness fixed at 90:

par(mfcol=c(1,4))
leg <- TRUE
for(cl in levels(data$Clim)) {
    pr <- predict(zi0modB, t="response",
                  newdata=data.frame(HF10=0:100/100,
                                     Clim=cl,
                                     Richness=90,
                                     LatDD=20,
                                     Bat_genus=NA, Species=NA, IDstudy=NA),
                  re.form=NULL)
    plot(0:100/100,pr,ylim=c(0,0.15),t="l",xlab="HF10",
         ylab="Predicted infection probability for average bat",
         main=paste(cl))
    abline(h=0)
    lt  <-  c(0,10,40,60)
    for(i in 1:4) {
        pr <- predict(zi0modB, t="response",
                      newdata=data.frame(HF10=0:100/100,
                                         Clim=cl,
                                         Richness=90,
                                         LatDD=lt[i],
                                         Bat_genus=NA, Species=NA, IDstudy=NA),
                      re.form=NULL)
        lines(0:100/100,pr,col=i+1)
    }
    if(leg) 
        legend("topleft",lty=1,col=c(2,3,1,4,5),
               legend=c(paste("abs(LatDD)=", c(0,10,20,40,60))))
    leg <- FALSE
}

### Check basic assumptions of our main model zi0Bmod and of the null model zi0mod:

zi0B <- zi0modB  ## we copy the model and assing class myglmmTMB, sucht that 
                 ## zi0modB will call the origninal glmmTMB simulator
class(zi0B) <- c("myglmmTMB","glmmTMB")
zi0Bsimres <- simulateResiduals(fittedModel = zi0B)
plot(zi0Bsimres)
testOutliers(zi0modB,type="bootstrap")


zi0 <- zi0mod  ## we copy the model and assing class myglmmTMB, sucht that 
                 ## zi0mod will call the origninal glmmTMB simulator
class(zi0) <- c("myglmmTMB","glmmTMB")
zi0simres <- simulateResiduals(fittedModel = zi0)
plot(zi0simres)
testOutliers(zi0mod,type="bootstrap")

## Check residuals for phylogenetic signal:

## First we have to make sure that the species names in our data are the same as the 
## species names in the tree.
## (In the analysis the tree is reduced to the species for which we have data.)

tree <- read.nexus("2015_ShiRabosky_tree.nex.gz")
translate <- c(
    Glossophaga_comissarisi="Glossophaga_commissarisi",
    Eptesicus_furnalis="Eptesicus_furinalis",
    Coelops_frithii_formosanus="Coelops_frithii",
    Barbastella_darjelingesis="Barbastella_darjelingensis",
    Scotophilus_dingani="Scotophilus_dinganii",
    Rhinolophus_hildebrandtii="Rhinolophus_hildebrandti",
    Rhinolophus_affinus="Rhinolophus_affinis",
    Rousettus_leschenaulti="Rousettus_leschenaultii",
    Pteronotus_parnelli="Pteronotus_parnellii",
    Hipposideros_commersonii="Hipposideros_commersoni",
    Myotis_capaccini="Myotis_capaccinii",
    Myotis_formosus_chofukusei="Myotis_formosus",
    Murina_usuriensis="Murina_ussuriensis",
    Pipistrellus_javaniscas="Pipistrellus_javanicus",
    Phyllostomus_discolour="Phyllostomus_discolor",
    Mormoops_megalohyla="Mormoops_megalophylla",
    Miniopterus_schreibersi="Miniopterus_schreibersii",
    Rhinolophus_darlingi_damarensis="Rhinolophus_darlingi",
    Eptesicus_hottentottus="Eptesicus_hottentotus",
    Nycticeinops_schlieffenii="Nycticeinops_schlieffeni",
    Neoromicia_nana="Neoromicia_nanus",
    Scotophilus_dingamii="Scotophilus_dinganii",
    Rhinolophus_rouxi="Rhinolophus_rouxii",
    Myotis_oxygenatus="Myotis_oxygnathus"
    )
sp <- gsub(" ","_",data$Species)
res <- residuals(simres)
resB <- residuals(zi0Bsimres)
res0 <- residuals(zi0simres)
nres <- numeric()
nresB <- numeric()
nres0 <- numeric()
for(i in 1:nrow(data)) {
    s <- sp[i]
    if(length(which(tree$tip.label==s))==1) {
        nres <- c(nres,res[i])
        nresB <- c(nresB,resB[i])
        nres0 <- c(nres0,res0[i])
        names(nres)[length(nres)] <- s
        names(nresB)[length(nresB)] <- s
        names(nres0)[length(nres0)] <- s
    } else {
        trs <- translate[s]
        if(!is.na(trs)) {
            nres <- c(nres,res[i])
            names(nres)[length(nres)] <- trs
            nresB <- c(nresB,resB[i])
            names(nresB)[length(nresB)] <- trs
            nres0 <- c(nres0,res0[i])
            names(nres0)[length(nres0)] <- trs
        }
    }
}
( tree <- keep.tip(tree,unique(names(nres))) )


## Phylip contrasts analysis accounting for within-species variation:
## (change path to phylip executables.)

Rcontrast(tree, nres,
          "/usr/opt/src/phylip-3.697/exe/")
Rcontrast(tree, nresB,
          "/usr/opt/src/phylip-3.697/exe/")
Rcontrast(tree, nres0,
          "/usr/opt/src/phylip-3.697/exe/")
