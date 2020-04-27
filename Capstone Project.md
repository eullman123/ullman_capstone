---
title: "Capstone 2020"
author: "Elijah Ullman"
date: "4.26.2020"
output: rmarkdown::github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
rmarkdown::render()
```

```{r include=FALSE}
library(tidyverse)
library(ggplot2)
library(nlfitr)
library(tidyr)
library(broom)
library(minpack.lm)
library(ez)
library(DescTools)

```

**1) Depression is a widespread mental illness affecting an estimated 17.3 million adult Americans without a unanimous consensus as to the cause It is characterized by loss of interest in previously enjoyed activities, decrease in mood and energy, and in the most severe cases, suicidal ideologies. The classical psychedelic drugs such as LSD and psilocybin share the N-methylated indoleamine core with serotonin. Serotonin and the receptor itself have been targets for treatment of depression even prior to the monoamine hypothesis via the tricyclic series of compounds and selective serotonin agonists like Losartin.  The metabotropic 5HT2 receptor family consists of  5HT2A-C and couples preferentially to Gq/11 following agonist binding. This canonical interaction recruits phospholipase C (PLC) which leads to a downstream activation of protein kinase C (PKC) and subsequent influx of calcium. The non-canonical GPCR pathway results in phosphorylation of the receptor by GPCR kinase (GRK) and subsequent recruitment of 𝛽-arrestin1/2 which leads to internalization of the receptor. 𝛽-arrestin1/2 also acts as a scaffold for recruitment and integration of the MAPK/ERK pathway for gene transcription. Biased agonism is when a ligand elicits more activity towards the non-canonical pathway relative to the cononical pathway.**

**2)It is not known whether biased agonism through the 𝛽-arrestin1/2 pathwayis required for the significant anti-depressant effects observed with the hallucinogenic 5HT2A receptor partial agonists.**

**3) While both LSD and its analog lisuride both act as partial agonists at the 5HT2A receptor, only LSD elicits hallucinatory activity, indicating a molecular bias towards hallucinatory properties of certain 5HT2A ligands. Although the binding affinity of psilocybin to the various serotonin receptors have been characterized, signaling bias towards the 𝛽-arrestin1/2 pathway has not been characterized with psilocybin. I hypothesize that the antidepressant activity of psilocybin is partially mediated through the 𝛽-arrestin1/2 biased pathway. Therefore, psilocybin should recruit greater amounts of 𝛽-arrestin1/2 than lizuride ** 

**4) The dependent variable is 𝛽-arrestin1/2 recruitment and the predictor variable is the drug evaluated, whether it is DMSO as a negative control, 2C-N which is an established 5HT2A biased agonist as a positive control, or lisuride which binds the 5HT2AR but does not produce anti-depressant effects. The dependent variable is measured and the independent variable is sorted. **

**5) Null Hypothesis: There will be no statistically significant difference in 𝛽-arrestin1/2 recruitment between lisuride and psilocybin. Alternate Hypothesis: The 𝛽-arrestin1/2 recruitment for psilocybin will be greater than in lisuride.**


**6)I will use One way ANOVA with a pairwise T test between the ratios of beta arrestin recruitment and IP3 signaling. The pairiwse T test is utilized as HEK293 are immortalized clones of one another and thus inherintly linked. One way ANOVA is used as I am comparing more than two groups but only one variable. The meaningful comparison that I am looking for is an increased ratio of biased agonism in psilocybin in comparison to lisuride. I normalize the ratios in comparison to serotonin as it is the endogenous ligand of the 5HT2AR. **

**7) I would choose HEK293 cells as they have been utilized to transfect 5HT2AR and for quantification of beta-arrestin recruitment. I utilized DMSO as a negative control to ensure that there is no aberrant signaling without a ligand. 2CN is used here as it is another partial agonist of the 5HT2AR and functions as a biased agonist, so I wanted another biased agonist to compare Psilocybin to. I chose to utilize max beta-arrestin and IP3 signaling rather than signaling @ EC50 as this would inform us about the efficacy of the drugs and normalize endogenous ligand so I could determine if the assayed drugs act as biased ligands. The independent replicates are the HEK cells in different passages while the technical replicates are the numbers of times the experiments are repeated**

```{r}
set.seed(12357)
doses <- c(1e-3, 3e-3,1e-2,3e-2,1e-1,3e-1,1,3,10)
#Beta-Arrestin Recruitment


ser <- simlogdr(doses,logk = -0.5,ylo=0,yhi = 5000, sd=200,reps=5,h=1)$data #5HT
psil <- simlogdr(doses,logk = -0.3,ylo=0,yhi = 5000, sd=200,reps=5,h=1)$data  #Psil
liz <- simlogdr(doses,logk = -1,ylo=0,yhi = 3000, sd=200,reps=5,h=1)$data  #Liz
tcn <- simlogdr(doses,logk = -0.1,ylo=0,yhi = 6000, sd=200,reps=5,h=1)$data  #2CN
dmso <- simlogdr(doses,logk = -1,ylo=0,yhi = 100, sd=20,reps=5,h=1)$data  # DMSO

add_name=function(df,name) {
  df %>% mutate(drug=rep(as.character(name), nrow(.)))
}
liz.df <- add_name(liz,"Lizuride")
psil.df <- add_name(psil,"Psilocybin")
DMSO.df <- add_name(dmso,"DMSO")
tcn.df <- add_name(tcn,"Two-CN")
ser.df <- add_name(ser,"Serotonin")

fullframe <- full_join(liz.df,psil.df) %>% full_join(.,DMSO.df) %>% full_join(.,tcn.df) %>% full_join(., ser.df)
fullframe

set.seed(1234)

ggplot(fullframe, aes(x,y,group=drug,color=drug))+
  geom_point()+
  geom_smooth(method="nlsLM",
              formula="y ~ ylo + (yhi-ylo)/(1+10^((logk-x)*h))",
              method.args=list(
                start=c(yhi=5000,ylo=0,logk=-1,h=1)),
              se=F)+
  ylab("Recruitment of B-Arrestin (A.U.)")+
  xlab("Concentration (A.U)")+
  ggtitle("Drug Induced B-Arrestin Recruitment")
```
![B-arrestin recruitment](https://user-images.githubusercontent.com/64387394/80413726-ea473680-889d-11ea-88e1-fe9b4d79eb1b.PNG)

```{r}
set.seed(123567)

simlogdr(doses,logk = -0.5,ylo=0,yhi = 300, sd=30,reps=5,h=1) #5HT
simlogdr(doses,logk = -0.3,ylo=0,yhi = 150, sd=30,reps=5,h=1) #Psil
simlogdr(doses,logk = -1,ylo=0,yhi = 150, sd=30,reps=5,h=1) #Liz
simlogdr(doses,logk = -0.1,ylo=0,yhi = 150, sd=30,reps=5,h=1) #2CN
simlogdr(doses,logk = -1,ylo=0,yhi = 10, sd=5,reps=5,h=1) # DMSO

ser2 <- simlogdr(doses,logk = -0.5,ylo=0,yhi = 300, sd=30,reps=5,h=1)$data #5HT
psil2 <- simlogdr(doses,logk = -0.3,ylo=0,yhi = 150, sd=30,reps=5,h=1)$data #Psil
liz2 <- simlogdr(doses,logk = -1,ylo=0,yhi = 150, sd=30,reps=5,h=1)$data #Liz
tcn2 <- simlogdr(doses,logk = -0.1,ylo=0,yhi = 150, sd=30,reps=5,h=1)$data #2CN
dmso2 <- simlogdr(doses,logk = -1,ylo=0,yhi = 10, sd=5,reps=5,h=1)$data # DMSO

add_name=function(df,name) {
  df %>% mutate(drug=rep(as.character(name), nrow(.)))
}

liz2.df <- add_name(liz2,"Lizuride")
psil2.df <- add_name(psil2,"Psilocybin")
DMSO2.df <- add_name(dmso2,"DMSO")
tcn2.df <- add_name(tcn2,"Two-CN")
ser2.df <- add_name(ser2,"Serotonin")

fullframe2 <- full_join(liz2.df,psil2.df) %>% full_join(.,DMSO2.df) %>% full_join(.,tcn2.df) %>% full_join(., ser2.df)
fullframe2

ggplot(fullframe2, aes(x,y,group=drug,color=drug))+
  geom_point()+
  geom_smooth(method="nlsLM",
              formula="y ~ ylo + (yhi-ylo)/(1+10^((logk-x)*h))",
              method.args=list(
                start=c(yhi=200,ylo=0,logk=-1,h=1)),
              se=F)+
  ylab("IP3 Signalling (A.U.)")+
  xlab("Concentration (A.U)")+
  ggtitle("Drug-Induced IP3 Signaling")

```
![IP3 recruitment](https://user-images.githubusercontent.com/64387394/80413724-ea473680-889d-11ea-9583-8024701e82ba.PNG)

```{r}
fitliz<- fitlogdr(x,y,liz,logk=-1,ylo=0,yhi=5000,h=1,weigh = F)
lizparams <- tidy(fitliz) %>% select(term,estimate, std.error) %>% mutate(drug=rep("liz",4));lizparams

fitpsil<- fitlogdr(x,y,psil,logk=-1,ylo=0,yhi=5000,h=1,weigh = F)
psilparams <- tidy(fitpsil) %>% select(term,estimate, std.error) %>% mutate(drug=rep("psil",4));psilparams

fitDMSO<- fitlogdr(x,y,dmso,logk=0.0003,ylo=0,yhi=100,h=1,weigh = F)
DMSOparams <- tidy(fitDMSO) %>% select(term,estimate, std.error) %>% mutate(drug=rep("DMSO",4));DMSOparams

fittcn<- fitlogdr(x,y,tcn,logk=0.6,ylo=0,yhi=5000,h=1,weigh = F)
tcnparams <- tidy(fittcn) %>% select(term,estimate, std.error) %>% mutate(drug=rep("TCN",4));tcnparams

fitser<- fitlogdr(x,y,ser,logk=0.6,ylo=0,yhi=5000,h=1,weigh = F)
serparams <- tidy(fitser) %>% select(term,estimate, std.error) %>% mutate(drug=rep("Ser",4));serparams

lisurideB <- lizparams %>% filter(term == "yhi") %>%  select(B = estimate, B.err = std.error) %>% mutate(drug = "liz")
PsilocybinB <- psilparams %>% filter(term == "yhi") %>%  select(B = estimate, B.err = std.error) %>% mutate(drug = "psil")
DMSOB <- DMSOparams %>% filter(term == "yhi") %>%  select(B = estimate, B.err = std.error) %>% mutate(drug = "dmso")
tcnB <- tcnparams %>% filter(term == "yhi") %>%  select(B = estimate, B.err = std.error) %>% mutate(drug = "tcn")
serB <- serparams %>% filter(term == "yhi") %>%  select(B = estimate, B.err = std.error) %>% mutate(drug = "ser")

lisurideB %>% rbind(., serB)
Bs <- rbind(lisurideB, PsilocybinB, DMSOB, tcnB, serB);Bs


fitliz2<- fitlogdr(x,y,liz2,logk=-1,ylo=0,yhi=200,h=1,weigh = F)
lizparams2 <- tidy(fitliz2) %>% select(term,estimate, std.error) %>% mutate(drug=rep("liz",4));lizparams2

fitpsil2<- fitlogdr(x,y,psil2,logk=-1,ylo=0,yhi=200,h=1,weigh = F)
psilparams2 <- tidy(fitpsil2) %>% select(term,estimate, std.error) %>% mutate(drug=rep("psil",4));psilparams2

fitDMSO2<- fitlogdr(x,y,dmso2,logk=0.0003,ylo=0,yhi=10,h=1,weigh = F)
DMSOparams2 <- tidy(fitDMSO2) %>% select(term,estimate, std.error) %>% mutate(drug=rep("DMSO",4));DMSOparams2

fittcn2<- fitlogdr(x,y,tcn2,logk=0.6,ylo=0,yhi=200,h=1,weigh = F)
tcnparams2 <- tidy(fittcn2) %>% select(term,estimate, std.error) %>% mutate(drug=rep("TCN",4));tcnparams2

fitser2<- fitlogdr(x,y,ser2,logk=0.6,ylo=0,yhi=200,h=1,weigh = F)
serparams2 <- tidy(fitser2) %>% select(term,estimate, std.error) %>% mutate(drug=rep("Ser",4));serparams2

lisurideIP <- lizparams2 %>% filter(term == "yhi") %>%  select(IP = estimate, IP.err = std.error) %>% mutate(drug = "liz")
PsilocybinIP <- psilparams2 %>% filter(term == "yhi") %>%  select(IP = estimate, IP.err = std.error) %>% mutate(drug = "psil")
DMSOIP <- DMSOparams2 %>% filter(term == "yhi") %>%  select(IP = estimate, IP.err = std.error) %>% mutate(drug = "dmso")
tcnIP <- tcnparams2 %>% filter(term == "yhi") %>%  select(IP = estimate, IP.err = std.error) %>% mutate(drug = "tcn")
serIP <- serparams2 %>% filter(term == "yhi") %>%  select(IP = estimate, IP.err = std.error) %>% mutate(drug = "ser")

IPs <- rbind(lisurideIP, PsilocybinIP, DMSOIP, tcnIP, serIP);IPs
IPs

Bs;IPs

df <- full_join(Bs, IPs) %>% mutate(Ratio = B/IP, Error = B.err/IP.err) %>% select(drug, Ratio, Error) %>% filter(drug != "dmso") 
df

ggplot(df, aes(drug, Ratio, fill = drug))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = Ratio-Error, ymax = Ratio+Error),
                width = 0.1)+
  ylab("Ratio of B-Arrestin Recruitmen/IP3 Signal")
```
![Barrestin to IP3](https://user-images.githubusercontent.com/64387394/80414119-86713d80-889e-11ea-9ac4-e7b13334d0aa.PNG)

```{r}

b = df$Ratio[4] #expected basal outcome value - set to 1
a = df$Ratio[3] #expected fold-to-basal effect of our positive control
f = df$Ratio[2]
g = df$Ratio[1]
sd = mean(df$Error) #expected standard deviation of Outcome variable
n = 3 # number of independent replicates per group

sims = 100 #number of Monte Carlo simulations to run. 
CRdataMaker <- function(n, b, a, f,g,sd) { 
  
  a1 <- rnorm(n, b, sd) #basal or negative ctrl
  a2 <- rnorm(n, (a), sd) #positive control or some other treatment
  a3 <- rnorm(n, (f), sd) #treatment effect
   a4 <- rnorm(n, (g), sd) #Comparison
    Outcome <- c(a1, a2, a3, a4)
    Predictor <- c(rep(c("Ser", "TCN", "Psil","Liz"), each = n))
    ID <- as.factor(rep(c(1:n), 4))
    df <-data.frame(ID, Predictor, Outcome)
    }

dat <- CRdataMaker(n,b,a,f,g,sd)
dat
ezANOVA(dat, wid = ID,
            dv = Outcome,
            within = Predictor,
            type = 2
            )


pt <- pairwise.t.test(dat$Outcome, dat$Predictor, p.adjust.method = "bonf", paired = T)
pt$p.value 

ggplot(dat, aes(Predictor, Outcome))+
  geom_jitter(width=0.1,size = 4, alpha=0.5)

```
![variance predictor](https://user-images.githubusercontent.com/64387394/80414122-86713d80-889e-11ea-9f7a-668603117076.PNG)
```{r}
pval <- replicate(
  sims, {
 
    sample.df <- CRdataMaker(n, b, a, f,g, sd)
    
    suppressMessages(sim.ezaov <- ezANOVA(
            data = sample.df, 
            wid = ID,
            dv = Outcome,
            within = Predictor,
            type = 2
            )
    )
  pval <- sim.ezaov$ANOVA[1,5]
    
    }
  )
pval
pwr.pct <- sum(pval<0.05)/sims*100
paste(pwr.pct, sep="", "% power.")
```

```{r pressure, echo=FALSE}
ggplot(data.frame(pval))+
  geom_histogram(aes(pval), fill="#d28e00")+ 
  labs(x="p-value")
```

```{r}
n = 12

pval <- replicate(
  sims, {
 
    sample.df <- CRdataMaker(n, b, a, f,g, sd)
    
    pt <- pairwise.t.test(sample.df$Outcome, sample.df$Predictor, p.adjust.method = "bonf", paired = T)

    pval <-  pt$p.value[1, 1]
    
    }
  )
pval
pwr.pct <- sum(pval<0.05)/sims*100
paste(pwr.pct, sep="", "% power.")
```
![p value](https://user-images.githubusercontent.com/64387394/80414115-85d8a700-889e-11ea-8cc5-01d90c74a47d.PNG)

