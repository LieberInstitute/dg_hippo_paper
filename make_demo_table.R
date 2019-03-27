library(SummarizedExperiment)
library(jaffelab)

## load data
load("count_data/astellas_dg_hg38_rseGene_n263.rda")
pd = colData(rse_gene)

## groups
pd$metricGroups = factor(pd$Dx,levels=c("Control", "Schizo", "Bipolar", "MDD"))

# make group indices
groupIndexes=splitit(pd$metricGroups)

# number
N = sapply(groupIndexes,length)
N

# sex
sex = sapply(groupIndexes, function(x) table(pd$Sex[x]))
sexF= signif(prop.table(sex,2),3)[1,]

# race
race = sapply(groupIndexes, function(x) table(pd$Race[x]))
raceCauc =  sapply(race, function(x) prop.table(x)["CAUC"])

### age
age = tapply(pd$Age, pd$metricGroups, mean)
agesd = tapply(pd$Age, pd$metricGroups, sd)

### RIN
rin = tapply(pd$RIN, pd$metricGroups, mean,na.rm=TRUE)
rinsd = tapply(pd$RIN, pd$metricGroups, sd,na.rm=TRUE)
paste0(round(rin,1), "(", round(rinsd,1), ")")

### mapping rate ##
mr = tapply(100*pd$overallMapRate, pd$metricGroups, mean)
mrsd = tapply(100*pd$overallMapRate, pd$metricGroups, sd)

## assign rate
ar = tapply(100*pd$totalAssignedGene, pd$metricGroups, mean)
arsd = tapply(100*pd$totalAssignedGene, pd$metricGroups, sd)

##############
## any significant?
summary(glm(factor(pd$Sex) ~ pd$metricGroups,family="binomial"))
summary(glm(factor(pd$Race == "CAUC") ~ pd$metricGroups,family="binomial"))

## any significant?
summary(lm(Age ~ metricGroups, data=pd))$coef
summary(lm(RIN ~ metricGroups, data=pd))$coef
summary(lm(overallMapRate ~ metricGroups, data=pd))$coef
summary(lm(totalAssignedGene ~ metricGroups, data=pd))$coef

## make table
tab = rbind(N, sexF*100,signif(raceCauc*100,3),
	paste0(signif(age,3), "(", signif(agesd,3), ")"),
	paste0(round(rin,1), "(", round(rinsd,1), ")"),
	paste0(round(mr,1), "(", round(mrsd,1), ")"),
	paste0(round(ar,1), "(", round(arsd,1), ")"))
rownames(tab) = c("N", "Sex(%F)","Race(%Cauc)",
	"Age:Mean(SD)", "RIN:Mean(SD)", "MapRate:Mean(SD)",
	"ExonRate:Mean(SD)")
colnames(tab)= c("CONT", "SCZD", "BPD", "MDD")

## add signif
tab[3,3:4] = paste0(tab[3,3:4], "*")
tab[4,2:3] = paste0(tab[4,2:3], "*")
tab[5,4] = paste0(tab[5,4], "*")
tab[6,2] = paste0(tab[6,2], "*")
tab[7,2:3] = paste0(tab[7,2:3], "*")

# write out
write.csv(tab, "tables/supp_table_demo.csv")
