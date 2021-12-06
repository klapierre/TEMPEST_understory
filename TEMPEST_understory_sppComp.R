################################################################################
##  TEMPEST_understory_sppComp.R
##
##  Author: Kimberly Komatsu
##  Date created: December 4, 2021
################################################################################

library(codyn)
library(vegan)
library(grid)
library(PerformanceAnalytics)
library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\SERC TEMPEST\\data') #kim's laptop
setwd('C:\\Users\\komatsuk\\Dropbox (Smithsonian)\\SERC TEMPEST\\data') #kim's desktop

##### functions and themes #####
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  



##### import cover data ####
data <- read.csv('TEMPEST_SERC_understory sppcomp_2021.csv')%>%
mutate(rep=paste(plot, subplot, sep='::'))

trt <- data%>%
  select(project_name, calendar_year, treatment_year, season, plot, subplot, rep)%>%
  unique()


##### species list #####
sppList <- data%>%
  filter(!(species %in% c('infrastructure cover', 'canopy tree basal cover', 'bare ground')))%>%
  select(project_name, calendar_year, treatment_year, plot, species)%>%
  unique()


#### calculate relative cover ####
cover <- data%>%
  filter(!(species %in% c('infrastructure cover', 'canopy tree basal cover', 'bare ground')))%>%
  filter(subplot!='additional')

totCover <- cover%>%
  group_by(rep)%>%
  summarise(tot_cover=sum(cover))%>%
  ungroup()

relCover <- cover%>%
  left_join(totCover)%>%
  mutate(rel_cover=cover/tot_cover)


##### calculate community metrics #####
subplotMetrics <- community_structure(df=relCover, time.var=NULL, abundance.var="rel_cover", replicate.var="rep")%>%
  left_join(trt)

plotMetrics <- sppList%>%
  group_by(project_name, calendar_year, treatment_year, plot)%>%
  summarise(richness=length(species))%>%
  ungroup()

#analysis
summary(aov(richness ~ plot, data=subplotMetrics)) #p=0.00275
summary(aov(Evar ~ plot, data=subplotMetrics)) #p=0.156


#figures
subplotRichnessFig <- ggplot(data=barGraphStats(data=subplotMetrics, variable="richness", byFactorNames=c("plot")), aes(x=plot, y=mean)) +
  geom_bar(stat='identity', fill='light gray', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Plot') + ylab('Subplot Species Richness')

subplotEvennessFig <- ggplot(data=barGraphStats(data=subplotMetrics, variable="Evar", byFactorNames=c("plot")), aes(x=plot, y=mean)) +
  geom_bar(stat='identity', fill='light gray', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Plot') + ylab('Subplot Evenness')

plotRichnessFig <- ggplot(data=plotMetrics, aes(x=plot, y=richness)) +
  geom_bar(stat='identity', fill='light gray', color='black') +
  xlab('Plot') + ylab('Plot Species Richness') +
  annotate("text", x= 1, y = 74, label= "70", size = 7)+
  annotate("text", x= 2, y = 77, label= "73", size = 7)+
  annotate("text", x= 3, y = 63, label= "59", size = 7)
  
pushViewport(viewport(layout=grid.layout(1,3)))
print(subplotRichnessFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(subplotEvennessFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(plotRichnessFig, vp=viewport(layout.pos.row=1, layout.pos.col=3))
#export at 1400x800


##### community composition #####
#NMDS
sppMatrix <- relCover%>%
  select(plot, rep, species, rel_cover)%>%
  spread(key=species, value=rel_cover, fill=0)

sppBC <- metaMDS(sppMatrix[,3:92])

plots <- 1:nrow(sppMatrix)
plotData <- sppMatrix[,1:2]
plot(sppBC$points,col=as.factor(plotData$plot))
ordiellipse(sppBC, groups = as.factor(plotData$plot), kind = "sd", display = "sites", label = T)


#PERMANOVA
permanovaData <- sppMatrix %>%
  select(-plot, -rep) 
Permanova <- adonis2(formula = permanovaData ~ plot, data = plotData, permutations = 999, method= "bray")
print(Permanova) #p = .001

#betadisper
Veg <- vegdist(permanovaData, method = "bray")
dispersion <- betadisper(Veg, plotData$plot)
permutest(dispersion, pairwise = TRUE, permutations = 999) 

#SIMPER
simper(permanovaData, group=plotData$plot)

#figure
ggplot(data=barGraphStats(data=subset(relCover, species %in% c('Berberis thunbergii', 'Botrychium dissectum', 'Polygonum virginianum', 'Carex sp', 'Elaeagnus umbellata', 'Epifagus virginiana', 'Galium circaezans', 'Lonicera japonica', 'Lindera benzoin', 'Ilex opaca', 'Mitchella repens', 'Parthenocissus quinquefolia', 'Rhus radicans', 'Rubus phoenicolasius', 'Symphyotrichum lateriflorus', 'Sceptridium biternatum')), variable="rel_cover", byFactorNames=c("plot", "species")), aes(x=plot, y=mean)) +
  geom_bar(stat='identity', fill='light gray', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Plot') + ylab('Relative Cover') +
  facet_wrap(~species, scales='free_y')
#export at 1200x1200




