################################################################################
##  TEMPEST_understory_sppComp.R
##
##  Author: Kimberly Komatsu
##  Date created: December 4, 2021
################################################################################

library(codyn)
library(grid)
library(PerformanceAnalytics)
library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\SERC TEMPEST\\data') #kim's laptop


##### import cover data ####
data <- read.csv('TEMPEST_SERC_understory sppcomp_2021.csv')


#### calculate relative cover ####
cover <- data%>%
  filter(!(species %in% c('infrastructure cover', 'canopy tree basal cover')))%>%
  mutate(rep=paste(plot, subplot, sep='::'))%>%
  filter(subplot!='additional')

totCover <- cover%>%
  group_by(rep)%>%
  summarise(tot_cover=sum(cover))%>%
  ungroup()

relCover <- cover%>%
  left_join(totCover)%>%
  mutate(rel_cover=cover/tot_cover)


commMetrics <- community_structure(df=relCover, time.var=NULL, abundance.var="rel_cover", replicate.var="rep")



