## this script loads the TMB (C++) models into the global
## workspace
setwd(here())
if('hmmTMB' %in% .packages()) detach(package:hmmTMB)
if('TMB' %in% .packages()) detach(package:TMB)
library(TMB)

TMB::runExample('simple')
obj.simple <- obj

library(dsem)
obj.dlm <- readRDS('models/dlm/obj.dlm.RDS')
obj.dlm$retape()
obj.causal <- readRDS('models/causal/obj.causal.RDS')
obj.causal$retape()


setwd(here('models/pollock/'))
dyn.load('pollock')
obj.pollock <- readRDS('obj.pollock.RDS')
obj.pollock$retape()
setwd(here())

library(sdmTMB)
obj.sdmTMB <- readRDS('models/sdmTMB/obj.sdmTMB.RDS')
obj.sdmTMB$retape()

setwd(here('models/swallows/'))
dyn.load('swallows')
obj.swallows <- readRDS('obj.swallows.RDS')
obj.swallows$retape()
setwd(here())

setwd(here('models/wildf/'))
dyn.load('wildf')
obj.wildf <- readRDS('obj.wildf.RDS')
obj.wildf$retape()
setwd(here())

library(glmmTMB)
setwd(here('models/salamanders'))
obj.salamanders <- readRDS('obj.salamanders.RDS')
obj.salamanders$retape()
setwd(here())

setwd(here('models/gp_pois_regr'))
obj.gp_pois_regr <- readRDS('obj.gp_pois_regr.RDS')
dyn.load('gp_pois_regr')
obj.gp_pois_regr$retape()
setwd(here())


library(hmmTMB)
setwd(here('models/petrel'))
obj.petrel <- readRDS('obj.petrel.RDS')
obj.petrel$retape()
setwd(here())

library(wham)
setwd(here('models/wham'))
obj.wham <- readRDS('obj.wham.RDS')
obj.wham$retape()
setwd(here())


library(stockassessment)
setwd(here('models/sam'))
obj.sam <- readRDS('obj.sam.RDS')
obj.sam$retape()
setwd(here())
