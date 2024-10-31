#GLM on all cages throughout expansion
source("/home/users/mcbitter/OrchardProject/Code/config.R")
source("/home/users/mcbitter/OrchardProject/Code/helper_functions.R")
source("/home/users/mcbitter/OrchardProject/Code/load_packages.R")
source("/home/users/mcbitter/OrchardProject/Code/plotting_functions.R")
source("/home/users/mcbitter/OrchardProject/Code/workflow_functions.R")
source("/home/users/mcbitter/OrchardProject/Code/general_cage_functions.R")
source("/home/users/mcbitter/OrchardProject/Code/IndoorCageExp/indoor.cage.functions.R")

setwd('~/dpetrov/MarkB/IndoorCageExp/')

vec = c(0:9)
load(paste0('./RData/indoor_cages_filteredx2.RData'))
df.glm = sites
samps = samps %>% mutate(tpt = as.numeric(as.character(tpt)))
afmat = cbind(samps, t(afmat))
eec = cbind(samps, t(eec))
afmat = as.data.frame(afmat %>%
    filter(phase %in% c('exp1', 'exp2')))
eec = as.data.frame(eec %>%
filter(phase %in% c('exp1', 'exp2')))
samps = afmat[,1:ncol(samps)]
afmat = afmat[,-c(1:ncol(samps))]
afmat = as.data.frame(t(afmat))
eec = eec[,-c(1:ncol(samps))]
eec = as.data.frame(t(eec))
res = as.data.frame(fit_GLM_ContinuousTime(afMatrix = afmat ,rdMatrix = eec, vec = vec,
                                               sampleData = samps, model.vars = 'tpt', poolCt=100, ncores = 16))
df.glm = as.data.frame(cbind(df.glm, res))
save(df.glm, file = paste0('./GLM/glm.FullExpansionx2.RData'))