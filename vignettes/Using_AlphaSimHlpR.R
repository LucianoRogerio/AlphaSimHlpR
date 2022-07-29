## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Load packages------------------------------------------------------------
# Make sure you have the right packages installed
neededPackages <- c("AlphaSimR", "dplyr", "tidyr", "plotrix", "lme4", "sommer", "optiSel")
for (p in neededPackages) if (!require(p, character.only=T)) install.packages(p)
suppressMessages(library(AlphaSimHlpR))

## ----Define Population--------------------------------------------------------
bsp <- specifyPopulation(ctrlFileName="../inst/PopulationCtrlFile_Small.txt")
bsp <- specifyPipeline(bsp, ctrlFileName="../inst/ControlFile_Small.txt")
bsp <- specifyCosts(bsp, ctrlFileName="../inst/CostsCtrlFile_Small.txt")
print(bsp)

nReplications <- 3
bsp$nCyclesToRun <- 6

## ----Replicate Scheme---------------------------------------------------------
replicRecords <- lapply(1:nReplications, runBreedingScheme, nCycles=bsp$nCyclesToRun, initializeFunc=initFuncADChk, productPipeline=prodPipeFncChk, populationImprovement=popImprov1Cyc, bsp)

## ----Calculate means----------------------------------------------------------
plotData <- plotRecords(replicRecords)
meanMeans <- tapply(plotData$genValMean, list(plotData$year, plotData$stage), mean)
meanMeans <- meanMeans[,c("F1", bsp$stageNames)]
stdErrMeans <- tapply(plotData$genValMean, list(plotData$year, plotData$stage), std.error)
stdErrMeans <- stdErrMeans[,c("F1", bsp$stageNames)]
print(meanMeans)
print(stdErrMeans)

