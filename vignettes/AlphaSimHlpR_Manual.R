## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  comment = "#>"
)

## ----Install AlphaSimHlpR, eval=FALSE-----------------------------------------
#  if (!require("devtools")) install.packages("devtools")
#  devtools::install_github("jeanlucj/AlphaSimHlpR")
#  suppressMessages(library(AlphaSimHlpR))

## ----Specify simulation, eval=FALSE-------------------------------------------
#  bsp <- specifyPopulation(ctrlFileName="PopulationCtrlFile.txt")
#  bsp <- specifyPipeline(bsp, ctrlFileName="ControlFile.txt")
#  bsp <- specifyCosts(bsp, ctrlFileName="CostsCtrlFile.txt")

## ----Toy example, eval=FALSE--------------------------------------------------
#  pathToSimDirectory <- path.package("AlphaSimHlpR")
#  setwd(pathToSimDirectory)
#  bsp <- specifyPopulation(ctrlFileName="PopulationCtrlFile_Small.txt")
#  bsp <- specifyPipeline(bsp, ctrlFileName="ControlFile_Small.txt")
#  bsp <- specifyCosts(bsp, ctrlFileName="CostsCtrlFile_Small.txt")

## ----One or more simulations, eval=FALSE--------------------------------------
#  simOut <- runBreedingScheme(nCycles=bsp$nCyclesToRun, initializeFunc=initializeScheme, productPipeline=productPipeline, populationImprovement=popImprov1Cyc, bsp=bsp)
#  
#  repSimOut <- lapply(1:nReplications, runBreedingScheme, nCycles=bsp$nCyclesToRun, initializeFunc=initializeScheme, productPipeline=productPipeline, populationImprovement=popImprov1Cyc, bsp)

## ----Plot replicated simulations, eval=FALSE----------------------------------
#  plotData <- plotRecords(repSimOut)

## ----optimizeByLOESS, eval=FALSE----------------------------------------------
#  optimizeByLOESS(batchSize=200, targetBudget=8000, percentRanges=percentRanges, startCycle=4, wgtPopImprov=0.5, tolerance=0.40, baseDir=baseDir, maxNumBatches=50, initializeFunc=initializeScheme, productPipeline=productPipeline, populationImprovement=popImprov1Cyc, bsp=bsp, randomSeed=rs, nCores=10)

