#' initFuncSimp function
#'
#' Deprecated function to initialize simulation of a breeding program
#'
#' @param bsp A list of breeding scheme parameters
#' @return The return of \code{initializeScheme}
#'
#' @details Deprecated in favor of the simply-named \code{initializeScheme}
#'
#' @export
initFuncSimp <- function(bsp){
  print("initFuncSimp deprecated. Please use initializeScheme")
  return(initializeScheme(bsp))
}

#' initFuncADChk function
#'
#' Deprecated function to initialize simulation of a breeding program
#'
#' @param bsp A list of breeding scheme parameters
#' @return The return of \code{initializeScheme}
#'
#' @details Deprecated in favor of the simply-named \code{initializeScheme}
#'
#' @export
initFuncADChk <- function(bsp){
  print("initFuncADChk deprecated. Please use initializeScheme")
  return(initializeScheme(bsp))
}

#' initializeScheme function
#'
#' function to initialize simulation of a breeding program. A single additive-dominance trait is simulated. Check are used in this scheme
#'
#' @param bsp A list of breeding scheme parameters.  See \code{specifyPipeline} and \code{specifyPopulation}
#' @param nThreadsForMacs uses the nThreads argument in \code{runMacs2}, parallelizes founder sim by chrom.
#' @return A list containing: 1. The simulation parameters in \code{SP}; 2. The initial records of the breeding program in \code{records}. See \code{fillPipeline} for details; 3. A completed \code{bsp} object
#'
#' @details Creates the founders and the initial records at the beginning of the simulation of a breeding program.
#'
#' @examples
#' bsp <- specifyPopulation(bsp)
#' bsp <- specifyPipeline()
#' initList <- initializeScheme(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#'
#' @export
initializeScheme <- function(bsp,nThreadsForMacs = NULL){
  # Create haplotypes for founder population of outbred individuals
  nF1 <- bsp$nCrosses * bsp$nProgeny + max(bsp$nChks)
  if (bsp$quickHaplo){
    founderHap <- quickHaplo(nInd=nF1, nChr=bsp$nChr, segSites=bsp$segSites)
  } else{
    founderHap <- runMacs2(nInd=nF1, nChr=bsp$nChr, segSites=bsp$segSites, Ne=bsp$effPopSize,nThreads = nThreadsForMacs)
  }

  # New global simulation parameters from founder haplotypes
  SP <- SimParam$new(founderHap)
  SP$restrSegSites(minQtlPerChr=1, minSnpPerChr=10, overlap=FALSE)
  # Additive, dominance, and epistatic trait architecture
  SP$addTraitADE(nQtlPerChr=bsp$nQTL, var=bsp$genVar, meanDD=bsp$meanDD, varDD=bsp$varDD, relAA=bsp$relAA, useVarA=FALSE)
  # Observed SNPs per chromosome
  SP$addSnpChip(bsp$nSNP)

  founders <- newPop(founderHap, simParam=SP)
  if (any(bsp$nChks > 0)){
    bsp$checks <- selectInd(founders, nInd=max(bsp$nChks), use="rand", simParam=SP)
    # remove checks from founders
    founders <- founders[-which(founders@id %in% bsp$checks@id)]
  } else bsp$checks <- NULL

  records <- fillPipeline(founders, bsp, SP)

  return(list(SP=SP, records=records, bsp=bsp))
}

#' fillPipeline function
#'
#' function to create initial records at the start of a simulation
#'
#' @param founders Pop-class object of the founders of the breeding program
#' @param bsp A list of product pipeline parameters. See \code{runBreedingScheme} for details
#' @return A \code{records} object. A list of lists containing nStages+1 lists. The first list contains one Pop-class of progeny per year of the scheme. The remaining lists contain one matrix per year that has individual id, mother, father, stage, phenotypes, and error variances. The individuals have been phenotyped using \code{setPheno}. The matrix may contain a mix of experimental and check phenotypes with different levels of replication
#'
#' @details This is a structure for a records object that will be used to simulate breeding schemes
#'
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' nF1 <- bsp$nCrosses * bsp$nProgeny
#' founderHap <- runMacs(nInd=nF1, nChr=bsp$nChr, segSites=bsp$segSites)
#' SP <- SimParam$new(founderHap)
#' SP$addTraitA(nQtlPerChr=bsp$nQTL, var=bsp$genVar)
#' SP$addSnpChip(bsp$nSNP)
#' founders <- newPop(founderHap, simParam=SP)
#' bsp <- c(bsp, checks=list(NULL))
#' records <- fillPipeline(founders, bsp, SP)
#'
#' @export
fillPipeline <- function(founders, bsp=NULL, SP){
  nF1 <- bsp$nCrosses * bsp$nProgeny
  records <- list(founders)
  for (year in 1 + -(bsp$nStages:1)){
    toAdd <- list()
    for (stage in 1:(year+bsp$nStages)){
      if (stage==1){ # Stage 1: F1 progeny population: random selection use pop
        # Select from the most recent F1s
        indToAdv <- nInd(records[[1]]) - nF1 + sort(sample(nF1, bsp$nEntries[stage]))
      } else{
        # Don't allow checks to be advanced: use 1:bsp$nEntries[stage-1]
        sourcePop <- last(records[[stage]])[1:bsp$nEntries[stage-1],]
        indToAdv <- order(sourcePop$pheno, decreasing=T)[1:bsp$nEntries[stage]]
        indToAdv <- sourcePop$id[sort(indToAdv)]
      }
      entries <- records[[1]][indToAdv]
      varE <- bsp$gxyVar + (bsp$gxlVar + bsp$gxyxlVar + bsp$errVars[stage] / bsp$nReps[stage]) / bsp$nLocs[stage]
      # reps=1 because varE is computed above
      entries <- setPheno(entries, varE=varE, reps=1, simParam=SP)
      phenoRec <- phenoRecFromPop(entries, bsp, stage)
      if(!is.null(bsp$checks) & bsp$nChks[stage] > 0){
        varE <- bsp$gxyVar + (bsp$gxlVar + bsp$gxyxlVar + bsp$errVars[stage] / bsp$chkReps[stage]) / bsp$nLocs[stage]
        chkPheno <- setPheno(bsp$checks[1:bsp$nChks[stage]], varE=varE, reps=1, simParam=SP)
        chkRec <- phenoRecFromPop(chkPheno, bsp, stage, checks=T)
        phenoRec <- bind_rows(phenoRec, chkRec)
      }
      toAdd <- c(toAdd, list(phenoRec))
    }#END stages

    # Make the next F1s with mild selection using gv
    lastGen <- nInd(records[[1]]) - nF1 + 1:nF1
    parents <- selectInd(records[[1]][lastGen], nInd=nF1/1.5, use="gv", simParam=SP)
    toAdd <- c(list(randCross(parents, nCrosses = bsp$nCrosses, nProgeny = bsp$nProgeny, ignoreSexes = T, simParam=SP)), toAdd)

    # Actually fill the records
    records[[1]] <- c(records[[1]], toAdd[[1]])
    for (i in 2:length(toAdd)){
      if (i > length(records)){
        records <- c(records, list(toAdd[i]))
      } else{
        records[[i]] <- c(records[[i]], toAdd[i])
      }
    }
  }#END years
  names(records) <- c("F1", bsp$stageNames)
  # stageOutputs relies on knowing the year from the previous year
  return(c(records, stageOutputs=list(tibble(year=-1))))
}

#' Run burn-in breeding scheme simulations
#'
#' Allows users to run simulation and then continue again later. Output is direct input for \code{runSchemesPostBurnIn}.
#' Runs potentially multiple replications and optionally in parallel.
#'
#' @param nReplications Integer number of replications of the specific breeding scheme to run
#' @param nSimCores Integer, number of cores to optionally execute replicate simulations in parallel
#' @param bsp  A list of breeding scheme parameters.
#' @param nBurnInCycles Integer number of cycles to as 'burn-in' using the \code{selCritPop} and \code{selCritPipe} settings.
#' @param iniFunc string, Function to initialize the breeding program.
#' @param productFunc string, Function to advance the product pipeline by one generation
#' @param popImprovFunc string, Function to improve the breeding population and select parents to initiate the next cycle of the breeding scheme
#' @param nBLASthreads number of cores for each worker to use for multi-thread BLAS. Will speed up, for example, genomic predictions when using selCritGRM. Careful to balance with other forms of parallel processing.
#' @param nThreadsMacs2 uses the nThreads argument in \code{runMacs2}, parallelizes founder sim by chrom.
#' @param selCritPop string, overrides the selCrit in \code{bsp} for the burn-in stage.
#' @param selCritPipe string, overrides the selCrit in \code{bsp} for the burn-in stage.
#' @return A \code{records} object containing the phenotypic records retained of the breeding scheme
#'
#' @details A wrapper to initiate the breeding program then iterate cycles of product pipeline and population improvement
#'
#' @export
runBurnInSchemes<-function(bsp,
                           nBurnInCycles,
                           selCritPop="selCritIID",
                           selCritPipe="selCritIID",
                           iniFunc="initializeScheme",
                           productFunc="productPipeline",
                           popImprovFunc="popImprov1Cyc",
                           TrainingPopSel=NULL,
                           nReplications=1,nSimCores=1,
                           nBLASthreads=NULL,nThreadsMacs2=NULL){

  require(furrr); plan(multisession, workers = nSimCores)
  require(progressr)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")
  with_progress({
    p <- progressor(steps = nReplications)

    simulations<-tibble(SimRep=1:nReplications) %>%
      mutate(burnInSim=future_map(SimRep,function(SimRep,...){
        if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }
        cat("******", SimRep, "\n")

        # This initiates the founding population
        bsp[["initializeFunc"]] <- get(iniFunc)
        bsp[["productPipeline"]] <- get(productFunc)
        bsp[["populationImprovement"]] <- get(popImprovFunc)
        if(!is.null(TrainingPopSel)) {bsp[["TrainingPopSel"]] <- get(TrainingPopSel)}
        if(bsp$parentsFlowering > 100 | bsp$parentsFlowering <= 0) stop("parent flowering ratio should be between 1 to 100")
        initList <- bsp$initializeFunc(bsp,nThreadsForMacs=nThreadsMacs2)
        SP <- initList$SP
        bsp <- initList$bsp
        records <- initList$records
        ## set the selection criteria for burn-in
        bsp[["selCritPipeAdv"]] <- get(selCritPipe)
        bsp[["selCritPopImprov"]] <- get(selCritPop)

        # Burn-in cycles
        cat("\n"); cat("Burn-in cycles"); cat("\n")
        for (cycle in 1:nBurnInCycles){
          cat(cycle, " ")
          records <- bsp$productPipeline(records, bsp, SP)
          records <- bsp$populationImprovement(records, bsp, SP)
        }
        p()
        return(list(records=records,
                    bsp=bsp,
                    SP=SP))
      },
      bsp=bsp,
      nBurnInCycles=nBurnInCycles,
      selCritPop=selCritPop,
      selCritPipe=selCritPipe,
      iniFunc=iniFunc,
      productFunc=productFunc,
      popImprovFunc=popImprovFunc,
      nBLASthreads=nBLASthreads,
      nThreadsMacs2=nThreadsMacs2,
      TrainingPopSel=TrainingPopSel,
      p=p))
  })
  plan(sequential)
  return(simulations)
}

#' Continue running breeding schemes post burn-in
#'
#' Allows users to continue a simulation, potentially initiated or 'burned-in' with \code{runBurnInSchemes}.
#' Allows user to optionally change the \code{bsp} and selection criterion.
#' Input \code{simulations} are the output of e.g. \code{runBurnInSchemes}.
#' Runs potentially multiple replications and optionally in parallel.
#'
#' @param simulations tibble, each row is a simulation, 2 columns: SimRep and burnInSim.
#' SimRep is an identifier. burnInSim is a list with 3 named elements:
#' "records", "bsp" and "SP"
#' @param nSimCores Integer, number of cores to optionally execute replicate simulations in parallel
#' @param newBSP optional, so you can specify a different bsp for post-burn in sims.
#' @param nPostBurnInCycles Integer number of cycles to run the \code{selCritPop} and \code{selCritPipe} settings.
#' @param productFunc string, Function to advance the product pipeline by one generation
#' @param popImprovFunc string, Function to improve the breeding population and select parents to initiate the next cycle of the breeding scheme
#' @param nBLASthreads number of cores for each worker to use for multi-thread BLAS. Will speed up, for example, genomic predictions when using selCritGRM. Careful to balance with other forms of parallel processing.
#' @param selCritPop string, overrides the selCrit in \code{bsp} for the post burn-in stage.
#' @param selCritPipe string, overrides the selCrit in \code{bsp} for the post burn-in stage.
#' @return A \code{records} object containing the phenotypic records retained of the breeding scheme
#'
#' @details A wrapper to initiate the breeding program then iterate cycles of product pipeline and population improvement
#'
#' @export
runSchemesPostBurnIn<-function(simulations,
                               newBSP=NULL, # so you can change the scheme entirely after burn-in
                               nPostBurnInCycles,
                               selCritPop="parentSelCritBLUP",
                               selCritPipe="productSelCritBLUP",
                               productFunc="productPipelinePostBurnIn",
                               popImprovFunc="popImprovByParentSel",
                               nSimCores=1,
                               nBLASthreads=NULL){

  require(furrr); plan(multisession, workers = nSimCores)
  require(progressr)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")
  with_progress({
    p <- progressor(steps = length(simulations$SimRep))

    simulations<-simulations %>%
      dplyr::mutate(SimOutput=future_map2(SimRep,burnInSim,function(SimRep,burnInSim,...){
        # debug
        # burnInSim<-simulations$burnInSim[[1]]
        # SimRep <- 1
        if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }
        cat("\n", "\n", "******", SimRep, "\n")

        # This CONTINUES where previous sims left off
        ## no initialize step
        ## Keep burn-in stage sim params "SP"
        SP<-burnInSim$SP
        ## specify a potentially new bsp object
        ## (keep checks stored in burn-in stage's bsp)
        if(!is.null(newBSP)){
          bsp<-newBSP; bsp$checks<-burnInSim$bsp$checks
          bsp[["burnInBSP"]]<-burnInSim$bsp
          # years are indexed starting at year==0,
          ## so for 10 burn-in cycles, max value should be 9, store for later
          bsp[["maxYearBurnInStage"]]<-max(burnInSim$records$stageOutputs$year)

          if(any(!bsp[["burnInBSP"]]$stageNames %in% bsp$stageNames)){
            rmstage <- bsp[["burnInBSP"]]$stageNames[!bsp[["burnInBSP"]]$stageNames %in% bsp$stageNames]
            burnInSim$records[rmstage] <- NULL
          }
        } else { bsp<-burnInSim$bsp }
        ## 'historical' records from burn-in
        records<-burnInSim$records
        ## override burn-in specified product and population improvement funcs
        bsp[["productPipeline"]] <- get(productFunc)
        bsp[["populationImprovement"]] <- get(popImprovFunc)
        bsp[["selCritPipeAdv"]] <- get(selCritPipe)
        bsp[["selCritPopImprov"]] <- get(selCritPop)
        if(!is.null(bsp$TrainPopSel)) {bsp[["TrainingPopSel"]] <- get(bsp$TrainPopSel)}
        if(bsp$parentsFlowering > 100 | bsp$parentsFlowering <= 0) stop("parent flowering ratio should be between 1 to 100")

        # Post burn-in cycles
        cat("Post burn-in cycles"); cat("\n")
        for (cycle in 1:nPostBurnInCycles){
          cat(cycle, " ")
          records <- bsp$productPipeline(records, bsp, SP)
          records <- bsp$populationImprovement(records, bsp, SP)
        }

        p()

        # Finalize the stageOutputs
        records <- AlphaSimHlpR:::lastCycStgOut(records, bsp, SP)

        return(list(records=records,
                    bsp=bsp,
                    SP=SP))
      },
      nPostBurnInCycles=nPostBurnInCycles,
      selCritPop=selCritPop,
      selCritPipe=selCritPipe,
      productFunc=productFunc,
      popImprovFunc=popImprovFunc,
      nBLASthreads=nBLASthreads,
      newBSP=newBSP,
      p=p))
  })
  plan(sequential)
  return(simulations)
}
