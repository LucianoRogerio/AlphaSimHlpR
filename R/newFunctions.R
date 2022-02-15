#' Run population improvement using parent selection
#'
#' Function to improve a simulated breeding population by one cycle.
#' This version is adjusted relative to the original \code{popImprov1Cyc} function,
#' which drew \code{candidates} from the full \code{records$F1} (excluding potentially
#' indivs only scored during the current year, if \code{useCurrentPhenoTrain=FALSE}).
#' My changes:
#' \itemize{
#'  \item \code{nTrainPopCycles}: draw additional training pop clones only from this number of recent cycles.
#'  \item \code{nYrsAsCandidates}: candidates for selection only from this number of recent years
#'  \item \code{maxTrainingPopSize}: From the lines in the most recent cycles (indicated by \code{nTrainPopCycles}),
#'  subsample this number of lines for training data. This is \emph{in addition to} the "check" (\code{bsp$checks@id})
#'  and the lines indicates as selection `candidates` according to the setting of `nYrsAsCandidates`.
#'  All "historical" data will always be used, but the number of maximum training lines will be held constant.
#'  Replaces the stage-specific `bsp$trainingPopCycles`, which will be unused in this pipeline, but not deleted from the package.
#'  \item The selection criteria \code{parentSelCritBLUP} and \code{parentSelCritGEBV} will allow both `candidates` and the additional training pop lines to be selected for crossing. Previous version only allowed selection of the `candidates`.
#'  }
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param SP The AlphaSimR SimParam object
#' @return A records object with a new F1 Pop-class object of progeny coming out of a population improvement scheme
#'
#' @details This function uses penotypic records coming out of the product pipeline to choose individuals as parents to initiate the next breeding cycle
#' @export
popImprovByParentSel <- function(records, bsp, SP){
  # Which phenotypes can be included for model training?
  ### Current year phenotypes?
  trainRec <- records
  if (!bsp$useCurrentPhenoTrain){
    for (stage in bsp$stageNames){
      trainRec[[stage]] <- trainRec[[stage]][-length(trainRec[[stage]])]
    }
  }

  # Which individuals can be selection candidates?
  ## only individuals that have been genotyped in the last "nYrsAsCandidates"
  if(bsp$stageToGenotype=="F1"){
    NrecentProgenySelCands <- (bsp$nProgeny * bsp$nCrosses)
    candidates <- records$F1@id %>% tail(., n = NrecentProgenySelCands)
    if (bsp$nYrsAsCandidates > 1) {
      for(i in bsp$stageNames[1:(bsp$nYrsAsCandidates-1)]) {
        candidates <- c(candidates, (records[[i]] %>% tail(., n = 1) %>% .[[1]] %$% unique(id) %>% setdiff(., bsp$checks@id)))
      }
    }
    candidates <- unique(candidates)
  } else {
    candidates <- records[[bsp$stageToGenotype]] %>%
      # select the "id" clones from the stage that the clones are genotyped
      tail(., n =1) %>% .[[1]] %$% unique(id) %>%
      # exclude checks
      setdiff(., bsp$checks@id)
    # get the progenitor candidates of the advanced trials
    if (bsp$nYrsAsCandidates > 1) {
      for(i in bsp$stageNames[((match(bsp$stageToGenotype,
                                      bsp$stageNames) + 1) :
                               bsp$nStages)[1:(bsp$nYrsAsCandidates-1)]]) {
        candidates <- c(candidates, (records[[i]] %>%
                                       # Selecting the last year of the stage "i"
                                       tail(., n = 1) %>%
                                       .[[1]] %$%
                                       # remove the duplicated "id" names
                                       unique(id) %>%
                                       # exclude checks
                                       setdiff(., bsp$checks@id)))
      }
    }
# Last change
#    candidates<-records[[bsp$stageToGenotype]] %>%
#      tail(.,n=bsp$nYrsAsCandidates) %>%
#      map_df(.,rbind) %$%
#      unique(id) %>%
#      # exclude checks
#      setdiff(.,bsp$checks@id)
  }

  # How many additional individuals to use as training?
  ## these are individuals with phenotypes
  ## but not in the list of selection candidates
  ## Drawn from the most recent cycles according to "nTrainPopCycles"
  ## Potentially subsampled according to "maxTrainingPopSize"
  ## RmStagePhen from the bsp object allows to remove trials with no accurate information from the phenotyped lines
  phenotypedLines<-trainRec[bsp$stageNames[!bsp$stageNames%in%bsp$RmStagePhen]] %>%
    map(.,~tail(.,n = bsp$nTrainPopCycles)) %>%
    map_df(.,rbind) %$%
    unique(id)

  phenotypedLines_notSelCands<-setdiff(phenotypedLines, c(candidates, bsp$checks@id)) %>% .[order(as.integer(.))]
  ## maxTPsize is lesser of specified 'maxTrainingPopSize' and actual number of phenotyped lines not considered selection candidates
  maxTPsize<-min(bsp$maxTrainingPopSize,length(phenotypedLines_notSelCands))
  ## Make sure checks ARE included

  # sample from the list of non-selection candidates that also are NOT checks
  if(!is.null(bsp$TrainingPopSel)){
    trainingpop<-bsp$TrainingPopSel(phenotypedLines_notSelCands)
  } else {
    trainingpop<-sample(phenotypedLines_notSelCands,
                        size = maxTPsize, replace = F) %>%
      .[order(as.integer(.))]
  }
      # include the checks if you have
  if(!is.null(bsp$checks)){
    trainingpop<-c(trainingpop,bsp$checks@id) %>%
      # vanity: order the ids
      .[order(as.integer(.))]
  }

  # require two inputs for downstream SelCrit
  ## only compatible SelCrit so far will therefore be "parentSelCritGEBV"
  ## "candidates" and "trainingpop": non-overlapping sets,
  ## available pheno records (in "trainRec") for any of the "candidates"
  ## will be automatically included in predictions
  crit <- bsp$selCritPopImprov(trainRec, candidates, trainingpop, bsp, SP)
  nParSam <- ifelse(bsp$NParentsFlw)

  # Not sure if useOptContrib will work "as is"
  if (bsp$useOptContrib){
    progeny <- optContrib(records, bsp, SP, crit)
  } else {
    # select the top nParents based
    selectedParentIDs<-names(crit[order(crit, decreasing=T)][1:bsp$nParents]) %>%
      sample(x = ., size = round(bsp$nParents*bsp$parentsFlowering/100, digits = 0), replace = FALSE) %>%
      .[order(as.integer(.))]
    # extract a pop-object of those parents
    parents <- records$F1[selectedParentIDs]
    # make crosses
    progeny <- randCross(parents, nCrosses=bsp$nCrosses, nProgeny=bsp$nProgeny, ignoreSexes=T, simParam=SP)
  }
  # not 100% sure, but seems to store the "year" in the @fixEff slot of "progeny"
  progeny@fixEff <- rep(as.integer(max(records$stageOutputs$year) + 1), bsp$nSeeds)
  # parentsUsed <- unique(c(progeny@mother, progeny@father))
  # stgCyc <- sapply(parentsUsed, AlphaSimHlpR:::whereIsID, records=records)
  # stgCyc <- table(stgCyc[1,], stgCyc[2,])
  # strtStgOut <- nrow(records$stageOutputs) - bsp$nStages - 1
  # for (i in 1:nrow(stgCyc)){
  #   stage <- as.integer(rownames(stgCyc)[i])
  #   records$stageOutputs$nContribToPar[[strtStgOut + stage]] <- tibble(cycle=as.integer(colnames(stgCyc)), nContribToPar=stgCyc[i,])
  # }
  GS <- tibble(Year = as.integer(max(records$stageOutputs$year, na.rm = TRUE)),
               cycle = as.integer(max(records$stageOutputs$cycle, na.rm = TRUE)),
               first = first(candidates),
               last = last(candidates),
               grmSize = length(union(candidates, trainingpop)),
               grmId = list(tibble(id = candidates) %>% mutate(pop = "c") %>%
                              bind_rows(tibble(id = trainingpop)) %>% mutate(pop = "t")),
               Ne = grmSize*(1/(1 + mean(diag(make_grm(records, union(candidates, trainingpop),
                                                       bsp, SP, grmType="add"))) - 1)),
               accAtSel=cor(gv(records$F1[setdiff(union(candidates, trainingpop), bsp$checks@id)]), crit),
               genoValMean = mean(gv(records$F1[setdiff(union(candidates, trainingpop), bsp$checks@id)])),
               genValSD = sd(gv(records$F1[setdiff(union(candidates, trainingpop), bsp$checks@id)])))
  records$F1 <- c(records$F1, progeny)
  records[["GS"]] <- rbind(records[["GS"]], GS)
  return(records)
}

#' Generate GEBV as genomic parent selection criteria
#'
#' \code{parentSelCritGEBV} will compute GEBV of all genotyped individuals using all available phenotypes.
#' Parents are selected by these criteria e.g. by the \code{popImprov1Cyc populationImprovement function}
#'
#' Setting up to distinguish between parent and cross selection and additive and non-additive predictions.
#' Uses \code{genomicMateSelectR} functions, which will need to be installed for these to work.
#'
#' Modified original \code{selCritGRM} function.
#' Uses alternative functions \code{make_grm} and \code{gebvPhenoEval}.
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param candidates Character vector of ids of the candidates to be parents, not necessarily phenotyped but must be predicted
#' @param trainingpop chr. vector of ids with phenotypes , but not necessarily in the list of selection candidates, but who should be included in the grm / training model.
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object (needed to pull SNPs)
#' @return Character vector of the ids of the selected individuals
#' @details Accesses all individuals in \code{records} to pick the highest ones
#'
#' @examples
#'
#' @export
parentSelCritGEBV <- function(records, candidates, trainingpop, bsp, SP){
  # first construct the GRM
  indivs2keep<-union(candidates,trainingpop)
  grm <- make_grm(records, indivs2keep, bsp, SP, grmType="add")
  phenoDF <- framePhenoRec(records, bsp)
  # Remove individuals with phenotypes but who do not have geno records
  phenoDF <- phenoDF[(phenoDF$id %in% rownames(grm) & !phenoDF$stage %in% bsp$RmStagePhen),]
  crit <- grmPhenoEval(phenoDF, grm)
  # exclude the checks from consideration as candidates
  crit <- crit[names(crit) %in% setdiff(indivs2keep,bsp$checks@id)]
  return(crit)
}



#' Function to make genomic relation matrices
#'
#' Function to make a genomic relationship matrix to be used to analyze the
#' phenotypic \code{records}. So far only works with \code{parentSelCritGEBV}
#' and \code{popImprovementByParentSel}.
#'
#' Setting up to distinguish between parent and cross selection and additive and non-additive predictions.
#' Uses \code{genomicMateSelectR} functions, which will need to be installed for these to work.
#'
#' This version uses the \code{genomicMateSelectR::kinship}.
#'
#' Includes a new \code{grmType} param, can produce additive and dominance matrices.
#' Same matrix as \code{A.mat} by default \code{grmType="add"}.
#'
#' @param records The breeding program \code{records} object. See
#'   \code{fillPipeline} for details
#' @param indivs2keep chr. vector of id's to be included in the grm. Must be present in \code{c(records$F1,bsp$checks)}.
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object. Needed to pull the SNP genotypes
#' @param grmType
#' @return A genomic relationship matrix
#' @details \code{records} maintains the phenotypic and genotypic records across
#'   years and stages. For GEBV analysis, you need the GRM of these individuals.
#'   \code{makeGRM} assumes the first phenotyping stage (records[[2]]) has all
#'   individuals that have been phenotyped. The GRM also includes the
#'   unphenotyped new F1 individuals in records[[1]]
#'
#' @examples

#' @export
make_grm <- function(records, indivs2keep, bsp, SP, grmType="add"){
  return(genomicMateSelectR::kinship(pullSnpGeno(c(records$F1,bsp$checks)[indivs2keep], simParam=SP), type = grmType)) }


#' Predict the GEBV of all genotyped individuals
#'
#' Function to predict the GEBV of all genotyped individuals in the population using all available phenotypes.
#' Takes a data.frame coming from \code{framePhenoRec} and
#' an additive kinship matrix (GRM) and analyze them with individuals .
#' Uses \code{\link[sommer]{mmer}} to implement a genomic BLUP model.
#'
#' @param phenoDF A data.frame of phenotypic observations. See \code{framePhenoRec} for details
#' @param grm A genomic relationship matrix
#' @return Named real vector of the BLUPs of all individuals in phenoDF (names are the individual ids), with appropriate weights by error variance of the observation
#' @details Given all the phenotypic records calculate the GEBV for each individual using all its records
#'
#' @examples
#'
#' @export
gebvPhenoEval <- function(phenoDF, grm){
    require(sommer)
    phenoDF$id <- factor(phenoDF$id, levels=rownames(grm)) # Enable prediction
    phenoDF$wgt <- 1/phenoDF$errVar # Make into weights
    fm <- mmer(fixed = pheno~1,
               random = ~vs(id, Gu=grm),
               weights=wgt,
               data=phenoDF,
               verbose=F,
               date.warning=F)
    gebv <- fm$U[[1]][[1]]
  # Ensure output has variation: needed for optimal contributions
  if (sd(gebv) == 0){
    namesGEBV <- names(gebv)
    gebv <- tapply(phenoDF$pheno, phenoDF$id, mean)
    names(gebv) <- namesGEBV
  }
  detach("package:sommer",unload = T); detach("package:MASS",unload = T)
  return(gebv)
}



#' Generate BLUPs as phenotype-only parent selection criteria
#'
#' \code{parentSelCritBLUP} will compute BLUPs of the \code{union(candidates,trainpop)} using all historical records but a controlled set of clones.
#' Only phenotyped individuals are predicted as no genomic covariances are used.
#' Parents are selected by these criteria e.g. by the \code{popImprovByParentSel populationImprovement function}
#'
#' Setting up to distinguish between parent and cross selection and additive and non-additive predictions.
#' Uses \code{genomicMateSelectR} functions, which will need to be installed for these to work.
#'
#' Modified original \code{selCritIID} function.
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param candidates Character vector of ids of the candidates to be parents, not necessarily phenotyped but must be predicted
#' @param trainingpop chr. vector of ids with phenotypes , but not necessarily in the list of selection candidates, but who should be included in the grm / training model.
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object (needed to pull SNPs)
#' @return Character vector of the ids of the selected individuals
#' @details Accesses all individuals in \code{records} to pick the highest ones
#'
#' @export
parentSelCritBLUP <- function(records, candidates, trainingpop, bsp, SP){
  indivs2keep<-union(candidates,trainingpop)
  phenoDF <- framePhenoRec(records, bsp)
  # Remove individuals not designated as candidates or trainingpop
  phenoDF <- phenoDF %>% filter(id %in% indivs2keep, !stage %in% bsp$RmStagePhen)
  crit <- iidPhenoEval(phenoDF)
  # exclude the checks from consideration as candidates
  crit <- crit[names(crit) %in% setdiff(indivs2keep,bsp$checks@id)]
  return(crit)
}

#' Generate predictions of cross means or usefulness as genomic mate selection criteria
#'
#' Mates are selected by these criteria e.g. by the \code{popImprovByMateSel populationImprovement function}
#' Additional arguments to add to the \code{bsp} object to accomplish this:
#'
#' \itemize{
#'  \item \code{modelType}:
#'  \item \code{crossSelCrit}:
#'  \item \code{propSel}:
#' }
#' Setting up to distinguish between parent and cross selection and additive and non-additive predictions.
#' Uses \code{genomicMateSelectR} functions, which will need to be installed for these to work.
#'
#' Modified original \code{selCritIID} function.
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param candidates Character vector of ids of the candidates to be parents, not necessarily phenotyped but must be predicted
#' @param trainingpop chr. vector of ids with phenotypes , but not necessarily in the list of selection candidates, but who should be included in the grm / training model.
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object (needed to pull SNPs)
#' @return Character vector of the ids of the selected individuals
#' @details Accesses all individuals in \code{records} to pick the highest ones
#'
#' @export
genomicMateSelCrit<-function(records, candidates, trainingpop, bsp, SP){

  # first construct the GRM(s)
  indivs2keep<-union(candidates,trainingpop)
  grms<-list(A=make_grm(records, indivs2keep, bsp, SP, grmType="add"))
  if(bsp$modelType=="DirDom"){
    grms[["D"]]<-make_grm(records, indivs2keep, bsp, SP, grmType="domGenotypic")
  }
  # phenotypes
  phenoDF <- framePhenoRec(records, bsp)
  # Remove individuals with phenotypes but who do not have geno records
  phenoDF <- phenoDF[phenoDF$id %in% rownames(grms$A),]
  # format the blups to please my own program (genomicMateSelectR)
  blups<-tibble(Trait="trait",
                TrainingData=list(phenoDF %>%
                                    rename(drgBLUP=pheno) %>%
                                    mutate(WT=1/errVar,
                                           GID=id)))
  # pull the dosage matrix, including checks
  dosages<-pullSnpGeno(c(records$F1,bsp$checks)[indivs2keep], simParam=SP)
  # run genomic predictions (to get marker effects and genomic BLUPs)
  gpreds<-runGenomicPredictions(modelType=bsp$modelType,selInd=FALSE, SIwts=NULL,
                                getMarkEffs=TRUE,
                                returnPEV=FALSE,
                                blups=blups,grms=grms,dosages=dosages,
                                ncores=1,nBLASthreads=nBLASthreads)
  # get the genetic map
  genmap<-getSnpMap(simParam = SP)
  m<-genmap$pos*100; # convert it to centimorgans
  names(m)<-genmap$id
  # construct the recombination frequency matrix
  recombFreqMat<-1-(2*genmap2recombfreq(m,nChr = bsp$nChr))
  # pull the haplotype matrix
  haploMat<-pullSnpHaplo(c(records$F1,bsp$checks)[indivs2keep], simParam=SP)
  # change haplotype tags in rownames of the haploMat to please genomicMateSelectR
  rownames(haploMat) %<>%
    gsub("_1","_HapA",.) %>%
    gsub("_2","_HapB",.)
  # set-up for each possible crossSelCrit
  if(bsp$crossSelCrit=="MeanBV"){
    predof<-"GEBV"; predTheMeans<-TRUE; predTheVars<-FALSE; }
  if(bsp$crossSelCrit=="MeanTGV"){
    predof<-"GETGV"; predTheMeans<-TRUE; predTheVars<-FALSE; }
  if(bsp$crossSelCrit=="UCparent"){
    predof<-"GEBV"; predTheMeans<-TRUE; predTheVars<-TRUE; }
  if(bsp$crossSelCrit=="UCvariety"){
    predof<-"GETGV"; predTheMeans<-TRUE; predTheVars<-TRUE; }
  # parents for which to predict crosses
  parents<-gpreds$gblups[[1]]  %>%
    filter(predOf==predof,
           # don't let the checks be considered as parents
           !GID %in% bsp$checks@id) %>%
    arrange(desc(trait)) %>%
    slice(1:bsp$nParents) %$%
    GID
  # crosses to predict (all pairwise of parents)
  CrossesToPredict<-crosses2predict(parents)
  # predict crosses
  crossPreds<-predictCrosses(modelType=bsp$modelType,
                             stdSelInt = intensity(bsp$propSel),
                             selInd=FALSE, SIwts=NULL,
                             CrossesToPredict=CrossesToPredict,
                             snpeffs=gpreds$genomicPredOut[[1]],
                             dosages=dosages,
                             haploMat=haploMat,recombFreqMat=recombFreqMat,
                             ncores=bsp$nCrossPredCores,nBLASthreads=nBLASthreads,
                             predTheMeans = predTheMeans,
                             predTheVars = predTheVars)
  # post prediction: extract the cross selection criterion
  if(bsp$crossSelCrit=="MeanBV"){
    crit<-crossPreds$tidyPreds[[1]] %>%
      filter(predOf=="BV") %>%
      select(sireID,damID,predMean) %>%
      rename(crossSelCrit=predMean) %>%
      arrange(desc(crossSelCrit))
  }
  if(bsp$crossSelCrit=="MeanTGV"){
    crit<-crossPreds$tidyPreds[[1]] %>%
      filter(predOf=="TGV") %>%
      select(sireID,damID,predMean) %>%
      rename(crossSelCrit=predMean) %>%
      arrange(desc(crossSelCrit))
  }
  if(bsp$crossSelCrit=="UCparent"){
    crit<-crossPreds$tidyPreds[[1]] %>%
      filter(predOf=="BV") %>%
      select(sireID,damID,predUsefulness) %>%
      rename(crossSelCrit=predUsefulness) %>%
      arrange(desc(crossSelCrit))
  }
  if(bsp$crossSelCrit=="UCvariety"){
    crit<-crossPreds$tidyPreds[[1]] %>%
      filter(predOf=="TGV") %>%
      select(sireID,damID,predUsefulness) %>%
      rename(crossSelCrit=predUsefulness) %>%
      arrange(desc(crossSelCrit))
  }

  return(crit)
}


#' Run population improvement using mate selection
#'
#' Function to improve a simulated breeding population by one cycle.
#' Also specify \code{selCritPop="genomicMateSelCrit"} and use the following new \code{bsp} arguments to control the predictions and selection:
#' \itemize{
#'  \item \code{modelType}: modelTypes: "A", "AD", "DirDom"
#'  \item \code{crossSelCrit}: "MeanBV" (modelTypes: "A", "AD", "DirDom"); "MeanTGV" (only with modelType="DirDom"); "UCparent" (modelTypes: "A", "AD", "DirDom"); "UCvariety" (modelTypes: "AD", "DirDom")
#'  \item \code{propSel}: 0-1, prop. of predicted crosses to select, used ONLY to calc the standardized selection intensity (i) and subsequently the usefulness criteria (\eqn{\mu + i \times \sigma}).
#' }
#'
#' My changes:
#' \itemize{
#'  \item \code{nTrainPopCycles}: draw training pop clones only from this number of recent cycles.
#'  \item \code{nYrsAsCandidates}: candidates for selection only from this number of recent years
#'  \item \code{maxTrainingPopSize}: From the lines in the most recent cycles (indicated by \code{nTrainPopCycles}),
#'  subsample this number of lines for training data. This is \emph{in addition to} the "check" (\code{bsp$checks@id})
#'  and the lines indicates as selection `candidates` according to the setting of `nYrsAsCandidates`.
#'  All "historical" data will always be used, but the number of maximum training lines will be held constant.
#'  Replaces the stage-specific `bsp$trainingPopCycles`, which will be unused in this pipeline, but not deleted from the package.
#'  }
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param SP The AlphaSimR SimParam object
#' @return A records object with a new F1 Pop-class object of progeny coming out of a population improvement scheme
#'
#' @details This function uses phenotypic records coming out of the product pipeline to choose individuals as parents to initiate the next breeding cycle
#' @export
popImprovByMateSel <- function(records, bsp, SP){
  require(genomicMateSelectR)

  # Which phenotypes can be included for model training?
  ### Current year phenotypes?
  trainRec <- records
  if (!bsp$useCurrentPhenoTrain){
    for (stage in bsp$stageNames){
      trainRec[[stage]] <- trainRec[[stage]][-length(trainRec[[stage]])]
    }
  }

  # Which individuals can be selection candidates?
  ## only individuals that have been genoytped in the last "nYrsAsCandidates"
  if(bsp$stageToGenotype=="F1"){
    NrecentProgenySelCands<-(bsp$nProgeny*bsp$nCrosses)*bsp$nYrsAsCandidates
    candidates<-records$F1@id %>% tail(.,n = NrecentProgenySelCands)
  } else {
    candidates<-records[[bsp$stageToGenotype]] %>%
      tail(.,n=bsp$nYrsAsCandidates) %>%
      map_df(.,bind_rows) %$%
      unique(id) %>%
      # exclude checks
      setdiff(.,bsp$checks@id)
  }

  # How many additional individuals to use as training?
  ## these are individuals with phenotypes
  ## but not in the list of selection candidates
  ## Drawn from the most recent cycles according to "nTrainPopCycles"
  ## Potentially subsampled according to "maxTrainingPopSize"
  phenotypedLines<-trainRec[bsp$stageNames] %>%
    map(.,~tail(.,n = bsp$nTrainPopCycles)) %>%
    map_df(.,bind_rows) %$%
    unique(id)

  phenotypedLines_notSelCands<-setdiff(phenotypedLines,candidates)
  ## maxTPsize is lesser of specified 'maxTrainingPopSize' and actual number of phenotyped lines not considered selection candidates
  maxTPsize<-min(bsp$maxTrainingPopSize,length(phenotypedLines_notSelCands))
  ## Make sure checks ARE included
  if(!is.null(bsp$checks)){
    # sample from the list of non-selection candidates that also are NOT checks
    trainingpop<-sample(setdiff(phenotypedLines_notSelCands,bsp$checks@id),
                        size = maxTPsize, replace = F) %>%
      # include the checks
      c(.,bsp$checks@id) %>%
      # vanity: order the ids
      .[order(as.integer(.))]
  } else {
    trainingpop<-sample(phenotypedLines_notSelCands,
                        size = maxTPsize, replace = F) %>%
      .[order(as.integer(.))]
  }

  # require two inputs for downstream SelCrit
  ## only compatible SelCrit so far will therefore be "parentSelCritGEBV"
  ## "candidates" and "trainingpop": non-overlapping sets,
  ## available pheno records (in "trainRec") for any of the "candidates"
  ## will be automatically included in predictions
  crit <- bsp$selCritPopImprov(trainRec, candidates, trainingpop, bsp, SP)

  # select the top nCrosses
  crossingPlan<-crit %>%
    slice_max(order_by = crossSelCrit,
              n = bsp$nCrosses,
              with_ties = F) %>%
    select(sireID,damID) %>%
    as.matrix

  # extract a pop-object of those parents
  parents <- records$F1[crossingPlan %>% as.vector %>% unique]
  # make crosses
  progeny <- makeCross(pop = parents,
                       crossPlan = crossingPlan,
                       nProgeny = bsp$nProgeny, simParam=SP)

  # not 100% sure, but seems to store the "year" in the @fixEff slot of "progeny"
  progeny@fixEff <- rep(as.integer(max(records$stageOutputs$year) + 1), bsp$nSeeds)
  #parentsUsed <- unique(c(progeny@mother, progeny@father))
  # stgCyc <- sapply(parentsUsed, AlphaSimHlpR:::whereIsID, records=records)
  # stgCyc <- table(stgCyc[1,], stgCyc[2,])
  # strtStgOut <- nrow(records$stageOutputs) - bsp$nStages - 1
  # for (i in 1:nrow(stgCyc)){
  #   stage <- as.integer(rownames(stgCyc)[i])
  #   records$stageOutputs$nContribToPar[[strtStgOut + stage]] <- tibble(cycle=as.integer(colnames(stgCyc)), nContribToPar=stgCyc[i,])
  # }
  records$F1 <- c(records$F1, progeny)
  return(records)
}
