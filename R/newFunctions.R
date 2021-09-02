#' Run population improvement using parent selection
#'
#' Function to improve a simulated breeding population by one cycle.
#' This version is adjusted relative to the original \code{popImprov1Cyc} function,
#' which drew \code{candidates} from the full \code{records$F1} (excluding potentially
#' indivs only scored during the current year, if \code{useCurrentPhenoTrain=FALSE}).
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
#' @details This function uses penotypic records coming out of the product pipeline to choose individuals as parents to initiate the next breeding cycle
#' @export
popImprovByParentSel <- function(records, bsp, SP){
  # Which phenotypes can be included for model training?
  ### Current year phenotypes?
  trainRec <- records
  if (!bsp$useCurrentPhenoTrain){
    for (stage in 1+1:bsp$nStages){
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
      map_df(.,rbind) %$%
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
    map_df(.,rbind) %$%
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

  # Not sure if useOptContrib will work "as is"
  if (bsp$useOptContrib){
    progeny <- optContrib(records, bsp, SP, crit)
  } else {
    # select the top nParents based
    selectedParentIDs<-names(crit[order(crit, decreasing=T)][1:bsp$nParents])
    # extract a pop-object of those parents
    parents <- records$F1[selectedParentIDs]
    # make crosses
    progeny <- randCross(parents, nCrosses=bsp$nCrosses, nProgeny=bsp$nProgeny, ignoreSexes=T, simParam=SP)
  }
  # not 100% sure, but seems to store the "year" in the @fixEff slot of "progeny"
  progeny@fixEff <- rep(as.integer(max(records$stageOutputs$year) + 1), bsp$nSeeds)
  parentsUsed <- unique(c(progeny@mother, progeny@father))
  stgCyc <- sapply(parentsUsed, AlphaSimHlpR:::whereIsID, records=records)
  stgCyc <- table(stgCyc[1,], stgCyc[2,])
  strtStgOut <- nrow(records$stageOutputs) - bsp$nStages - 1
  for (i in 1:nrow(stgCyc)){
    stage <- as.integer(rownames(stgCyc)[i])
    records$stageOutputs$nContribToPar[[strtStgOut + stage]] <- tibble(cycle=as.integer(colnames(stgCyc)), nContribToPar=stgCyc[i,])
  }
  records$F1 <- c(records$F1, progeny)
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
  phenoDF <- phenoDF[phenoDF$id %in% rownames(grm),]
  crit <- gebvPhenoEval(phenoDF, grm)
  crit <- crit[candidates]
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
  phenoDF <- phenoDF %>% filter(id %in% indivs2keep)
  crit <- iidPhenoEval(phenoDF)
  crit <- crit[names(crit) %in% candidates]
  return(crit)
}


