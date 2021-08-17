#' parentSelCritGEBV function
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
#' @param candidates Character vector of ids of the candidates to be parents
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object (needed to pull SNPs)
#' @return Character vector of the ids of the selected individuals
#' @details Accesses all individuals in \code{records} to pick the highest ones
#' @examples
#' candidates <- records[[1]][[1]]@id
#' parents <- records[[1]][[1]][selCritGRM(records, candidates, bsp, SP)]
#'
#' @export
parentSelCritGEBV <- function(records, candidates, bsp, SP){
  grm <- make_grm(records, bsp, SP, grmType="add")
  if (!any(candidates %in% rownames(grm))){
    crit <- runif(length(candidates))
  } else {
    phenoDF <- framePhenoRec(records, bsp)
    # Remove individuals with phenotypes but who do not have geno records
    phenoDF <- phenoDF[phenoDF$id %in% rownames(grm),]
    crit <- gebvPhenoEval(phenoDF, grm)
    crit <- crit[candidates]
  }
  names(crit) <- candidates
  return(crit)
}


#' Function to make genomic relation matrices
#'
#' Function to make a genomic relationship matrix to be used to analyze the
#' phenotypic \code{records}
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
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object. Needed to pull the SNP genotypes
#' @param grmType
#' @return A genomic relationship matrix
#' @details \code{records} maintains the phenotypic and genotypic records across
#'   years and stages. For GEBV analysis, you need the GRM of these individuals.
#'   \code{makeGRM} assumes the first phenotyping stage (records[[2]]) has all
#'   individuals that have been phenotyped. The GRM also includes the
#'   unphenotyped new F1 individuals in records[[1]]
#' @export
make_grm <- function(records, bsp, SP, grmType="add"){

  allPop <- records$F1
  # Only make GRM of individuals that are specified in the TP
  allID <- NULL
  tpc <- bsp$trainingPopCycles[1]
  if (tpc){
    lastGen <- max(records$F1@fixEff)
    tpc <- min(lastGen, tpc)
    allID <- records$F1@id[records$F1@fixEff %in% lastGen + (1 - tpc):0]
  }
  for (stageNum in 1 + 1:bsp$nStages){
    tpc <- bsp$trainingPopCycles[stageNum]
    if (tpc){
      lastGen <- length(records[[stageNum]])
      tpc <- min(lastGen, tpc)
      for (cyc in lastGen + (1 - tpc):0){
        allID <- c(allID, records[[stageNum]][[cyc]]$id)
      }
    }
  }
  allID <- unique(allID)
  if (!is.null(bsp$checks)) allID <- setdiff(allID, bsp$checks@id)
  allID <- allID[order(as.integer(allID))]
  allPop <- allPop[allID]

  if (!is.null(bsp$checks)){
    putInChks <- setdiff(bsp$checks@id, allPop@id)
    if (length(putInChks > 0)) allPop <- c(allPop, bsp$checks[putInChks])
  }
  return(genomicMateSelectR::kinship(pullSnpGeno(allPop, simParam=SP), type = grmType))
}

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
#' phenoDF <- framePhenoRec(records, bsp)
#' grm <- makeGRM(records, bsp, SP)
#' grmBLUPs <- gebvPhenoEval(phenoDF, grm)
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
  return(gebv)
}

