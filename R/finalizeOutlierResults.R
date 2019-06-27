source('./R/calculateOutlier.R')
source('./R/calculateExpressionFilter.R')
source('./R/calculatePercentageFilter.R')

prepareOutlierMatrixForPromoters <- function(aliquotsPerGroup.byHistotype, idConversionTable, promoterExpression, totalExpression, outlierList) {
  sample.all <- colnames(promoterExpression)[-(1:ncol(idConversionTable))]
  expressedPromoters.all <- sapply(names(aliquotsPerGroup.byHistotype), function(currentGroup) {as.character(unique(totalExpression$promoterId[totalExpression[[paste0(currentGroup, '.isPromoterExpressed')]]]))})
  names(expressedPromoters.all) <- names(aliquotsPerGroup.byHistotype)
  
  promoterId.all <- unique(as.character(idConversionTable$promoterId))
  idTable <- idConversionTable %>% dplyr::select(one_of('promoterId', 'geneId')) %>% distinct()
  outlierMatrix.final <- matrix(data = rep(NaN, length(promoterId.all) * length(sample.all)), nrow = length(promoterId.all), dimnames = list(promoterId.all, sample.all))
  
  print('Prepare promoter outlier matrix...')
  for (currentGroup in names(aliquotsPerGroup.byHistotype)) {
    print(currentGroup)
    
    aliquots <- aliquotsPerGroup.byHistotype[[currentGroup]]
    expressedPromoters.logical <- promoterId.all %in% expressedPromoters.all[[currentGroup]]
    promoterIds <- as.vector(promoterId.all[expressedPromoters.logical])
    
    outlierList.currentGroup <- outlierList[[currentGroup]]
    
    # Create outlier matrix
    outlierMatrix <- t(sapply(outlierList.currentGroup, function(o) {o$outlierState}))
    dimnames(outlierMatrix) <- list(promoterIds, aliquots)
    outlierMatrix.final[dimnames(outlierMatrix)[[1]], dimnames(outlierMatrix)[[2]]] <- outlierMatrix
  }
  
  outlierMatrix.final <- cbind(idTable, outlierMatrix.final)
  return(outlierMatrix.final)
}   
#########################################
#########################################

preparePercentageFilterMatrixForPromoters <- function(aliquotsPerGroup.byHistotype, idConversionTable, promoterExpression, percentageFilterList) {
  sample.all <- colnames(promoterExpression)[-(1:ncol(idConversionTable))]
  promoterId.all <- unique(as.character(idConversionTable$promoterId))
  idTable <- idConversionTable %>% dplyr::select(one_of('promoterId', 'geneId')) %>% distinct()
  
  print('Prepare promoter percentage filter matrix...')
  filterMatrix.final <- matrix(data = rep(NaN, length(promoterId.all) * length(sample.all)), nrow = length(promoterId.all), dimnames = list(promoterId.all, sample.all))
  for (currentGroup in names(aliquotsPerGroup.byHistotype)) {
    print(currentGroup)
    # Load filter matrix
    filterMatrix <- percentageFilterList[[currentGroup]]
    filterMatrix.final[dimnames(filterMatrix)[[1]], dimnames(filterMatrix)[[2]]] <- filterMatrix
  }
  
  filterMatrix.final <- cbind(idTable, filterMatrix.final)
  return(filterMatrix.final)
}

###############################################################
### Outlier analysis - Element Table to Gene Centric Table  ###
###############################################################

filterPromoterOutlierTable <- function(outlier, expressionFilter, percentageFilter, promoterRanges.all, exonReducedRanges) {
  print('Filter promoter outlier based on expression and percentage filters...')
  # print('Find expressed promoter index')
  isExpressed <- !(rowSums(is.na(outlier[,-c(1,2)])) == (dim(outlier)[2] - 2))
  
  # print('Remove internal promoters')
  promoterRanges.all <- promoterRanges.all[match(outlier$promoterId, promoterRanges.all$promoterId)]
  promoterOverlaps <- findOverlaps(promoterRanges.all, exonReducedRanges)
  nonInternalPromoters.logical <- exonReducedRanges$MaxMergedExonRank[subjectHits(promoterOverlaps)] == 1
  
  # print('Use expression filter')
  outlier[expressionFilter == FALSE] <- 0
  print('Use percentage filter')
  outlier[percentageFilter == 0] <- 0
  print('Get only the expressed noninternal promoter outliers')
  outlier <- outlier[isExpressed & nonInternalPromoters.logical,]
  
  outlier.promoter <- as.data.frame(matrix(NaN, nrow = dim(expressionFilter)[1], ncol = dim(expressionFilter)[2]))
  colnames(outlier.promoter) <- colnames(expressionFilter)
  outlier.promoter[, 1:2] <- expressionFilter[, 1:2]
  outlier.promoter[match(as.character(outlier$promoterId), as.character(outlier.promoter$promoterId)), 3:(dim(outlier.promoter)[2])] <- outlier[, 3:(dim(outlier)[2])]
  return(outlier.promoter)
}

prepareOutlierResultsForPromoter <- function(metadata.pcawg, idConversionTable, promoterExpression_raw, promoterExpression_rel, totalExpression, promoterRanges, exonReducedRanges) {
  aliquotsPerGroup.byHistotype <- split(metadata.pcawg$aliquot_id, metadata.pcawg$histology_abbreviation)
  
  expressionFilter <- prepareExpressionFilter(aliquotsPerGroup.byHistotype, idConversionTable, promoterExpression_raw)
  outlierList <- prepareOutlierListForPromoter(aliquotsPerGroup.byHistotype, totalExpression, promoterExpression = promoterExpression_rel)
  percentageFilterList <- preparePercentageFilterForPromoter(aliquotsPerGroup.byHistotype, totalExpression, promoterExpression = promoterExpression_rel, outlierList)
  exonReducedRanges <- prepareReducedExonRanges(txdb)

  outlierOriginal <- prepareOutlierMatrixForPromoters(aliquotsPerGroup.byHistotype, idConversionTable, promoterExpression = promoterExpression_raw, totalExpression, outlierList)
  rm(outlierList)
  gc(verbose = FALSE)
  percentageFilter <- preparePercentageFilterMatrixForPromoters(aliquotsPerGroup.byHistotype, idConversionTable, promoterExpression = promoterExpression_raw, percentageFilterList)
  rm(percentageFilterList)
  gc(verbose = FALSE)
  
  promoterOutliers.filtered <- filterPromoterOutlierTable(outlier = outlierOriginal, expressionFilter, percentageFilter, promoterRanges.all = promoterRanges, exonReducedRanges)
  return(promoterOutliers.filtered)
}

getOutlierState <- function(os) {
    state <- NULL
    n <- length(os)
    if (all(is.na(os))) {
        state <- NaN
    } else {
        os <- os[!is.na(os)]
        if (all(os == 0)) {
            state <- 0
        } else {
            state <- 1
        }
    }
    return(rep(state, n))
}

prepareOutlierResultsForGene <- function(outlier.promoter) {
  outlier.gene <- outlier.promoter %>% group_by(geneId) %>% mutate_each(funs(getOutlierState), -one_of('promoterId', 'geneId')) %>% ungroup() %>% dplyr::select(-promoterId) %>% distinct()
  outlier.gene <- as.data.frame(outlier.gene)
  return(outlier.gene)
}
