
calculatePercentageFilter <- function(currentGroup, aliquots.currentGroup, elementIds.currentGroup,  outlierExpression, outlierList) {
  referenceMean <- sapply(outlierList, function(o) {o$referenceMean})
  # Create filter matrix
  filterMatrix <- matrix(data = rep(FALSE, (ncol(outlierExpression) * nrow(outlierExpression))), nrow = nrow(outlierExpression))
  for (j in 1:nrow(outlierExpression)) {
      filterMatrix[j, which(abs((unlist(outlierExpression[j,])) - referenceMean[j]) > 0.2)] <- TRUE
  }
  dimnames(filterMatrix) <- list(elementIds.currentGroup, aliquots.currentGroup)
  return(filterMatrix)
}

preparePercentageFilterForPromoter <- function(aliquotsPerGroup.byHistotype, totalExpression, promoterExpression, outlierList) {
  expressedPromoters.all <- sapply(names(aliquotsPerGroup.byHistotype), function(currentGroup) {as.character(unique(totalExpression$promoterId[totalExpression[[paste0(currentGroup,'.isPromoterExpressed')]]]))})
  names(expressedPromoters.all) <- names(aliquotsPerGroup.byHistotype)
  
  promoterExpression <- promoterExpression[, -(1:2)] %>% group_by(promoterId) %>% distinct() %>% ungroup()
  promoterIds.all <- promoterExpression$promoterId
  
  percentageFilterList <- list() 
  for (currentGroup in names(aliquotsPerGroup.byHistotype)) {
    print(paste0('Calculating percentage filter for: ', currentGroup))
    
    aliquots.currentGroup <- aliquotsPerGroup.byHistotype[[currentGroup]]
    expressedPromoters.logical <- promoterExpression$promoterId %in% expressedPromoters.all[[currentGroup]]
    promoterIds.currentGroup <- as.vector(promoterIds.all[expressedPromoters.logical])
    promoterActivity <- promoterExpression[expressedPromoters.logical, match(aliquots.currentGroup, colnames(promoterExpression))]
    # promoterActivity[is.na(promoterActivity)] <- 0
    
    percentageFilterList[[currentGroup]] <- calculatePercentageFilter(currentGroup, aliquots.currentGroup, 
                                                                      elementIds.currentGroup = promoterIds.currentGroup, 
                                                                      outlierExpression = promoterActivity, outlierList[[currentGroup]])
  }
  return(percentageFilterList)
}
