
calculateOutlier.Unimodal <- function(x, outlierDetectionMethod = 'sd', sdWindow = 2) {
    # The distribution is unimodal so do MLE to find unimodal fit
    unimodalFit <- fitdistr(x[!is.nan(x)], 'normal')
    referenceMean <- unimodalFit$estimate[['mean']]
    referenceSd <- unimodalFit$estimate[['sd']]

    # Since unimodal fit, set bimodal checks to FALSE
    check.mean <- FALSE
    check.area <- FALSE

    outlierState <- rep(0, length(x))
    if (outlierDetectionMethod == 'sd') {
        # If the data is less than the lower threshold set state to -1
        outlierState[x < (referenceMean - (sdWindow * referenceSd))] <- -1
        # If the data is greater than the lower threshold set state to 1
        outlierState[x > (referenceMean + (sdWindow * referenceSd))] <- 1
        # If there is missing data place NaN in its place
        outlierState[is.nan(x)] <- NaN
    } else {
        print('For future work & outlier analysis techniques')
    }

    result <- list(outlierState = outlierState, check.mean = check.mean, check.area = check.area, referenceMean = referenceMean, referenceSd = referenceSd)
    return(result)
}

calculateOutlierForCurrentGroup <- function(currentGroup, outlierExpression) {
  outlierResult <- lapply(1:nrow(outlierExpression), function(i) {calculateOutlier.Unimodal(unlist(outlierExpression[i,]), outlierDetectionMethod = 'sd', sdWindow = 2)})
  return(outlierResult)
}

prepareOutlierListForPromoter <- function(aliquotsPerGroup.byHistotype, totalExpression, promoterExpression) {
  expressedPromoters.all <- sapply(names(aliquotsPerGroup.byHistotype), function(currentGroup) {as.character(unique(totalExpression$promoterId[totalExpression[[paste0(currentGroup,'.isPromoterExpressed')]]]))})
  names(expressedPromoters.all) <- names(aliquotsPerGroup.byHistotype)
  
  promoterExpression <- promoterExpression[, -(1:2)] %>% group_by(promoterId) %>% distinct() %>% ungroup()
  
  outlierList <- list()
  for (currentGroup in names(aliquotsPerGroup.byHistotype)) {
    print(paste0('Calculating outliers for: ', currentGroup))
    
    aliquots.currentGroup <- aliquotsPerGroup.byHistotype[[currentGroup]]
    expressedPromoters.logical <- promoterExpression$promoterId %in% expressedPromoters.all[[currentGroup]]
    promoterActivity <- promoterExpression[expressedPromoters.logical, match(aliquots.currentGroup, colnames(promoterExpression))]
    promoterActivity[is.na(promoterActivity)] <- 0
    
    outlierList[[currentGroup]] <- calculateOutlierForCurrentGroup(currentGroup, outlierExpression = promoterActivity)
  }
  return(outlierList)
}
