
calculateMeanElementExpression <- function(idConversionTable, elementExpression, aliquotsPerGroup.byHistotype) {
	meanExpression <- idConversionTable
	for (currentGroup in names(aliquotsPerGroup.byHistotype)) {
		meanExpression[paste0(currentGroup, '.expression')] <- rowMeans(elementExpression[, match(aliquotsPerGroup.byHistotype[[currentGroup]], colnames(elementExpression))])
	}
	return(meanExpression)
}

calculateMeanExpressionFilter <- function(meanExpression, idConversionTable) {
	meanExpressionFilter <- meanExpression
	meanExpressionFilter[, -(1:ncol(idConversionTable))] <- meanExpressionFilter[, -(1:ncol(idConversionTable))] > 1
  return(meanExpressionFilter)
}

calculateElementExpressionFilter <- function(elementExpression, idConversionTable) {
	elementExpressionFilter <- elementExpression
  elementExpressionFilter[, -(1:ncol(idConversionTable))] <- elementExpressionFilter[, -(1:ncol(idConversionTable))] > 1
  return(elementExpressionFilter)
}

calculateElementExpressionMeanFilter <- function(aliquotsPerGroup.byHistotype, elementExpressionFilter, meanExpressionFilter, idConversionTable) {
	aliquot.all <- names(elementExpressionFilter)[-(1:ncol(idConversionTable))]
	aliquotToGroup.all <- sapply(aliquot.all, function(aliquot) {names(which(sapply(aliquotsPerGroup.byHistotype, function(currentGroupAliquots) {aliquot %in% currentGroupAliquots})))})
	elementExpressionMeanFilter <- meanExpressionFilter[paste0(aliquotToGroup.all, '.expression')]
	names(elementExpressionMeanFilter) <- names(aliquotToGroup.all)
	elementExpressionMeanFilter <- dplyr::bind_cols(idConversionTable, elementExpressionMeanFilter)
	return(elementExpressionMeanFilter)
}

prepareExpressionFilter <- function(aliquotsPerGroup.byHistotype, idConversionTable, elementExpression) {
  print('Prepare expression filter...')
  meanExpression <- calculateMeanElementExpression(idConversionTable, elementExpression, aliquotsPerGroup.byHistotype)
  meanExpressionFilter <- calculateMeanExpressionFilter(meanExpression, idConversionTable)
  elementExpressionFilter <- calculateElementExpressionFilter(elementExpression, idConversionTable)
  elementExpressionMeanFilter <- calculateElementExpressionMeanFilter(aliquotsPerGroup.byHistotype, elementExpressionFilter, meanExpressionFilter, idConversionTable)

	expressionFilter <- tbl_df(as.data.frame(as.matrix(elementExpressionFilter[, -(1:ncol(idConversionTable))]) | as.matrix(elementExpressionMeanFilter[, -(1:ncol(idConversionTable))])))
  expressionFilter <- dplyr::bind_cols(idConversionTable, expressionFilter)
  
  expressionFilter <- expressionFilter[, 3:(dim(expressionFilter)[2])] %>% group_by(promoterId) %>% distinct() %>% ungroup()
  expressionFilter <- as.data.frame(expressionFilter)
  
  return(expressionFilter)
}
