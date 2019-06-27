
calculateGeneExpressionFromTranscript <- function(idConversionTable, tx.isoform.pcawg, metadata.pcawg) {
	tx.isoform.pcawg.annotated <- tbl_df(cbind(idConversionTable, tx.isoform.pcawg))
	geneExpression <- tx.isoform.pcawg.annotated %>% group_by(geneId) %>% mutate_each(funs(sum(., na.rm = TRUE)), -one_of(colnames(idConversionTable))) %>% ungroup
	geneExpression <- as.data.frame(geneExpression)
	return(geneExpression)
}

calculateRawPromoterActivity <- function(idConversionTable, tx.isoform.pcawg, metadata.pcawg) {
	tx.isoform.pcawg.annotated <- tbl_df(cbind(idConversionTable, tx.isoform.pcawg))
	promoterExpression_raw <- tx.isoform.pcawg.annotated %>% group_by(promoterId) %>% mutate_each(funs(sum(., na.rm = TRUE)), -one_of('transcriptId', 'tssId', 'promoterId', 'geneId')) %>% ungroup
	promoterExpression_raw <- as.data.frame(promoterExpression_raw)
	return(promoterExpression_raw)
}

calculateRelativePromoterActivity <- function(idConversionTable, geneExpression, promoterExpression_raw, metadata.pcawg) {
	promoterExpression_rel <- promoterExpression_raw[, -(1:4)] / geneExpression[, -(1:4)]
	promoterExpression_rel <- bind_cols(tbl_df(idConversionTable), promoterExpression_rel)
	promoterExpression_rel <- as.data.frame(promoterExpression_rel)
	return(promoterExpression_rel)
}

calculateTotalGeneExpressionByGroup <- function(geneExpression, metadata.pcawg, projectConfiguration) {
	total.gene <- t(apply(geneExpression, 1, tapply, metadata.pcawg$histology_abbreviation, sum, na.rm = TRUE))
	saveRDS(total.gene, file = file.path(projectConfiguration[['projectBase']], projectConfiguration[['analysisFolder']],
                                                        paste0('totalGeneExpression.', projectConfiguration[['dataSource']], '.', projectConfiguration[['annotationVersion']], '.rds')))
}

calculateTotalPromoterExpressionByGroup <- function(promoterExpression_raw, metadata.pcawg, projectConfiguration) {
    total.promoter <- t(apply(promoterExpression_raw, 1, tapply, metadata.pcawg$histology_abbreviation, sum, na.rm = TRUE))
    saveRDS(total.promoter, file = file.path(projectConfiguration[['projectBase']], projectConfiguration[['analysisFolder']],
                                                        paste0('totalPromoterExpression.', projectConfiguration[['dataSource']], '.', projectConfiguration[['annotationVersion']], '.rds')))
}

calculateTotalExpressionTableForPromoter <- function(idConversionTable, tx.isoform.pcawg, geneExpression, promoterExpression_raw, metadata.pcawg) {
	
  total.transcript <- t(apply(tx.isoform.pcawg, 1, tapply, metadata.pcawg$histology_abbreviation, sum, na.rm = TRUE))
  total.transcript.pcawg <- rowSums(tx.isoform.pcawg, na.rm = TRUE)

  total.gene <- t(apply(geneExpression[, -(1:4)], 1, tapply, metadata.pcawg$histology_abbreviation, sum, na.rm = TRUE))
  total.gene.pcawg <- rowSums(geneExpression[, -(1:4)], na.rm = TRUE)

  total.promoter <- t(apply(promoterExpression_raw[, -(1:4)], 1, tapply, metadata.pcawg$histology_abbreviation, sum, na.rm = TRUE))
	total.promoter.pcawg <- rowSums(promoterExpression_raw[, -(1:4)], na.rm = TRUE)
	
	totalExpression <- cbind(idConversionTable, total.transcript.pcawg, total.transcript, total.promoter.pcawg > 0, total.promoter > 0, total.gene.pcawg > 0, total.gene > 0)
	columnNames.custom <- c(colnames(idConversionTable), 'pcawg.expression', paste0(colnames(total.transcript), '.expression'),
                                                        'pcawg.isPromoterExpressed', paste0(colnames(total.promoter), '.isPromoterExpressed'),
                                                        'pcawg.isGeneExpressed', paste0(colnames(total.gene), '.isGeneExpressed'))
	colnames(totalExpression) <- columnNames.custom
	return(totalExpression)
}
