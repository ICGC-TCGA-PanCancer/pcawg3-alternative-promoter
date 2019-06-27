
prepareTranscriptRanges <- function(txdb, idConversionTable) {
  tx.granges <- transcripts(txdb, columns = c('tx_id', 'tx_name', 'gene_id'))
  seqlevelsStyle(tx.granges) <- 'UCSC'
  tx.granges <- tx.granges[match(idConversionTable$transcriptId, tx.granges$tx_name)]
  tx.granges$gene_id <- as.character(tx.granges$gene_id)
  tx.granges$promoterId <- idConversionTable$promoterId
  tx.granges$tssOrder <- findTSSOrder(tx.granges, promoterIdList = tx.granges$promoterId)
  return(tx.granges)
}

getExpressedPromoters <- function(promoterCoordinates, totalExpression) {
  currentGroup <- 'pcawg'
  expressedPromoters.pid <- as.character(unique(totalExpression$promoterId[totalExpression[[paste0(currentGroup, '.isPromoterExpressed')]]]))
  expressedPromoters.index <- match(expressedPromoters.pid, promoterCoordinates$promoterId)
  promoterCoordinates.expressed <- promoterCoordinates[expressedPromoters.index]
  return(promoterCoordinates.expressed)
}

preparePromoterCoordinates <- function(txdb, idConversionTable, tx.isoform.pcawg, totalExpression) {
  tx.granges <- prepareTranscriptRanges(txdb, idConversionTable)
  majorTSS.pcawg <- prepareMajorTssPerPromoter(aliquotsPerGroup.byHistotype, tx.isoform.pcawg, tx.granges)
  tssCoordinates.all <- promoters(tx.granges, upstream = 0, downstream = 1)
  names(tssCoordinates.all) <- tssCoordinates.all$tx_name
  promoterRanges <- tssCoordinates.all[majorTSS.pcawg]
  names(promoterRanges) <- promoterRanges$promoterId
  return(promoterRanges)
}
