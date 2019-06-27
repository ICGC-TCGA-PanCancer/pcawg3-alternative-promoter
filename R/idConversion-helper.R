tssGrouper <- function(geneId, geneToExon) {
    geneInfo <- geneToExon[[geneId]]

    if (runValue(strand(geneInfo) == '+')[[1]]) {
        startPos <- sapply(start(geneInfo), head, 1)
    } else {
        startPos <- sapply(end(geneInfo), head, 1)
    }

    transcriptIds <- lapply(unique(startPos), function(pos) {names(startPos[startPos == pos])})
    geneIds <- rep(geneId, length(transcriptIds))

    tssInformationDF <- data.table(geneId = geneIds, transcriptId = transcriptIds)

    return(tssInformationDF)
}

promoterGrouper <- function(geneId, geneToExon) {
    geneInfo <- geneToExon[[geneId]]
    firstExonRanges <- Reduce(c, sapply(geneInfo, head, 1))
    names(firstExonRanges) <- names(geneInfo)

    promoterRanges <- reduce(firstExonRanges, with.revmap = TRUE, min.gapwidth = 0L)
    revmap <- mcols(promoterRanges)$revmap

    geneIds <- rep(geneId, length(promoterRanges))
    transcriptIds <- lapply(revmap, function(idx) {names(firstExonRanges[idx])})

    promoterInformationDF <- data.table(geneId = geneIds, transcriptId = transcriptIds)

    return(promoterInformationDF)
}

preparePromoterIdConversion <- function(txdb, mc.cores = 1, force = FALSE) {
  tx.ranges <- transcripts(txdb, columns = c("tx_name", "gene_id"))
  tx.all <- mcols(tx.ranges)$tx_name
  txToGene <- as.character(mcols(tx.ranges)$gene_id)

  gene.ranges <- genes(txdb, columns = c("tx_name"))
  geneToTx <- as.vector(mcols(gene.ranges)$tx_name)

  txByGenes <- transcriptsBy(txdb, by = 'gene')
  exonsByTx <- exonsBy(txdb, by = 'tx')
  names(exonsByTx) <- tx.all

  geneList <- names(txByGenes)

  if (mc.cores == 1) {
    geneToExon <- lapply(geneToTx, function(txs) {exonsByTx[txs]})
  } else {
    require(parallel)
    geneToExon <- mclapply(geneToTx, function(txs) {exonsByTx[txs]}, mc.cores = mc.cores)
  }
  
  if (mc.cores == 1) {
    tssList <- lapply(geneList, tssGrouper, geneToExon)
  } else {
    require(parallel)
    tssList <- mclapply(geneList, tssGrouper, geneToExon, mc.cores = mc.cores)
  }
  tssDF <- Reduce(rbind, tssList)
  tssIds <- paste0(rep('tss.', dim(tssDF)[1]), 1:dim(tssDF)[1])
  tssDF <- data.frame(tssId = tssIds, tssDF)

  if (mc.cores == 1) {
    promoterList <- lapply(geneList, promoterGrouper, geneToExon)
  } else {
    require(parallel)
    promoterList <- mclapply(geneList, promoterGrouper, geneToExon, mc.cores = mc.cores)
  }
  promoterDF <- Reduce(rbind, promoterList)
  promoterIds <- paste0(rep('prmtr.', dim(promoterDF)[1]), 1:dim(promoterDF)[1])
  promoterDF <- data.frame(promoterId = promoterIds, promoterDF)
    
	c.tssId <- unlist(apply(tssDF, 1, function(t) {rep(t$tssId, length(t$transcriptId))}))
	c.geneId <- unlist(apply(tssDF, 1, function(t) {rep(t$geneId, length(t$transcriptId))}))
	c.transcriptId <- unlist(tssDF$transcriptId)
	c.promoterId2 <- unlist(apply(promoterDF, 1, function(p) {rep(p$promoterId, length(p$transcriptId))}))
	names(c.promoterId2) <- unlist(promoterDF$transcriptId)
	c.promoterId <- c.promoterId2[c.transcriptId]

	idConversionTable <- data.frame(transcriptId = c.transcriptId, tssId = c.tssId, promoterId = c.promoterId, geneId = c.geneId, row.names = c.transcriptId, stringsAsFactors = FALSE)
  return(idConversionTable)
}
