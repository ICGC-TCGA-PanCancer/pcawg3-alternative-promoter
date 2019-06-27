
processAliquotId <- function(id) {
  if (substring(id, 1, 1) == 'X') {
    return (substring(id, 2, nchar(id)))
  } else {
    return (id)
  }
}

prepareTranscriptExpressionTable <- function(txdb, idConversionTable, metadata.pcawg, transcriptExpressionPath) {
    tx.granges <- transcripts(txdb, columns = c('tx_id', 'tx_name', 'gene_id'))
    tmp <- read.table(file = transcriptExpressionPath, header = TRUE, row.names = c('Feature'))
    aliquot_ids <- sapply(gsub('.', '-', names(tmp), fixed = TRUE), processAliquotId)
    colnames(tmp) <- aliquot_ids
    tmp <- tmp[rownames(tmp) %in% tx.granges$tx_name, ]

    tx.isoform.pcawg <- data.frame(matrix(NA, ncol = dim(tmp)[2], nrow = nrow(idConversionTable)))
    rownames(tx.isoform.pcawg) <- idConversionTable$transcriptId
    colnames(tx.isoform.pcawg) <- colnames(tmp)
    tx.isoform.pcawg[rownames(tmp), colnames(tmp)] <- tmp

    return(tx.isoform.pcawg)
}

getMostActive <- function(activityScores, minActivity, elementOrder) {
    if (all(is.na(activityScores))) {
        res <- elementOrder == 1
    } else {
        res <- activityScores == max(activityScores, na.rm = TRUE)  & activityScores >= minActivity & !duplicated(activityScores)
    }
    return (res)
}

findTSSOrder <- function(granges, promoterIdList) {
    granges$order <- 1:length(granges)
    granges$promoterIdList <- promoterIdList
    granges <- tbl_df(as.data.frame(granges))
    granges <- arrange(granges, seqnames, gene_id, promoterIdList, start, end, tx_id)
    granges <- group_by(granges, promoterIdList) %>%
        mutate(elementOrder = ifelse(strand == '+', rank(start, ties.method = 'first'), rank(-end, ties.method = 'first'))) %>%
        ungroup() %>%
        arrange(order)
    return(granges$elementOrder)
}

findMajorTSS <- function(granges, promoterIdList, elementOrder, activity, minActivity = 0) {
    granges$order <- 1:length(granges)
    granges$promoterIdList <- promoterIdList
    granges$activityScores <- activity
    granges$elementOrder <- elementOrder
    granges <- tbl_df(as.data.frame(granges))
    granges <- group_by(granges, promoterIdList) %>% arrange(elementOrder) %>% mutate(majorTSS = getMostActive(activityScores, minActivity, elementOrder)) %>% ungroup() %>% arrange(order)
    return(granges$majorTSS)
}
