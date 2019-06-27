getVotes.promoter <- function(tx.isoform.pcawg, tx.granges) {
    votes <- apply(tx.isoform.pcawg, 2, function(activity) {findMajorTSS(granges = tx.granges,
                                                                       promoterIdList = tx.granges$promoterId,
                                                                       elementOrder = tx.granges$tssOrder,
                                                                       activity = activity,
                                                                       minActivity = 0)})
    return(votes)
}

getMajorTSS <- function(votes, tx.granges) {
    votes.combined <- rowSums(votes, na.rm = TRUE)
    majorTSS.logical <- findMajorTSS(granges = tx.granges, promoterIdList = tx.granges$promoterId, elementOrder = tx.granges$tssOrder, activity = votes.combined, minActivity = 0)
    majorTSS.transcript <- tx.granges$tx_name[majorTSS.logical]
    return(majorTSS.transcript)
}

prepareMajorTssPerPromoter <- function(aliquotsPerGroup.byHistotype, tx.isoform.pcawg, tx.granges) {
  votes.pcawg <- getVotes.promoter(tx.isoform.pcawg, tx.granges)
  majorTSS.pcawg <- getMajorTSS(votes.pcawg, tx.granges)
  return(majorTSS.pcawg)
}

