
generateExonTables <- function(txdb) {
  exonsRanges = exons(txdb, columns = c('EXONID', 'EXONRANK', 'TXID', 'GENEID', 'TXNAME')) 
  seqlevelsStyle(exonsRanges) <- 'UCSC'
  
  exonsReducedRanges = reduce(exonsRanges, with.revmap = T)
  exonsReducedRanges$MergedExonId = 1:length(exonsReducedRanges)  # add mergedExon ID
  
  ## Add mergedExon ID to exonsRanges
  revmapTMP = exonsReducedRanges$revmap
  names(revmapTMP) <- exonsReducedRanges$MergedExonId
  revmapTMP = unlist(revmapTMP)
  exonsRanges$MergedExonId = as.integer(names(sort(revmapTMP)))
  
  exonsTable = tbl_df(as.data.frame(exonsRanges))  # tbl_df for manipulation
  
  ## (1) Annotate exonsRanges, for each potential first exon, annotate by transcript and gene
  exonsByTxidTable = AnnotationDbi::select(txdb, keys = unique(as.character(unlist(exonsRanges$TXID))), keytype = 'TXID', 
                                           columns = c('EXONID', 'EXONRANK'))
  exonsByTxidTable <- exonsByTxidTable[, c('EXONID', 'EXONRANK')]
  
  exonsByTxidTable = mutate(group_by(exonsByTxidTable, EXONID), MinExonRank = min(EXONRANK))
  exonsByTxidTable = mutate(group_by(exonsByTxidTable, EXONID), MaxExonRank = max(EXONRANK))
  
  ## ad txLength here? select which Tx?
  exonsUniqueByTxid = group_by(exonsByTxidTable, EXONID) %>%
    filter(EXONRANK == MinExonRank) %>%
    distinct() %>%
    dplyr::select(EXONID, MinExonRank, MaxExonRank)
  exonsTable=inner_join(exonsTable, exonsUniqueByTxid, by = 'EXONID')
  
  ## ----- (2) Annotate exonsReducedRanges -----
  ##extract per mergedExon information from exonsTable
  exonsTableTMP <- dplyr::select(exonsTable, MergedExonId, MinExonRank, MaxExonRank) %>%
    group_by(MergedExonId) %>%
    mutate(MinMergedExonRank = min(MinExonRank),
           MaxMergedExonRank = max(MaxExonRank),
           SumMergedExons = n()) %>%
    filter(MinExonRank == min(MinExonRank)) %>%
    filter(MaxExonRank == max(MaxExonRank)) %>%
    distinct() %>%
    dplyr::select(MergedExonId, MinMergedExonRank, MaxMergedExonRank, SumMergedExons)

  exonsReducedTable <- tbl_df(as.data.frame(exonsReducedRanges))  # tbl_df for manipulation
  exonsReducedTable <- inner_join(exonsReducedTable, exonsTableTMP, by = 'MergedExonId')
  
  exonsReducedRanges=GRanges(exonsReducedTable$seqnames,
                             ranges = IRanges(start = exonsReducedTable$start, exonsReducedTable$end),
                             strand = exonsReducedTable$strand,
                             data.frame(dplyr::select(exonsReducedTable,
                                                      MergedExonId,
                                                      MinMergedExonRank,
                                                      MaxMergedExonRank,
                                                      SumMergedExons)))
  return(list(exonsReducedRanges = exonsReducedRanges, exonsTable = exonsTable))
}


prepareReducedExonRanges <- function(txdb) {
  exonTables.all <- generateExonTables(txdb)
	exonReducedRanges <- exonTables.all[['exonsReducedRanges']]
  return(exonReducedRanges)	
}

