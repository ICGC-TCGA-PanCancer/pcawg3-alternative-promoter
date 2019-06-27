library(GenomicFeatures)
library(dplyr)
library(GenomicAlignments)
library(data.table)
library(graphics)
library(MASS)

source('./R/geneCentricAnalysis-helper.R')
source('./R/idConversion-helper.R')
source('./R/splitReadRatios-helper.R')
source('./R/calculateVotes.R')
source('./R/calculateElementActivity.R')
source('./R/calculateElementCoordinates.R')
source('./R/finalizeOutlierResults.R')

#####

# The following code generates the promoter coordinates, promoter activity
# estimates and the promoter/gene outlier tables for the PCAWG project using the
# Gencode v19 annotation and the PCAWG transcript expression estimates. The
# following 3 files are required to run the alternative promoter analysis: 
# 1. PCAWG Metadata file: This file has been published as part of the PCAWG data and can be obtained from
# https://dcc.icgc.org/api/v1/download?fn=/PCAWG/transcriptome/metadata/rnaseq.extended.metadata.aliquot_id.V4.tsv.gz
# 2. PCAWG Transcript expression table file: This file has been published as
# part of the PCAWG data and can be obtained from
# https://dcc.icgc.org/api/v1/download?fn=/PCAWG/transcriptome/transcript_expression/pcawg.rnaseq.transcript.expr.fpkm.tsv.gz
# 3. Gencode v19 txdb file: We created a txdb object using the GTF file for
# Gencode v19 (https://www.gencodegenes.org/human/release_19.html) 
# The txdb file can be downloaded from
# https://drive.google.com/open?id=1Uhp65KqCPupRUWwrPgyxTCFuSO5Yul2d

path_to_gencode_v19_txdb <- './temp-data/gencode.v19.annotation.sqlite'
path_to_pcawg_metadata <- './data/rnaseq.extended.metadata.aliquot_id.V4.tsv'
path_to_pcawg_tx_expression <- './data/fixed_fpkm.tsv'

print('### Load PCAWG metadata...')
metadata.pcawg <- read.delim(path_to_pcawg_metadata, header = TRUE, stringsAsFactors = FALSE)
metadata.pcawg <- metadata.pcawg[metadata.pcawg$wgs_white_black_gray == 'Whitelist',]
metadata.pcawg <- metadata.pcawg[order(metadata.pcawg$project_code),]

print('### Load Gencode v19 txdb object...')
txdb <- loadDb(path_to_gencode_v19_txdb)

print('### Prepare the id mapping between transcripts, tss, promoters and genes...')
idConversionTable <- preparePromoterIdConversion(txdb)

print('### Prepare transcript expression table...')
transcriptExpression <- prepareTranscriptExpressionTable(txdb, idConversionTable, metadata.pcawg, transcriptExpressionPath = path_to_pcawg_tx_expression)

print('### Prepare gene expression tables...')
geneExpression <- calculateGeneExpressionFromTranscript(idConversionTable, transcriptExpression, metadata.pcawg)

print('### Prepare raw promoter expression tables...')
promoterExpression_raw <- calculateRawPromoterActivity(idConversionTable, transcriptExpression, metadata.pcawg)

print('### Prepare relative promoter expression tables...')
promoterExpression_rel <- calculateRelativePromoterActivity(idConversionTable, geneExpression, promoterExpression_raw, metadata.pcawg)

print('### Prepare total promoter expression tables...')
totalExpression <- calculateTotalExpressionTableForPromoter(idConversionTable, transcriptExpression, geneExpression, promoterExpression_raw, metadata.pcawg)

print('### Determine tss coordinate for each promoter...')
promoterRanges <- preparePromoterCoordinates(txdb, idConversionTable, transcriptExpression, totalExpression)
promoterRanges.expressed <- getExpressedPromoters(promoterRanges, totalExpression)

print('### Prepare promoter centric outlier tables...')
promoterOutliers <- prepareOutlierResultsForPromoter(metadata.pcawg, idConversionTable, promoterExpression_raw, promoterExpression_rel, totalExpression, 
                                                     promoterRanges, exonReducedRanges)

print('### Prepare gene centric outlier tables...')
geneOutliers <- prepareOutlierResultsForGene(promoterOutliers)
