

## pcawg3-alternative-promoter

This pipeline is developed as part of PCAWG project to estimate promoter
activity from RNA-Seq data using transcript expression estimates. The 
promoter and gene centric tables can be created for alternative promoter
phenotype using this pipeline.

### Required files

The following 3 files are required to run the alternative promoter analysis: 

1. PCAWG Metadata file: This file has been published as part of the PCAWG data and can be obtained from:
https://dcc.icgc.org/api/v1/download?fn=/PCAWG/transcriptome/metadata/rnaseq.extended.metadata.aliquot_id.V4.tsv.gz
2. PCAWG Transcript expression table file: This file has been published as part of the PCAWG data and can be obtained from:
https://dcc.icgc.org/api/v1/download?fn=/PCAWG/transcriptome/transcript_expression/pcawg.rnaseq.transcript.expr.fpkm.tsv.gz
3. Gencode v19 txdb file: We created a txdb object using the GTF file for Gencode v19 (https://www.gencodegenes.org/human/release_19.html). The txdb file can be downloaded from: https://drive.google.com/open?id=1Uhp65KqCPupRUWwrPgyxTCFuSO5Yul2d

### How to run

The main file can be found at "./R/geneCentricAnalysis-promoter.R". 
The user needs to change the paths for the input files. 
The script will calculate the promoter activities, coordinates, promoter and gene centric tables. 

## Citing this pipeline

If you use this pipeline, please cite: Calabrese, C., Davidson, N.R., Fonseca, N.A., He, Y., 
Kahles, A., Lehmann, K.-V., Liu, F., Shiraishi, Y., Soulette, C.M., Urban, L., et al. (2018). 
"Genomic basis for RNA alterations revealed by whole-genome analyses of 27 cancer types" bioRxiv, 183889.

## Contributors

This pipeline is developed and maintained by Deniz Demircioglu and Jonathan
Göke.
