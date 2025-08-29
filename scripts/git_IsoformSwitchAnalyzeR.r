library(IsoformSwitchAnalyzeR)
library(BSgenome.Mmusculus.UCSC.mm39)
library(plyranges)
full_gtf = rtracklayer::import("/Users/lipeihang/novel_transcript_extract/novel_Taf1_mouse.gtf")
aSwitchList <- importGTF(pathToGTF = "/Users/lipeihang/novel_transcript_extract/novel_Taf1_mouse.gtf")

# Add ORF analysis
aSwitchList <- analyzeORF(aSwitchList,
                         genomeObject = BSgenome.Mmusculus.UCSC.mm39, # mouse genome
                         orfMethod = 'mostUpstream',
                         minORFlength = 100)



coding_isoforms = aSwitchList$orfAnalysis %>% filter(PTC == FALSE) %>% pull(isoform_id)

extractSequence(aSwitchList,onlySwitchingGenes = FALSE, pathToOutput = "/Users/lipeihang/Desktop/AA_output")

#re-import gtf having real novel Taf1 transcripts with known Taf1
full_gtf = rtracklayer::import("/Users/lipeihang/Desktop/novel_transcript_extract_try/PTC_NCBITaf1.gtf")
aSwitchList <- importGTF(pathToGTF = "/Users/lipeihang/Desktop/novel_transcript_extract_try/PTC_NCBITaf1.gtf")

# Add ORF analysis
aSwitchList <- analyzeORF(aSwitchList,
                          genomeObject = BSgenome.Mmusculus.UCSC.mm39, # mouse genome
                          orfMethod = 'mostUpstream',
                          minORFlength = 100) 

#Extract AA dierecty
extractSequence(aSwitchList,onlySwitchingGenes = FALSE, pathToOutput = "/Users/lipeihang/Desktop/AA_output/NoPTC_Novel_known/")

# Pfam predicts protein domain
aSwitchList <- analyzePFAM(
  switchAnalyzeRlist = aSwitchList,
  pathToPFAMresultFile = "/Users/lipeihang/Desktop/AA_output/NoPTC_Novel_known/novel_isoform_pfam.txt", 
  showProgress = TRUE
)

# IUPred2A predicts protein IDR
aSwitchList <- analyzeIUPred2A(
  switchAnalyzeRlist = aSwitchList,
  pathToIUPred2AresultFile = "/Users/lipeihang/Desktop/AA_output/NoPTC_Novel_known/isoformSwitchAnalyzeR_isoform_AA.result",
  showProgress = TRUE
)

# Predict signal peptide
aSwitchList <- analyzeSignalP(
  switchAnalyzeRlist = aSwitchList,
  pathToSignalPresultFile = "/Users/lipeihang/Downloads/SignalP_output_all_results/SignalP_sum.txt",
  minSignalPeptideProbability = 0.5,
  quiet=FALSE
)

switchPlotTranscript(aSwitchList, isoform_id = "MSTRG.1.1")
switchPlotTranscript(aSwitchList)

outputPlots = TRUE
