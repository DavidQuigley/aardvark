genelist=c("ATM","ATR","BRCA1","BRCA2","BAP1","BARD1","BRIP1","CDK12","CHEK1","CHEK2","FANCA","FANCC","FANCL","MRE11","PALB2","RAD50","RAD51C","RAD51D")
transcript_ids = rep(NA, length(genelist))

ensembl = biomaRt::useDataset("hsapiens_gene_ensembl", mart = useMart("ensembl") )
transcript_cache_hg38 = list()
sequence_cache_hg38 = list()
for( i in 1:length( genelist ) ){
    print( genelist[i] )
    transcript_ids[i] = hg38_ensembl_coding_transcripts$transcript_id[hg38_ensembl_coding_transcripts$symbol==genelist[i] ]
    transcript = aardvark::TranscriptData( ensembl, EnsDb.Hsapiens.v86, transcript_ids[i] )
    ts = biomaRt::getSequence( id=transcript_ids[i],
                               type = "ensembl_transcript_id",
                               seqType = "transcript_exon_intron",
                               mart = ensembl)
    transcript_cache_hg38 = c( transcript_cache_hg38, list( transcript ) )
    sequence_cache_hg38 = c( sequence_cache_hg38, list( ts ) )
}
names(transcript_cache_hg38) = transcript_ids
#names(sequence_cache_hg38) = transcript_ids
save( transcript_cache_hg38, sequence_cache_hg38, exon_cache_hg38, file='/notebook/code/aardvark/data/ensembl_cache.RData' )


#hg38_ensembl_coding_transcripts = read.table('/notebook/code/aardvark/exec/hg38_ensembl_coding_transcripts.txt', header=TRUE, sep='\t')
#save( hg38_ensembl_coding_transcripts, file="/notebook/code/aardvark/exec/hg38_ensembl_coding_transcripts.RData")

