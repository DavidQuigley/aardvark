ts_start = Sys.time()
suppressPackageStartupMessages( require( optparse, quietly=TRUE, warn.conflicts=FALSE ) )

option_list <- list(
    make_option(c("-s", "--sample_id"), type = "character", help = "Sample identifer string, used as prefix for output files. [Required]"),
    make_option(c("-i", "--fn_bam"), type = "character", help = "Path to BAM formatted file with sequence data to read [Required]"),
    make_option(c("-o", "--dir_out"), type = "character", help = "Path to output folder to write results [Required]"),

    make_option(c("-t", "--transcript_id"), type = "character", help = "Ensembl transcript identifier for target gene (format ENST###########) [Required]"),
    make_option(c("-c", "--chromosome"), type = "character", help = "Target variant chromosome to find in fn_bam [Required]"),
    make_option(c("-p", "--position"), type = "numeric", help = "Target variant position on chromosome in fn_bam [Required]"),
    make_option(c("-r", "--REF"), type = "character", help = "Target variant reference nucleotide sequence in fn_bam [Required]"),
    make_option(c("-a", "--ALT"), type = "character", help = "Target variant alternative nucleotide sequence in fn_bam [Required]"),
    make_option(c("-w", "--window_size"), type = "numeric", help = "Search space in nucleotides centered on pathogenic mutation; default 1000", default=1000),
    make_option(c("-g", "--genome_draft"), type = "numeric", help = "Human genome draft for gene models, one of {19,38}; default: 38", default=38),

    make_option(c("-n", "--near_bound"), type = "numeric", help = "Size of nucleotide window adjacent to target variant for locally remapping read ends; default 6", default=6),
    make_option(c("-l", "--allow_insertions_in_realign"), type = "logical", action="store_true", help = "Allow insertions in local realignment of soft-clipped read ends; default TRUE", default=TRUE),
    make_option(c("-m", "--min_nt_for_distant_realign"), type = "numeric", help = "Minimum number of perfectly aligned contiguous nucleotides to accept distant realignment of a soft-clipped read end; default: 15", default=15),
    make_option(c("-e", "--min_percent_realigned"), type = "numeric", help = "Minimum percent of a perfectly aligned soft-clipped read to accept realignment, default 0.9 (90%)", default=0.9),

    make_option(c("-q", "--min_nt_qual"), type = "numeric", help = "Minimum quality score for realigned nucleotides, default: 20", default=20),
    make_option(c("-x", "--proxy_http"), type = "character", help = "Address of proxy server to use for HTTP requests, only relevant if behind a proxy server, default empty"),
    make_option(c("-y", "--proxy_https"), type = "character", help = "Address of proxy server to use for HTTPS requests, only relevant if behind a proxy server, default empty"),
    make_option(c("-f", "--write_rdata"), type = "logical", action="store_true", help = "Write RData file with analysis results to output folder, default FALSE", default=FALSE),
    make_option(c("-v", "--verbose"), type = "logical", action="store_true", help = "Write progress to stdout, default TRUE", default=TRUE)
)

opt = parse_args( OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]") )

if( sum( names(opt)=="fn_bam")==0 ){
    stop( "invalid parameter(s) for input. Input BAM file required in parameter fn_bam")
}
if( !file.exists( opt$fn_bam ) ){
    stop( paste0( "Error: input BAM file ", opt$fn_bam," not found") )
}
if( sum( names(opt)=="dir_out")==0 ){
    stop( "invalid parameter(s) for input. Output folder required in parameter dir_out")
}
if( ! dir.exists( opt$dir_out ) ){
    stop( paste0( "Error: output directory ", opt$dir_out," not found") )
}
if( sum( names(opt)=="sample_id")==0 ){
    stop( "invalid parameter(s) for input. Sample identifier required in parameter sample_id")
}
if( opt$genome_draft != 38 & opt$genome_draft != 19 ){
    stop( "human genome draft in parameter genome_draft must be one of 38, 19")
}
if( sum( names(opt)=="proxy_http")==1 & sum( names(opt)=="proxy_https")==1 ){
    cat(paste("MESSAGE: Set https proxy to", opt$proxy_https, "\n") )
    Sys.setenv( http_proxy = opt$proxy_http )
    Sys.setenv( https_proxy = opt$proxy_https )
}

if( opt$REF == "" & opt$ALT == "" ){
    stop( paste0( "invalid parameter(s) for input: Reference and Alternate sequence cannot both be empty" ) )
}
if( opt$REF == opt$ALT ){
    stop( paste0( "invalid parameter(s) for input: Reference and Alternate sequence cannot be identical." ) )
}

if( opt$verbose ){  cat("MESSAGE: Loading required packages\n") }
suppressPackageStartupMessages( require( GenomicAlignments, quietly=TRUE, warn.conflicts=FALSE) )
suppressPackageStartupMessages( require( Rsamtools, quietly=TRUE, warn.conflicts=FALSE ) )
suppressPackageStartupMessages( require( Biostrings, quietly=TRUE, warn.conflicts=FALSE ) )
suppressPackageStartupMessages( require( stringr, quietly=TRUE, warn.conflicts=FALSE ) )
suppressPackageStartupMessages( require( BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE, warn.conflicts=FALSE ) )
suppressPackageStartupMessages( require( BSgenome.Hsapiens.UCSC.hg19, quietly=TRUE, warn.conflicts=FALSE ) )
suppressPackageStartupMessages( require( biomaRt, quietly=TRUE, warn.conflicts=FALSE ) )
suppressPackageStartupMessages( require( ensembldb, quietly=TRUE, warn.conflicts=FALSE ) )
suppressPackageStartupMessages( require( VariantAnnotation, quietly=TRUE, warn.conflicts=FALSE ) )
suppressPackageStartupMessages( require( aardvark, quietly=TRUE, warn.conflicts=FALSE ) )

if( opt$verbose ){  cat( paste0( "MESSAGE: Loading Ensembl data and fetching transcript ",
                                   opt$transcript_id, "\n" ) ) }

if( opt$genome_draft== 38 ){
    suppressPackageStartupMessages( require( EnsDb.Hsapiens.v86, quietly=TRUE, warn.conflicts=FALSE ) )
    suppressPackageStartupMessages( require( BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE, warn.conflicts=FALSE ) )
    ensembl = biomaRt::useDataset("hsapiens_gene_ensembl", mart = biomaRt::useMart("ensembl") )
}else{
    suppressPackageStartupMessages( require( BSgenome.Hsapiens.UCSC.hg19, quietly=TRUE, warn.conflicts=FALSE ) )
    suppressPackageStartupMessages( require( EnsDb.Hsapiens.v75, quietly=TRUE, warn.conflicts=FALSE ) )
    ensembl <-useMart(biomart="ENSEMBL_MART_ENSEMBL",
                      host="https://grch37.ensembl.org",
                      path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
}

if( opt$genome_draft == 38 ){
    transcript = aardvark::TranscriptData( ensembl, EnsDb.Hsapiens.v86, opt$transcript_id )
    Hsapiens_version = BSgenome.Hsapiens.UCSC.hg38
}else if( opt$genome_draft == 19 ){
    transcript = aardvark::TranscriptData( ensembl, EnsDb.Hsapiens.v75, opt$transcript_id )
    Hsapiens_version = BSgenome.Hsapiens.UCSC.hg19
}

window_start = opt$position - round(opt$window_size/2)
window_end =  opt$position + round(opt$window_size/2)
alignment_helper = aardvark::AlignmentWindow(Hsapiens_version,
                                             opt$chromosome,
                                             window_start = window_start,
                                             window_end = window_end)

pathogenic_mutation = aardvark::Mutation(chrom = opt$chromosome,
                                       pos = opt$position,
                                       seq_ref = opt$REF,
                                       seq_alt = opt$ALT,
                                       transcript = transcript)

if( opt$verbose ){  cat( paste0( "MESSAGE: Reading BAM ", opt$fn_bam, '\n' ) ) }
if( opt$verbose ){  cat( paste0( "MESSAGE: Loading reads in range ",
                                 opt$chrom, ":", window_start, "-", window_end, "\n" ) ) }

bb = aardvark::BamData( fn_bam = opt$fn_bam,
                        chrom = opt$chromosome,
                        start = window_start,
                        end = window_end)
reads = vector( mode = "list", length = bb$N )
gr_pathogenic = genomicRangesFromMutation( pathogenic_mutation )

if( opt$verbose ){ cat( paste("MESSAGE: Processing", bb$N, "reads\n") ) }
for( ctr in 1:bb$N ){
    if( opt$verbose & ctr %% 100 == 0 ){ cat( paste( "MESSAGE: Processed read", ctr, "of", bb$N, "\n" ) ) }
    read = aardvark::read_from_BamData( bb, ctr )
    read = aardvark::realign_read( read,
                                   alignment_helper,
                                   pathogenic_mutation = pathogenic_mutation,
                                   gr_pathogenic = gr_pathogenic,
                                   allow_insertions_in_realign = opt$allow_insertions_in_realign,
                                   min_nt_for_distant_realign = opt$min_nt_for_distant_realign,
                                   min_percent_realigned = opt$min_percent_realigned,
                                   near_bound = opt$near_bound,
                                   min_nt_qual = opt$min_nt_qual)
    reads[[ ctr ]] = aardvark::assess_reversion(read,
                                                transcript,
                                                pathogenic = pathogenic_mutation,
                                                align_window = alignment_helper,
                                                min_nt_qual = opt$min_nt_qual,
                                                gr_pathogenic = gr_pathogenic,
                                                gr_exclude = alignment_helper$homopolymer_regions)
}
if( opt$verbose ){ cat( "MESSAGE: Completed realignment, writing output files\n" ) }

read_summary = aardvark::summarize_candidates( reads, transcript, pathogenic_mutation = pathogenic_mutation )

ss = paste0( opt$sample_id, "_", opt$chromosome, "_", opt$position )
fn_out_Rdata = paste0( opt$dir_out, "/", ss, "_AARDVARK.RData")
fn_out_txt = paste0( opt$dir_out, "/", ss, "_AARDVARK_reversion_summary.txt")
fn_out_reads = paste0( opt$dir_out, "/", ss, "_AARDVARK_reversion_summary_with_reads.txt")
fn_out_settings = paste0( opt$dir_out, "/", ss, "_AARDVARK_run_settings.txt")

# Write out a summary of these results
aardvark::write_reversion_summary( read_summary, pathogenic_mutation, fn_out_txt, opt$sample_id, opt$transcript_id )
if( opt$verbose ){ cat( paste("MESSAGE: Saved output reversion summary file", fn_out_txt, "\n") )  }

aardvark::write_read_summary( read_summary, pathogenic_mutation, fn_out_reads, opt$sample_id, opt$transcript_id )
if( opt$verbose ){ cat( paste("MESSAGE: Saved output reversion summary file with reads:", fn_out_reads, "\n") ) }

# Write out the parameters used to generate these results
df_param = data.frame( param=c(), value=c(), stringsAsFactors = FALSE )
df_param = rbind( df_param, c( "time_start", format( ts_start, "%F %T" ) ) )
df_param = rbind( df_param, c( "time_end", format( Sys.time(), "%F %T" ) ) )
df_param = rbind( df_param, c( "aardvark_version", as.character(packageVersion( "aardvark" ) ) ) )
df_param = rbind( df_param, c( "sample_id", opt$sample_id ) )
df_param = rbind( df_param, c( "input_bam", opt$fn_bam ) )
df_param = rbind( df_param, c( "output_dir", opt$dir_out ) )
df_param = rbind( df_param, c( "transcript_identifier", opt$transcript_id ) )
df_param = rbind( df_param, c( "chromosome", opt$chromosome ) )
df_param = rbind( df_param, c( "variant_position", opt$position ) )
df_param = rbind( df_param, c( "variant_REF", opt$REF ) )
df_param = rbind( df_param, c( "variant_ALT", opt$ALT ) )
df_param = rbind( df_param, c( "window_size", opt$window_size ) )
df_param = rbind( df_param, c( "near_bound", opt$near_bound ) )
df_param = rbind( df_param, c( "allow_insertions_in_realign", opt$allow_insertions_in_realign )  )
df_param = rbind( df_param, c( "min_nt_for_distant_realign", opt$min_nt_for_distant_realign ) )
df_param = rbind( df_param, c( "min_percent_realigned", opt$min_percent_realigned ) )
df_param = rbind( df_param, c( "min_nt_qual", opt$min_nt_qual ) )
df_param = rbind( df_param, c( "human_genome_draft", opt$genome_draft ) )
df_param = rbind( df_param, c( "n_reads_unreverted", read_summary$n_unreverted ) )
df_param = rbind( df_param, c( "n_reads_uninformative", read_summary$n_uninformative ) )
df_param = rbind( df_param, c( "n_reads_considered", bb$N ) )
df_param = rbind( df_param, c( "n_reads_excluded", read_summary$n_excluded ) )
dimnames(df_param)[[2]] = c("parameter", "value" )
write.table( df_param, file = fn_out_settings, sep='\t', row.names = FALSE, quote = FALSE )
if( opt$verbose ){ cat( paste("MESSAGE: Saved run settings", fn_out_settings,"\n" ) ) }

# Write out the AARDVARK data for optional further analysis
if( opt$write_rdata ){
    if( opt$verbose ){ cat( paste("MESSAGE: Writing AARDVARK R objects in Rdata file", fn_out_Rdata, "\n" ) ) }
    save( reads, read_summary, file = fn_out_Rdata )
    if( opt$verbose ){ cat( paste("MESSAGE: Saved AARDVARK R object file\n" ) ) }
}

if( opt$verbose ){ cat( "MESSAGE: AARDVARK run completed.\n") }

