suppressPackageStartupMessages( require( optparse, quietly=TRUE, warn.conflicts=FALSE ) )
suppressPackageStartupMessages( require( EnsDb.Hsapiens.v86) )
suppressPackageStartupMessages( require( ensembldb ) )
suppressPackageStartupMessages( require( biomaRt ) )
suppressPackageStartupMessages( require( aardvark, quietly=TRUE, warn.conflicts=FALSE ) )

option_list <- list(
    make_option(c("-f", "--fn_summary"), type = "character", help = "Path to AARDVARK reversion summary file [required]"),
    make_option(c("-s", "--pos_start"), type = "numeric", help = "genomic position to start the plot [required]"),
    make_option(c("-e", "--pos_end"), type = "numeric", help = "genomic position to end the plot [required]"),
    make_option(c("-t", "--fig_height"), type = "numeric", help = "height in inches for PDF, default 8", default=8),
    make_option(c("-w", "--fig_width"), type = "numeric", help = "width in inches for PDF, default 10", default=10),
    make_option(c("-o", "--fn_out"), type = "character", help = "Path to PDF to create [required]"),
    make_option(c("-g", "--genome_draft"), type = "numeric", help = "Human genome draft for gene models, one of {19,38}; default: 38", default=38),
    make_option(c("-x", "--proxy_http"), type = "character", help = "Address of proxy server to use for HTTP requests, only relevant if behind a proxy server, default empty"),
    make_option(c("-y", "--proxy_https"), type = "character", help = "Address of proxy server to use for HTTPS requests, only relevant if behind a proxy server, default empty"),
    make_option(c("-l", "--exclude_blacklist"), type = "logical", action="store_true", help = "Exclude reversions that overlap blacklisted regions, default TRUE", default=TRUE),
    make_option(c("-m", "--min_freq"), type = "numeric", help = "minimum number of times to see a reversion allele to plot it, default: 1", default=1)
)

opt = parse_args( OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]") )

if( sum( names(opt)=="fn_summary")==0 ){
    stop( "invalid parameter(s) for input. Input summary file required in parameter fn_summary")
}
if( sum( names(opt)=="pos_start")==0 ){
    stop( "invalid parameter(s) for input. Start position required in parameter pos_start")
}
if( sum( names(opt)=="pos_end")==0 ){
    stop( "invalid parameter(s) for input. End position required in parameter pos_end")
}
if( !file.exists( opt$fn_summary ) ){
    stop( paste0( "Error: input summary file ", opt$fn_summary," not found") )
}
if( opt$genome_draft != 38 & opt$genome_draft != 19 ){
    stop( "human genome draft in parameter genome_draft must be one of 38, 19")
}
if( sum( names(opt)=="proxy_http")==1 & sum( names(opt)=="proxy_https")==1 ){
    cat(paste("MESSAGE: Set https proxy to", opt$proxy_https, "\n") )
    Sys.setenv( http_proxy = opt$proxy_http )
    Sys.setenv( https_proxy = opt$proxy_https )
}

bm <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
sums = read.table( opt$fn_summary, sep='\t', header=TRUE)
pdf( opt$fn_out, height=opt$fig_height, width=opt$fig_width)
if( opt$genome_draft==38 ){
    suppressPackageStartupMessages( require( BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE, warn.conflicts=FALSE ) )
    plot_reversion_summary( sums,
                            genome_version=opt$genome_draft,
                            BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
                            biomart_object=bm,
                            pos_start=opt$pos_start,
                            pos_end = opt$pos_end,
                            exclude_blacklist=opt$exclude_blacklist,
                            min_freq=opt$min_freq)
}else{
    suppressPackageStartupMessages( require( BSgenome.Hsapiens.UCSC.hg19, quietly=TRUE, warn.conflicts=FALSE ) )
    plot_reversion_summary( sums,
                            genome_version=opt$genome_draft,
                            BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                            biomart_object = bm,
                            pos_start = opt$pos_start,
                            pos_end = opt$pos_end,
                            exclude_blacklist = opt$exclude_blacklist,
                            min_freq = opt$min_freq)
}
dev.off()
