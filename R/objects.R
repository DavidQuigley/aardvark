
# they're really constructors that return lists, not objects

#' AlignmentWindows hold a region's sequence as both a DNAString and a vector of characters
#'
#'
#' @param Hsapiens_version A \pkg{BSgenome} object used to retreive sequence. See \pkg{BSgenome::available.genomes} for choices.
#' @param chrom A valid chromosome in Hsapiens_version.
#' @param window_start The genomic position on chrom at which to start the region
#' @param window_end The genomic position on chrom at which to end the region
#' @param min_length_for_homopolymer The minimum number of repeated nucleotides to declare a region a homopolymer. Defaults to 5.
#' Homopolymers are repeats of a single nucleotide within the reference sequence.
#' Any single base alteration called by the sequencer that fall within a homopolymer site are likely
#' to be artifacts, and can be discarded.
#' @returns A list with elements start, end, seq, vec. seq is a Biostrings object;
#' @export
AlignmentWindow = function( Hsapiens_version, chrom, window_start, window_end, min_length_for_homopolymer=5 ){

    chrom_ucsc = chrom
    if( length( grep("chr", chrom))==0 ){
        chrom_ucsc = paste0("chr", chrom)
    }
    seq = Biostrings::getSeq( Hsapiens_version, chrom_ucsc, start = window_start, end = window_end)

    r = rle( strsplit( as.character(seq), "" )[[1]] )
    idx_hom = which( r$lengths>= min_length_for_homopolymer)
    idx_start = window_start
    starts = c()
    ends = c()
    for(i in 1:length(r$lengths)){
        if( r$lengths[i]>= min_length_for_homopolymer ){
            starts = c( starts, idx_start )
            ends = c( ends, idx_start + r$lengths[i] - 1 )
        }
        idx_start = idx_start + r$lengths[i]
    }
    homopolymer_regions = GenomicRanges::makeGRangesFromDataFrame(
        data.frame( chrom = rep( chrom, length(starts)),
                    start=starts,
                    end=ends,
                    stringsAsFactors=FALSE ) )

    list( start=window_start,
          end=window_end,
          seq = seq,
          chrom = chrom,
          homopolymer_regions=homopolymer_regions,
          vec = strsplit( as.character(seq), "")[[1]] )
}


#' helper function to return a slice of an AlignmentWindow character vector between pos_start and pos_end
#'
#' If pos_start or pos_end is out of bounds for AW, message the problem and return NULL
#' @param AW Alignment Window object to examine
#' @param pos_start first position to return
#' @param pos_end last position to return
#' @returns vector of characters within AW$vec
AW_vec = function( AW, pos_start, pos_end){
    if( pos_start < AW$start ){
        message( paste("request for position",pos_start,"out of range, AlignmentWindow range starts at",AW$start) )
        NULL
    }else if( pos_end > AW$end ){
        message( paste("request for position", pos_end, "out of range, AlignmentWindow range ends at",AW$end) )
        NULL
    }else{
        idx_start = pos_start - AW$start + 1
        idx_end = pos_end - AW$start + 1
        AW$vec[ (pos_start - AW$start + 1) : ( pos_end - AW$start + 1 ) ]
    }
}

#' helper function to return a slice of an AlignmentWindow DNAString between pos_start and pos_end
#'
#' @param AW Alignment Window object to examine
#' @param pos_start first position to return
#' @param pos_end last position to return
#' @returns DNAString subsetted by pos_start, pos_end
AW_seq = function( AW, pos_start, pos_end ){
    if( pos_start < AW$start ){
        stop( paste("request for position",pos_start,"out of range, AlignmentWindow range starts at",AW$start) )
    }
    if( pos_end > AW$end ){
        stop( paste("request for position", pos_end, "out of range, AlignmentWindow range ends at",AW$end) )
    }
    AW$seq[ (pos_start - AW$start + 1) : ( pos_end - AW$start + 1 ) ]
}

#' helper function to return a specific nucleotide in AlignmentWindow DNAString at pos
#'
#' @param AW Alignment Window object to examine
#' @param pos position to return
#' @returns character at pos
AW_get_nucleotide = function( AW, pos ){
    if( pos < AW$start ){
        stop( paste("request for position", pos, "out of range, AlignmentWindow range starts at",AW$start) )
    }
    if( pos > AW$end ){
        stop( paste("request for position", pos, "out of range, AlignmentWindow range ends at",AW$end) )
    }
    as.character( AW$seq[ (pos - AW$start + 1) ] )
}

#' convert UTF format to integers, used for sanger conversion
#' @param utf8 input string
#' @return converted to integer
convert_UTF = function( utf8 ){
    utf8ToInt( as.character( utf8 ) ) - 33 # sanger conversion
}

#' Load a BAM file into memory
#'
#' This function is a thin wrapper around Rsamtools::scanBam to return reads in a bam file.
#' You can call scanBam directly if you want different choices for the flags. The flag settings
#' are:
#'
#' \itemize{
#'   \item isDuplicate=FALSE
#'   \item isPaired=TRUE
#'   \item isUnmappedQuery=FALSE
#'   \item isSecondaryAlignment=FALSE
#'   \item isSupplementaryAlignment=FALSE
#'   \item isNotPassingQualityControls=FALSE
#' }
#'
#' @param fn_bam path to Bam file
#' @param chrom chromosome to read
#' @param start read start position
#' @param end read end position
#' @returns list of reads (data frame) and qualities (list)
#' @export
BamData = function( fn_bam, chrom, start, end ){
    if( ! file.exists( fn_bam ) ){
        stop( paste( "cannot open file:", fn_bam ) )
    }
    flag1 = scanBamFlag( isDuplicate=FALSE,
                         isPaired=TRUE,
                         isUnmappedQuery=FALSE,
                         isSecondaryAlignment=FALSE,
                         isSupplementaryAlignment=FALSE,
                         isNotPassingQualityControls=FALSE)
    bam_which = GenomicRanges::makeGRangesFromDataFrame( data.frame( chrom=chrom, start=start, stop=end) )
    bam_fields=c("qname", "rname", "pos", "strand", "qual", "cigar", "mapq", "seq")
    bam_param = ScanBamParam(flag=flag1,
                             what=bam_fields,
                             which=bam_which)
    bb = scanBam(fn_bam,  param=bam_param)
    qualities=lapply( bb[[1]]$qual, convert_UTF)
    reads = data.frame(
        qname = bb[[1]]$qname,
        cigar = bb[[1]]$cigar,
        chrom = as.character( bb[[1]]$rname ),
        pos = bb[[1]]$pos,
        seq = as.character( bb[[1]]$seq ),
        stringsAsFactors=FALSE)
    keep = which( reads$pos >= start & reads$pos <= end )
    reads = reads[keep,]
    qualities = qualities[ keep ]
    list( reads=reads, qualities=qualities, N = dim(reads)[1] )
}


#' Instantiate a Read from a list of reads
#'
#' This function consumes a BamData object and the index of a read within that object and
#' generates a single aardvark::Read object.
#'
#' @param bamdata BamData object
#' @param read_index the index of the read to return.
#' @returns An aardvark::read object
#' @export
read_from_BamData = function( bamdata, read_index ){
    # wrapper to make it easier for end user to iterate through ScanBam lists
    aardvark::Read( qname = bamdata$reads$qname[ read_index ],
                    cigar = bamdata$reads$cigar[ read_index ],
                    chrom = bamdata$reads$chrom[ read_index ],
                    pos = bamdata$reads$pos[ read_index ],
                    seq = DNAString( bamdata$reads$seq[ read_index ] ),
                    qual = bamdata$qualities[[ read_index ]] )
}

#' Generate an object that holds information about a read
#'
#' The Read object can be decorated with additional information by
#' generated through local realignment.
#'
#' @param qname read name, usually assigned by the sequencer
#' @param cigar CIGAR string
#' @param chrom chromosome
#' @param pos read start position, first non-clipped nucleotide
#' @param seq DNAString object that represents the nucleotide sequence as
#' @param qual vector of quality scores for nucleotides in seq; must have same length as seq.
#' @returns a list with elements that describe the read
#' @seealso \code{\link{read_from_BamData}}
#' @seealso \code{\link{assess_reversion}}
#' @seealso \code{\link{realign_read}}
#' @export
Read = function( qname, cigar, chrom, pos, seq, qual ){
    rd = list( qname = qname,
               cigar = cigar,
               cigar_original = cigar,
               chrom = chrom,
               pos = pos,
               pos_original = pos,
               seq = seq,
               qualities = qual,
               cigar_rle = GenomicAlignments::cigarToRleList( cigar )[[1]],
               positions = rep(0, length( seq ) ),
               has_del =                       FALSE,
               has_realigned_softclipped=      FALSE,
               has_realigned_softclipped_left= FALSE,
               has_realigned_softclipped_right=FALSE,
               has_realigned_read_end =        FALSE,
               has_insertion =                 FALSE)

    rd$cigar_ranges = data.frame( cigarRangesAlongQuerySpace( rd$cigar )[[1]] )
    rd$cigar_ranges$cigar_code_original = runValue( rd$cigar_rle )
    rd$cigar_ranges$cigar_code =          rd$cigar_ranges$cigar_code_original
    rd$cigar_ranges$width =               rep(NA, dim(rd$cigar_ranges)[1] )
    rd$cigar_ranges$ref_start =           rep(NA, dim(rd$cigar_ranges)[1] )
    rd$cigar_ranges$ref_end =             rep(NA, dim(rd$cigar_ranges)[1] )

    if( length(rd$qualities) != length( seq )){
        stop( paste("read qual vector length", length(qual), "different from read seq length", length(seq)))
    }

    # set ref_start, ref_end, width for M, D, and I cigar_ranges elements
    # populate cigar_ranges data frame with genome positions for cigar_ranges elements
    cumulative_adjustment = 0

    # insertions do not advance us in the ref_start or ref_end.
    # the value recorded for both ref_start and ref_end is the insertion point, which we need to know later
    # when we translate the insertion into an amino acid sequence.
    insertion_adjustment = 0
    if( dim(rd$cigar_ranges)[1] > 0 ){
        for(ii in 1:dim(rd$cigar_ranges)[1] ){
            if( rd$cigar_ranges$cigar_code[ii] == "S" ){
                cumulative_adjustment = cumulative_adjustment - rd$cigar_ranges$end[ii]
                rd$cigar_ranges$width[ii] = rd$cigar_ranges$end[ii] - rd$cigar_ranges$start[ii] + 1
            }
            if( rd$cigar_ranges$cigar_code[ii] == "M" | rd$cigar_ranges$cigar_code[ii]=="D"){
                rd$cigar_ranges$ref_start[ii] = rd$pos - 1 + rd$cigar_ranges$start[ii] + cumulative_adjustment
                rd$cigar_ranges$ref_end[ii] = rd$cigar_ranges$ref_start[ii] + rd$cigar_rle@lengths[ii] - 1
                rd$cigar_ranges$width[ii] = rd$cigar_ranges$ref_end[ii] - rd$cigar_ranges$ref_start[ii] + 1
            }
            if( rd$cigar_ranges$cigar_code[ii]=="D" ){
                cumulative_adjustment = cumulative_adjustment + rd$cigar_ranges$width[ii]
                rd$has_del=TRUE
            }
            if( rd$cigar_ranges$cigar_code[ii] == "I" ){
                rd$cigar_ranges$ref_start[ii] = rd$pos + rd$cigar_ranges$start[ii] - 1
                rd$cigar_ranges$ref_end[ii] = rd$cigar_ranges$ref_start[ii]
                rd$cigar_ranges$width[ii] =  rd$cigar_rle@lengths[ii]
                cumulative_adjustment = cumulative_adjustment - rd$cigar_ranges$width[ii]
                rd$has_insertion = TRUE
            }
        }

        # check that cumulative adjustment is consistent with the sequence length
        length_calc_from_cigar = 0
        for(ii in 1:dim(rd$cigar_ranges)[1] ){
            if( rd$cigar_ranges$cigar_code[ii] == "M" | rd$cigar_ranges$cigar_code[ii] == "I" | rd$cigar_ranges$cigar_code[ii] == "S" ){
                length_calc_from_cigar = length_calc_from_cigar + rd$cigar_ranges$width[ii]
            }
        }
        if( length_calc_from_cigar != length( seq ) ){
            stop( paste0( "Read ", qname, " sequence length inconsistent with cigar annotation" ) )
        }

    }
    vec_qual = rd$qualities
    for(i in 1:dim(rd$cigar_ranges)[1] ){
        start_end_range = rd$cigar_ranges$start[i] : rd$cigar_ranges$end[i]
        if( rd$cigar_ranges$cigar_code[i] == "M" ){
            rd$positions[ start_end_range ] = rd$cigar_ranges$ref_start[i] : rd$cigar_ranges$ref_end[i]
        }
        if( rd$cigar_ranges$cigar_code[i] == "M" | rd$cigar_ranges$cigar_code[i] == "S" ){
            rd$qualities[ start_end_range ] = vec_qual[ rd$cigar_ranges$start[i] : rd$cigar_ranges$end[i]  ]
        }
    }

    rd
}


#' Generate an object that holds genome location information about transcript
#'
#' This function uses biomaRt to fetch relevant gene model data about a transcript.
#'
#' @param biomart_object instantated biomart object
#' @param ensdb_object instantiated ensemblDB object
#' @param transcript_id ENBEMBL identity of a transcript to look up using biomaRt
#' @return a list with elements that describe the transcript: transcript sequence as a DNAString object, the strand (+ or -) and a data frame of the nucleotides listing pos, seq (nucleotide), exon_id, strand, is_splice (whether the position is a splice junction)
#' @export
TranscriptData = function( biomart_object, ensdb_object, transcript_id ){
    transcript = biomaRt::getSequence( id=transcript_id,
                              type = "ensembl_transcript_id",
                              seqType = "transcript_exon_intron",
                              mart = biomart_object)
    if( dim(transcript)[1]==0 ){
        stop( paste( "transcript", transcript_id, "not found in ensembl database. Do not include versions (e.g. use ENST00000380152 not ENST00000380152.5)" ) )
    }

    tr_seq = strsplit( transcript$transcript_exon_intron, "" )[[1]]
    tr_exons = data.frame( ensembldb::cdsBy( ensdb_object, filter = ~ tx_id == transcript_id) )
    transcript_loc = biomaRt::getBM(attributes=c('chromosome_name','transcript_start','transcript_end','strand'),
                           filters = 'ensembl_transcript_id',
                           values = transcript_id,
                           mart = biomart_object)
    strand = "+"
    tr_pos = transcript_loc$transcript_start : transcript_loc$transcript_end
    if( transcript_loc$strand != 1 ){
        strand = "-"
        tr_seq = reverse_complement(tr_seq)
    }

    tr_pos = transcript_loc$transcript_start : transcript_loc$transcript_end
    tr_exon_id = rep("", length(tr_pos))
    tr_splice_junction = rep(FALSE, length(tr_pos))
    for(i in 1:dim(tr_exons)[1]){
        tr_exon_id[ tr_pos>=tr_exons$start[i] & tr_pos<=tr_exons$end[i] ] = tr_exons$exon_id[i]
    }
    tr_splice_junction[ (tr_exons$start - transcript_loc$transcript_start + 1 - 1) ] = TRUE
    tr_splice_junction[ (tr_exons$start - transcript_loc$transcript_start + 1 - 2) ] = TRUE
    tr_splice_junction[ (tr_exons$end - transcript_loc$transcript_start + 1 + 1) ] = TRUE
    tr_splice_junction[ (tr_exons$end - transcript_loc$transcript_start + 1 + 2) ] = TRUE

    list( sequence = DNAString( transcript$transcript_exon_intron ),
          strand=strand,
          nucleotides = data.frame( pos=tr_pos,
                                    seq=tr_seq,
                                    exon_id=tr_exon_id,
                                    strand=rep(strand, length(tr_seq)),
                                    is_splice=tr_splice_junction,
                                    stringsAsFactors=FALSE),
          transcript_id = transcript_id
    )
}

#' Generate an object that holds information about a mutation.
#'
#' seq_ref and seq_alt are strings that should hold either a dash or a nucleotide.
#' To indicate a two nucleotide deletion, where the reference sequence that was deleted is AC,
#' pass "AC" for seq_ref and "--" for seq_alt To indicate a one nucleotide insertion, where
#' the mutation inserts a T, pass "-" for seq_ref and "T" for seq_alt To indicate a single
#' nucleotide variation where the reference sequence was a "C" and the mutation converts the C to
#' an, pass "C" for seq_ref and "A" for seq_alt
#'
#' @param chrom chromosome
#' @param pos start position of the mutation in on chromsosome chrom
#' @param seq_ref character string indicating the reference sequence
#' @param seq_alt character string indicating the mutation sequence
#' @param transcript an aardvark::TranscriptData object used to assess the impact of the mutation
#' @returns a list with elements that describe the consequences of the mutation
#' @export
Mutation = function( chrom, pos, seq_ref, seq_alt, transcript ){

    # VCF format would encode a two nt deletion as REF=CTT ALT=C pos=32339657
    # we store that as REF=TT ALT=-- pos 32339658, with incremented position
    # if first nt is not the same letter, this is a mismatch, don't change anything

    # currently assumes first
    vec_ref = strsplit( seq_ref, "")[[1]]
    vec_alt = strsplit( seq_alt, "")[[1]]
    n_ref = length( vec_ref )
    n_alt = length( vec_alt )
    is_splice_variant = FALSE
    total_frameshift = 0
    pos_original = pos

    if( vec_ref[1] == vec_alt[1] ){
        if( n_ref < n_alt ){
            vec_ref = c( vec_ref, rep("-", n_alt-n_ref) ) # insertion, pad ref
        }else if( n_ref > n_alt ){
            vec_alt = c( vec_alt, rep("-", n_ref-n_alt) ) # deletion, pad alt
        }
        # trim leading nucleotide and advance to position where different
        # this is only relevant if variant is an indel; if missense, should use pos_original
        vec_ref = vec_ref[2:length(vec_ref)]
        vec_alt = vec_alt[2:length(vec_alt)]
        pos = pos + 1
    }
    seq_ref = paste( vec_ref, sep="", collapse="" )
    seq_alt = paste( vec_alt, sep="", collapse="" )
    cigar = data.frame()
    rr = data.frame( pos = rep(0, length(seq_alt)),
                     ref = strsplit( seq_ref, "")[[1]],
                     mut = strsplit( seq_alt, "")[[1]] )

    cigar_code = rep(NA, dim(rr)[1])
    cur_pos = pos

    for(i in 1 : dim(rr)[1] ){
        if( rr$ref[i] == "-" & rr$mut[i] != "-" ){
            cigar_code[i] = "I"
            mutation_class = "insertion"
        }else if( rr$ref[i] != "-" & rr$mut[i] == "-" ){
            cigar_code[i] = "D"
            rr$pos[i] = cur_pos
            cur_pos = cur_pos+1
            mutation_class = "deletion"
        }else if( rr$ref[i] != "-" & rr$mut[i] != "-" & rr$ref[i] != rr$mut[i] ){
            cigar_code[i] = "M"
            rr$pos[i] = cur_pos
            cur_pos = cur_pos+1
            mutation_class = "mismatch"
        }
    }
    rr$cigar_code = cigar_code

    if( mutation_class=="mismatch" ){
        pos = pos_original
    }
    cigar_ranges = data.frame(
        start=c(1),
        end=c(dim(rr)[1]),
        width=c(dim(rr)[1]),
        cigar_code_original = cigar_code[1],
        cigar_code = cigar_code[1],
        ref_start = rr$pos[1],
        ref_end = rr$pos[ dim(rr)[1] ],
        stringsAsFactors = FALSE)

    # assess whether this mutation impacts a splice site
    idx_in_ts = which( transcript$nucleotides$pos >= cigar_ranges$ref_start[1] &
                       transcript$nucleotides$pos <= cigar_ranges$ref_end[1] )
    is_splice_variant = sum( transcript$nucleotides$is_splice[ idx_in_ts ] ) > 0

    if( mutation_class=="deletion"){
        # if deletion affects splice site it doesn't cause a frameshift
        total_frameshift = -1 * sum( transcript$nucleotides$exon_id[ idx_in_ts ] != "" )
    }else if( mutation_class=="insertion"){
        total_frameshift = sum( rr$cigar_code=="I" )
    }
    list( chrom = chrom,
          pos = pos,
          cigar_rle = rle( cigar_code ),
          location = rr,
          mutation_class = mutation_class,
          cigar_ranges = cigar_ranges,
          is_splice_variant = is_splice_variant,
          total_frameshift = total_frameshift )
}


#' Helper function to convert a mutation to a GenomicRanges object
#'
#' @param germline_mutation mutation object
#' @returns a GenomicRanges object
#' @export
genomicRangesFromMutation = function( germline_mutation ){
    GenomicRanges::makeGRangesFromDataFrame( data.frame(chrom=germline_mutation$chrom,
                                         start=germline_mutation$pos,
                                         end=germline_mutation$pos + dim(germline_mutation$location)[1] - 1,
                                         stringsAsFactors = FALSE ),
                              ignore.strand = TRUE )
}
