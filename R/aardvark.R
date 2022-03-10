# AARDVARK
# An Automated Reversion Detector for Variants Affecting Resistance Kinetics
# Install from source:
# install.packages("/data1/projects/software/aardvark_0.2.0.tar.gz", repos=NULL, type="source")


#' Align a nucleotide sequence in a region of the reference
#'
#' performs a Biostrings::pairwise alignment of query_str and ref_seq
#' and identifies the longest perfet alignment between the two sequences
#'
#' @param ref_pos chromosome position in genome at which original alignment started
#' @param ref_seq DNAString of nucleotides in reference sequence
#' @param query_str character string of nucleotides to align in ref_seq
#' @return list with length of longest alignment, start and end in reference and sequence space, and the alignment object
realign_sequence_in_region = function( ref_pos, ref_seq, query_str ){

    pwa = pairwiseAlignment(DNAString( query_str ),
                            ref_seq,
                            gapOpening=5, gapExtension=5, type="global-local")

    vec_full = strsplit(toString( aligned(pwa) ), "")[[1]]  # aligned including leading/trailing
    vec_ref = strsplit(as.character(ref_seq), "")[[1]]      # aligned including leading/trailing

    # find the largest continuous match in res_vec, meaning not a dash
    idx_start_longest=1
    idx_end_longest=1
    idx_start_cur=1
    idx_end_cur=1
    in_match=FALSE
    for(i in 2:length(vec_full)){
        if( in_match & ( vec_full[i] != vec_ref[i] | vec_full[i] == "-" ) ){
            in_match = FALSE
            if( idx_end_cur - idx_start_cur > idx_end_longest - idx_start_longest ){
                idx_start_longest = idx_start_cur
                idx_end_longest = idx_end_cur
            }
        }else{
            if( vec_full[i] == vec_ref[i] | vec_full[i] == "-" ){
                if( !in_match ){
                    in_match = TRUE
                    idx_start_cur = i
                }
                idx_end_cur = i
            }
        }
    }

    # calculate the deletion that was created between the newly aligned
    # region and the previously aligned portion of the read, which starts at ref_pos
    list( ref_start = ref_pos + idx_start_longest - 1,
          ref_end = ref_pos + idx_end_longest - 1,
          idx_start=idx_start_longest,
          idx_end=idx_end_longest,
          longest_alignment = idx_end_longest - idx_start_longest , # DEBUG removed + 1 because it was returning an impossible result of 16 for ctr 3693
          pwa=pwa)
}

#' attempt to realign a soft-clipped region on the left side of an alignment
#'
#' @param read aardvark::Read object on which to attempt the realignment
#' @param min_percent_realigned minimum percent of the soft-clipped read segment to accept a new alignment
#' @param align_window aardvark::AlignmentWindow object that holds reference sequence search space for realignment
#' @return A modified aardvark::read object
realign_left_softclip = function( read, min_percent_realigned, align_window ){

    # should use the same code for large and small segments, just vary the amount of search space.
    # small segments are searched immediately adjacent, large are searched farther away.
    clipped_width = read$cigar_ranges$end[1] - read$cigar_ranges$start[1] + 1
    clipped_seq = str_sub( as.character( read$seq ), 1, clipped_width )
    align_result = realign_sequence_in_region( ref_pos = align_window$start, align_window$seq, clipped_seq )

    if( align_result$longest_alignment >= clipped_width * min_percent_realigned ){
        existing_row = read$cigar_ranges[1,]
        existing_row$cigar_code = "M"
        existing_row$ref_start = align_result$ref_start
        existing_row$ref_end = existing_row$ref_start + clipped_width - 1
        existing_row$width = clipped_width
        new_row = data.frame( start = clipped_width+1,
                              end = clipped_width+1,
                              width = read$pos - align_result$ref_end - 1,
                              cigar_code_original="S",
                              cigar_code="D",
                              ref_start= align_result$ref_end+1,
                              ref_end = read$pos-1, stringsAsFactors = FALSE )
        read$cigar_ranges = replace_row( read$cigar_ranges, 1, rbind(existing_row, new_row) )
        read$pos = existing_row$ref_start
        read$has_del=TRUE
        read$has_realigned_softclipped=TRUE
        read$has_realigned_softclipped_left=TRUE
    }
    read
}

#'  attempt to realign a soft-clipped region on the right side of an alignment
#'
#' @param read aardvark::Read object on which to attempt the realignment
#' @param min_percent_realigned minimum percent of the soft-clipped read segment to accept a new alignment
#' @param align_window aardvark::AlignmentWindow object that holds reference sequence search space for realignment
#' @return A modified aardvark::ead object
realign_right_softclip = function( read, min_percent_realigned, align_window){

    ii = dim(read$cigar_ranges)[1]
    clipped_width = read$cigar_ranges$end[ii] - read$cigar_ranges$start[ii] + 1
    clipped_seq = str_sub( as.character( read$seq ), read$cigar_ranges$start[ii], read$cigar_ranges$end[ii] )
    align_result = realign_sequence_in_region(ref_pos = align_window$start, align_window$seq, clipped_seq)

    if( align_result$longest_alignment >= clipped_width * min_percent_realigned ){
        # for soft clipping, read position is reported as last non-clipped locus
        # newly detected deletion starts at ref_end and moves to pos, the first aligned position.
        read$cigar_ranges$cigar_code[ii] = "D" # was soft-clipped, reclassifying as deletion
        read$cigar_ranges$ref_start[ii] = read$pos + read$cigar_ranges$start[ii] - 1
        read$cigar_ranges$ref_end[ii] = align_result$ref_start - 1  # immediately before start of newly aligned read
        read$cigar_ranges$width[ii] = align_result$ref_start - read$cigar_ranges$ref_start[ii]
        read$cigar_ranges$end[ii] = read$cigar_ranges$start[ii]
        read$cigar_ranges = rbind(
            read$cigar_ranges,
            data.frame( start = read$cigar_ranges$start[ii],
                        end = read$cigar_ranges$start[ii] + align_result$longest_alignment - 1,
                        width=str_length(clipped_seq),
                        cigar_code_original="S",
                        cigar_code="M",
                        ref_start= align_result$ref_start,
                        ref_end=align_result$ref_end-1), stringsAsFactors = FALSE ) # DEBUG added -1

        read$has_del=TRUE
        read$has_realigned_softclipped=TRUE
        read$has_realigned_softclipped_right=TRUE
    }
    read
}


#' Attempt to merge the softclipped portion of a read segment into the segment for realignment
#'
#' This merge is only saved as the final aligned read if the new alignment is an improvement.
#'
#' @param read aardvark::Read object containing segments to test
#' @param idx_from index of cigar_ranges that is the source of clipping
#' @param idx_into index of cigar_ranges that is the target of merge
#' @return Returns the updated read object
merge_softclipped_into_read = function( read, idx_from, idx_into ){
    if( idx_from > dim(read$cigar_ranges)[1] | idx_into > dim(read$cigar_ranges)[1] ){
        stop( "bounds error: merge_softclipped_into_read() was passed impossible index")
    }
    clipped_width = read$cigar_ranges$end[idx_from] - read$cigar_ranges$start[idx_from] + 1
    if( idx_from < idx_into ){
        read$cigar_ranges$start[ idx_into ] = read$cigar_ranges$start[idx_into] - clipped_width
        read$cigar_ranges$ref_start[ idx_into ] = read$cigar_ranges$ref_start[ idx_into ] - clipped_width
        read$positions[ read$cigar_ranges$start[idx_from] : read$cigar_ranges$end[idx_from] ] =
            ( read$cigar_ranges$ref_start[idx_into] ) : ( read$cigar_ranges$ref_start[idx_into] + clipped_width - 1 )
        read$pos = read$positions[ 1 ]
    }else{
        read$cigar_ranges$end[idx_into] = read$cigar_ranges$end[idx_into] + clipped_width
        read$cigar_ranges$ref_end[ idx_into ] = read$cigar_ranges$ref_end[ idx_into ]
        read$positions[ read$cigar_ranges$start[idx_from] : read$cigar_ranges$end[idx_from] ] =
            ( read$cigar_ranges$ref_end[idx_into] + 1 ) : ( read$cigar_ranges$ref_end[idx_into] + clipped_width )
    }
    read$cigar_ranges$width[idx_into] = read$cigar_ranges$width[idx_into] + clipped_width
    read$cigar_ranges = remove_row( read$cigar_ranges, idx_from )
    read
}


#' attempt to improve the alignment of a read through local realignment
#'
#' @param read aardvark::Read object under consideration
#' @param align_window aardvark::AlignmentWindow object that holds reference sequence search space for realignment
#' @param gr_pathogenic GenomicRanges spanning the pathogenic mutation
#' @param allow_insertions_in_realign whether to allow insertions in the realignment
#' @param min_nt_for_distant_realign minimum size for a local realignment that moves a soft-clipped segment to a distant locus
#' @param min_percent_realigned minimum percentage of the total length of a soft-clipped segment that must be aligned perfectly
#' @param near_bound number of nucleotides adjacent to a germline mutation to consider off of the end of a read for local realignment
#' @param min_nt_qual minimal nucleotide quality score to count when evaluating perfect matches in a candidate alignment
#' @return Returns the updated read object
#' @export
realign_read = function( read,
                         align_window,
                         gr_pathogenic,
                         allow_insertions_in_realign = TRUE,
                         min_nt_for_distant_realign = 15,
                         min_percent_realigned = 0.9,
                         near_bound = 6,
                         min_nt_qual = 20){

    clipped_width = 0
    clipped_left = read$cigar_ranges$cigar_code[1] == "S"
    clipped_right = read$cigar_ranges$cigar_code[ dim( read$cigar_ranges)[1] ] == "S"

    if( clipped_left ){
        clipped_width = read$cigar_ranges$end[1] - read$cigar_ranges$start[1] + 1
        if( clipped_width >= min_nt_for_distant_realign ){
            read = realign_left_softclip( read, min_percent_realigned, align_window )
        }else{
            read_merged = merge_softclipped_into_read( read=read, idx_from=1, idx_into=2 )
            read_merged = locally_realign_read( read_merged, gr_pathogenic,
                                                align_window, near_bound, min_nt_qual,
                                                allow_insertions_in_realign, "left" )
            if( read_merged$has_realigned_softclipped_left ){
                read = read_merged
            }
        }
    }

    if( clipped_right ){
        clipped_width = read$cigar_ranges$end[ dim( read$cigar_ranges)[1] ] -
                        read$cigar_ranges$start[ dim( read$cigar_ranges)[1] ] + 1
        if( clipped_width >= min_nt_for_distant_realign ){
            read = realign_right_softclip( read, min_percent_realigned, align_window )
        }else{
            read_merged = merge_softclipped_into_read( read = read,
                                                       idx_from = dim( read$cigar_ranges)[1],
                                                       idx_into=(dim( read$cigar_ranges)[1]) - 1 )
            read_merged = locally_realign_read( read_merged, gr_pathogenic,
                                                align_window, near_bound, min_nt_qual,
                                                allow_insertions_in_realign, "right" )
            if( read_merged$has_realigned_softclipped_right ){
                read = read_merged
            }
        }
    }
    if( !clipped_left & !clipped_right ){
        # case where read is not clipped, but aligned with mismatches
        read = locally_realign_read( read, gr_pathogenic, align_window, near_bound, min_nt_qual, allow_insertions_in_realign, "both" )
    }

    read$cigar_ranges = read$cigar_ranges[order( read$cigar_ranges$start, read$cigar_ranges$end), ]
    idx_d = which( read$cigar_ranges$cigar_code == "D" )
    if( length(idx_d) > 0 ){
        read$del_lengths = read$cigar_ranges$width[ idx_d ]
    }

    # rewrite positions to reflect any changes made during realignment
    for(i in 1:dim(read$cigar_ranges)[1] ){
        start_end_range = read$cigar_ranges$start[i] : read$cigar_ranges$end[i]
        if( read$cigar_ranges$cigar_code[i] == "M" ){
            read$positions[ start_end_range ] = read$cigar_ranges$ref_start[i] : read$cigar_ranges$ref_end[i]
        }
        if( read$cigar_ranges$cigar_code[i] == "M" | read$cigar_ranges$cigar_code[i] == "S" ){
            read$qualities[ start_end_range ] = read$qualities[ read$cigar_ranges$start[i] : read$cigar_ranges$end[i]  ]
        }
    }
    read$encoding = encode_variants( read )
    read
}


#' Update read properties to match a new local realignment
#'
#' there were two obvious choices for the values to write into ref_start and ref_end in a cigar string:
#' a) cigar ref_start and ref_end refer to the section that was originally there
#' b) cigar ref_start and ref_end refer to the new genomic location
#' I chose option A. Option A means that a 7 base deletion looks like:
#' start end width cigar_code_original cigar_code ref_start  ref_end
#'     1  21    21                   M          M  32339609 32339629
#'    22  21     7                   D          D  32339630 32339636
#'    22  42    21                   M          M  32339637 32339657
#'
#' the end of the deletion on the second line is rolled back to the start (21) because it doesn't take
#' up any read space. The ref_start and ref_end of the deletion indicate which nucleotides were
#' deleted (32339630 through 32339636, seven nucleotides using the UCSC 0-1 addressing).
#'
#' For a 3 base insertion:
#' start end width cigar_code_original cigar_code ref_start  ref_end
#'     1   5     5                   M          M  32339768 32339772
#'     6   8     3                   I          I  32339773 32339773
#'     9 151   143                   M          M  32339773 32339915
#'
#' The start and end of the deletion on the second line are 6-8 because the three inserted bases
#' take up three spaces in the read. The ref_start and reF_end meanwhile are not updated because
#' the insertion does not exist in reference space, so the ref_start for the insertion on line 2
#' and the Match segment on line 3 are the same. The ref_start value indicates the position at
#' which the insertion took place; we use this information in translate_cigar.
#'
#'
#' @param read aardvark::Read to correct
#' @param rr dataframe that represents the new alignment values
#' @return returns the updated read
rebuild_read_from_realignment = function( read, rr ){
    cur_code=""
    read_position_start=0
    ref_position_start=0
    read_position_cur=0
    ref_position_cur=0

    i_start= min( which( rr$read != "-" ) )
    rr = rr[i_start:dim(rr)[1],]
    cigar_new = data.frame()
    cigar_code = rep(NA, dim(rr)[1])
    for(i in 1 : dim(rr)[1] ){
        if( rr$ref[i] == "-" & rr$read[i] != "-" ){
            cigar_code[i] = "I"
        }else if( rr$ref[i] != "-" & rr$read[i] == "-" ){
            cigar_code[i] = "D"
        }else if( rr$ref[i] != "-" & rr$read[i] != "-" ){
            cigar_code[i] = "M"
        }
    }
    cigar_rle = rle( cigar_code )
    starts = c(); ends = c(); widths = c(); ref_starts = c(); ref_ends = c()
    vals = cigar_rle$values
    lens = cigar_rle$lengths
    cur_pos_read = 1
    cur_pos_ref = 1 # tracked separately because insertions advance us in read but not on reference
    ref_offset = 0
    for(i in 1:length(vals)){
        if( vals[i] == "M"){
            starts = c(starts, cur_pos_read)
            ends = c( ends, cur_pos_read + lens[i] - 1 )
            widths = c( widths, lens[i] )
            ref_starts = c( ref_starts, rr$pos[ cur_pos_ref ] + ref_offset )
            ref_ends = c( ref_ends, rr$pos[ cur_pos_ref + lens[i] - 1 ] + ref_offset )
            cur_pos_read = cur_pos_read + lens[i]
            cur_pos_ref = cur_pos_ref + lens[i]
        }else if( vals[i] == "I"){
            starts = c(starts, cur_pos_read)
            ends = c( ends, cur_pos_read + lens[i] - 1 )
            widths = c( widths, lens[i] )
            ref_starts = c(ref_starts, rr$pos[ cur_pos_ref ] + ref_offset )
            ref_ends = c(ref_ends, rr$pos[ cur_pos_ref ] + ref_offset )
            cur_pos_read = cur_pos_read + lens[i]
            read$has_insertion = TRUE
            read$has_indel = TRUE
        }else if( vals[i] == "D" ){
            starts = c(starts, cur_pos_read)
            ends = c( ends, cur_pos_read )
            widths = c(widths, lens[i] )
            ref_starts = c(ref_starts, rr$pos[ cur_pos_ref ] + ref_offset )
            ref_ends = c(ref_ends, rr$pos[ cur_pos_ref + lens[i] - 1 ] + ref_offset )
            ref_offset = ref_offset + lens[i] # cumulative because can have multiple deletions
            read$has_del = TRUE
            read$has_indel = TRUE
            # don't update cur_pos_ref or cur_pos_read
        }
    }
    cigar_new = data.frame(
        start=starts,
        end=ends,
        cigar_code_original=vals,
        cigar_code=vals,
        width=widths,
        ref_start=ref_starts,
        ref_end=ref_ends, stringsAsFactors = FALSE
    )
    read$cigar_ranges = cigar_new
    read$pos = min( read$cigar_ranges$ref_start[ read$cigar_ranges$cigar_code  != "S" ] )
    read
}


#' Local realignment of read ends to improve alignment near a germline mutation
#'
#' Primarily effective when the read ends at or just over the germline mutation
#' if both left and right side are softclipped, need to consider them separately;
#' if neither is softclipped can consider both sides.
#'
#' @param read Read to modify
#' @param gr_mut GenomicRanges location of mutation
#' @param align_window AlignmentWindow object that holds reference sequence search space for realignment
#' @param near_bound number of nucleotides adjacent to a germline mutation to consider off of the end of a read for local realignment
#' @param min_nt_qual minimal nucleotide quality score to count when evaluating perfect matches in a candidate alignment
#' @param allow_insertions_in_realign whether to allow insertions in the realignment
#' @param realign_side which side to attempt realignment, one of "left", "right", "both"
#' @return returns the updated read
locally_realign_read = function( read,
                                 gr_mut,
                                 align_window,
                                 near_bound = 10,
                                 min_nt_qual = 20,
                                 allow_insertions_in_realign = FALSE,
                                 realign_side = "both"){


    # alignment_idx may be a subset of the whole read e.g. if 7S137M7S,
    # and align left, alignment_idx_right would be 8, 144
    # set mutation_search_left or mutation_search_right to reflect region near end
    # of read where mutation might be lurking
    if( ! realign_side=="left" & ! realign_side=="right" & ! realign_side=="both" ){
        stop( "locally_realign_read must be passed either left or right or both" )
    }
    cigar_preserve = data.frame(start=c(),end=c(),width=c(),cigar_code_original=c(),
                                cigar_code=c(),ref_start=c(),ref_end=c(),stringsAsFactors = FALSE)
    n_cigar=dim(read$cigar_ranges)[1]
    if( realign_side == "left" ){
        alignment_idx_left = read$cigar_ranges$start[1]
        alignment_idx_right = read$cigar_ranges$end[1]
        mutation_search_left = read$positions[ alignment_idx_left ] - near_bound
        mutation_search_right = read$positions[ alignment_idx_right ]
        if( n_cigar > 1 ){
            cigar_preserve = read$cigar_ranges[ 2:n_cigar,]
        }
    }else if( realign_side == "right" ){
        alignment_idx_left = read$cigar_ranges$start[ dim(read$cigar_ranges)[1] ]
        alignment_idx_right = read$cigar_ranges$end[ dim(read$cigar_ranges)[1] ]
        mutation_search_left = read$positions[ alignment_idx_left ]
        mutation_search_right = read$positions[ alignment_idx_right ]  + near_bound
        if( n_cigar > 1 ){
            cigar_preserve = read$cigar_ranges[ 1:(n_cigar-1),]
        }
    }else{
        alignment_idx_left = read$cigar_ranges$start[1]
        alignment_idx_right = read$cigar_ranges$end[ dim(read$cigar_ranges)[1] ]
        mutation_search_left = read$positions[ alignment_idx_left ] - near_bound
        mutation_search_right = read$positions[ alignment_idx_right ]  + near_bound
    }
    alignment_idx = alignment_idx_left : alignment_idx_right    # piece of read we'll attempt to realign

    # if mutation overlaps read
    if( ( end(gr_mut) <= mutation_search_right  & end(gr_mut) >= mutation_search_left ) ||
        ( start(gr_mut) >= mutation_search_left & start(gr_mut) <= mutation_search_right ) ){
        vec_read = strsplit( as.character( read$seq ), "")[[1]]
        rr_realign = data.frame(
            pos = read$positions[ alignment_idx ],
            read = vec_read[ alignment_idx ],
            qual = read$qualities[ alignment_idx ],
            ref = rep("-", length( alignment_idx ) ),
            stringsAsFactors = FALSE )
        ref_positions = min(rr_realign$pos[ rr_realign$pos>0] ) : max( rr_realign$pos ) # insertions have position = 0
        ref_in_range = AW_vec( align_window,
                               pos_start = ref_positions[1],
                               pos_end = ref_positions[ length( ref_positions ) ])
        m = match.idx( rr_realign$pos, ref_positions )
        rr_realign$ref[ m$idx.A ] = ref_in_range[ m$idx.B ]
        n_match_pre = sum( rr_realign$read==rr_realign$ref | rr_realign$qual<min_nt_qual )
        idx_left = (1:near_bound)[ rr_realign$qual[1:near_bound] >= min_nt_qual ]
        idx_right = ( dim(rr_realign)[1] - near_bound ) : dim(rr_realign)[1]
        idx_right = idx_right[ rr_realign$qual[ idx_right ] >= min_nt_qual ]

        has_mismatch_at_left =  sum( rr_realign$read[idx_left]  != rr_realign$ref[idx_left] ) > 0
        has_mismatch_at_right = sum( rr_realign$read[idx_right] != rr_realign$ref[idx_right] )  > 0

        if( has_mismatch_at_left | has_mismatch_at_right ){
            # dels extend the reference space covered by the read
            #sum_of_deletions = sum( read$cigar_ranges$width[ read$cigar_ranges$cigar_code=="D"] )
            sum_of_deletions=0 # modified because we're restricting this to just one range with type M
            # extend range for reference to be larger than read by near_bound, only in direction of
            ref_start = rr_realign$pos[1] - sum_of_deletions
            ref_end = rr_realign$pos[ dim(rr_realign)[1] ] + sum_of_deletions
            if( realign_side=="left" | realign_side == "both" ){
                ref_start = ref_start - ( near_bound * 3 ) # size of gap may be larger than interior bound for mutation
            }
            if( realign_side=="right" | realign_side == "both" ){
                ref_end = ref_end + ( near_bound * 3 )
            }
            seq_ref = AW_seq( align_window, pos_start = ref_start, pos_end = ref_end)
            pwa = pairwiseAlignment( DNAString( paste( rr_realign$read, collapse="" ) ),
                                     seq_ref, gapOpening=2, gapExtension=1,
                                     patternQuality = PhredQuality( as.integer( rr_realign$qual  ) ),
                                     subjectQuality = PhredQuality( as.integer( rep(40, length(seq_ref) ) ) ),
                                     type="global-local")

            vec_full = strsplit(toString( aligned( pwa ) ), "")[[1]]      # aligned including leading/trailing
            pos_start = ref_start + min( which( vec_full != "-" ) ) - 1 # first aligned pos with ACGT

            vec_read = strsplit( as.character( pattern(pwa) ), "")[[1]] # aligned (trimmed)
            vec_ref = strsplit( as.character( subject(pwa) ), "")[[1]]  # reference (trimmed)
            rr_combined = data.frame(
                pos = pos_start : ( pos_start + length(vec_ref) - 1 ),
                read = vec_read,
                ref = vec_ref,
                qual = rep(0, length(vec_read)), stringsAsFactors = FALSE )
            rr_combined$qual[ which( rr_combined$read != "-" ) ] = rr_realign$qual
            # low-quality nucleotides that are called inserts rewritten as a mismatch
            # to avoid preserving insert when we call rebuild_read_from_realignment()
            idx_lowqual_inserts = which( rr_combined$read != "-" & rr_combined$ref=="-" & rr_combined$qual < min_nt_qual)
            if( length( idx_lowqual_inserts>0 )){
                for(i in 1:length( idx_lowqual_inserts )){
                    rr_combined$ref[ idx_lowqual_inserts[i] ] = AW_get_nucleotide(align_window,
                                                                                  rr_combined$pos[ idx_lowqual_inserts[i]] )
                }
            }
            n_new_inserts = sum( rr_combined$read != "-" & rr_combined$ref=="-" & rr_combined$qual >= min_nt_qual)

            n_match_post = sum( rr_combined$read == rr_combined$ref |
                                    (rr_combined$read != "-" & rr_combined$qual < min_nt_qual) )
            if( n_match_post > n_match_pre &
                (n_new_inserts==0 | allow_insertions_in_realign) ){
                read = rebuild_read_from_realignment( read, rr_combined )
                read$cigar_ranges = rbind( read$cigar_ranges, cigar_preserve )
                read$cigar_ranges = read$cigar_ranges[order( read$cigar_ranges$start ), ]
                read$has_realigned_read_end = TRUE
                read$has_realigned_softclipped = TRUE
                if( realign_side == "right" ){
                    read$has_realigned_softclipped_right = TRUE
                }
                if( realign_side == "left" ){
                    read$has_realigned_softclipped_left = TRUE
                }
            }
        }
    }
    read
}

#' Classify a read to assess the evidence that it contains a reversion
#'
#' Assesses the evidence that alterations present in a read will result
#' in a functional trancript when the gene is translated. Sequence alterations are used
#' to predict the amino acid sequence that would be generated.
#'
#' @param read aardvark::Read object to evaluate
#' @param transcript aardvark::TranscriptData object for a gene's transcript
#' @param pathogenic aardvark::Mutation object describing pathogenic mutation
#' @param min_nt_qual Minimum nucleotide quality score to write read sequence as part of translated sequence.
#' Reference sequence will be used for nucleotides in the read that are below this quality score.
#' @param gr_pathogenic GenomicRanges object spanning the pathogenic mutation
#' @param gr_exclude GenomicRanges of regions to exclude from consideration, such as regions
#' where homopolymers are present. An AlignmentWindow object automatically calculates a
#' GenomicRangeList of homopolymer regions. Providing this list, or another GenomicRangesList
#' value to gr_exclude can reduce this source of technical artifacts. Single base deletions
#' that fall in gr_exclude will be marked read_harbors_variant_in_excluded_region. Larger
#' alterations are not excluded because they are less likely to be an artifact.
#'
#' @return aardvark::Read object, with updated values for evidence.
#' The AlignmentWindow object automatically calculates a GenomicRangeList of homopolymer regions. Providing this list, or another
#' GenomicRangesList value to gr_exclude can reduce this source of technical artifacts.
#'
#' Possible values for the evidence column:
#' \itemize{
#'   \item read_not_informative
#'   \item read_harbors_variant_in_excluded_region (e.g. a possible homopolymer artifact)
#'   \item read_harbors_pathogenic_variant_no_reversion
#'   \item reversion_read_includes_pathogenic_variant
#'   \item reversion_read_deletion_spans_pathogenic_variant
#'   \item reversion_read_does_not_include_pathogenic_variant
#' }
#' @export
assess_reversion = function( read, transcript, pathogenic, min_nt_qual=20, gr_pathogenic=NULL, gr_exclude=NULL ){

    read$includes_pathogenic_mutation = FALSE
    gene_in_frame = TRUE
    distance_to_germline = NA
    read$evidence = "read_not_informative"
    read$total_frameshift = 0
    read$pathogenic_is_deleted = FALSE
    single_del_in_excluded = FALSE
    read_pos_start = read$pos
    read_pos_end =  max(read$positions)

    # TODO: if the read is 151M and downstream of + (or upstream of - strand) then it cannot cause a
    # reversion, and we can skip translation and converting the quality entirely.
    if( is.null( gr_pathogenic ) ){
        # for high-througput processing, pass this creating the same GR object repeatedly
        gr_pathogenic = aardvark::genomicRangesFromMutation( pathogenic )
    }
    for( i in 1:dim( read$cigar_ranges )[1] ){

        if( read$cigar_ranges$cigar_code[i] == "D" ){
            del_start = max( c( read$cigar_ranges$ref_start[i], transcript$nucleotides$pos[1] ) )
            del_end = min( c( read$cigar_ranges$ref_end[i], transcript$nucleotides$pos[ length( transcript$nucleotides$pos ) ] ) )

            if( !is.null( gr_exclude ) ){
                if( del_start == del_end ){
                    gr_del=GenomicRanges::GRanges(seqnames=read$chrom, IRanges::IRanges(del_start, del_end, names = "del") )
                    if( length( findOverlaps( gr_exclude, gr_del ) ) > 0 ){
                        single_del_in_excluded = TRUE # single base del has del_end == del_start
                    }
                }
            }

            idx_in_ts = which( transcript$nucleotides$pos >= del_start &
                               transcript$nucleotides$pos <= del_end )
            n_coding_nt_deleted = sum( transcript$nucleotides$exon_id[ idx_in_ts ] != "" )
            read$total_frameshift = read$total_frameshift - n_coding_nt_deleted

            if( del_start <= start( gr_pathogenic ) & del_end >= end( gr_pathogenic ) ){
                read$pathogenic_is_deleted = TRUE
            }else if( read$cigar_ranges$width[i]<0) {
                if( del_end <= start( gr_pathogenic ) & del_start >= end( gr_pathogenic ) ){
                    read$pathogenic_is_deleted = TRUE
                }
            }
            if( pathogenic$mutation_class=="deletion" ){
                if( del_start == pathogenic$pos[1] &
                    del_end == pathogenic$location$pos[ length(pathogenic$location$pos ) ] ){
                    # read with reversion that also includes pathogenic mutation is the strongest evidence
                    read$includes_pathogenic_mutation = TRUE
                }
            }
        }else if( read$cigar_ranges$cigar_code[i] == "I"){

            # if the insertion start point is in or adjacent to an exon, framshift
            idx_in_ts = which( transcript$nucleotides$pos == read$cigar_ranges$ref_start[i] )
            if( transcript$nucleotides$exon_id[idx_in_ts] != "" |
                transcript$nucleotides$exon_id[idx_in_ts-1] != "" |
                transcript$nucleotides$exon_id[idx_in_ts+1] != ""){
                read$total_frameshift = read$total_frameshift + read$cigar_ranges$width[i]
            }
            if( pathogenic$mutation_class=="insertion" ){
                if( read$cigar_ranges$ref_start[i] == pathogenic$pos[1] &
                    read$cigar_ranges$ref_end[i] == pathogenic$pos[ length( pathogenic$pos ) ] ){
                    # read with reversion that also includes pathogenic mutation is the strongest evidence
                    read$includes_pathogenic_mutation = TRUE
                }
            }
        }
    }

    if( pathogenic$mutation_class=="insertion" | pathogenic$mutation_class=="deletion" ){
        if( single_del_in_excluded ){
            read$evidence = "read_harbors_variant_in_excluded_region"
        }else{
            if( read$includes_pathogenic_mutation & read$total_frameshift  %% 3 == 0 ){
                read$evidence = "reversion_read_includes_pathogenic_variant"
            }else if( read$total_frameshift != 0 ){
                AA_sequence = aardvark::translate_cigar( transcript, read, pathogenic, min_nt_qual )
                transcript_is_inframe = stringr::str_count( as.character( AA_sequence ), stringr::fixed("*") ) == 1
                if( transcript_is_inframe ){
                    read$evidence = "reversion_read_does_not_include_pathogenic_variant"
                }
            }
        }
    }else if( pathogenic$mutation_class == "missense" ){
        # missense mutations can only be reverted by deleting them or by a missense that
        # directly alters the mutation to eliminate a stop codon.
        AA_sequence = aardvark::translate_cigar( transcript, read, pathogenic, min_nt_qual)
        transcript_is_inframe = stringr::str_count( as.character( AA_sequence ), stringr::fixed("*") ) == 1
        if( transcript_is_inframe ){
            read$evidence = "reversion_read_includes_pathogenic_variant"
        }
    }

    if( read$evidence=="read_not_informative" & read$includes_pathogenic_mutation ){
        read$evidence = "read_harbors_pathogenic_variant_no_reversion"
    }

    read
}

#' Generate the amino acid sequence by applying the alterations present in a read to a given transcript
#'
#' If pathogenic is not passed, the native sequence modified by the value of read is returned.
#' If both read and pathogenic are passed, the impact of both alterations can affect the amino
#' acid sequence. Alterations obsered in the read supercede alterations in pathogenic if the read overlaps
#' the pathogenic alteration.
#' @param transcript an aardvark::TranscriptData object for affected gene
#' @param read aardvark::Read object to evaluate
#' @param pathogenic aardvark::Mutation object describing a pathogenic mutation (optional)
#' @param min_nt_qual Minimum nucleotide quality score to write read sequence as part of translated sequence.
#' Reference sequence will be used for nucleotides in the read that are below this quality score.
#' @return A Biostrings::AAString object representing the amino acid sequence of complete transcript; * is STOP.
#' @export
translate_cigar = function( transcript, read, pathogenic=NULL, min_nt_qual ){
    # create a local copy of nucleotides from transcript
    # edit to reflect the read
    # translate into predicted amino acid sequence

    ts = transcript$nucleotides
    ts$keep = rep(TRUE, dim(ts)[1]) # cheaper than resizing the data frame

    # see the documentation for rebuild_read_from_realignment
    # after I process a cigar_range with an insertion, the genomic positions in any subsequent cigar_range values
    # are now out of date so the wrong nucleotides could get overwritten.

    for( i in 1:dim(read$cigar_ranges)[1] ){
        if( read$cigar_ranges$cigar_code[i] == "D"){
            idx_remove = which( ts$pos[ts$keep] >= read$cigar_ranges$ref_start[i]  &
                                ts$pos[ts$keep] <= read$cigar_ranges$ref_end[i]  )
            ts$keep[idx_remove] = FALSE
            #idx_keep = setdiff( 1:dim(ts)[1], idx_remove )
            #ts = ts[idx_keep,]
        }else if( read$cigar_ranges$cigar_code[i] == "I" ){

            # insertions actually have to modify ts itself, can't just mark ts$keep
            seq_to_insert = strsplit( as.character( read$seq[ read$cigar_ranges$start[i] : read$cigar_ranges$end[i] ] ), "" )[[1]]
            idx_split = which( ts$pos == read$cigar_ranges$ref_start[i]  )
            ts_pre = ts[1:(idx_split-1),]
            ts_post = ts[ idx_split : dim(ts)[1],]
            ts_to_insert = data.frame( pos = rep(0, length(seq_to_insert) ),
                                       seq=seq_to_insert,
                                       exon_id = rep( ts$exon_id[idx_split], length(seq_to_insert)),
                                       strand = rep( ts$strand[idx_split], length(seq_to_insert)),
                                       is_splice = rep( ts$is_splice[ idx_split], length(seq_to_insert)),
                                       keep = rep( TRUE, length(seq_to_insert) ),
                                       stringsAsFactors=FALSE)

            ts = rbind( ts_pre, ts_to_insert, ts_post )

        }else if( read$cigar_ranges$cigar_code[i] == "M" ){
            # M can contain mismatches, update the sequence of the transcript accordingly.

            # find indexes in readspace (e.g. 1:151) of this cigar_range, limited to high quality
            idx_in_readspace = read$cigar_ranges$start[i] : read$cigar_ranges$end[i]
            has_hiqual = read$qualities[ idx_in_readspace ] >= min_nt_qual
            idx_in_readspace = idx_in_readspace[ has_hiqual ]

            # find genomic positions for indexes in readspace and match that to indexes in transcript
            pos_to_update = read$positions[ idx_in_readspace ]
            idx_in_ts = match.idx( ts$pos, pos_to_update )$idx.A
            seq_to_update = strsplit( as.character( read$seq[ idx_in_readspace ] ), "" )[[1]]
            ts$seq[ idx_in_ts ] = seq_to_update
        }
    }

    if( ! is.null(pathogenic) ){
        # if the read doesn't include the pathogenic mutation, write it into the reference sequence.
        # analysis of whether the overall transcript is a reversion is then contingent on the assumption
        # that the alterations present in this read are in phase with the pathogenic mutation. The fact
        # that the alterations are not on the same read as the pathogneic mutation is noted in the evidence
        # value for the read, so the user can make their own judgements
        read_start = min( read$positions[ read$positions>0 ] )
        read_end = max( read$positions )
        mut_start = pathogenic$cigar_ranges$ref_start[1]
        mut_end = pathogenic$cigar_ranges$ref_end[1]
        if( !( ( mut_end <= read_end  & mut_end >= read_start ) ||
               ( mut_start >= read_start & mut_start <= read_end ) ) ){
            # if read ref_start, ref_end doesn't overlap pathogenic
            if( pathogenic$cigar_ranges$cigar_code[1] == "D"){
                idx_remove = which( ts$pos >= pathogenic$cigar_ranges$ref_start[1] &
                                    ts$pos <= pathogenic$cigar_ranges$ref_end[1] )
                #idx_keep = setdiff( 1:dim(ts)[1], idx_remove )
                #ts = ts[idx_keep,]
                ts$keep[idx_remove] = FALSE
            }else if( pathogenic$cigar_ranges$cigar_code[1] == "I" ){
                seq_to_insert = strsplit( as.character( pathogenic$seq ), "" )[[1]]
                idx_split = which( ts$pos == mut_start )
                ts_pre = ts[1:(idx_split-1),]
                ts_post = ts[ idx_split : dim(ts)[1],]
                ts_post$pos = ts_post$pos + length( seq_to_insert )
                ts_to_insert = data.frame( pos = mut_start : mut_end,
                                           seq = seq_to_insert,
                                           exon_id = rep( ts$exon_id[idx_split], length(seq_to_insert)),
                                           strand = rep( ts$strand[idx_split], length(seq_to_insert)),
                                           is_splice = rep( ts$is_splice[ idx_split], length(seq_to_insert)),
                                           stringsAsFactors=FALSE)
                ts = rbind( ts_pre, ts_to_insert, ts_post )
            }else{
                idx_update = which( ts$pos == pathogenic$location$pos[pathogenic$location$cigar_code=="M"] )
                ts$seq[ idx_update ] = pathogenic$location$mut
            }
        }
    }

    # translate exons
    # If indel(s) generated a predicted transcript that is out of frame, avoid "last base was ignored" warning
    # by clipping off the last incomplete set of trailing nucleotides. This isn't strictly accurate since in principle
    # translation could continue out of frame until it hits a stop codon, but in practice any realistic out of frame
    # mutation is going to hit a stop codon far before the end of a transcript.
    final_seq = ts$seq[ ts$keep & ts$exon_id != "" ]
    if( length(final_seq) %% 3 != 0 ){
        final_seq = final_seq[ 1 : ( length(final_seq) - length(final_seq) %% 3 ) ]
    }
    final_seq = DNAString( paste( final_seq, collapse="") )
    if( transcript$strand == "-"){
        # minus strand genes (e.g. BRCA1) are reported by ensembl on their native strand
        # we standardize nucleotides to positive strand in TranscriptData objects to match
        # the orientation of reads from alignment. If gene is negative strand,
        # reverse complement sequence back before translation
        Biostrings::translate( Biostrings::reverseComplement( final_seq ) )
    }else{
        Biostrings::translate( DNAString( paste( final_seq, collapse="") ) )
    }
}

#' Create a string represetation of all cigar codes that distinguish this from perfect match
#'
#' @param rd aardvark::Read object
#' @return String representation of the cigar codes (not CIGAR formatted)
encode_variants = function( rd ){
    code = ""
    if( rd$has_del | rd$has_insertion | rd$has_realigned_softclipped ){
        for(i in 1:dim( rd$cigar_ranges)[1] ){
            cigar_code=rd$cigar_ranges$cigar_code[i]
            if( cigar_code=="D" | cigar_code=="I" ){
                cc=paste( cigar_code,
                          rd$cigar_ranges$width[i],
                          rd$cigar_ranges$ref_start[i],
                          rd$cigar_ranges$ref_end[i], sep=":")
                if( code == "" ){
                    code=cc
                }else{
                    code = paste( code, cc, sep="_" )
                }
            }
        }
    }
    code
}


#' Report candidate reversions from a list of Reads
#'
#' @param reads A list of aardvark::Read objects
#' @param transcript An aardvark::TranscriptData object describing the genomic region being analyzed
#'
#' @return a list with the summary dataframe, cigar ranges for each allele, read names for each allee, and number of unreverted reads.
#' Possible values for the evidence column in the summary dataframe:
#' \itemize{
#'   \item read_not_informative
#'   \item read_harbors_variant_in_excluded_region (e.g. a possible homopolymer artifact)
#'   \item read_harbors_pathogenic_variant_no_reversion
#'   \item reversion_read_includes_pathogenic_variant
#'   \item reversion_read_deletion_spans_pathogenic_variant
#'   \item reversion_read_does_not_include_pathogenic_variant
#' }
#' @export
summarize_candidates = function( reads, transcript ){
    pos = c()
    evidence = c()
    realigned = c()
    n_obs = c()
    list_cigar_ranges = list()
    list_qnames = c()
    hsh_cigar = hsh_new()
    hsh_qnames = hsh_new()
    encodings = c()
    n_unreverted = 0
    n_uninformative = 0
    N = length(reads)
    for( ctr in 1:N ){
        rr = reads[[ ctr ]]
        if( rr$evidence=="reversion_read_includes_pathogenic_variant" |
            rr$evidence=="reversion_read_deletion_spans_pathogenic_variant" |
            rr$evidence=="reversion_read_does_not_include_pathogenic_variant" ){
            if( hsh_in( hsh_cigar, rr$encoding ) ){
                idx = hsh_get( hsh_cigar, rr$encoding )
                list_qnames[[ idx ]] = c( list_qnames[[ idx ]], rr$qname )
                n_obs[idx] = n_obs[ idx ] + 1
            }else{
                encodings = c( encodings, rr$encoding )
                hsh_set( hsh_cigar, rr$encoding, length(encodings) )
                if( transcript$nucleotides$strand[1] == "+" ){
                    pos_first = as.integer( strsplit( strsplit( rr$encoding, "_")[[1]][1], ":")[[1]][3] )
                }else{
                    s = strsplit( rr$encoding, "_")[[1]]
                    s = s[ length(s) ]
                    pos_first = as.integer( strsplit( s, ":")[[1]][3] )
                }
                list_qnames = c( list_qnames, list( c( rr$qname ) ) )
                pos = c( pos, pos_first )
                evidence = c( evidence, rr$evidence )
                n_obs = c(n_obs, 1 )
                list_cigar_ranges = c( list_cigar_ranges, list( rr$cigar_ranges ) )
            }
        }else if( rr$evidence=="read_harbors_pathogenic_variant_no_reversion" ){
            n_unreverted = n_unreverted+1
        }else if( rr$evidence == "read_not_informative" ){
            n_uninformative = n_uninformative + 1
        }
    }
    df_sum = data.frame( N=n_obs,
                         evidence,
                         pos = pos,
                         row.names = encodings,
                         stringsAsFactors = FALSE)
    list( summary = df_sum,
          cigar_ranges = list_cigar_ranges,
          qnames = list_qnames,
          n_unreverted = n_unreverted,
          n_uninformative = n_uninformative)
}
