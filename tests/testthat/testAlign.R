context("alignment code")

#create data for lazy loading
#ensembl = useDataset("hsapiens_gene_ensembl", mart = useMart("ensembl") )
#transcript_BRCA2 = TranscriptData( ensembl, EnsDb.Hsapiens.v86, "ENST00000380152" )
#transcript_BRCA1 = TranscriptData( ensembl, EnsDb.Hsapiens.v86, "ENST00000357654" )
#transcript_ATM = TranscriptData( ensembl, EnsDb.Hsapiens.v86, "ENST00000452508" )
#save( transcript_BRCA2, file='/notebook/code/aardvark/data/transcript_BRCA2.RData', compress="bzip2")
#save( transcript_BRCA1, file='/notebook/code/aardvark/data/transcript_BRCA1.RData', compress="bzip2")
#save( transcript_ATM, file='/notebook/code/aardvark/data/transcript_ATM.RData', compress="bzip2")
test_that("alignment works", {
    pkg = "aardvark"


################################################################################################################
# Testing code for annotate_read()
################################################################################################################


GM = aardvark::Mutation( chrom="chr13", pos=32339657, seq_ref = "CTT", seq_alt = "C", transcript=transcript_BRCA2)
gr_pathogenic = aardvark::genomicRangesFromMutation(GM)
AW = aardvark::AlignmentWindow( BSgenome.Hsapiens.UCSC.hg38, chrom="chr13",
                                GM$pos - 3000, GM$pos + 3000,
                                min_length_for_homopolymer = 5   )


# Test homopolymer detection in the region chr13:
expect_equal( AW$start, 32336658 )
expect_equal( AW$end, 32342658 )
expect_equal( length( AW$homopolymer_regions ), 62)
expect_equal( start(AW$homopolymer_regions)[1], 32336786)
expect_equal( end(AW$homopolymer_regions)[1], 32336790)

# soft clipped split read deletion
qq=rep(37, 151)
qq[70]=25
qq[145]=25
rd = aardvark::Read( qname="A00887:299:HWFYGDSXY:2:2674:25211:28682",
                     cigar="72S79M",
                     chrom="chr13",
                     pos=32340564,
                     seq=DNAString("CAGCCTTAGCTTTTTACACAAGTTGTAGTAGAAAAACTTCTGTGAGTCAGACTTCATTACTTGAAGCAAAAAAAAGTTCCTTACACAAAGTTAAGGGAGTGTTAGAGGAATTTGATTTAATCAGAACTGAGCATAGTCTTCACTATTCACC"),
                     qual=qq )

rd=aardvark::realign_read( rd, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)

expect_equal( rd$cigar_ranges$ref_start[1], 32339355)
expect_equal( rd$cigar_ranges$ref_end[1],   32339426)
expect_equal( rd$cigar_ranges$ref_start[2], 32339427)
expect_equal( rd$cigar_ranges$ref_end[2],   32340563)
expect_equal( rd$cigar_ranges$ref_start[3], 32340564)
expect_equal( rd$cigar_ranges$ref_end[3],   32340642)
expect_equal( rd$qualities[70], 25)
expect_equal( rd$qualities[145], 25)





# chr13:32,339,355-32,339,426
# chr13:32,340,564-32,340,642
# with 1137 (divis by 3)

# this is a classic reversion allele with two deletions
read = aardvark::Read( qname="A00887:299:HWFYGDSXY:2:2235:15501:11929",
                       cigar="21M7D21M2D109M",
                       chrom="chr13",
                       pos=32339609,
                       seq=DNAString("CATTCTGATGAGGTATATAATGATATCTCTCAAAAAATAAACGATTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCACAAACTGTAA"),
                       qual=rep(37, 151) )
rd=aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)

# to test with BLAT: chr13:32339609-32339768
expect_equal( rd$cigar_ranges$ref_start[1], 32339609 ); expect_equal( rd$cigar_ranges$ref_end[1], 32339629 )
expect_equal( rd$cigar_ranges$ref_start[2], 32339630 ); expect_equal( rd$cigar_ranges$ref_end[2], 32339636 )
expect_equal( rd$cigar_ranges$ref_start[3], 32339637 ); expect_equal( rd$cigar_ranges$ref_end[3], 32339657 )
expect_equal( rd$cigar_ranges$ref_start[4], 32339658 ); expect_equal( rd$cigar_ranges$ref_end[4], 32339659 )
expect_equal( rd$cigar_ranges$ref_start[5], 32339660 ); expect_equal( rd$cigar_ranges$ref_end[5], 32339768 )
expect_equal( rd$positions[1],   32339609)
expect_equal( rd$positions[21],  32339629)
expect_equal( rd$positions[22],  32339637)
expect_equal( rd$positions[22],  32339637)
expect_equal( rd$positions[42],  32339657)
expect_equal( rd$positions[43],  32339660)
expect_equal( rd$positions[151], 32339768)

# three nt insertion
read = aardvark::Read(qname="A00887:299:HWFYGDSXY:2:1533:16161:23719",
                      cigar="5M3I143M",
                      chrom="chr13",
                      pos=32339768,
                      seq=DNAString("AATGATTCAGATATTTGCGTTGAGGAACTTGTGACTAGCTCTTCACCCTGCAAAAATAAAAATGCAGCCATTAAATTGTCCATATCTAATAGTAATAATTTTGAGGTAGGGCCACCTGCATTTAGGATAGCCAGTGGTAAAATCGTTTGTG"),
                      qual=rep(37, 151) )
rd=aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)

expect_equal( rd$cigar_ranges$ref_start[1], 32339768 ); expect_equal( rd$cigar_ranges$ref_end[1], 32339772 )
expect_equal( rd$cigar_ranges$ref_start[2], 32339773 ); expect_equal( rd$cigar_ranges$ref_end[2], 32339773 )
expect_equal( rd$cigar_ranges$ref_start[3], 32339773 ); expect_equal( rd$cigar_ranges$ref_end[3], 32339915 )
expect_equal( rd$positions[1], 32339768)
expect_equal( rd$positions[9], 32339773)
expect_equal( rd$positions[6], 0 )

# TODO: test right side soft clipped reads

################################################################################################################
# Testing code for aardvark::realign_read_end_at_germline_deletion()
################################################################################################################

# four nt with perfect matches
# AATAAACTTGATTCTGG  ref
# AATAAAC--GATTCTGG  germline
#      AAACGATTCTGG  read
#    AAAC--GATTCTGG  outcome
read = aardvark::Read(
            qname="test",
            chrom="chr13",
            pos=32339656,
            seq=DNAString("AAACGATTCTGG"),
            qual=c(37,37,37,37,37,37,37,37,37,37,37,37),
            cigar="12M")
rd=aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)
expect_equal( dim( rd$cigar_ranges)[1],3 )
expect_equal( rd$cigar_ranges$start[1], 1 )
expect_equal( rd$cigar_ranges$end[1], 4 )
expect_equal( rd$cigar_ranges$ref_start[1], 32339654 )
expect_equal( rd$cigar_ranges$ref_end[1], 32339657 )
expect_equal( rd$cigar_ranges$start[2], 5 )
expect_equal( rd$cigar_ranges$end[2], 5 )
expect_equal( rd$cigar_ranges$ref_start[2], 32339658 )
expect_equal( rd$cigar_ranges$ref_end[2], 32339659 )
expect_equal( rd$cigar_ranges$start[3], 5 )
expect_equal( rd$cigar_ranges$end[3], 12 )
expect_equal( rd$cigar_ranges$ref_start[3], 32339660 )
expect_equal( rd$cigar_ranges$ref_end[3], 32339667 )

# four nt with mismatch that has disqualifying low quality
# AATAAACTTGATTCTGG  ref
# AATAAAC--GATTCTGG  germline
#      AAACGATTCTGG  read
#    AAGC--GATTCTGG  outcome
read = aardvark::Read(
    qname="test",
    chrom="chr13",
    pos=32339656,
    seq=DNAString("AAACGATTCTGG"),
    qual=c(37,37,11,37,37,37,37,37,37,37,37,37),
    cigar="12M")

rd=aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)

expect_equal( dim( rd$cigar_ranges)[1], 3 )
expect_equal( rd$cigar_ranges$start[1], 1 )
expect_equal( rd$cigar_ranges$end[1], 4 )
expect_equal( rd$cigar_ranges$ref_start[1], 32339654 )
expect_equal( rd$cigar_ranges$ref_end[1], 32339657 )

expect_equal( rd$cigar_ranges$start[2], 5 )
expect_equal( rd$cigar_ranges$end[2], 5 )
expect_equal( rd$cigar_ranges$ref_start[2], 32339658 )
expect_equal( rd$cigar_ranges$ref_end[2], 32339659 )

expect_equal( rd$cigar_ranges$start[3], 5 )
expect_equal( rd$cigar_ranges$end[3], 12 )
expect_equal( rd$cigar_ranges$ref_start[3], 32339660 )
expect_equal( rd$cigar_ranges$ref_end[3], 32339667 )

# two nt perfect match overhang
# AATAAACTTGATTCTGGTA  ref
# AATAAAC--GATTCTGGTA  germline
#        ACGATTCTGGTA  read
#      AC--GATTCTGGTA  post

read = aardvark::Read(
    qname="test",
    chrom="chr13",
    pos=32339658,
    seq=DNAString("ACGATTCTGGTA"),
    qual=c(37,37,37,37,37,37,37,37,37,37,37,37),
    cigar="12M")
rd=aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)

expect_equal( dim( rd$cigar_ranges)[1],3 )
expect_equal( rd$cigar_ranges$start[1], 1 ); expect_equal( rd$cigar_ranges$end[1], 2 ); expect_equal( rd$cigar_ranges$ref_start[1], 32339656 ); expect_equal( rd$cigar_ranges$ref_end[1], 32339657 );
expect_equal( rd$cigar_ranges$start[2], 3 ); expect_equal( rd$cigar_ranges$end[2], 3 ); expect_equal( rd$cigar_ranges$ref_start[2], 32339658 ); expect_equal( rd$cigar_ranges$ref_end[2], 32339659 )
expect_equal( rd$cigar_ranges$start[3], 3 ); expect_equal( rd$cigar_ranges$end[3], 12 ); expect_equal( rd$cigar_ranges$ref_start[3], 32339660 ); expect_equal( rd$cigar_ranges$ref_end[3], 32339669 )


# does not intersect germline locus, should't do anything
# AATAAACTTGATTCTGGTATT  ref
# AATAAAC--GATTCTGGTATT  germline
#          GATTCTGGTATT  read
#          GATTCTGGTATT  outcome

read = aardvark::Read(
    qname="test",
    chrom="chr13",
    pos=32339660,
    seq=DNAString("GATTCTGGTATT"),
    qual=c(37,37,37,37,37,37,37,37,37,37,37,37),
    cigar="12M")
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)

expect_equal( dim( rd$cigar_ranges)[1],1 )
expect_equal( rd$cigar_ranges$start[1], 1 )
expect_equal( rd$cigar_ranges$end[1], 12 )
expect_equal( rd$cigar_ranges$ref_start[1], 32339660 )
expect_equal( rd$cigar_ranges$ref_end[1], 32339671 )

# intersects germline locus with WT sequence, should't do anything
# AATAAACTTGATTCTGGTA  ref
# AATAAAC--GATTCTGGTA  germline
#        TTGATTCTGGTA  read
#        TTGATTCTGGTA  outcome
read = aardvark::Read(
    qname="test",
    chrom="chr13",
    pos=32339658,
    seq=DNAString("TTGATTCTGGTA"),
    qual=c(37,37,37,37,37,37,37,37,37,37,37,37),
    cigar="12M")
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)

expect_equal( dim( rd$cigar_ranges)[1],1 )
expect_equal( rd$cigar_ranges$start[1], 1 )
expect_equal( rd$cigar_ranges$end[1], 12 )
expect_equal( rd$cigar_ranges$ref_start[1], 32339658 )
expect_equal( rd$cigar_ranges$ref_end[1], 32339669 )

# intersects germline locus with wrong sequence (leading GG), should't do anything
# AATAAACTTGATTCTGGTA  ref
# AATAAAC--GATTCTGGTA  germline
#        GGGATTCTGGTA  read
read = aardvark::Read(
    qname="test",
    chrom="chr13",
    pos=32339658,
    seq=DNAString("GGGATTCTGGTA"),
    qual=c(37,37,37,37,37,37,37,37,37,37,37,37),
    cigar="12M")
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)

expect_equal( dim( rd$cigar_ranges)[1],1 )
expect_equal( rd$cigar_ranges$start[1], 1 )
expect_equal( rd$cigar_ranges$end[1], 12 )
expect_equal( rd$cigar_ranges$ref_start[1], 32339658 )
expect_equal( rd$cigar_ranges$ref_end[1], 32339669 )

# intersects germline locus perfect match with adjacent deletion
# AATAAACTTGATTCTGGTA  ref
# AATAAAC--GATTCTGGTA  germline
#    TA--ACGATTCTGGTA  read
#   TAA-C--GATTCTGGTA  outcome

read = aardvark::Read(
    qname="test",
    chrom="chr13",
    pos=32339654,
    seq=DNAString("TAACGATTCTGGTA"),
    qual=c(37,37,37,37,37,37,37,37,37,37,37,37,37,37),
    cigar="2M2D12M")
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)

expect_equal( dim( rd$cigar_ranges)[1], 5 )

expect_equal( rd$cigar_ranges$start[1], 1 )
expect_equal( rd$cigar_ranges$end[1], 3 )
expect_equal( rd$cigar_ranges$ref_start[1], 32339653 )
expect_equal( rd$cigar_ranges$ref_end[1], 32339655 )
expect_equal( rd$cigar_ranges$start[2], 4 )
expect_equal( rd$cigar_ranges$end[2], 4 )
expect_equal( rd$cigar_ranges$ref_start[2], 32339656 )
expect_equal( rd$cigar_ranges$ref_end[2], 32339656 )
expect_equal( rd$cigar_ranges$start[3], 4)
expect_equal( rd$cigar_ranges$end[3], 4 )
expect_equal( rd$cigar_ranges$ref_start[3], 32339657 )
expect_equal( rd$cigar_ranges$ref_end[3], 32339657 );
expect_equal( rd$cigar_ranges$start[4], 5 )
expect_equal( rd$cigar_ranges$end[4], 5 )
expect_equal( rd$cigar_ranges$ref_start[4], 32339658 )
expect_equal( rd$cigar_ranges$ref_end[4], 32339659 )
expect_equal( rd$cigar_ranges$start[5], 5 )
expect_equal( rd$cigar_ranges$end[5], 14 )
expect_equal( rd$cigar_ranges$ref_start[5], 32339660 )
expect_equal( rd$cigar_ranges$ref_end[5], 32339669 )


# right side 3 nt with perfect matches
# AAAATAAACTTGATTCTGGTA reference
# AAAATAAAC--GATTCTGGTA germline
# AAAATAAACGAT read
# AAAATAAAC--GAT read
read = aardvark::Read(
    qname="test",
    chrom="chr13",
    pos=32339649,
    seq=DNAString("AAAATAAACGAT"),
    qual=c(37,37,37,37,37,37,37,37,37,37,37,37),
    cigar="12M")
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)

expect_equal( dim( rd$cigar_ranges)[1],3 )
expect_equal( rd$cigar_ranges$start[1], 1 )
expect_equal( rd$cigar_ranges$end[1], 9 )
expect_equal( rd$cigar_ranges$ref_start[1], 32339649 )
expect_equal( rd$cigar_ranges$ref_end[1], 32339657 )
expect_equal( rd$cigar_ranges$start[2], 10 )
expect_equal( rd$cigar_ranges$end[2], 10 )
expect_equal( rd$cigar_ranges$ref_start[2], 32339658 )
expect_equal( rd$cigar_ranges$ref_end[2], 32339659 )
expect_equal( rd$cigar_ranges$start[3], 10 )
expect_equal( rd$cigar_ranges$end[3], 12 )
expect_equal( rd$cigar_ranges$ref_start[3], 32339660 )
expect_equal( rd$cigar_ranges$ref_end[3], 32339662 )

# right side intersects germline with wrong sequence called at high confidence
# AAAATAAACTTGAT reference
# AAAATAAAC--GAT germline
# AAAATAAACCAT read
# AAAATAAAC[CA]T best

read = aardvark::Read(
    qname="test",
    chrom="chr13",
    pos=32339649,
    seq=DNAString("AAAATAAACCAT"),
    qual=c(37,37,37,37,37,37,37,37,37,37,37,37),
    cigar="12M")
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)

expect_equal( dim( rd$cigar_ranges)[1], 3 )
expect_equal( rd$cigar_ranges$start[1], 1 )
expect_equal( rd$cigar_ranges$end[1], 9 )
expect_equal( rd$cigar_ranges$ref_start[1], 32339649 )
expect_equal( rd$cigar_ranges$ref_end[1], 32339657 )
expect_equal( rd$cigar_ranges$start[2], 10)
expect_equal( rd$cigar_ranges$end[2], 11 )
expect_equal( rd$cigar_ranges$ref_start[2], 32339658 )
expect_equal( rd$cigar_ranges$ref_end[2], 32339658 )
expect_equal( rd$cigar_ranges$cigar_code[2]=="I", TRUE)
expect_equal( rd$cigar_ranges$start[3], 12 )
expect_equal( rd$cigar_ranges$end[3], 12 )
expect_equal( rd$cigar_ranges$ref_start[3], 32339658 )
expect_equal( rd$cigar_ranges$ref_end[3], 32339658 )



# right side intersects germline with wrong sequence called at low confidence
# AAAATAAACTTGAT reference
# AAAATAAAC--GAT germline
# AAAATAAACGAA read
# AAAATAAAC--GAA read
read = aardvark::Read(
    qname="test",
    chrom="chr13",
    pos=32339649,
    seq=DNAString("AAAATAAACGAA"),
    qual=c(37,37,37,37,37,37,37,37,37,37,37,11),
    cigar="12M")
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=FALSE)

expect_equal( dim( rd$cigar_ranges)[1],3 )
expect_equal( rd$cigar_ranges$start[1], 1 )
expect_equal( rd$cigar_ranges$end[1], 9 )
expect_equal( rd$cigar_ranges$ref_start[1], 32339649 )
expect_equal( rd$cigar_ranges$ref_end[1], 32339657 )
expect_equal( rd$cigar_ranges$start[2], 10 )
expect_equal( rd$cigar_ranges$end[2], 10 )
expect_equal( rd$cigar_ranges$ref_start[2], 32339658 )
expect_equal( rd$cigar_ranges$ref_end[2], 32339659 )
expect_equal( rd$cigar_ranges$start[3], 10 )
expect_equal( rd$cigar_ranges$end[3], 12 )
expect_equal( rd$cigar_ranges$ref_start[3], 32339660 )
expect_equal( rd$cigar_ranges$ref_end[3], 32339662 )

# mutation that is too far away for realignment, note different pos
# AAAATAAACTTGAT reference
# AAAATAAAC--GAT germline
# AAAATAAACGAA read
# AAAATAAAC--GAA read
read = aardvark::Read(
    qname="test",
    chrom="chr13",
    pos=31339649,
    seq=DNAString("AAAATAAACGAA"),
    qual=c(37,37,37,37,37,37,37,37,37,37,37,11),
    cigar="12M")
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)
expect_equal( dim( rd$cigar_ranges)[1], 1)
expect_equal( rd$cigar_ranges$cigar_code[1], "M")
expect_equal( rd$cigar_ranges$width[1], 12)

# germline but soft-clipped edge
# AATAAACTTGATTCTGGTA  ref
# AATAAAC--GATTCTGGTA  germline
#        ACGATTCTGGTA  read
#      AC--GATTCTGGTA  outcome
read = aardvark::Read(
    qname="test",
    chrom="chr13",
    pos=32339660,
    seq=DNAString("ACGATTCTGGTA"),
    qual=c(37,37,37,37,37,37,37,37,37,37,37,37),
    cigar="2S10M")
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)

expect_equal( dim( rd$cigar_ranges)[1],3 )
expect_equal( rd$cigar_ranges$start[1], 1 )
expect_equal( rd$cigar_ranges$end[1], 2 )
expect_equal( rd$cigar_ranges$ref_start[1], 32339656 )
expect_equal( rd$cigar_ranges$ref_end[1], 32339657 )

expect_equal( rd$cigar_ranges$start[2], 3 )
expect_equal( rd$cigar_ranges$end[2], 3 )
expect_equal( rd$cigar_ranges$ref_start[2], 32339658 )
expect_equal( rd$cigar_ranges$ref_end[2], 32339659 )

expect_equal( rd$cigar_ranges$start[3], 3 )
expect_equal( rd$cigar_ranges$end[3], 12 )
expect_equal( rd$cigar_ranges$ref_start[3], 32339660 )
expect_equal( rd$cigar_ranges$ref_end[3], 32339669 )

# Test code in assess_reversion to check for homopolymers
# this deletion is in a homopolymer and should be marked as blacklisted
rd = aardvark::Read( qname="A00887:299:HWFYGDSXY:2:1309:14787:31704",
                     cigar="110M1D41M",
                     chrom="chr13",
                     pos=32339590,
                     seq = DNAString("GTCTAACAGCTATTCCTACCATTCTGATGAGGTATATAATGATTCAGGATATCTCTCAAAAAATAAACTTGATTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAG"),
                     qual = rep(37, 151) )

rd = aardvark::assess_reversion(rd, transcript_BRCA2, pathogenic=GM, align_window=AW, min_nt_qual=20, gr_pathogenic = gr_pathogenic, gr_exclude = AW$homopolymer_regions)
expect_equal( rd$evidence, "read_harbors_variant_in_excluded_region" )

# this is a classic reversion allele with two deletions
# there are multiple reasonable alignments
# REF CATTCTGATGAGGTATATAATGATTCAGGATATCTCTCAAAAAATAAACTTGATTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCAC
#      H  S  D  E  V  Y  N  D  S  G  Y  L  S  K  N  K  R  F  W  Y
# ALT CATTCTGATGAGGTATATAATGAT       ATCTCTCAAAAAATAAAC  GATTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCACAAACTGTAA
read = aardvark::Read( qname="A00887:299:HWFYGDSXY:2:2235:15501:11929",
                       cigar="24M4D2M2D125M",
                       chrom="chr13",
                       pos=32339609,
                       seq=DNAString("CATTCTGATGAGGTATATAATGATATCTCTCAAAAAATAAACGATTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCACAAACTGTAA"),
                       qual=rep(37, 151) )
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)
rd = aardvark::assess_reversion(rd, transcript_BRCA2, GM, AW, min_nt_qual=20, gr_pathogenic = gr_pathogenic, gr_exclude = AW$homopolymer_regions)
expect_equal( rd$evidence, "reversion_read_includes_pathogenic_variant" )

# This looks like it might be a reversion allele since there's a deletion that would put the ORF back in-frame,
# but it puts the sequence back in frame AFTER the pathogenic stop codon so it's not a reversion
# REF CATTCTGATGAGGTATATAATGATTCAGGATATCTCTCAAAAAATAAACTTGATTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCAC
#      H  S  D  E  V  Y  N  D  S  G  Y  L  S  K  N  K  R  F  W  Y
# MUT CATTCTGATGAGGTATATAATGATTCAGGATATCTCTCAAAAAATAAAC--GATTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCAC
#      H  S  D  E  V  Y  N  D  S  G  Y  L  S  K  N  K    R  F  W  Y  *
# RD                                                     GATTCTGGTATTGAGCCAGTATT-AAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCACA
read = aardvark::Read( qname="A00887:299:HWFYGDSXY:2:2235:15501:11929",
                       cigar="23M1D77M",
                       chrom="chr13",
                       pos=32339660,
                       seq=DNAString("GATTCTGGTATTGAGCCAGTATTAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCACA"),
                       qual=rep(37, 100) )
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)
rd = aardvark::assess_reversion(rd, transcript_BRCA2, GM, AW,  min_nt_qual=20, gr_pathogenic = gr_pathogenic, gr_exclude = AW$homopolymer_regions)
expect_equal( rd$evidence, "read_not_informative" )


# This is a reversion allele, with a deletion that does not include the pathogenic mutation
# REF CATTCTGATGAGGTATATAATGATTCAGGATATCTCTCAAAAAATAAACTTGATTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCAC
#      H  S  D  E  V  Y  N  D  S  G  Y  L  S  K  N  K  R  F  W  Y
# MUT CATTCTGATGAGGTATATAATGATTCAGGATATCTCTCAAAAAATAAAC--GATTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCAC
#      H  S  D  E  V  Y  N  D  S  G  Y  L  S  K  N  K    R  F  W  Y  *
# RD                                                     GA-TCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCACA
#      H  S  D  E  V  Y  N  D  S  G  Y  L  S  K  N  K    R   S  G  I  E  P  V  L  K  N  V  E  D  Q  K  N  T  S  F  S  K  V  I  S  N  V  K  D  A  N  A  Y  P
read = aardvark::Read( qname="A00887:299:HWFYGDSXY:2:2235:15501:11929",
                       cigar="2M1D98M",
                       chrom="chr13",
                       pos=32339660,
                       seq=DNAString("GATCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCACA"),
                       qual=rep(37, 100) )
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)
rd = aardvark::assess_reversion(rd, transcript_BRCA2, GM, AW, min_nt_qual=20, gr_pathogenic = gr_pathogenic, gr_exclude = AW$homopolymer_regions)
expect_equal( rd$evidence, "reversion_read_does_not_include_pathogenic_variant" )


# intersects germline locus perfect match
# AATAAACTTGATTCTGGTA  ref
# AATAAAC--GATTCTGGTA  germline
#    AAAC--GATTCTGGTA  read


read = aardvark::Read(
    qname="test",
    chrom="chr13",
    pos=32339654,
    seq=DNAString("AAACGATTCTGGTA"),
    qual=c(37,37,37,37,37,37,37,37,37,37,37,37,37,37),
    cigar="4M2D10M")
rd = aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)
rd = aardvark::assess_reversion(rd, transcript_BRCA2, GM, AW, min_nt_qual=20, gr_pathogenic = gr_pathogenic, gr_exclude = AW$homopolymer_regions)
expect_equal( rd$evidence, "read_harbors_pathogenic_variant_no_reversion")
expect_equal( rd$pathogenic_is_deleted, FALSE)


read = aardvark::Read( qname="A00887:299:HWFYGDSXY:2:2235:15501:11929",
                       cigar="73M1D28M",
                       chrom="chr11",
                       pos=108288903,
                       seq=DNAString("ATTTTGGAAGTTCACTGGTCTATGAACAAAACTTTTTAAAACGATGACTGTATTTTTTCCCTTAACTCTGTTAGGATTTGGATCCTGCTCCTAATCCACCA"),
                       qual=rep(37, 101) )
GM = aardvark::Mutation( chrom="chr11", pos=108288975, seq_ref = "AG", seq_alt = "A", transcript=transcript_ATM)
gr_pathogenic = aardvark::genomicRangesFromMutation(GM)
AW = aardvark::AlignmentWindow( BSgenome.Hsapiens.UCSC.hg38, chrom="chr11",
                                GM$pos - 300, GM$pos + 300,
                                min_length_for_homopolymer = 5   )
rd=aardvark::realign_read( read, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)
rd = assess_reversion( rd, transcript_ATM, GM, AW)

# Test code in assess_reversion to check for homopolymers
# this insertion is in a homopolymer and should be marked as blacklisted
rd = aardvark::Read( qname="A00887:299:HWFYGDSXY:2:1309:14787:31704",
                     cigar="44M1D28M1I27M",
                     chrom="chr13",
                     pos=32379813,
                     seq = DNAString("ACAGAATTTATCATCTTGCAACTTCAAAATCTAAAAGTAAATCTAAAGAGCTAACATACAGTTAGCAGCGACAAAAAAAAACTCAGTATCAACAACTACC"),
                     qual = rep(37, 100) )
GM = aardvark::Mutation( chrom="chr13", pos=32379813, seq_ref = "TG", seq_alt = "T", transcript=transcript_BRCA2)
gr_pathogenic = aardvark::genomicRangesFromMutation(GM)
AW = aardvark::AlignmentWindow( BSgenome.Hsapiens.UCSC.hg38, chrom="chr13",
                                GM$pos - 300, GM$pos + 300,
                                min_length_for_homopolymer = 5   )
rd = aardvark::assess_reversion(rd, transcript_BRCA2, pathogenic=GM, align_window=AW, min_nt_qual=20, gr_pathogenic = gr_pathogenic, gr_exclude = AW$homopolymer_regions)
expect_equal( rd$evidence, "read_harbors_variant_in_excluded_region" )



GM = aardvark::Mutation( chrom="chr17", pos=43093898, seq_ref = "CTT", seq_alt = "C", transcript=transcript_BRCA1)
gr_pathogenic = aardvark::genomicRangesFromMutation(GM)
AW = aardvark::AlignmentWindow( BSgenome.Hsapiens.UCSC.hg38, chrom="chr17",
                                GM$pos - 300, GM$pos + 300,
                                min_length_for_homopolymer = 5   )
rd = aardvark::Read( qname="HWI-D00108:1315:HFG52BCX3:1:1104:17461:26128",
                     cigar="79M21S",
                     chrom="chr17",
                     pos=43093800,
                     seq=DNAString("TTCTTTTTCGAGTGATTCTATTGGGTTAGGATTTTTCTCATTCTGAATAGAATCACCTTTTGTTTTATTCTCATGACCATTCTGCTCCGTTTGGTTAGTT"),
                     qual=rep(37, 100))
rd=aardvark::realign_read( rd, AW, gr_pathogenic, allow_insertions_in_realign=TRUE)
expect_equal( rd$cigar_ranges$ref_start[1], 43093800)
expect_equal( rd$cigar_ranges$ref_end[1],   43093878)
expect_equal( rd$cigar_ranges$ref_start[2], 43093879)
expect_equal( rd$cigar_ranges$ref_end[2],   43093905) # in the actual read this is 43093904
expect_equal( rd$cigar_ranges$ref_start[3], 43093906) # in the actual read this is 43093905
expect_equal( rd$cigar_ranges$ref_end[3],   43093925)
expect_equal( rd$cigar_ranges$width[2], 27)
rd = aardvark::assess_reversion(rd, transcript_BRCA1, pathogenic=GM, align_window=AW, min_nt_qual=20, gr_pathogenic = gr_pathogenic, gr_exclude = AW$homopolymer_regions)
expect_equal( rd$evidence, "reversion_read_deletion_spans_pathogenic_variant")

} )
