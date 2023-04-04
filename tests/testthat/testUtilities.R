context("utility code")

test_that("alignment works", {
    pkg = "aardvark"

    ####################################################################################
    # Test read_constructor
    ####################################################################################

    qname="test"
    pos=32339656
    positions=32339656:(32339656+11)
    seq=DNAString("AAACGATTCTGG")
    qual=c(37,37,37,37,37,37,37,37,37,37,37,37)
    cigar = "6S6M"
    chrom="chr13"
    rd=Read( qname, cigar, chrom, pos, seq, qual )
    expect_equal( dim( rd$cigar_ranges )[1], 2 )
    expect_equal( rd$cigar_ranges$start[1], 1 )
    expect_equal( rd$cigar_ranges$start[2], 7 )
    expect_equal( rd$cigar_ranges$end[1], 6 )
    expect_equal( rd$cigar_ranges$end[2], 12 )
    expect_equal( rd$cigar_ranges$width[1], 6 )
    expect_equal( rd$cigar_ranges$width[2], 6 )
    expect_equal( is.na(rd$cigar_ranges$ref_start[1]), TRUE ) # because softclipped
    expect_equal( is.na(rd$cigar_ranges$ref_end[1]), TRUE ) # because softclipped

    expect_equal( rd$cigar_ranges$ref_start[2], 32339656 )
    expect_equal( rd$cigar_ranges$ref_end[2], 32339661 )


    # this is a classic reversion allele with two deletions
    rd = Read( qname="A00887:299:HWFYGDSXY:2:2235:15501:11929",
                           cigar="21M7D21M2D109M",
                           chrom="chr13",
                           pos=32339609,
                           seq=DNAString("CATTCTGATGAGGTATATAATGATATCTCTCAAAAAATAAACGATTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCACAAACTGTAA"),
                           qual=rep(37, 151) )

    expect_equal( dim(rd$cigar_ranges)[1], 5 )
    expect_equal( rd$cigar_ranges$start[1], 1 )
    expect_equal( rd$cigar_ranges$start[2], 22 )
    expect_equal( rd$cigar_ranges$start[3], 22 )
    expect_equal( rd$cigar_ranges$start[4], 43 )
    expect_equal( rd$cigar_ranges$start[5], 43 )
    expect_equal( rd$cigar_ranges$end[1], 21 )
    expect_equal( rd$cigar_ranges$end[2], 21 )
    expect_equal( rd$cigar_ranges$end[3], 42 )
    expect_equal( rd$cigar_ranges$end[4], 42 )
    expect_equal( rd$cigar_ranges$end[5], 151 )
    expect_equal( rd$cigar_ranges$width[1], 21 )
    expect_equal( rd$cigar_ranges$width[2], 7 )
    expect_equal( rd$cigar_ranges$width[3], 21 )
    expect_equal( rd$cigar_ranges$width[4], 2 )
    expect_equal( rd$cigar_ranges$width[5], 109 )

    ####################################################################################
    # Test replace_row and remove_row
    ####################################################################################

    df = data.frame( a=1:3, b=c("first","second","third"), stringsAsFactors=FALSE)
    dft = remove_row( df, 1 )
    expect_equal( dim(dft)[1], 2 );
    expect_equal( dft$a[1], 2 )

    dft = remove_row( df, 2 )
    expect_equal( dim(dft)[1], 2 );
    expect_equal( dft$a[1], 1 );
    expect_equal( dft$a[2], 3 )

    dft = remove_row( df, 3 )
    expect_equal( dim(dft)[1], 2 );
    expect_equal( dft$a[1], 1 );
    expect_equal( dft$a[2], 2 )

    df_rep = data.frame( a=4:5, b=c("fourth","fifth"), stringsAsFactors=FALSE)
    dft = replace_row( df, 1, df_rep)
    expect_equal( dim(dft)[1], 4 );
    expect_equal( dft$a[1], 4 );
    expect_equal( dft$a[3], 2 )

    dft = replace_row( df, 3, df_rep)
    expect_equal( dim(dft)[1], 4 );
    expect_equal( dft$a[1], 1 );
    expect_equal( dft$a[3], 4 )

    dft = replace_row( df, 2, df_rep)
    expect_equal( dim(dft)[1], 4 );
    expect_equal( dft$a[1], 1 );
    expect_equal( dft$a[2], 4 );
    expect_equal( dft$a[3], 5 ) ;
    expect_equal( dft$a[4], 3 )

    x=reverse_complement( c("A","C","G","T", "T") )
    expect_equal( x[1], "A" )
    expect_equal( x[2], "A" )
    expect_equal( x[3], "C" )
    expect_equal( x[4], "G" )
    expect_equal( x[5], "T" )

    s_ref = "TTCTT"
    s_alt = "T"
    RA = convert.VCF.REFALT.to.dash.format( s_ref, s_alt )
    expect_equal( length( RA ), 2 )
    expect_equal( RA$REF, "TTCTT" )
    expect_equal( RA$ALT, "T----" )
    s_ref = "T"
    s_alt = "TTCTT"
    RA = convert.VCF.REFALT.to.dash.format( s_ref, s_alt )
    expect_equal( RA$REF, "T----" )
    expect_equal( RA$ALT, "TTCTT" )
    s_ref = "ACT"
    s_alt = "GAG"
    RA = convert.VCF.REFALT.to.dash.format( s_ref, s_alt )
    expect_equal( RA$REF, "ACT" )
    expect_equal( RA$ALT, "GAG" )

    # test that read_to_genome_sequence generates string representations
    # this is a classic reversion allele with two deletions
    rd_test = Read( qname="A00887:299:HWFYGDSXY:2:2235:15501:11929",
               cigar="21M7D21M2D109M",
               chrom="chr13",
               pos=32339609,
               seq=DNAString("CATTCTGATGAGGTATATAATGATATCTCTCAAAAAATAAACGATTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCACAAACTGTAA"),
               qual=rep(37, 151) )
    AW_test = aardvark::AlignmentWindow( BSgenome.Hsapiens.UCSC.hg38, chrom="chr13",
                                    32339609 - 300, 32339609 + 300,
                                    min_length_for_homopolymer = 5   )

    nt = read_to_genome_sequence( rd_test, AW_test )
    expect_equal( dim(nt)[1], 160 )
    expect_equal( nt$nt[21], "T" )
    expect_equal( nt$nt[22], "-" )
    expect_equal( nt$nt[28], "-" )
    expect_equal( nt$nt[29], "G" )
    expect_equal( nt$nt[48], "A" )
    expect_equal( nt$nt[49], "C" )
    expect_equal( nt$nt[50], "-" )
    expect_equal( nt$nt[51], "-" )
    expect_equal( nt$nt[ nt$pos==32339658],  "-" )
    expect_equal( nt$nt[ nt$pos==32339659],  "-" )

} )
