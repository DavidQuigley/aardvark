context("mutation code")
test_that("mutation object works", {
    pkg = "aardvark"

    chromosome="chr13"
    pos = 32339657
    REF="CTT"
    ALT="C"
    m = Mutation( chromosome, pos, REF, ALT )
    expect_equal( m$chrom, "chr13" )
    expect_equal( m$pos, pos+1)
    expect_equal( dim(m$cigar_ranges)[1], 1 )
    expect_equal( m$cigar_ranges$cigar_code[1], "D" )
    expect_equal( m$cigar_ranges$ref_start[1], pos+1 )
    expect_equal( m$cigar_ranges$ref_end[1], pos+2 )
    expect_equal( m$total_frameshift, -2 )

    pos = 32339657
    REF="C"
    ALT="CTT"
    m = Mutation( chromosome, pos, REF, ALT )
    expect_equal( m$chrom, "chr13" )
    expect_equal( m$pos, pos+1)
    expect_equal( dim(m$cigar_ranges)[1], 1 )
    expect_equal( m$cigar_ranges$cigar_code[1], "I" )
    expect_equal( m$cigar_ranges$ref_start[1], 0 )
    expect_equal( m$cigar_ranges$ref_end[1], 0 )
    expect_equal( m$total_frameshift, 2 )

    pos = 32339657
    REF="C"
    ALT="G"
    m = Mutation( chromosome, pos, REF, ALT )
    expect_equal( m$chrom, "chr13" )
    expect_equal( m$pos, pos)
    expect_equal( dim(m$cigar_ranges)[1], 1 )
    expect_equal( m$cigar_ranges$cigar_code[1], "M" )
    expect_equal( m$cigar_ranges$ref_start[1], 32339657 )
    expect_equal( m$cigar_ranges$ref_end[1], 32339657 )
    expect_equal( m$total_frameshift, 0 )

} )
