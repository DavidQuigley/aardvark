context("translation code")

test_that("translation works", {
    pkg = "aardvark"
    min_nt_qual = 20
    pathogenic = aardvark::Mutation(chrom="chr13", pos = 32339657,
                                  seq_ref = "CTT",
                                  seq_alt = "C")

    # read overlaps mutation and reports no mutation, so should be eliminated
    pos=32339628
    seq=DNAString("ATGATTCAGGATATCTCTCAAAAAATAAACTTGATTCTGGTATTGAGCCAGTATTGAAGA")
    qual=rep(37, length(seq) )
    cigar = "60M"
    chrom="chr13"
    read=Read( "test", cigar, chrom, pos, seq, qual )
    aa=translate_cigar( transcript_BRCA2, read, pathogenic, min_nt_qual )
    expect_equal( str_count(as.character(aa), stringr::fixed("*")), 1 )
    expect_equal( length(aa), str_locate(as.character(aa), stringr::fixed("*") )[1] )

    # read that does not overlap pathogenic mutation, expect frameshifts
    pos=32338587
    seq=DNAString("CTACTAAAACGGAGCAAAATATAAAAGATTTTGAGACTTCTGATACATTTTTTCAGACTG")
    qual=rep(37, length(seq) )
    cigar = "60M"
    chrom="chr13"
    read=Read( "test", cigar, chrom, pos, seq, qual )
    aa=translate_cigar( transcript_BRCA2, read, pathogenic, min_nt_qual )
    expect_equal( str_count(as.character(aa), stringr::fixed("*")), 167 )
    expect_equal( str_locate(as.character(aa), stringr::fixed("*") )[1], 1772 )

    # single base deletion reversion outside of pathogenic mutation
    #
    #
    # REF AATGATTCAGGATATCTCTCAAAAAATAAACTTGATTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGA
    #      N  D  S  G  Y  L  S  K  N  K  L  D  S  G  I  E  P  V  L  K  N  V  E
    # GL  AATGATTCAGGATATCTCTCAAAAAATAAACGA--TTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGA
    #      N  D  S  G  Y  L  S  K  N  K  R    F  W  Y  *  A
    # REV ATGA-TTCAGGATATCTCTCAAAAAATAAACGA--TTCTGGTATTGAGCCAGTATTGAAGAATGTTGAAGA
    #      M  I  Q  D  I  S  Q  K  I  N  D    S  G  I  E  P  V  L  K  N  V  E
    pos=32339628
    seq=DNAString("TGATTCAGGATATCTCTCAA")
    qual=rep(37, length(seq) )
    cigar = "1D20M"
    chrom="chr13"
    read=Read( "test", cigar, chrom, pos, seq, qual )
    aa=translate_cigar( transcript_BRCA2, read, pathogenic, min_nt_qual )
    expect_equal( str_count(as.character(aa),  stringr::fixed("*") ), 1 )
    expect_equal( length(aa), str_locate(as.character(aa), stringr::fixed("*") )[1] )

    # one base insertion outside of pathogenic mutation that should not revert
    # REF ATGATTCAGGATATCTCTCA
    #     N  D  S  G  Y  L  S
    # ALT GATGATTCAGGATATCTCTC
    #     R  *  I

    pos=32339628
    seq=DNAString("GATGATTCAGGATATCTCTC")
    qual=rep(37, length(seq) )
    cigar = "1I19M"
    chrom="chr13"
    read=Read( "test", cigar, chrom, pos, seq, qual )
    aa=translate_cigar( transcript_BRCA2, read, pathogenic, min_nt_qual )
    expect_equal( str_count(as.character(aa), stringr::fixed("*")), 86 )
    expect_equal( str_locate(as.character(aa), stringr::fixed("*") )[1], 1759 )


    #2I19M implies 21 total nucleotides, while sequence is only 19 nucleotides
    expect_error( Read( "test", "2I17M", "chr13", 32339628, DNAString("GAATGATTCAGGATATCTCT"), qual )  )

    # two base insertion outside of pathogenic mutation that should  revert
    pos=32339628
    seq=DNAString("GAATGATTCAGGATATCTCT")
    qual=rep(37, length(seq) )
    cigar = "2I18M"
    chrom="chr13"
    read=Read( "test", cigar, chrom, pos, seq, qual )
    aa=translate_cigar( transcript_BRCA2, read, pathogenic, min_nt_qual )
    expect_equal( str_count(as.character(aa), stringr::fixed("*")), 1 )
    expect_equal( length(aa), str_locate(as.character(aa), stringr::fixed("*") )[1] )


    # test translation with perfect match that doesn't intersect pathogenic read and
    # and the pathogenic allele present, should introduce stop codon
    #ATGCCTATTGGATCCAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAGCA
    pos=32316461
    seq=DNAString("ATGCCTATTGGATCCAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAGCA")
    qual=rep(37, length(seq) )
    cigar = "66M"
    chrom="chr13"
    rd=Read( "test", cigar, chrom, pos, seq, qual )
    aa=translate_cigar( transcript_BRCA2, rd, pathogenic, min_nt_qual )
    expect_equal( length(aa), 3418 )
    expect_equal( as.character( aa[1:5] ), "MPIGS" )
    expect_equal( str_count( as.character(aa), stringr::fixed("*") ), 167 )


    # test translation with perfect match that also includes a "N" with low quality,
    # and the pathogenic allele present, should introduce stop codon
    # NTGCCTATTGGATCCAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAGCA
    pos=32316461
    seq=DNAString("NTGCCTATTGGATCCAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAGCA")
    qual=rep(37, length(seq) )
    qual[1] = 2 # first nucleotide is low quality, should be ignored
    cigar = "66M"
    chrom="chr13"
    rd=Read( "test", cigar, chrom, pos, seq, qual )
    aa=translate_cigar( transcript_BRCA2, rd, pathogenic, min_nt_qual )
    expect_equal( length(aa), 3418 )
    expect_equal( as.character( aa[1:5] ), "MPIGS" )
    expect_equal( str_count( as.character(aa), stringr::fixed("*") ), 167 )


    # test translation with no pathogenic read, should introduce stop codon
    pos=32316461
    # REF ATGCCTATTGGATCCAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAAC
    #      M  P  I  G  S  K  E  R  P  T  F  F  E  I  F  K  T  R  C  N
    # ALT ATGCCTATTGGAT--AAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAGCAGG
    #      M  P  I  G    *
    seq=DNAString("ATGCCTATTGGATAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAGCAGG")
    qual=rep(37, length(seq) )
    cigar = "13M2D53M"
    chrom="chr13"
    rd=Read( "test", cigar, chrom, pos, seq, qual )
    aa=translate_cigar( transcript_BRCA2, rd, min_nt_qual=20 )
    expect_equal( length(aa), 3418 )
    expect_equal( as.character( aa[1:5] ), "MPIG*" )

    # REF ATGCCTATTGGATCCAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAAC
    #      M  P  I  G  S  K  E  R  P  T  F  F  E  I  F  K  T  R  C  N
    # ALT ATGCCTATTGGATAGCCAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAG
    #      M  P  I  G  *

    #                           AG
    seq=DNAString("ATGCCTATTGGATAGCCAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAG")
    pos=32316461
    qual=rep(37, length(seq) )
    cigar = "13M2I51M"
    chrom="chr13"
    rd = Read( "test", cigar, chrom, pos, seq, qual )
    aa=translate_cigar( transcript_BRCA2, rd, pathogenic, min_nt_qual ) # < there's a bug here?
    expect_equal( length(aa), 3419 )
    expect_equal( as.character( aa[1:5] ), "MPIG*" )

    # test negative strand
    pathogenic = aardvark::Mutation(chrom="chr17",
                                  pos = 43093367,
                                  seq_ref = "CAA",
                                  seq_alt="C")

    # remove deletion, should produce original amino acid sequence
    qname="test"
    pos=43093338
    seq=DNAString("TTTTTCTTCTCTTGGAAGGCTAGGATTGACAAATTCTTTAAGTTCACTGGTATTTGAACA")
    qual=rep(37, length(seq) )
    cigar = "60M"
    chrom="chr17"
    rd=Read( qname, cigar, chrom, pos, seq, qual )
    aa=translate_cigar( transcript_BRCA1, rd, pathogenic, min_nt_qual = 20 )
    expect_equal( as.character( head(aa )), "MDLSAL")
    expect_equal( str_count(as.character(aa), stringr::fixed("*")), 1 )
    expect_equal( length(aa), str_locate(as.character(aa), stringr::fixed("*") )[1] )


    # test that pathogenic mutations that change sequence actually affect translation
    # The mutation should put a TAA onto amino acid 1767
    pathogenic = aardvark::Mutation(chrom="chr13",
                                    pos = 32339654,
                                    seq_ref = "A",
                                    seq_alt="T")
    # read overlaps mutation and reports no mutation, we have changed the pathogenic mutation
    seq=DNAString("AAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATGCAAATGCATACCCAC")
    qual=rep(37, length(seq) )
    rd=Read( "test", "60M", "chr13", 32339700, seq, qual )
    aa=translate_cigar( transcript_BRCA2, rd, pathogenic, min_nt_qual = 20 )
    expect_equal( str_count(as.character(aa), stringr::fixed("*")), 2 )
    expect_equal( "*", as.character(aa[1767]) )
    expect_equal( "*",  as.character(aa[ length(aa) ] ) )


} )
