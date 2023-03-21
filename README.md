# AARDVARK

An Automated Reversion Detector for Variants Affecting Resistance Kinetics

<img src="https://quigleylab.ucsf.edu/sites/g/files/tkssra5646/f/wysiwyg/AARDVARK.jpg" style="height: 169px; width:362px;"><br />

AARDVARK is an R package that identifies reversion mutations in DNA sequence data.

AARDVARK was developed in the [Quigley lab](https://quigleylab.ucsf.edu) at UCSF.

## Quick Start: within R

**Summary**

1) Define a pathogenic mutation
2) Tell AARDVARK a region in the genome to use for realignments
3) Define a read object to realign in that region
4) Realign the read and assess whether it creates a reversion

In real use cases, all of the reads in a region are automatically read from a BAM file and the pathogenic mutation can be read from a VCF file of candidate locations. In this use case, most of the reads will not harbor a reversion. AARDVARK ships with wrappers that allow users to perform searches from the command line without writing any R code.

**Example**

```
library( aardvark )
library( BSgenome.Hsapiens.UCSC.hg38 )

# define the location and nature of the pathogenic mutation

pathogenic_mut = aardvark::Mutation( chrom="chr13", 
                         pos=32339657, 
                         seq_ref = "CTT", 
                         seq_alt = "C", 
                         transcript=aardvark::transcript_BRCA2)

# tell AARDVARK where to locally realign
AW = aardvark::AlignmentWindow( Hsapiens_version = BSgenome.Hsapiens.UCSC.hg38, 
                                chrom="chr13",
                                window_start = pathogenic_mut$pos - 3000, 
                                window_end = pathogenic_mut$pos + 3000)

# create a read object (normally automated with the read_from_BamData() function )
read_nt = "CAGCCTTAGCTTTTTACACAAGTTGTAGTAGAAAAACTTCTGTGAGTCAGACTTCATTACTTGAAGCAAAAAAAAGTTCCTTACACAAAGTTAAGGGAGTGTTAGAGGAATTTGATTTAATCAGAACTGAGCATAGTCTTCACTATTCACC"
read_qualities = rep(37, 151)
read_original = aardvark::Read( qname="A00887:299:HWFYGDSXY:2:2674:25211:28682",
                     cigar = "72S79M",
                     chrom = "chr13",
                     pos = 32340564,
                     seq = DNAString( read_nt ),
                     qual = read_qualities )

# Realign the read and show that while the original read has a 72bp soft clip,
# the corrected read has a large deletion flanked by two perfect matches.
# Translation of the original read alignment predicts a first stop codon at p.1772
# Translation of the realigned read predicts stop codon at p.3040, producing a reversion.

print( read_original$cigar_ranges )
read_realigned = aardvark::realign_read( read = read_original, 
                                         align_window = AW, 
                                         gr_pathogenic = aardvark::genomicRangesFromMutation(GM) )

print( read_realigned$cigar_ranges )
aa = translate_cigar( transcript_BRCA2, read_original, pathogenic_mut, min_nt_qual=20 )
print( which( strsplit(toString( aa ), "")[[1]] == "*" ) )

aa = translate_cigar( transcript_BRCA2, read_realigned, pathogenic_mut, min_nt_qual=20 )
print( which( strsplit(toString( aa ), "")[[1]] == "*" ) )

# Assess the predicted consequences of the realigned read and show the read
# is predicted to produce a reversion that spans the pathogenic variant.
read_realigned = aardvark::assess_reversion(read = read_realigned,
                           transcript = aardvark::transcript_BRCA2,
                           pathogenic = pathogenic_mut,
                           align_window = AW,
                           gr_pathogenic =  aardvark::genomicRangesFromMutation(GM) )
                                                    
                                                    
read_summary = aardvark::summarize_candidates( list( read_realigned ), 
                                               transcript=aardvark::transcript_BRCA2 )
print( read_summary$summary )
```

## Quick Start: command line

This script reads in the reads from an aligned and indexed bam file bam_input.bam and realign any reads that fall within a 3000 base window centered around any variant described in the VCF file. The BAM must be aligned and indexed, and the reference must match the genome draft you specify. Using this approach, any number of alleles could be tested.

**IMPORTANT**: This does *not* create a new BAM file or in any way modify the BAM you are using; the realigned reads are used to generate a report that will be written to the output folder.

```
Rscript realign_BAM_region_from_VCF.R \
  --sample_id test \
  --fn_bam bam_input.bam \
  --window_size 3000 \
  --genome_draft 38 \
  --fn_vcf /path/to/vcf/with/candidate/pathogenic/alleles.vcf \
  --dir_out /path/to/output

```

