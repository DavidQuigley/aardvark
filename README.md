# AARDVARK

**An Automated Reversion Detector for Variants Affecting Resistance Kinetics**

<img src="https://quigleylab.ucsf.edu/sites/g/files/tkssra5646/f/wysiwyg/AARDVARK.jpg" style="height: 169px; width:362px;"><br />

AARDVARK is an R package that identifies reversion mutations in DNA sequence data. For 
motivation, you could read our 2017 paper [Quigley et al. Cancer Discovery 2017](https://pubmed.ncbi.nlm.nih.gov/28450426/) where we demonstrated that a common form of PARP inhibitor resistance, called reversion mutations, can be detected in advanced prostate cancer by liquid biopsy.

AARDVARK was developed in the [Quigley lab](https://quigleylab.ucsf.edu) at UCSF.

## Installation and dependencies

AARDVARK is an R package and can be installed in the usual way from source:

```
install.packages("aardvark_0.2.5.tar.gz", repos=NULL, type="source")
```

The AARDVARK package is dependent on the following R packages, all of which are available through 
CRAN or Bioconductor:

GenomicRanges, VariantAnnotation, GenomicAlignments, Rsamtools, Biostrings, stringr, 
IRanges, BSgenome.Hsapiens.UCSC.hg38, BSgenome.Hsapiens.UCSC.hg19, biomaRt, ensembldb, 
testthat, EnsDb.Hsapiens.v86, Gviz

## Quick Start: within R

**Summary**

1) Define a pathogenic mutation (e.g. g.32339657 CTT > C)
2) Define region in the genome to use for realignments
3) Define a read in that region to test
4) Locally realign the read and predict whether that read creates a reversion

In typical use cases, you define define one or more pathogenic mutations in a VCF 
file and all sequencing reads in a window around that region are automatically read from a BAM file. 

AARDVARK ships with command-line scripts that allow users to use AARDVARK without writing any R code.

**Example**

```
library( aardvark )
library( BSgenome.Hsapiens.UCSC.hg38 )

# ------------------------------------------------------------------------------
# define the location and nature of the pathogenic mutation and
# tell AARDVARK where to locally realign
# ------------------------------------------------------------------------------

pathogenic_mut = aardvark::Mutation( chrom="chr13", 
                         pos=32339657, 
                         seq_ref = "CTT", 
                         seq_alt = "C", 
                         transcript=aardvark::transcript_BRCA2)

AW = aardvark::AlignmentWindow( Hsapiens_version = BSgenome.Hsapiens.UCSC.hg38, 
                                chrom="chr13",
                                window_start = pathogenic_mut$pos - 3000, 
                                window_end = pathogenic_mut$pos + 3000)

# ------------------------------------------------------------------------------
# create a read object (can be automated with the read_from_BamData() function )
# Realign the read and show that while the original read has a 72bp soft clip,
# the corrected read has a large deletion flanked by two perfect matches.
# ------------------------------------------------------------------------------

read_nt = "CAGCCTTAGCTTTTTACACAAGTTGTAGTAGAAAAACTTCTGTGAGTCAGACTTCATTACTTGAAGCAAAAAAAAGTTCCTTACACAAAGTTAAGGGAGTGTTAGAGGAATTTGATTTAATCAGAACTGAGCATAGTCTTCACTATTCACC"
read_original = aardvark::Read( qname="A00887:299:HWFYGDSXY:2:2674:25211:28682",
                     cigar = "72S79M",
                     chrom = "chr13",
                     pos = 32340564,
                     seq = DNAString( read_nt ),
                     qual = rep(37, 151) )
                     
print( read_original$cigar_ranges )
read_realigned = aardvark::realign_read( read = read_original, 
                                         align_window = AW, 
                                         pathogenic_mutation = pathogenic_mut )

print( read_realigned$cigar_ranges )

# ------------------------------------------------------------------------------
# Assess the predicted consequences of the realigned read and show the read
# is predicted to produce a reversion that spans the pathogenic variant.
# ------------------------------------------------------------------------------

read_realigned = aardvark::assess_reversion(read = read_realigned,
                           transcript = aardvark::transcript_BRCA2,
                           pathogenic = pathogenic_mut,
                           align_window = AW,  
                           gr_pathogenic =  aardvark::genomicRangesFromMutation(pathogenic_mut) )
                                                    
                                                    
read_summary = aardvark::summarize_candidates( list( read_realigned ), 
                                               transcript=aardvark::transcript_BRCA2,
                                               pathogenic_mutation=pathogenic_mut)
print( read_summary$summary )
```

The read summary will be a data frame containing:

alias|N|reversion|evidence|pos|chrom|transcript_id|pathogenic_mutation
--|--|--|--|--|--|--|--
D:1137:32339427:32340563|1|D:1137:32339427:32340563|reversion_read_deletion_spans_pathogenic_variant|32339427|chr13|ENST00000380152|D:2:32339658:32339659

## Quick Start: command line

This script loads in the reads from an aligned and indexed bam file bam_input.bam and realigns any reads that fall within a 3000 base window centered around any variant described in the VCF file. The BAM must be aligned and indexed, and the reference must match the genome draft you specify. Using this approach, any number of alleles could be tested.

**IMPORTANT**: Running this program does *not* create a new BAM file or in any way modify the BAM you are using; the realigned reads are used to generate a report that will be written to the output folder.

```
Rscript realign_BAM_region_from_VCF.R \
  --sample_id test \
  --fn_bam bam_input.bam \
  --window_size 3000 \
  --genome_draft 38 \
  --fn_vcf /path/to/vcf/with/candidate/pathogenic/alleles.vcf \
  --dir_out /path/to/output

```

### Where to find the command line scripts

The realign_BAM_region_from_VCF.R file can be found in the */exec* folder under wherever R installs aardvark. To find it on your installation, use the built-in *.libPaths()* function in R. 

On my current build *.libPaths()* returns  
*/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library*  
so on my installation the *realign_BAM_region_from_VCF.R* script is at  
*/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/aardvark/exec/realign_BAM_region_from_VCF.R*

### Example data

To try out AARDVARK on a real BAM file, you can [download a small segment of an indexed BAM file and paired VCF](https://quigleylab.ucsf.edu/publications/2023/aardvark_example.tar.gz) that is suitable for AARDVARK. Note that this BAM file only includes a tiny piece of sequence surrounding the pathogenic mutation on *BRCA2* and cannot be used to reconstruct anything else about the genome from this person.
