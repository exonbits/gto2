# FASTQ Tools

The toolkit has a set of tools dedicated to manipulating FASTQ files. Some of these tools allow the data conversion to/from different formats, i. e., there are tools designed to convert a FASTQ file into a sequence or a FASTA/Multi-FASTA format, or converting DNA in some of those formats to FASTQ.

There are also tools for data manipulation in this format, which are designed to exclude 'N', remove low quality scored reads, following different metrics and randomize DNA sequences. Succeeding the manipulation, it is also possible to perform analyses over these files, simulations and mutations. The current available tools for FASTQ format analysis and manipulation include:

- **gto2_fq_to_fa**: it converts a FASTQ file format to a pseudo FASTA file.
- **gto2_fq_to_mfa**: it converts a FASTQ file format to a pseudo Multi-FASTA file.
- **gto2_fq_exclude_n**: it discards the FASTQ reads with the minimum number of "N" symbols.
- **gto2_fq_extract_quality_scores**: it extracts all the quality-scores from FASTQ reads.
- **gto2_fq_info**: it analyses the basic information of FASTQ file format.
- **gto2_fq_maximum_read_size**: it filters the FASTQ reads with the length higher than the value defined.
- **gto2_fq_minimum_quality_score**: it discards reads with average quality-score below of the defined.
- **gto2_fq_minimum_read_size**: it filters the FASTQ reads with the length smaller than the value defined.
- **gto2_fq_rand_extra_chars**: it substitues in the FASTQ files, the DNA sequence the outside ACGT chars by random ACGT symbols.
- **gto2_fq_from_seq**: it converts a genomic sequence to pseudo FASTQ file format.
- **gto2_fq_mutate**: it creates a synthetic mutation of a FASTQ file given specific rates of mutations, deletions and additions.
- **gto2_fq_split**: it splits Paired End files according to the direction of the strand ('/1' or '/2').
- **gto2_fq_pack**: it packages each FASTQ read in a single line.
- **gto2_fq_unpack**: it unpacks the FASTQ reads packaged using the **gto2_fastq_pack** tool.
- **gto2_fq_quality_score_info**: it analyses the quality-scores of a FASTQ file.
- **gto2_fq_quality_score_min**: it analyses the minimal quality-scores of a FASTQ file.
- **gto2_fq_quality_score_max**: it analyses the maximal quality-scores of a FASTQ file.
- **gto2_fq_cut**: it cuts read sequences in a FASTQ file. 
- **gto2_fq_minimum_local_quality_score_forward**: it filters the reads considering the quality score average of a defined window size of bases.
- **gto2_fq_minimum_local_quality_score_reverse**: it filters the reverse reads, considering the average window size score defined by the bases.
- **gto2_fq_xs**: it is a skilled FASTQ read simulation tool, flexible, portable and tunable in terms of sequence complexity.
- **gto2_fq_complement**: it replaces the ACGT bases with their complements in a FASTQ file format.
- **gto2_fq_reverse**: it reverses the ACGT bases order for each read in a FASTQ file format.
- **gto2_fq_variation_map**: it identifies the variation that occours in the sequences relative to the reads or a set of reads.
- **gto2_fq_variation_filter**: it filters and segments the regions of singularity from the output of **gto2_fq_variation_map**.
- **gto2_fq_variation_visual**: it depites the regions of singularity using the output from **gto2_fq_variation_filter** into an SVG image.
- **gto2_fq_metagenomics**: it measures similarity between any FASTQ file, independently from the size, against any multi-FASTA database.


## gto2_fq_to_fa

The **gto2_fq_to_fa** converts a FASTQ file format to a pseudo FASTA file. However, it does not align the sequence. Also, it extracts the sequence and adds a pseudo header.

For help type:

```sh
./gto2_fq_to_fa -h
```

In the following subsections, we explain the input and output paramters.

### Input parameters {-}

The **gto2_fq_to_fa** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_to_fa [options] [[--] args]
   or: ./gto2_fq_to_fa [options]

It converts a FASTQ file format to a pseudo FASTA file.
It does NOT align the sequence.
It extracts the sequence and adds a pseudo header.

    -h, --help            show this help message and exit

Basic options
    < input.fastq         Input FASTQ file format (stdin)
    > output.fasta        Output FASTA file format (stdout)

Example: ./gto2_fq_to_fa < input.fastq > output.fasta
```

An example of such an input file is:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGG
+SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
+SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_to_fa** program a FASTA file.
Using the input above, an output example for this is the following:

```sh
> Computed with Fastq2Fasta
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGG
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
```


## gto2_fq_to_mfa
to do

## gto2_fq_exclude_n
to do

## gto2_fq_extract_quality_scores
to do

## gto2_fq_info
to do

## gto2_fq_maximum_read_size
to do

## gto2_fq_minimum_quality_score
to do

## gto2_fq_minimum_read_size
to do

## gto2_fq_rand_extra_chars
to do

## gto2_fq_from_seq
to do

## gto2_fq_mutate
to do

## gto2_fq_split
to do

## gto2_fq_pack
to do

## gto2_fq_unpack
to do

## gto2_fq_quality_score_info
to do

## gto2_fq_quality_score_min
to do

## gto2_fq_quality_score_max
to do

## gto2_fq_cut
to do

## gto2_fq_minimum_local_quality_score_forward
to do

## gto2_fq_minimum_local_quality_score_reverse
to do

## gto2_fq_xs
to do

## gto2_fq_complement
to do

## gto2_fq_reverse
to do

## gto2_fq_variation_map
to do

## gto2_fq_variation_filter
to do

## gto2_fq_variation_visual
to do

## gto2_fq_metagenomics
to do
