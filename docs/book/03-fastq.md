# FASTQ Tools

The toolkit has a set of tools dedicated to manipulating FASTQ files. Some of these tools allow the data conversion to/from different formats, i. e., there are tools designed to convert a FASTQ file into a sequence or a FASTA/Multi-FASTA format, or converting DNA in some of those formats to FASTQ.

There are also tools for data manipulation in this format, which are designed to exclude 'N', remove low quality scored reads, following different metrics and randomize DNA sequences. Succeeding the manipulation, it is also possible to perform analyses over these files, simulations and mutations. The current available tools for FASTQ format analysis and manipulation include:

- **gto2_fq_to_fa**: to convert a FASTQ file format to a pseudo FASTA file.
- **gto2_fq_to_mfa**: to convert a FASTQ file format to a pseudo Multi-FASTA file.
- **gto2_fq_exclude_n**: to discard the FASTQ reads with the minimum number of "N" symbols.
- **gto2_fq_extract_quality_scores**: to extract all the quality-scores from FASTQ reads.
- **gto2_fq_info**: to analyse the basic information of FASTQ file format.
- **gto2_fq_maximum_read_size**: to filter the FASTQ reads with the length higher than the value defined.
- **gto2_fq_minimum_quality_score**: to discard reads with average quality-score below of the defined.
- **gto2_fq_minimum_read_size**: to filter the FASTQ reads with the length smaller than the value defined.
- **gto2_fq_rand_extra_chars**: to substitue in the FASTQ files, the DNA sequence the outside ACGT chars by random ACGT symbols.
- **gto2_fq_from_seq**: to convert a genomic sequence to pseudo FASTQ file format.
- **gto2_fq_mutate**: to create a synthetic mutation of a FASTQ file given specific rates of mutations, deletions and additions.
- **gto2_fq_split**: to split Paired End files according to the direction of the strand ('/1' or '/2').
- **gto2_fq_pack**: to package each FASTQ read in a single line.
- **gto2_fq_unpack**: to unpack the FASTQ reads packaged using the **gto2_fq_pack** tool.
- **gto2_fq_quality_score_info**: to analyse the quality-scores of a FASTQ file.
- **gto2_fq_quality_score_min**: to analyse the minimal quality-scores of a FASTQ file.
- **gto2_fq_quality_score_max**: to analyse the maximal quality-scores of a FASTQ file.
- **gto2_fq_cut**: to cut read sequences in a FASTQ file. 
- **gto2_fq_minimum_local_quality_score_forward**: to filter the reads considering the quality score average of a defined window size of bases.
- **gto2_fq_minimum_local_quality_score_reverse**: to filter the reverse reads, considering the average window size score defined by the bases.
- **gto2_fq_xs**: a skilled FASTQ read simulation tool, flexible, portable and tunable in terms of sequence complexity.
- **gto2_fq_complement**: to replace the ACGT bases with their complements in a FASTQ file format.
- **gto2_fq_reverse**: to reverse the ACGT bases order for each read in a FASTQ file format.
- **gto2_fq_variation_map**: to identify the variation that occours in the sequences relative to the reads or a set of reads.
- **gto2_fq_variation_filter**: to filter and segments the regions of singularity from the output of **gto2_fq_variation_map**.
- **gto2_fq_variation_visual**: to depict the regions of singularity using the output from **gto2_fq_variation_filter** into an SVG image.
- **gto2_fq_metagenomics**: to measure the similarity between any FASTQ file, independently from the size, against any multi-FASTA database.


## Program gto2_fq_to_fa

The **gto2_fq_to_fa** converts a FASTQ file format to a pseudo FASTA file. However, this tool does not align the sequence, instead, it extracts the sequence and adds a pseudo-header.

For help type:

```sh
./gto2_fq_to_fa -h
```

In the following subsections, we explain the input and output parameters.

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

The output of the **gto2_fq_to_fa** program is a FASTA file.
Using the input above, an output example of this is the following:

```sh
> Computed with Fastq2Fasta
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGG
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
```

## Program gto2_fq_to_mfa

The **gto2_fq_to_mfa** converts a FASTQ file format to a pseudo Multi-FASTA file. However, this tool does not align the sequence, instead, it extracts the sequence and adds a pseudo header.

For help type:

```sh
./gto2_fq_to_mfa -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_to_mfa** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_to_mfa [options] [[--] args]
   or: ./gto2_fq_to_mfa [options]

It converts a FASTQ file format to a pseudo Multi-FASTA file.
It does NOT align the sequence.
It extracts the sequence and adds each header in a Multi-FASTA format.

    -h, --help            show this help message and exit

Basic options
    < input.fastq         Input FASTQ file format (stdin)
    > output.mfasta       Output Multi-FASTA file format (stdout)

Example: ./gto2_fq_to_mfa < input.fastq > output.mfasta
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

The output of the **gto2_fq_to_mfa** program is a Multi-FASTA file.
Using the input above, an output example of this is the following:

```sh
>SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGG
>SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
```

## Program gto2_fq_exclude_n

The **gto2_fq_exclude_n** discards the FASTQ reads with the minimum number of ''N'' symbols, and it will erase the second header (after +), if presented.

For help type:

```sh
./gto2_fq_exclude_n -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_exclude_n** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_exclude_n [options] [[--] args]
   or: ./gto2_fq_exclude_n [options]

It discards the FASTQ reads with the minimum number of "N" 
symbols. 
If present, it will erase the second header (after +).

    -h, --help            show this help message and exit

Basic options
    -m, --max=<int>       The maximum of of "N" symbols in 
    					  the read
    < input.fastq         Input FASTQ file format (stdin)
    > output.fastq        Output FASTQ file format (stdout)

Example: ./gto2_fq_exclude_n -m <max> < input.fastq > 
output.fastq

Console output example :
<FASTQ non-filtered reads>
Total reads    : value
Filtered reads : value
```

An example of such an input file is:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GNNTGATGGCCGCTGCCGATGGCGNANAATCCCACCAANATACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_exclude_n** program is a set of all the filtered FASTQ reads, followed by the execution report.
The execution report only appears in the console.

Using the input above with the max value as 5, an output example for this is the following:

```sh
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
Total reads    : 2
Filtered reads : 1
```

## Program gto2_fq_extract_quality_scores

The **gto2_fq_extract_quality_scores** extracts all the quality-scores from FASTQ reads.

For help type:

```sh
./gto2_fq_extract_quality_scores -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_extract_quality_scores** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_extract_quality_scores [options] [[--]args]
   or: ./gto2_fq_extract_quality_scores [options]

It extracts all the quality-scores from FASTQ reads.

    -h, --help            show this help message and exit

Basic options
    < input.fastq         Input FASTQ file format (stdin)
    > output.fastq        Output FASTQ file format (stdout)

Example: ./gto2_fq_extract_quality_scores < input.fastq >
output.fastq

Console output example:
<FASTQ quality scores>
Total reads          : value
Total Quality-Scores : value
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

The output of the **gto2_fq_extract_quality_scores** program is a set of all the quality scores from the FASTQ reads, followed by the execution report.
The execution report only appears in the console. Using the input above, an output example of this is the following:

```sh
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
Total reads          : 2
Total Quality-Scores : 144
```

## Program gto2_fq_info

The **gto2_fq_info** analyses the basic information of FASTQ file format.

For help type:

```sh
./gto2_fq_info -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_info** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_info [options] [[--] args]
   or: ./gto2_fq_info [options]

It analyses the basic information of FASTQ file format.

    -h, --help            show this help message and exit

Basic options
    < input.fastq         Input FASTQ file format (stdin)
    > output              Output read information (stdout)

Example: ./gto2_fq_info < input.fastq > output

Output example:
Total reads     : value
Max read length : value
Min read length : value
Min QS value    : value
Max QS value    : value
QS range        : value
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

The output of the **gto2_fq_info** program is a set of information related to the file read. Using the input above, an output example of this is the following:

```sh
Total reads     : 2
Max read length : 72
Min read length : 72
Min QS value    : 41
Max QS value    : 73
QS range        : 33
```

## Program gto2_fq_maximum_read_size

The **gto2_fq_maximum_read_size** filters the FASTQ reads with the length higher than the value defined.

For help type:

```sh
./gto2_fq_maximum_read_size -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_maximum_read_size** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_maximum_read_size [options] [[--] args]
   or: ./gto2_fq_maximum_read_size [options]

It filters the FASTQ reads with the length higher than the 
value defined. 
If present, it will erase the second header (after +).

    -h, --help            show this help message and exit

Basic options
    -s, --size=<int>      The maximum read length
    < input.fastq         Input FASTQ file format (stdin)
    > output.fastq        Output FASTQ file format (stdout)

Example: ./gto2_fq_maximum_read_size -s <size> < input.fastq 
> output.fastq

Console output example :
<FASTQ non-filtered reads>
Total reads    : value
Filtered reads : value
```

An example of such an input file is:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=59
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_maximum_read_size** program is a set of all the filtered FASTQ reads, followed by the execution report.
The execution report only appears in the console.

Using the input above with the size values as 59, an output example for this is the following:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=59
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDII
Total reads    : 2
Filtered reads : 1
```

## Program gto2_fq_minimum_quality_score

The **gto2_fq_minimum_quality_score** discards reads with average quality-score below of the defined.

For help type:

```sh
./gto2_fq_minimum_quality_score -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_minimum_quality_score** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_minimum_quality_score [options] [[--] args]
   or: ./gto2_fq_minimum_quality_score [options]

It discards reads with average quality-score below value.

    -h, --help            show this help message and exit

Basic options
    -m, --min=<int>       The minimum average quality-score 
    					  (Value 25 or 30 is commonly used)
    < input.fastq         Input FASTQ file format (stdin)
    > output.fastq        Output FASTQ file format (stdout)

Example: ./gto2_fq_minimum_quality_score -m <min> < 
input.fastq > output.fastq

Console output example:
<FASTQ non-filtered reads>
Total reads    : value
Filtered reads : value
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
54599<>77977==6=?I6IBI::33344235521677999>>><<<@@A@BBCDGGBFF
```

### Output {-}

The output of the **gto2_fq_minimum_quality_score** program is a set of all the filtered FASTQ reads, followed by the execution report. Using the input above with the minimum averge value as 30, an output example of this is the following:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
Total reads    : 2
Filtered reads : 1
```

## Program gto2_fq_minimum_read_size

The **gto2_fq_minimum_read_size** filters the FASTQ reads with the length smaller than the value defined.

For help type:

```sh
./gto2_fq_minimum_read_size -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_minimum_read_size** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_minimum_read_size [options] [[--] args]
   or: ./gto2_fq_minimum_read_size [options]

It filters the FASTQ reads with the length smaller than the 
value defined. 
If present, it will erase the second header (after +).

    -h, --help            show this help message and exit

Basic options
    -s, --size=<int>      The minimum read length
    < input.fastq         Input FASTQ file format (stdin)
    > output.fastq        Output FASTQ file format (stdout)

Example: ./gto2_fq_minimum_read_size -s <size> < input.fastq
> output.fastq

Console output example:
<FASTQ non-filtered reads>
Total reads    : value
Filtered reads : value
```

An example of such an input file is:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=50
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_minimum_read_size** program is a set of all the filtered FASTQ reads, followed by the execution report. The execution report only appears in the console. Using the input above with the size values as 55, an output example of this is the following:

```sh
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
Total reads    : 2
Filtered reads : 1
```

## Program gto2_fq_rand_extra_chars

The **gto2_fq_rand_extra_chars** substitutes the outside ACGT chars by random ACGT symbols in the DNA sequence of FASTQ files.

For help type:

```sh
./gto2_fq_rand_extra_chars -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_rand_extra_chars** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_rand_extra_chars [options] [[--] args]
   or: ./gto2_fq_rand_extra_chars [options]

It substitues in the FASTQ files, the DNA sequence the 
outside ACGT chars by random ACGT symbols.

    -h, --help            show this help message and exit

Basic options
    < input.fastq         Input FASTQ file format (stdin)
    > output.fastq        Output FASTQ file format (stdout)

Example: ./gto2_fq_rand_extra_chars < input.fastq > 
output.fastq
```

An example of such an input file is:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GNNTGATGGCCGCTGCCGATGGCGNANAATCCCACCAANATACCCTTAACAACTTAAGGG
+SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
+SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_rand_extra_chars** program is a FASTQ file. Using the input above, an output example of this is the following:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GTGTGATGGCCGCTGCCGATGGCGCATAATCCCACCAACATACCCTTAACAACTTAAGGG
+SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGATTGTATTTTAAGAATCTGCAGCAGAAGTCGATGATAATACGCG
+SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

## Program gto2_fq_from_seq

The **gto2_fq_from_seq** converts a genomic sequence to pseudo FASTQ file format.

For help type:

```sh
./gto2_fq_from_seq -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_from_seq** program needs two streams for the computation, namely the input and output standard. The input stream is a sequence group file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_from_seq [options] [[--] args]
   or: ./gto2_fq_from_seq [options]

It converts a genomic sequence to pseudo FASTQ file format.

    -h, --help            show this help message and exit

Basic options
    < input.seq           Input sequence file (stdin)
    > output.fastq        Output FASTQ file format (stdout)

Optional options
    -n, --name=<str>      The read's header
    -l, --lineSize=<int>  The maximum of chars for line

Example: ./gto2_fq_from_seq -l <lineSize> -n <name> < 
input.seq > output.fastq
```

An example of such an input file is:

```sh
ACAAGACGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCCTGGAGGGTCCACCGCT
GCCCTGCTGCCATTGTCCCCGGCCCCACCTAAGGAAAAGCAGCCTCCTGACTTTCCTCGC
TTGGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAAGTGGTTTGAGTGGACCTCCG
GGCCCCTCATAGGAGAGGAAGCTCGGGAGGTGGCCAGGCGGCAGGAAGCAGGCCAGTGCC
GCGAATCCGCGCGCCGGGACAGAATCTCCTGCAAAGCCCTGCAGGAACTTCTTCTGGAAG
```

### Output {-}

The output of the **gto2_fq_from_seq** program is a pseudo FASTQ file. An example, using the size line as 60 and the read's header as ''SeqToFastq'', for the input, is:

```sh
@SeqToFastq1
ACAAGACGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCCTGGAGGGTCCACCGCT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@SeqToFastq2
GCCCTGCTGCCATTGTCCCCGGCCCCACCTAAGGAAAAGCAGCCTCCTGACTTTCCTCGC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@SeqToFastq3
TTGGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAAGTGGTTTGAGTGGACCTCCG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@SeqToFastq4
GGCCCCTCATAGGAGAGGAAGCTCGGGAGGTGGCCAGGCGGCAGGAAGCAGGCCAGTGCC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@SeqToFastq5
GCGAATCCGCGCGCCGGGACAGAATCTCCTGCAAAGCCCTGCAGGAACTTCTTCTGGAAG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```

## Program gto2_fq_mutate

The **gto2_fq_mutate** creates a synthetic mutation of a FASTQ file given specific rates of mutations, deletions and additions. All these paramenters are defined by the user, and their are optional.

For help type:

```sh
./gto2_fq_mutate -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_mutate** program needs two streams for the computation, namely the input and output standard. However, optional settings can be supplied too, such as the starting point to the random generator, and the edition, deletion and insertion rates. Also, the user can choose to use the ACGTN alphabet in the synthetic mutation. The input stream is a FASTQ File.

The attribution is given according to:

```sh
Usage: ./gto2_fq_mutate [options] [[--] args]
   or: ./gto2_fq_mutate [options]

Creates a synthetic mutation of a FASTQ file given specific 
rates of mutations, deletions and additions.

    -h, --help      show this help message and exit

Basic options
    < input.fasta   Input FASTQ file format (stdin)
    > output.fasta  Output FASTQ file format (stdout)

Optional
    -s              Starting point to the random generator
    -m 				Defines the mutation rate (default 0.0)
    -d 				Defines the deletion rate (default 0.0)
    -i 				Defines the insertion rate (default 0.0)
    -a 				When active, the application uses the 
    				ACGTN alphabet

Example: ./gto2_fq_mutate -s <seed> -m <mutation rate> 
-d <deletion rate> -i <insertion rate> -a 
< input.fastq > output.fastq

```

An example of such an input file is:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_mutate** program is a FASTQ file whith the synthetic mutation of input file. Using the input above with the seed value as 1 and the mutation rate as 0.5, an output example of this is the following:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GGACTTTGAGGTGTGGCGATAGACTGAAAACACTTCAGGGTAAAATCACTCGCAAAAGTG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGAGCCTTTACCGTAGGGGTGTAAGATTTTATACAAAAAGTCCAGGTCAAGAGGAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

## Program gto2_fq_split

The **gto2_fq_split** splits Paired End files according to the direction of the strand ('/1' or '/2'). It writes by default singleton reads as forward stands. 

For help type:

```sh
./gto2_fq_split -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_split** program needs a stream for the computation, namely the input standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_split [options] [[--] args]
   or: ./gto2_fq_split [options]

It writes by default singleton reads as forward stands.

    -h, --help            show this help message and exit

Basic options
    -f, --forward=<str>   Output forward file
    -r, --reverse=<str>   Output reverse file
    < input.fastq         Input FASTQ file format (stdin)
    > output         	  Output read information (stdout)

Example: ./gto2_fq_split -f <output_forward.fastq> 
-r <output_reverse.fastq> < input.fastq > output

Output example :
Total reads      : value
Singleton reads  : value
Forward reads    : value
Reverse reads    : value
```

An example of such an input file is:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60 1
GNNTGATGGCCGCTGCCGATGGCGNANAATCCCACCAANATACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60 2
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_split** program is a set of information related to the file read.
Using the input above, an output example of this is the following:

```sh
Total reads     : 2
Singleton reads : 0
Forward reads   : 65536
Reverse reads   : 1
```

Also, this program generates two FASTQ files, with the reverse and forward reads.

An example of the forward reads, for the input, is: 

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60 1
GNNTGATGGCCGCTGCCGATGGCGNANAATCCCACCAANATACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
```

## Program gto2_fq_pack

The **gto2_fq_pack** packages each FASTQ read in a single line. It can show the read score first or the dna sequence, depending on the execution mode. 

For help type:

```sh
./gto2_fq_pack -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_pack** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_pack [options] [[--] args]
   or: ./gto2_fq_pack [options]

It packages each FASTQ read in a single line.

    -h, --help            show this help message and exit

Basic options
    < input.fastq         Input FASTQ file format (stdin)
    > output.fastqpack    Output packaged FASTQ file format 
    					  (stdout)

Optional
    -s, --scores          When active, the application show 
    					  the scores first

Example: ./gto2_fq_pack -s < input.fastq > output.fastqpack
```

An example of such an input file is:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GNNTGATGGCCGCTGCCGATGGCGNANAATCCCACCAANATACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_pack** program is a packaged FASTQ file.
Using the input above, an output example of this is the following:

```sh
GNNTGATGGCCGCTGCCGATGGCGNANAATCCCACCAANATACCCTTAACAACTTAAGGG
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60+	0
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60+	1
```

Another example for the same input, but using the scores first (flag ''s''), is:

```sh
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
GNNTGATGGCCGCTGCCGATGGCGNANAATCCCACCAANATACCCTTAACAACTTAAGGG
SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60+	0
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60+	1
```

## Program gto2_fq_unpack

The **gto2_fq_unpack** unpacks the FASTQ reads packaged using the **gto2_fq_pack** tool.

For help type:

```sh
./gto2_fq_unpack -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_unpack** program needs two streams for the computation, namely the input and output standard. The input stream is a packaged FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_unpack [options] [[--] args]
   or: ./gto2_fq_unpack [options]

It unpacks the FASTQ reads packaged using the gto2_fq_pack 
tool.

    -h, --help            show this help message and exit

Basic options
    < input.fastq         Input FASTQ file format (stdin)
    > output.fastq        Output FASTQ file format (stdout)

Optional
    -s, --scores          When active, the application show 
    					  the scores first
    
Example: ./gto2_fq_unpack -s < input.fastqpack > out.fastq
```

An example of such an input file is:

```sh
GNNTGATGGCCGCTGCCGATGGCGNANAATCCCACCAANATACCCTTAACAACTTAAGGG
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60+	0
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60+ 1
```

### Output {-}

The output of the **gto2_fq_unpack** program is a FASTQ file.
Using the input above, an output example of this is the following:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GNNTGATGGCCGCTGCCGATGGCGNANAATCCCACCAANATACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

## Program gto2_fq_quality_score_info

The **gto2_fq_quality_score_info** analyses the quality-scores of a FASTQ file.

For help type:

```sh
./gto2_fq_quality_score_info -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_quality_score_info** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_quality_score_info [options] [[--] args]
   or: ./gto2_fq_quality_score_info [options]

It analyses the quality-scores of a FASTQ file.

    -h, --help            show this help message and exit

Basic options
    < input.fastq         Input FASTQ file format (stdin)
    > output              Output read information (stdout)
    
Optional
    -m, --max=<int>       The lenght of the maximum window

Example: ./gto2_fq_quality_score_info -m <max> < input.fastq
> output

Output example :
Total reads     : value
Max read length : value
Min read length : value
Min QS value    : value
Max QS value    : value
QS range        : value
```

An example of such an input file is:

```sh
@111 071112_SLXA-EAS1_s_7:5:1:817:345 length=60 1
GNNTGATGGCCGCTGCCGATGGCGNANAATCCCACCAANATACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@222 071112_SLXA-EAS1_s_7:5:1:801:338 length=60 2
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_quality_score_info** program is a set of information related to the file read. Using the input above with the max window value as 30, an output example of this is the following:

```sh
Total reads     : 2
Max read length : 60
Min read length : 60
Min QS value    : 54
Max QS value    : 73
QS range        : 20
 1  ...  24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 
--+ ... +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
73  ...  73 73 73 73 73 73 73 65 73 62 65 69 70 73 73 73 73 
                               *     *  *  *  *             
                               *     *  *  *  *             
                               *     *  *  *  *             
                               *     *  *  *                
                               *     *  *                   
                               *     *  *                   
                               *     *  *                   
                               *     *  *                   
                                     *                      
                                     *                      
                                     *   
```

## Program gto2_fq_quality_score_min

The **gto2_fq_quality_score_min** analyses the minimal quality-scores of a FASTQ file.

For help type:

```sh
./gto2_fq_quality_score_min -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_quality_score_min** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_quality_score_min [options] [[--] args]
   or: ./gto2_fq_quality_score_min [options]

It analyses the minimal quality-scores of a FASTQ file.

    -h, --help        show this help message and exit

Basic options
    < input.fastq     Input FASTQ file format (stdin)
    > output          Output read information (stdout)

Optional
    -m, --max=<int>   The maximum window length (default 40)

Example: ./gto2_fq_quality_score_min -m <max> < input.fastq 
> output
```

An example of such an input file is:

```sh
@111 071112_SLXA-EAS1_s_7:5:1:817:345 length=60 1
GNNTGATGGCCGCTGCCGATGGCGNANAATCCCACCAANATACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@222 071112_SLXA-EAS1_s_7:5:1:801:338 length=60 2
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_quality_score_min** program is a set of information related to the file read, considering the minimum quality scores. Using the input above with the max window value as 20, an output example of this is the following:

```sh
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
73 73 73 73 73 73 73 73 73 73 73 73 73 73 73 73 73 73 73 73
```

## Program gto2_fq_quality_score_max

The **gto2_fq_quality_score_max** analyses the maximal quality-scores of a FASTQ file.

For help type:

```sh
./gto2_fq_quality_score_max -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_quality_score_max** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_quality_score_max [options] [[--] args]
   or: ./gto2_fq_quality_score_max [options]

It analyses the maximal quality-scores of a FASTQ file.

    -h, --help            show this help message and exit

Basic options
    < input.fastq         Input FASTQ file format (stdin)
    > output              Output read information (stdout)

Optional
    -m, --max=<int>       The maximum window length (default 40)

Example: ./gto2_fq_quality_score_max -m <max> < input.fastq > output
```

An example of such an input file is:

```sh
@111 071112_SLXA-EAS1_s_7:5:1:817:345 length=60 1
GNNTGATGGCCGCTGCCGATGGCGNANAATCCCACCAANATACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@222 071112_SLXA-EAS1_s_7:5:1:801:338 length=60 2
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_quality_score_max** program is a set of information related to the file read, considering the maximal quality scores. Using the input above with the max window value as 20, an output example of this is the following:

```sh
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
73 73 73 73 73 73 73 73 73 73 73 73 73 73 73 73 73 73 73 73
```

## Program gto2_fq_cut
to do

## Program gto2_fq_minimum_local_quality_score_forward
to do

## Program gto2_fq_minimum_local_quality_score_reverse
to do

## Program gto2_fq_xs
to do

## Program gto2_fq_complement
to do

## Program gto2_fq_reverse
to do

## Program gto2_fq_variation_map
to do

## Program gto2_fq_variation_filter
to do

## Program gto2_fq_variation_visual
to do

## Program gto2_fq_metagenomics
to do
