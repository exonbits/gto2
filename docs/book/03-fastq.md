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
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
+
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

It converts a FASTQ file format to a pseudo Multi-FASTA 
file. It does NOT align the sequence. It extracts the 
sequence and adds each header in a Multi-FASTA format.

    -h, --help       show this help message and exit

Basic options
    < input.fastq    Input FASTQ file format (stdin)
    > output.mfasta  Output Multi-FASTA file format (stdout)

Example: ./gto2_fq_to_mfa < input.fastq > output.mfasta
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
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
+
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
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
+
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
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
+
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
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
NTTCAGGGATACGACGNTTGTATTTTAAGAATCTGNAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_rand_extra_chars** program is a FASTQ file. Using the input above, an output example of this is the following:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GTGTGATGGCCGCTGCCGATGGCGCATAATCCCACCAACATACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGATTGTATTTTAAGAATCTGCAGCAGAAGTCGATGATAATACGCG
+
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

Example: ./gto2_fq_from_seq -l <lineSize> -n <name> 
< input.seq > output.fastq
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

    -h, --help        show this help message and exit

Basic options
    < input.fastq     Input FASTQ file format (stdin)
    > output          Output read information (stdout)

Optional
    -m, --max=<int>   The maximum window length (default 40)

Example: ./gto2_fq_quality_score_max -m <max> 
< input.fastq > output
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

The **gto2_fq_cut** cuts read sequences in a FASTQ file. It requires that the initial and end positions for the cut.

For help type:

```sh
./gto2_fq_cut -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_cut** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_cut [options] [[--] args]
   or: ./gto2_fq_cut [options]

It cuts read sequences in a FASTQ file.

    -h, --help            show this help message and exit

Basic options
    -i, --initial=<int>   Starting position to the cut
    -e, --end=<int>       Ending position to the cut
    < input.fastq         Input FASTQ file format (stdin)
    > output.fastq        Output FASTQ file format (stdout)

Example: ./gto2_fq_cut -i <initial> -e <end> < input.fastq 
> output.fastq
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

The output of the **gto2_fq_cut** program is a FASTQ file cut. Using the initial value as 10 and the end value as 30, an example of this input, is the following:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
CGCTGCCGATGGCGTCAAATC
+
IIIIIIIIIIIIIIIIIIII9
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
ACGACGTTTGTATTTTAAGAA
+
IIIIIIIIIIIIIIIIIIIII
```

## Program gto2_fq_minimum_local_quality_score_forward

The **gto2_fq_minimum_local_quality_score_forward** filters the reads considering the quality score average of a defined window size of bases.

For help type:

```sh
./gto2_fq_minimum_local_quality_score_forward -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_minimum_local_quality_score_forward** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_minimum_local_quality_score_forward [options] [[--] args]
   or: ./gto2_fq_minimum_local_quality_score_forward [options]

It filters the reads considering the quality score average 
of a defined window size of bases.

    -h, --help        show this help message and exit

Basic options
    -k                The window size of bases (default 5)
    -w                The minimum average of quality score 
    				  (default 25)
    -m                The minimum value of the quality score 
    				  (default 33)
    < input.fastq     Input FASTQ file format (stdin)
    > output.fastq    Output FASTQ file format (stdout)

Example: ./gto2_fq_minimum_local_quality_score_forward 
-k <windowsize> -w <minavg> -m <minqs> 
< input.fastq > output.fastq

Console output example:
Minimum QS       : value
<FASTQ output>
Total reads      : value
Trimmed reads    : value
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
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIII-I)8IIIIIIIIIIIIIIIIII
```

### Output {-}

The output of the **gto2_fq_minimum_local_quality_score_forward** program is a FASTQ file with the reads filtered following a quality score average of a defined window of bases. The execution report only appears in the console. Using the input above with the default values, an output example of this is the following:

```sh
Minimum QS     : 33
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBII
Total reads    : 2
Trimmed reads  : 1
```

## Program gto2_fq_minimum_local_quality_score_reverse

The **gto2_fq_minimum_local_quality_score_reverse** filters the reverse reads, considering the quality score average of a defined window size of bases.

For help type:

```sh
./gto2_fq_minimum_local_quality_score_reverse -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_minimum_local_quality_score_reverse** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_minimum_local_quality_score_reverse [options] [[--] args]
   or: ./gto2_fq_minimum_local_quality_score_reverse [options]

It filters the reverse reads, considering the quality score 
average of a defined 
window size of bases.

    -h, --help        show this help message and exit

Basic options
    -k                The window size of bases (default 5)
    -w                The minimum average of quality score 
    				  (default 25)
    -m                The minimum value of the quality score
    				  (default 33)
    < input.fastq     Input FASTQ file format (stdin)
    > output.fastq    Output FASTQ file format (stdout)

Example: ./gto2_fq_minimum_local_quality_score_reverse 
-k <windowsize> -w <minavg> -m <minqs> 
< input.fastq > output.fastq

Console output example:
Minimum QS       : value
<FASTQ output>
Total reads      : value
Trimmed reads    : value
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
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIII-I)8IIIIIIIIIIIIIIIIII
```

### Output {-}

The output of the **gto2_fq_minimum_local_quality_score_reverse** program is a FASTQ file with the reads filtered following a quality score average of a defined window of bases. The execution report only appears in the console. Using the input above with the default values, an output example of this is the following:

```sh
Minimum QS: 33
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
GTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIII
Total reads    : 2
Trimmed reads  : 1
```

## Program gto2_fq_xs

The **gto2_fq_xs** is a skilled FASTQ read simulation tool, flexible, portable (does not need a reference sequence) and tunable in terms of sequence complexity. XS handles Ion Torrent, Roche-454, Illumina and ABI-SOLiD simulation sequencing types. It has several running modes, depending on the time and memory available, and is aimed at testing computing infrastructures, namely cloud computing of large-scale projects, and testing FASTQ compression algorithms. Moreover, XS offers the possibility of simulating the three main FASTQ components individually (headers, DNA sequences and quality-scores). Quality-scores can be simulated using uniform and Gaussian distributions.

For help type:

```sh
./gto2_fq_xs -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_xs** program needs a FASTQ file to compute.

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

The output of the **gto2_fq_xs** program is a FASTQ file. Using the input above using the common usage with 5 reads (-n 5), an output example of this is the following:

```sh

```

## Program gto2_fq_complement

The **gto2_fq_complement** replaces the ACGT bases with their complements in a FASTQ file format.

For help type:

```sh
./gto2_fq_complement -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_complement** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_complement [options] [[--] args]
   or: ./gto2_fq_complement [options]

It replaces the ACGT bases with their complements in a FASTQ
file format.

    -h, --help            Show this help message and exit

Basic options
    < input.fastq         Input FASTQ file (stdin)
    > output.fastq        Output FASTQ file (stdout)

Example: ./gto2_fq_complement < input.fastq > output.fastq
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

The output of the **gto2_fq_complement** program is the FASTQ file with the ACGT base complements. Using the input above, an output example of this is the following:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=60
CCCACTACCGGCGACGGCTACCGCAGTTTAGGGTGGTTCAATGGGAATTGTTGAATTCCC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=60
CAAGTCCCTATGCTGCAAACATAAAATTCTTAGACTTCGTCTTCAGCTACTATTATGCGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

## Program gto2_fq_reverse

The **gto2_fq_reverse** reverses the ACGT bases order for each read in a FASTQ file format.

For help type:

```sh
./gto2_fq_reverse -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_reverse** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ./gto2_fq_reverse [options] [[--] args]
   or: ./gto2_fq_reverse [options]

It reverses the ACGT bases order for each read in a FASTQ 
file.

    -h, --help            Show this help message and exit

Basic options
    < input.fastq         Input FASTQ file (stdin)
    > output.fastq        Output FASTQ file (stdout)

Example: ./gto2_fq_reverse < input.fastq > output.fastq
```

An example of such an input file is:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGI
```

### Output {-}

The output of the **gto2_fq_reverse** program is the FASTQ file complement with the flag ''(Reversed)'' added in the header.
Using the input above, an output example of this is the following:

```sh
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 (Reversed)
GGGAATTCAACAATTCCCATTGAACCACCCTAAACTGCGGTAGCCGTCGCCGGTAGTGGG
+
IIIDIIIIIIIIIIIIIIIIIIIICI9GI9IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 (Reversed)
GCGCATAATAGTAGCTGAAGACGAAGTCTAAGAATTTTATGTTTGCAGCATAGGGACTTG
+
IGIIIIIIIIIIIIIIIIIIIIIIIBI6IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

## Program gto2_fq_variation_map

The **gto2_fq_variation_map** identifies the variation that occours in the sequences relative to the reads or a set of reads.

For help type:

```sh
./gto2_fq_variation_map -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_variation_map** program needs FASTQ, FASTA or SEQ files to be used as reference and target files.

The attribution is given according to:

```sh
Usage: ./gto2_fq_variation_map <OPTIONS>... [FILE]:<...> 
[FILE]:<...>

The gto2_fq_variation_map is a tool to map relative 
singularity regions. The (probabilistic) Bloom filter is 
automatically set.   
                                                         
  -v                       verbose mode,                 
  -a                       about CHESTER,                
  -s <value>               bloom size,                   
  -i                       use inversions,               
  -p                       show positions/words,         
  -k <value>               k-mer size (up to 30),        
                                                         
  [rFile1]:<rFile2>:<...>  reference file(s),            
  [tFile1]:<tFile2>:<...>  target file(s).               
                                                         
The reference files may be FASTA, FASTQ or DNA SEQ,      
while the target files may be FASTA OR DNA SEQ.          
Report bugs to <{pratas,raquelsilva,ap,pjf}@ua.pt>.
```

An example of a reference file (Multi-FASTA format) is:

```sh
>AB000264 |acc=AB000264|descr=Homo sapiens mRNA 
ACAAGACGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCCTGGAGGGTCCACCGCT
CCCTGCTGCCATTGTCCCCGGCCCCACCTAAGGAAAAGCAGCCTCCTGACTTTCCTCGCT
TGGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAAGTGGTTTGAGTGGACCTCCGG
GCCCCTCATAGGAGAGGAAGCTCGGGAGGTGGCCAGGCGGCAGGAAGCAGGCCAGTGCCG
CGAATCCGCGCGCCGGGACAGAATCTCCTGCAAAGCCCTGCAGGAACTTCTTCTGGAAGA
CCTTCTCCACCCCCCCAGCTAAAACCTCACCCATGAATGCTCACGCAAGTTTAATTACAG
ACCTGAA
>AB000263 |acc=AB000263|descr=Homo sapiens mRNA 
ACAAGATGCCATTGTCCCCCGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCACCG
CTGCCCTGCCCCTGGAGGGTGGCCCCACCGGCCGAGACAGCGAGCATATGCAGGAAGCGG
CAGGAATAAGGAAAAGCAGCCTCCTGACTTTCCTCGCTTGGTGGTTTGAGTGGACCTCCC
AGGCCAGTGCCGGGCCCCTCATAGGAGAGGAAGCTCGGGAGGTGGCCAGGCGGCAGGAAG
GCGCACCCCCCCAGCAATCCGCGCGCCGGGACAGAATGCCCTGCAGGAACTTCTTCTGGA
AGACCTTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAGTTTAATTACA
GACCTGAA
```

An example for the target file (FASTQ format) is:

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

The output of the **gto2_fq_variation_map** program is a text file identifying the relative regions. Using the inputs above, an output example of this is the following:

```sh
111111111111111111111111111110000000000000000000000000000000
111111111111111111111111111110000000000000000000000000000000
000000000000
```

## Program gto2_fq_variation_filter

The **gto2_fq_variation_filter** filters and segments the regions of singularity from the output of **gto2_fq_variation_map**.

For help type:

```sh
./gto2_fq_variation_filter -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_variation_filter** program needs a the output of **gto2_fq_variation_map** to compute.

The attribution is given according to:

```sh
Usage: ./gto2_fq_variation_filter <OPTIONS>... [FILE]:<...>      
The gto2_fq_variation_filter is a tool to filter maps.  
                                                     
  -v                       verbose mode,             
  -a                       about CHESTER,            
  -t <value>               threshold [0.0;1.0],      
  -w <value>               window size,              
  -u <value>               sub-sampling,             
                                                     
  [tFile1]:<tFile2>:<...>  target file(s).           
                                                     
The target files may be generated by gto2_fq_variation_map.    
Report bugs to <{pratas,raquelsilva,ap,pjf}@ua.pt>. 
```

An example of such an input file is:

```sh
111111111111111111111111111110000000000000000000000000000000
111111111111111111111111111110000000000000000000000000000000
000000000000
```

### Output {-}

The output of the **gto2_fq_variation_filter** program is a text file with the coordenates of the segmented regions. Using the inputs above, an output example of this is the following:

```sh
#132#132
30:60
90:130
```

## Program gto2_fq_variation_visual

The **gto2_fq_variation_visual** depites the regions of singularity using the output from **gto2_fq_variation_filter** into an SVG image.

For help type:

```sh
./gto2_fq_variation_visual -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_variation_visual** program needs a the output of **gto2_fq_variation_filter** to compute.

The attribution is given according to:

```sh
Usage: ./gto2_fq_variation_visual <OPTIONS>... [FILE]:<...>
The gto2_fq_variation_visual is a tool to visualize relative
singularity regions.
                                                     
  -v                       verbose mode,             
  -a                       about CHESTER,            
  -e <value>               enlarge painted regions,  
                                                     
  [tFile1]:<tFile2>:<...>  target file(s).           
                                                     
Report bugs to <{pratas,raquelsilva,ap,pjf}@ua.pt>. 
```

An example of such an input file is:

```sh
#132#132
30:60
90:130
```

### Output {-}

The output of the **gto2_fq_variation_visual** program is a SVG plot with the maps. In the following Figure, is represented the plot using the input above.

<div class="figure">
<img src="images/gto2FqVariationVisual.png" alt="Execution plot of the variation visual tool using the previous input." width="100%" />
<p class="caption">(\#fig:filters)Execution plot of the variation visual tool using the previous input.</p>
</div>

## Program gto2_fq_metagenomics

The **gto2_fq_metagenomics** is an ultra-fast method to infer metagenomic composition of sequenced reads relative to a database. gto2_fq_metagenomics measures similarity between any FASTQ file (or FASTA), independently from the size, against any multi-FASTA database, such as the entire set of complete genomes from the NCBI. gto2_fq_metagenomics supports single reads, paired-end reads, and compositions of both. It has been tested in many plataforms, such as Illumina MySeq, HiSeq, Novaseq, IonTorrent.

gto2_fq_metagenomics is efficient to detect the presence and authenticate a given species in the FASTQ reads. The core of the method is based on relative data compression. gto2_fq_metagenomics uses variable multi-threading, without multiplying the memory for each thread, being able to run efficiently in a common laptop.

For help type:

```sh
./gto2_fq_metagenomics -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_fq_metagenomics** program needs a FASTQ file to compute.

The attribution is given according to:

```sh
NAME                                                                     
      The gto2_fq_metagenomics is a tool to infer 
      metagenomic composition.            
                                                                         
SYNOPSIS                                                                 
      gto2_fq_metagenomics [OPTION]... [FILE1]:[FILE2]:... 
      [FILE]                      
                                                                         
SAMPLE                                                                   
      gto2_fq_metagenomics -v -F -l 47 -Z -y pro.com 
      reads1.fq:reads2.fq DB.fa         
                                                                         
DESCRIPTION                                                              
      It infers metagenomic sample composition of sequenced 
      reads. The core of the method uses a cooperation 
      between multiple context and tolerant context models 
      with several depths. The reference sequences must be 
      in a multi-FASTA format. The sequenced reads must be 
      trimmed and in FASTQ format.           
                                                                         
      Non-mandatory arguments:                                           
                                                                         
      -h                   give this help,                               
      -F                   force mode (overwrites top file),             
      -V                   display version number,                       
      -v                   verbose mode (more information),              
      -Z                   database local similarity,                    
      -s                   show compression levels,                      
                                                                         
      -l <level>           compression level [1;47],                    
      -p <sample>          subsampling (default: 1),                    
      -t <top>             top of similarity (default: 20),              
      -n <nThreads>        number of threads (default: 2),              
                                                                         
      -x <FILE>            similarity top filename,                      
      -y <FILE>            profile filename (-Z must be on).             
                                                                         
      Mandatory arguments:                                               
                                                                         
      [FILE1]:[FILE2]:...  metagenomic filename (FASTQ),                 
                           Use ":" for splitting files.                
                                                                         
      [FILE]               database filename (Multi-FASTA).                
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

The output of the **gto2_fq_metagenomics** program is a CSV file (top.csv) with the highest probability of being contained in the samples. An example for this CSV file is the following:

```sh
2  66725   12.263   NC_037703.1_Saccharomycodes_ludwigii...
1  66725   12.263   NC_037703.1_Saccharomycodes_ludwigii...
3  107123  11.492   NC_012621.1_Nakaseomyces_bacillispor...
4  107123  11.492   NC_012621.1_Nakaseomyces_bacillispor...
5  16592   11.153   NC_024030.1_Equus_przewalskii_mitoch...
6  14583   10.851   NC_021120.1_Bursaphelenchus_mucronat...
7  162504  10.607   NC_018415.1_Candidatus_Carsonella_ru...
8  10315   10.586   NC_016117.1_Mnemiopsis_leidyi_mitoch...
9  162589  10.550   NC_018414.1_Candidatus_Carsonella_ru...
10 166163  10.476   NC_018416.1_Candidatus_Carsonella_ru...
```
