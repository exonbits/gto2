# Amino Acid Tools

A more specific subset of tools is the Amino Acid Sequence tools, designed to manipulate amino acid sequences. The main features of those tools are grouping sequences, for instance by their properties, such as electric charge (positive and negative), uncharged side chains, hydrophobic side chains and special cases. It is also possible generating pseudo-DNA with characteristics passed by amino acid sequences, or for data compression, using cooperation between multiple contexts and substitutional tolerant context models. The current available amino acid sequence tools, for analysis and manipulation, are:

- **gto2_aa_to_group**: to convert an amino acid sequence to a group sequence.
- **gto2_aa_to_pseudo_dna**: to convert an amino acid (protein) sequence to a pseudo DNA sequence.
- **gto2_aa_compressor**: a new lossless compressor to compress efficiently amino acid sequences (proteins).
- **gto2_aa_from_fa**: to convert DNA sequences in FASTA or Multi-FASTA file format to an amino acid sequence.
- **gto2_aa_from_fq**: to convert DNA sequences in the FASTQ file format to an amino acid sequence.
- **gto2_aa_from_seq**: to convert DNA sequences to an amino acid sequence.

## Program gto2_aa_to_group

The **gto2_aa_to_group** converts an amino acid sequence to a group 
sequence.

For help type:

```sh
./gto2_aa_to_group -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_aa_to_group** program needs two streams for the computation, namely the input and output standard. The input stream is an amino acid sequence. 
The attribution is given according to:

```sh
Usage: ./gto2_aa_to_group [options] [[--] args]
   or: ./gto2_aa_to_group [options]

It converts a amino acid sequence to a group sequence.

    -h, --help        show this help message and exit

Basic options
    < input.prot      Input amino acid sequence file (stdin)
    > output.group    Output group sequence file (stdout)

Example: ./gto2_aa_to_group < input.prot > output.group
Table:
Prot	Group
R		P
H		P  Electric charged side chains: POSITIVE
K		P
-		-
D		N
E		N  Electric charged side chains: NEGATIVE
-		-
S		U
T		U
N		U  Amino acids with electric UNCHARGED side chains
Q		U
-		-
C		S
U		S
G		S  Special cases
P		S
-		-
A		H
V		H
I		H
L		H
M		H  Amino acids with hydrophobic side chains
F		H
Y		H
W		H
-		-
*		*  Others
X		X  Unknown
```

This tool can be used to group amino acids by properties, such as electric charge (positive
and negative), uncharged side chains, hydrophobic side chains and special cases.
An example of such an input file is:

```sh
IPFLLKKQFALADKLVLSKLRQLLGGRIKMMPCGGAKLEPAIGLFFHAIGINIKLGYGMT
ETTATVSCWHDFQFNPNSIGTLMPKAEVKIGENNEILVRGGMVMKGYYKKPEETAQAFTE
DGFLKTGDAGEFDEQGNLFITDRIKELMKTSNGKYIAPQYIESKIGKDKFIEQIAIIADA
KKYVSALIVPCFDSLEEYAKQLNIKYHDRLELLKNSDILKMFE
```

### Output {-}

The output of the **gto2_aa_to_group** program is a group sequence.
Using the input above, an output example of this is the following:

```sh
HSHHHPPUHHHHNPHHHUPHPUHHSSPHPHHSSSSHPHNSHHSHHHPHHSHUHPHSHSHU
NUUHUHUSHPNHUHUSUUHSUHHSPHNHPHSNUUNHHHPSSHHHPSHHPPSNNUHUHHUN
NSHHPUSNHSNHNNUSUHHHUNPHPNHHPUUUSPHHHSUHHNUPHSPNPHHNUHHHHHNH
PPHHUHHHHSSHNUHNNHHPUHUHPHPNPHNHHPUUNHHPHHN
```

## Program gto2_aa_to_pseudo_dna

The **gto2_aa_to_pseudo_dna** converts an amino acid (protein) sequence to a pseudo DNA sequence.

For help type:

```sh
./gto2_aa_to_pseudo_dna -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_aa_to_pseudo_dna** program needs two streams for the computation, namely the input and output standard. The input stream is an amino acid sequence.
The attribution is given according to:

```sh
Usage: ./gto2_aa_to_pseudo_dna [options] [[--] args]
   or: ./gto2_aa_to_pseudo_dna [options]

It converts a protein sequence to a pseudo DNA sequence.

    -h, --help        show this help message and exit

Basic options
    < input.prot      Input amino acid sequence file (stdin)
    > output.dna      Output DNA sequence file (stdout)

Example: ./gto2_aa_to_pseudo_dna < input.prot > output.dna
Table:
Prot	DNA
A		GCA
C		TGC
D		GAC
E		GAG
F		TTT
G		GGC
H		CAT
I		ATC
K		AAA
L		CTG
M		ATG
N		AAC
P		CCG
Q		CAG
R		CGT
S		TCT
T		ACG
V		GTA
W		TGG
Y		TAC
*		TAG
X		GGG
```

It can be used to generate pseudo-DNA with characteristics passed by amino acid (protein) sequences. An example of such an input file is:

```sh
IPFLLKKQFALADKLVLSKLRQLLGGRIKMMPCGGAKLEPAIGLFFHAIGINIKLGYGMT
ETTATVSCWHDFQFNPNSIGTLMPKAEVKIGENNEILVRGGMVMKGYYKKPEETAQAFTE
DGFLKTGDAGEFDEQGNLFITDRIKELMKTSNGKYIAPQYIESKIGKDKFIEQIAIIADA
KKYVSALIVPCFDSLEEYAKQLNIKYHDRLELLKNSDILKMFE
```

### Output {-}

The output of the **gto2_aa_to_pseudo_dna** program is a DNA sequence. Using the input above, an output example of this is the following:

```sh
ATCCCGTTTCTGCTGAAAAAACAGTTTGCACTGGCAGACAAACTGGTACTGTCTAAACTG
CGTCAGCTGCTGGGCGGCCGTATCAAAATGATGCCGTGCGGCGGCGCAAAACTGGAGCCG
GCAATCGGCCTGTTTTTTCATGCAATCGGCATCAACATCAAACTGGGCTACGGCATGACG
GAGACGACGGCAACGGTATCTTGCTGGCATGACTTTCAGTTTAACCCGAACTCTATCGGC
ACGCTGATGCCGAAAGCAGAGGTAAAAATCGGCGAGAACAACGAGATCCTGGTACGTGGC
GGCATGGTAATGAAAGGCTACTACAAAAAACCGGAGGAGACGGCACAGGCATTTACGGAG
GACGGCTTTCTGAAAACGGGCGACGCAGGCGAGTTTGACGAGCAGGGCAACCTGTTTATC
ACGGACCGTATCAAAGAGCTGATGAAAACGTCTAACGGCAAATACATCGCACCGCAGTAC
ATCGAGTCTAAAATCGGCAAAGACAAATTTATCGAGCAGATCGCAATCATCGCAGACGCA
AAAAAATACGTATCTGCACTGATCGTACCGTGCTTTGACTCTCTGGAGGAGTACGCAAAA
CAGCTGAACATCAAATACCATGACCGTCTGGAGCTGCTGAAAAACTCTGACATCCTGAAA
ATGTTTGAG
```

## Program gto2_aa_compressor

The **gto2_aa_compressor** is a new lossless compressor to compress efficiently amino acid sequences (proteins). It uses a cooperation between multiple context and substitutional tolerant context models. The cooperation between models is balanced with weights that benefit the models with better performance according to a forgetting function specific for each model.

For help type:

```sh
./gto2_aa_compressor -h
```

The **gto2_aa_compressor** program needs a file with amino acid sequences to compress. In the following example, it will be downloaded nine amino acid sequences and compress and decompress one of the smallest (HI). Finally, it compares if the uncompressed sequence is equal to the original.

```sh
wget http://sweet.ua.pt/pratas/datasets/AminoAcidsCorpus.zip
unzip AminoAcidsCorpus.zip
cp AminoAcidsCorpus/HI .
./gto2_aa_compressor -v -l 2 HI
./gto2_aa_decompressor -v HI.co
cmp HI HI.de
```

## Program gto2_aa_from_fa

The **gto2_aa_from_fasta** converts DNA sequences in FASTA or Multi-FASTA file format to an amino acid sequence.

For help type:

```sh
./gto2_aa_from_fasta -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_aa_from_fasta** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTA or Multi-FASTA file.

The attribution is given according to:

```sh
Usage: ./gto2_aa_from_fasta [options] [[--] args]
   or: ./gto2_aa_from_fasta [options]

It converts FASTA or Multi-FASTA file format to an amino 
acid sequence (translation).

    -h, --help        Show this help message and exit

Basic options
    < input.mfasta    Input FASTA or Multi-FASTA file format 
    				  (stdin)
    > output.prot     Output amino acid sequence file 
    				  (stdout)

Optional
    -f    	          Translation codon frame (1, 2 or 3)

Example: ./gto2_aa_from_fasta < input.mfasta > output.prot
```

An example of such an input file is:

```sh
>AB000264 |acc=AB000264|descr=Homo sapiens mRNA 
ACAAGACGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCCTGGAGGGTCCACCGCT
GCCCTGCTGCCATTGTCCCCGGCCCCACCTAAGGAAAAGCAGCCTCCTGACTTTCCTCGC
TTGGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAAGTGGTTTGAGTGGACCTCCG
GGCCCCTCATAGGAGAGGAAGCTCGGGAGGTGGCCAGGCGGCAGGAAGCAGGCCAGTGCC
GCGAATCCGCGCGCCGGGACAGAATCTCCTGCAAAGCCCTGCAGGAACTTCTTCTGGAAG
ACCTTCTCCACCCCCCCAGCTAAAACCTCACCCATGAATGCTCACGCAAGTTTAATTACA
GACCTGAA
>AB000263 |acc=AB000263|descr=Homo sapiens mRNA 
ACAAGATGCCATTGTCCCCCGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCACCG
CTGCCCTGCCCCTGGAGGGTGGCCCCACCGGCCGAGACAGCGAGCATATGCAGGAAGCGG
CAGGAATAAGGAAAAGCAGCCTCCTGACTTTCCTCGCTTGGTGGTTTGAGTGGACCTCCC
AGGCCAGTGCCGGGCCCCTCATAGGAGAGGAAGCTCGGGAGGTGGCCAGGCGGCAGGAAG
GCGCACCCCCCCAGCAATCCGCGCGCCGGGACAGAATGCCCTGCAGGAACTTCTTCTGGA
AGACCTTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAGTTTAATTACA
GACCTGAA
```

### Output {-}

The output of the **gto2_aa_from_fasta** program is an amino acid sequence. Using the input above, an output example of this is the following:

```sh
TRRPPAAAALRGHGPGGSTAALLPLSPAPPKEKQPPDFPRLGRDSEHMQEAAGSGLSGPP
GPS-ERKLGRWPGGRKQASAANPRAGTESPAKPCRNFFWKTFSTPPAKTSPMNAHASLIT
DLTRCHCPPASCCCCSPGPRPPLPCPWRVAPPAETASICRKRQE-GKAAS-LSSLGGLSG
PPRPVPGPS-ERKLGRWPGGRKAHPPSNPRAGTECPAGTSSGRPSPPANKTSPMNAHASL
ITDL
```

## Program gto2_aa_from_fq

The **gto2_aa_from_fastq** converts DNA sequences in the FASTQ file format to an amino acid sequence.

For help type:

```sh
./gto2_aa_from_fastq -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_aa_from_fastq** program needs two streams for the computation, namely the input and output standard. The input stream is a FASTQ file.

The attribution is given according to:

```sh
Usage: ../../bin/gto2_aa_from_fastq [options] [[--] args]
   or: ../../bin/gto2_aa_from_fastq [options]

It converts FASTQ file format to an amino acid sequence 
(translation).

    -h, --help      Show this help message and exit

Basic options
    < input.fastq   Input FASTQ file format (stdin)
    > output.prot   Output amino acid sequence file (stdout)

Optional
    -f              Translation codon frame (1, 2 or 3)

Example: ./gto2_aa_from_fastq < input.fastq > output.prot
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

The output of the **gto2_aa_from_fastq** program is an amino acid sequence. Using the input above, an output example of this is the following:

```sh
G-WPLPMASNPTKLPLTT-GFSNRVQGYDVCILRI-SRSR--YASFYH
```


## Program gto2_aa_from_seq

The **gto2_aa_from_seq** converts DNA sequence to an amino acid sequence.

For help type:

```sh
./gto2_aa_from_seq -h
```

In the following subsections, we explain the input and output parameters.

### Input parameters {-}

The **gto2_aa_from_seq** program needs two streams for the computation, namely the input and output standard. The input stream is a DNA sequence.

The attribution is given according to:

```sh
Usage: ./gto2_aa_from_seq [options] [[--] args]
   or: ./gto2_aa_from_seq [options]

It converts DNA sequence to an amino acid sequence 
(translation).

    -h, --help      Show this help message and exit

Basic options
    < input.seq     Input sequence file (stdin)
    > output.prot   Output amino acid sequence file (stdout)

Optional
    -f     			Translation codon frame (1, 2 or 3)

Example: ./gto2_aa_from_seq < input.seq > output.prot
```

An example of such an input file is:

```sh
ACAAGACGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCCTGGAGGGTCCACCGCT
GCCCTGCTGCCATTGTCCCCGGCCCCACCTAAGGAAAAGCAGCCTCCTGACTTTCCTCGC
TTGGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAAGTGGTTTGAGTGGACCTCCG
GGCCCCTCATAGGAGAGGAAGCTCGGGAGGTGGCCAGGCGGCAGGAAGCAGGCCAGTGCC
GCGAATCCGCGCGCCGGGACAGAATCTCCTGCAAAGCCCTGCAGGAACTTCTTCTGGAAG
ACCTTCTCCACCCCCCCAGCTAAAACCTCACCCATGAATGCTCACGCAAGTTTAATTACA
GACCTGAAACAAGATGCCATTGTCCCCCGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCAC
GGCCACCGCTGCCCTGCCCCTGGAGGGTGGCCCCACCGGCCGAGACAGCGAGCATATGCA
GGAAGCGGCAGGAATAAGGAAAAGCAGCCTCCTGACTTTCCTCGCTTGGTGGTTTGAGTG
GACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGAGAGGAAGCTCGGGAGGTGGCCAGGCG
GCAGGAAGGCGCACCCCCCCAGCAATCCGCGCGCCGGGACAGAATGCCCTGCAGGAACTT
CTTCTGGAAGACCTTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAGTT
TAATTACAGACCTGAA
```

### Output {-}

The output of the **gto2_aa_from_seq** program is an amino acid sequence. Using the input above, an output example of this is the following:

```sh
TRRPPAAAALRGHGPGGSTAALLPLSPAPPKEKQPPDFPRLGRDSEHMQEAAGSGLSGPP
GPS-ERKLGRWPGGRKQASAANPRAGTESPAKPCRNFFWKTFSTPPAKTSPMNAHASLIT
DLKQDAIVPRPPAAAALRGHGHRCPAPGGWPHRPRQRAYAGSGRNKEKQPPDFPRLVV-V
DLPGQCRAPHRRGSSGGGQAAGRRTPPAIRAPGQNALQELLLEDLLLLQIKPHP-MLTQV
-LQT-
```
