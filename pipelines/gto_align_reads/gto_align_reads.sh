#!/bin/bash
#
TAG="X";
RUN="0";
THREADS="4";
REFERENCE="";
CONSENSUS=0;
FORWARD="blood1.fq";
REVERSE="blood2.fq";
#
################################################################################
#
SHOW_MENU () {
  echo " ------------------------------------------------------- ";
  echo " ";
  echo " gto_align_reads.sh >> FASTQ alignment pipeline with     ";
  echo "                       trim, sort, duplications removal, ";
  echo "                       indexes and consensus.            ";
  echo " ";
  echo " Program options --------------------------------------- ";
  echo " ";
  echo " -h, --help                     Show this, ";
  echo " -i, --install                  Installation (w/ conda), ";
  echo " ";
  echo " -c, --consensus                Create consensus seq,    ";
  echo " ";
  echo " -n <STR>, --name <STR>         TAG name for uniqueness, ";
  echo " -t <INT>, --threads <INT>      Number of threads.       ";
  echo " ";
  echo " Input options ----------------------------------------- ";
  echo " ";
  echo " -r <FILE>, --reference <FILE>  FASTA Reference filename,";
  echo " ";
  echo " -1 <FILE>, --forward <FILE>    FASTQ forward filename,  ";
  echo " -2 <FILE>, --reverse <FILE>    FASTQ reverse filename,  ";
  echo " ";
  echo " Example ----------------------------------------------- ";
  echo " ";
  echo " gto_align_reads.sh --threads 8 --reference B19.fa \\    ";
  echo " --forward FW.fq --reverse RV.fq --consensus             ";
  echo " ";
  echo " ------------------------------------------------------- ";
  }
#
################################################################################
#
CHECK_INPUT () {
  FILE=$1
  if [ -f "$FILE" ];
    then
    echo "Input filename: $FILE"
    else
    echo -e "\e[31mERROR: input file not found ($FILE)!\e[0m";
    SHOW_MENU;
    exit;
    fi
  }
#
################################################################################
#
PROGRAM_EXISTS () {
  if ! [ -x "$(command -v $1)" ];
    then
    echo -e "\e[41mERROR\e[49m: $1 is not installed." >&2;
    echo -e "\e[42mTIP\e[49m: Try: gto_align_reads.sh --install" >&2;
    exit 1;
    fi
  }
#
################################################################################
#
if [[ "$#" -lt 1 ]];
  then
  HELP=1;
  fi
#
POSITIONAL=();
#
while [[ $# -gt 0 ]]
  do
  i="$1";
  case $i in
    -h|--help|?)
      HELP=1;
      shift
    ;;
    -i|--install|--compile)
      INSTALL=1;
      shift
    ;;
    -n|--name|--tag)
      TAG="$2";
      shift 2;
    ;;
    -t|--threads)
      THREADS="$2";
      shift 2;
    ;;
    -r|--reference)
      REFERENCE="$2";
      RUN=1;
      shift 2;
    ;;
    -c|--consensus|--cons)
      CONSENSUS=1;
      shift
    ;;
    -1|--forward)
      FORWARD="$2";
      shift 2;
    ;;
    -2|--reverse)
      REVERSE="$2";
      shift 2;
    ;;
    -U|--unpaired)
      UNPAIRED="$2";
      shift 2;
    ;;
    -*) # unknown option with small
    echo "Invalid arg ($1)!";
    echo "For help, try: gto_align_reads.sh -h"
    exit 1;
    ;;
  esac
  done
#
set -- "${POSITIONAL[@]}" # restore positional parameters
#
################################################################################
#
if [[ "$HELP" -eq "1" ]];
  then
  SHOW_MENU;
  exit;
  fi
#
################################################################################
#
if [[ "$INSTALL" -eq "1" ]];
  then
  conda install -c bioconda trimmomatic --yes
  conda install -c bioconda samtools --yes
  conda install -c bioconda bowtie2 --yes
  conda install -c bioconda bedtools --yes
  conda install -c bioconda bedops --yes
  conda install -c bioconda tabix --yes
  conda install -c bioconda bcftools --yes
  exit;
  fi
#
################################################################################
#
if [[ "$RUN" -eq "1" ]];
  then
  #
  echo "Using $THREADS threads ...";
  #
  PROGRAM_EXISTS "trimmomatic";
  PROGRAM_EXISTS "bowtie2";
  PROGRAM_EXISTS "samtools";
  PROGRAM_EXISTS "bcftools";
  PROGRAM_EXISTS "bedtools";
  PROGRAM_EXISTS "tabix";
  #
  CHECK_INPUT "$REFERENCE";
  #
  CHECK_INPUT "$FORWARD";
  CHECK_INPUT "$REVERSE";
  #
  rm -f o_fw_pr.fq.gz o_fw_unpr.fq.gz o_rv_pr.fq.gz o_rv_unpr.fq.gz;
  #
  echo ">PrefixPE/1" > adapters.fa
  echo "TACACTCTTTCCCTACACGACGCTCTTCCGATCT" >> adapters.fa
  echo ">PrefixPE/2" >> adapters.fa
  echo "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" >> adapters.fa
  echo ">PE1" >> adapters.fa
  echo "TACACTCTTTCCCTACACGACGCTCTTCCGATCT" >> adapters.fa
  echo ">PE1_rc" >> adapters.fa
  echo "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA" >> adapters.fa
  echo ">PE2" >> adapters.fa
  echo "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" >> adapters.fa
  echo ">PE2_rc" >> adapters.fa
  echo "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" >> adapters.fa
  #
  trimmomatic PE -threads $THREADS -phred33 $FORWARD $REVERSE o_fw_pr.fq o_fw_unpr.fq o_rv_pr.fq o_rv_unpr.fq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
  #
  rm -f index_file_$TAG* aligned_sorted-$TAG.bam aligned_sorted-$TAG.bam.bai aligned_sorted-$TAG.bam
  #
  bowtie2-build $REFERENCE index_file_$TAG
  bowtie2 -a --threads $THREADS --very-sensitive -x index_file_$TAG -1 o_fw_pr.fq -2 o_rv_pr.fq -U "o_fw_unpr.fq,o_rv_unpr.fq" > aligned-$TAG.sam
  #
  samtools sort --threads $THREADS aligned-$TAG.sam > aligned_sorted-$TAG.bam
  samtools sort --threads $THREADS -n aligned_sorted-$TAG.bam > aligned_sorted_sorted-$TAG.bam
  samtools fixmate --threads $THREADS -m aligned_sorted_sorted-$TAG.bam aligned_sorted_sorted-$TAG-fixmate.bam
  samtools sort --threads $THREADS -o aligned_sorted_sorted-$TAG-fixmate-sort.bam aligned_sorted_sorted-$TAG-fixmate.bam
  samtools markdup --threads $THREADS -r aligned_sorted_sorted-$TAG-fixmate-sort.bam aligned_sorted-$TAG.bam
  samtools index -@ $THREADS aligned_sorted-$TAG.bam aligned_sorted-$TAG.bam.bai
  samtools flagstat -@ $THREADS aligned_sorted-$TAG.bam > Alignments-Stats-$TAG.txt
  #
  rm -f *.bt2 aligned-$TAG.sam aligned_sorted_sorted-$TAG.bam aligned_sorted_sorted-$TAG-fixmate.bam aligned_sorted_sorted-$TAG-fixmate-sort.bam
  #
  if [[ "$CONSENSUS" -eq "1" ]];
    then
    #
    Reference="$REFERENCE"; 
    Alignments="aligned_sorted-$TAG.bam";    # EXAMPLE: ttv_aligned_sorted-heart.bam
    Organ="organ";                            # Example: heart
    Label="$TAG";                             # Example: TTV
    #
    echo "Using Reference   : $Reference";
    echo "Using Alignments  : $Alignments";
    echo "Using Organ       : $Organ";
    echo "Using TAG         : $Label";
    #
    rm -f $Label-$Organ-calls.vcf.gz $Label-$Organ-calls.vcf.gz.csi $Label-$Organ-calls.norm.bcf $Label-$Organ-calls.norm.flt-indels.bcf $Label-$Organ-calls.norm.flt-indels.vcf.gz $Label-$Organ-calls.norm.flt-indels.vcf.gz.csi $Label-$Organ-calls.norm.vcf.gz;
    #
    bedtools genomecov -ibam $Alignments -bga > $Label-coverage-$Organ.bed
    awk '$4 < 1' $Label-coverage-$Organ.bed > $Label-zero-coverage-$Organ.bed
    #bedtools maskfasta -fi $Reference -bed $Label-$Organ-zero_coverage.bed -fo MASKED-$Label-consensus-$Organ.fa;
    samtools faidx $Reference # -P 9.9e-1                                         #here!
    bcftools mpileup -Ou -f $Reference $Alignments | bcftools call --ploidy 1 -P 9.9e-1 -mv -Oz -o $Label-$Organ-calls.vcf.gz
    bcftools index $Label-$Organ-calls.vcf.gz
    bcftools norm -f $Reference $Label-$Organ-calls.vcf.gz -Oz -o $Label-$Organ-calls.norm.vcf.gz
    bcftools filter --IndelGap 5 $Label-$Organ-calls.norm.vcf.gz -Oz -o $Label-$Organ-calls.norm.flt-indels.vcf.gz
    zcat $Label-$Organ-calls.norm.flt-indels.vcf.gz |vcf2bed --snvs > $Label-calls-$Organ.bed
    tabix -f $Label-$Organ-calls.norm.flt-indels.vcf.gz
    bcftools consensus -m $Label-zero-coverage-$Organ.bed -f $Reference $Label-$Organ-calls.norm.flt-indels.vcf.gz > $Label-consensus-$Organ.fa
    tail -n +2 $Label-consensus-$Organ.fa > $Label-$Organ-TMP_FILE.xki
    echo ">$Label consensus ($Organ)" > $Label-consensus-$Organ.fa
    cat $Label-$Organ-TMP_FILE.xki >> $Label-consensus-$Organ.fa
    rm -f $Label-$Organ-TMP_FILE.xki;
    #
    rm -f $Label-$Organ-calls.vcf.gz $Label-$Organ-calls.vcf.gz.csi $Label-$Organ-calls.norm.bcf $Label-$Organ-calls.norm.flt-indels.bcf $Label-$Organ-calls.norm.flt-indels.vcf.gz $Label-$Organ-calls.norm.flt-indels.vcf.gz.csi $Label-$Organ-calls.norm.vcf.gz;
    #
    fi
  #
  #
  echo "Reference  : $REFERENCE (index included)";
  echo "Alignments : aligned_sorted-$TAG.bam (index included)";
  if [[ "$CONSENSUS" -eq "1" ]];
    then
    echo "Consensus  : $Label-consensus-$Organ.fa";
  fi
  #
  fi
#
################################################################################
# 
