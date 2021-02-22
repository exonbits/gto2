#!/bin/bash
#
# $1 -> ref.fa
# $2 -> reads.fq
#
rm -rf AlienTrimmer_0.4.0.tar.gz AlienTrimmer_0.4.0/
wget http://ftp.pasteur.fr/pub/gensoft/projects/AlienTrimmer/AlienTrimmer_0.4.0.tar.gz
tar -vzxf AlienTrimmer_0.4.0.tar.gz
cd AlienTrimmer_0.4.0/src/
chmod +x JarMaker.sh
./JarMaker.sh
cp AlienTrimmer.jar ../../
cd ../../
#
java -jar AlienTrimmer.jar -i $2 -o fil_$2 -c 0 -l 30 -p 80
#
rm -f index_file* aligned-$2.sam aligned_sorted-$2.bam
bowtie2-build $1 index_file
#
bowtie2 -a --threads 4 -x index_file -U fil_$2 > aligned-$2.sam
samtools sort aligned-$2.sam > aligned_sorted-$2.bam
rm -f aligned-$2.sam 
#
samtools index aligned_sorted-$2.bam aligned_sorted-$2.bam.bai
