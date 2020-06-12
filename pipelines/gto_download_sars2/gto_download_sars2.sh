#!/bin/bash
rm -f SARS-CoV-2.fa;
echo "Downloading SARS-CoV-2 ...";
efetch -db nucleotide -format fasta -id "MT007544,MT126808,NC_045512,MT121215,MT135042,MT135041,MT135044,MT135043,MT123290,MT123293,MT123291,MT123292,MT093631,MT049951,MT039873,MT019529,MT019532,MT019531,MT019533,MT019530,MN996529,MN996530,MN996531,MN996528,MN996527,MN988669,MN988668,MN938384,MN975262,MN908947,MT050493,MT012098,MT066156,MT072688,MT039890,MT093571,MT192759,MT066176,MT066175,MT192765,MT188339,MT188340,MT188341,MT184907,MT184913,MT184909,MT184908,MT184911,MT184910,MT184912,MT163716,MT163719,MT163718,MT163717,MT159721,MT159711,MT159710,MT159708,MT159712,MT159707,MT159715,MT159716,MT159722,MT159714,MT159713,MT159706,MT159705,MT159719,MT159709,MT159717,MT159720,MT159718,MT152824,MT118835,MT106052,MT106053,MT106054,MT044258,MT044257,MT039887,MT039888,MT027063,MT027064,MT027062,MT020880,MT020881,MN994468,MN997409,MN994467,MN988713,MN985325,MT192772,MT192773" >> SARS-CoV-2.fa;
