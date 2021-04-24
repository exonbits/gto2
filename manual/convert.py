#!/usr/bin/python

import sys

f = open(sys.argv[1], "r")
data=f.read()
data=data.replace("\\begin{lstlisting}", "\n```sh") 
data=data.replace("\\end{lstlisting}", "```\n") 
data=data.replace("\\subsection*{Output}", "### Output {-}\n") 
data=data.replace("\\subsection*{Input parameters}", "### Input parameters {-}\n") 
data=data.replace("\\char`","") 
data=data.replace("gto_","gto2_")
data=data.replace("gto2_fastq","gto2_fq")
data=data.replace("gto2_fasta","gto2_fa")
data=data.replace("\nUsing the input above, an output example for this is the following:", "Using the input above, an output example of this is the following:")
data=data.replace("paramters", "parameters")
data=data.replace("length=72", "length=60")
data=data.replace("an output example for this is the following:","an output example of this is the following:")
data=data.replace("\\\\","\n")
data=data.replace("\n\n\n","\n\n")
newData=""
for x in data.split(" "):
	tmp = x
	if "texttt{" in x:
		tmp = tmp.replace("\\texttt{","**")
		tmp = tmp.replace("}","**")
	newData += tmp + " "

print(newData)