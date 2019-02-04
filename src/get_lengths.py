#Take fasta file of longest isoforms and calculate length
from Bio import SeqIO
import sys

genes = SeqIO.to_dict(SeqIO.parse(sys.argv[1],"fasta"))
outfile = sys.argv[2]
with open(outfile,'w') as o:
	o.write('Gene\tLength\n')

for rec in genes:
	gene = genes[rec].description.split(" ")[0]
	length = len(genes[rec])
	with open(outfile,'a') as o:
		o.write(gene+'\t'+str(length)+'\n')