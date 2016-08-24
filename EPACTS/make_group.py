#!/usr/local/bin/python
#make groups(group by gene) from annovar annotated .vcf file

from sys import argv
import vcf
import re

def usage():
	print "python make_group.py -i input-vcf-file -o group_file.grp"
	exit()

def take_arg(argv):
	infile = ""
	outfile = ""

	i = 1
	while i < len(argv):
		if argv[i] == "-i":
			infile = argv[i+1]
		if argv[i] == "-o":
			outfile = argv[i+1]
		i += 2
	if infile == "":
		usage()
	if outfile == "":
		usage()
	return infile,outfile


#make a record of gene and variants in it.

def gene_marker(gene,variants):
	pass

def get_gene(INFO):
	return INFO['Gene.refGene']

#construct a maker
#format:[CHROM]:[POS]_[REF]/[ALT]
def get_var(record):
	chrom = record.CHROM
	pos = record.POS
	ref = record.REF
	alt = ','.join(map(str,record.ALT))
	return chrom.strip('chr') + ':' + str(pos) + '_' + str(ref) + '/' + str(alt)

def make_group():
	infile,outfile = take_arg(argv)
	fin = open(infile,'rb')
	fout = open(outfile,'wb')
	vcf_reader = vcf.Reader(fin)
	gene_groups = {}
	i = 0
	for record in vcf_reader:
		i += 1
		gene = get_gene(record.INFO)
		gene = ','.join(gene)
		if gene in gene_groups:
			gene_groups[gene].append(get_var(record))
		else:
			gene_groups[gene] = [get_var(record)]
		if i%1000 == 0:
			print i,"record processed"
	for k,v in gene_groups.items():
		
		fout.write(k + "\t" + "\t".join(v) + "\n")
	
	fin.close()
	fout.close()



def main():
	make_group()


if __name__=="__main__":
	main()
