#!/usr/local/bin/python
from sys import argv
import gzip

def usage():
	print 'This Script used to veiw variants from certain genes from .vcf file.'
	print 'python select_by_gene.py -g gene1 gene2 gene3 ... geneN -p pedfile.ped -i input.vcf -o outfile.txt'
	exit()

def gene_in_line(genes,l):
	for gene in genes:
		if '='+gene+';' in l:
			return gene
	return False

def has_allele(genotype):
	genotype = genotype.split(':')[0]
	a,b = genotype.split('/')
	if (a != '0' and a != '.') or (b != '0' and b != '.'):
		print a,b,genotype
		return True
	else:
		return False


def make_info(string):
	subs = string.split(';')
	res = {}
	for i in subs:
		pair = i.split('=')
		if len(pair) != 2:
			res[pair[0]] = None
		else:
			res[pair[0]] = pair[1]
	return res

def find_maf(string):
	INFO = make_info(string)
	#mafs:esp6500siv2_all=.;1000g2014oct_all=.;1000g2014oct_afr=.;1000g2014oct_eas=.;1000g2014oct_eur=.;ExAC_ALL=.;ExAC_AFR=.;ExAC_AMR=.;ExAC_EAS=.;ExAC_FIN=.;ExAC_NFE=.;ExAC_OTH=.;ExAC_SAS=.;
	#snp:avsnp144=.
	res = []
	mafs = ['esp6500siv2_all','1000g2014oct_all','1000g2014oct_afr','1000g2014oct_eas','1000g2014oct_eur','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS']
	for maf in mafs:
		if INFO[maf] != '.':
			res.append(maf + ':' + INFO[maf])
	res.append('  avsnp144=' + INFO['avsnp144'])
	return ';'.join(res)


def make_record(last_gene,variants,fout,indiv,header):
	total_case = 48
	total_control = 34+251
	if last_gene == None:
		return
	print last_gene
	fout.write("Report on " + last_gene + "\n")
	total = len(indiv.keys())
	indi_has_the_gene = []
	for variant in variants:
		fout.write(variant)
		case = 0
		control = 0
		unknow = 0
		info = variant.split('\t')
		indi_has_the_var = []
		for genotype in info[9:]: #individual start at 9th
			ind = header[info.index(genotype)] #get individual name from header
			if ind in indiv:
				if indiv[ind] == '1': #control
					if has_allele(genotype):
						indi_has_the_var.append(ind)
						control += 1
				elif indiv[ind] == '2': #case
					if has_allele(genotype):
						indi_has_the_var.append(ind)
						case += 1
				elif indiv[ind] == '0':
					if has_allele(genotype):
						indi_has_the_var.append(ind)
						unknow += 1
				else:
					print "ERROR! this affect state should not appear"
		fout.write(str(case) + " Cases has this variant," + "Frequency in case:"+ str(float(case)/total_case)+ "\n")
		fout.write(str(control) + " Controls has this variant" + "Frequency in control:"+ str(float(control)/total_control)+ "\n")
		fout.write("Individuals has this mutation: "+';'.join(indi_has_the_var) + "\n")
		mafs = find_maf(info[7])
		fout.write("Reported MAF=" + mafs + '\n')
		update_indi_gene(indi_has_the_var,indi_has_the_gene)
	fout.write("\nIndividuals has mutation in this gene: " + ';'.join(indi_has_the_gene)+'\n')
	case,control = get_case_control(indi_has_the_gene,indiv)
	fout.write(str(case) + " Cases in this gene,"+ str(control) + " Controls in this gene\n")
	fout.write("\n")

def get_case_control(indi_has_the_gene,indiv):
	case,control = 0,0
	for indi in indi_has_the_gene:
		if indiv[indi] == '1':
			control += 1
		else:
			case += 1
	return case,control

def update_indi_gene(indi_has_the_var,indi_has_the_gene):
	for indi in indi_has_the_var:
		if indi not in indi_has_the_gene:
			indi_has_the_gene.append(indi)
	indi_has_the_var = []

def grab_gene(infile,genes,outfile,indiv):
	fout = open(outfile,'w')
	fin = gzip.open(infile, 'rb')
	last_gene = None
	variants = []
	header = ""
	for l in fin:
		if l.startswith('##'):
			#fout.write(l)
			continue
		if l.startswith('#'):
			header = l.strip().split('\t')
			fout.write(l)
		gene = gene_in_line(genes,l)
		if gene != False:
			if gene == last_gene: #Add variants into previous gene
				variants.append(l)
			else: #A new gene located
				make_record(last_gene,variants,fout,indiv,header)
				last_gene = gene
				variants = [l]
	make_record(last_gene,variants,fout,indiv,header)

	fin.close()
	fout.close()

def get_argvs(argv):
	genes = []
	infile = ""
	outfile = ""
	pedfile = ""
	if len(argv) <= 1:
		usage()
	i = 1
	while i < len(argv):
		if argv[i] == "-i":
			infile = argv[i+1]
			i += 2
		elif argv[i] == "-o":
			outfile = argv[i+1]
			i += 2
		elif argv[i] == "-p":
			pedfile = argv[i+1]
			i += 2
		elif argv[i] == '-g':
			while True:
				#print i,len(argv)
				#print argv[i+1]
				genes.append(argv[i+1])
				i += 1
				if i+1 > len(argv)-1 or argv[i+1][0] == '-':
					break
			i += 1
		else:
			print "UNKNOWN ARGV", argv[i]
			usage()
	if len(genes) == 0:
		print "No genes provided"
		exit()
	else:
		print "GENES are:",
		for i in genes:
			print i,
		print 
	if infile == "":
		print "Input vcf file is needed"
		exit()
	else:
		print "Input VCF file is:",infile
	if outfile == "":
		print "Output text filename is needed"
		exit()
	else:
		print "Output filename is:",outfile
	if pedfile == "":
		print "PED file contains individuals and their affected state is needed"
		exit()
	else:
		print "PED file is:",pedfile
	
	return genes,infile,outfile,pedfile

#Get affected situation of each individuals from the ped file
#Result in a dictionary. Individual id as key and affected state as value
def indi_stat(pedfile):
	res = {}
	with open(pedfile,'r') as ped:
		ped.readline() #skip head
		for l in ped:
			#1th is id, 5th is disease
			l = l.strip().split('\t')
			print l
			res[l[1]] = l[5]
	ped.close()
	return res

def main():
	genes,infile,outfile,pedfile = get_argvs(argv)
	indiv = indi_stat(pedfile)
	grab_gene(infile,genes,outfile,indiv)

if __name__=="__main__":
	main()

