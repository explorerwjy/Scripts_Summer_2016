#!/usr/local/bin
#Count each Individual's Variant number from a .vcf file

from sys import argv

#input_name: .vcf file name
#outout_name: output text file include counts on each kind of variant
def input_output_name(name):
	input_name = name
	output_name = input_name.split('.')[0] + "_count.txt"
	return input_name,output_name

def make_index_header(l):
	info = l.strip().split('\t')
	individuals = info[9:]
	dic_indi = {}
	for indi in individuals:
		dic_indi[indi] = {}
	
	"""
	res = {}
	for i in len(info):
		res[i] = info[i]
	"""
	# res: headers
	# individuals as key, empty list as value to hold different types of variants
	#return res,dic_indi
	return info,dic_indi

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
					
def has_allele(genotype):
	genotype = genotype.split(':')[0]
	a,b = genotype.split('/')
	if (a != '0' and a != '.') or (b != '0' and b != '.'):
		return True
	else:
		return False

# update individual dictionary with info of variants
def update(header,info,individuals):
	#different filter may use
	#Filter info[6]
	#INFO info[7] Func.refGene=exonic
	INFO = make_info(info[7])
	if INFO['Func.refGene'] != "exonic":
		return
	variant_type = INFO['ExonicFunc.refGene']
	for i in xrange(9,len(info)):
		genotype = info[i]
		if has_allele(genotype):
			indi = header[i]
			if variant_type not in individuals[indi]:
				individuals[indi][variant_type] = 1
			else:
				individuals[indi][variant_type] += 1
	return

# write the result
def write_res(index_header,individuals,fout):
	#sequence should be same as index_header
	for indi in index_header[9:]:
		variant_pair = sorted(individuals[indi].items(), key = lambda k:k[0])
		variant_pair = [k[0] + ':' + str(k[1]) for k in variant_pair]
		#print variant_pair
		variant = "\t".join(variant_pair)
		fout.write(indi+'\t'+variant+'\n')
	return

#check each individuals' variant type of each variant
def var_check(input_name,output_name):
	fin = open(input_name,'rb')
	fout = open(output_name,'wb')
	

	for l in fin:
		if l[0:2] == "##":
			continue #skip meta
		elif l[0] == "#":
			#headers
			index_header,individuals = make_index_header(l) 
		else: #variants
			info = l.strip().split('\t')
			update(index_header,info,individuals)
	fin.close()
	#print individuals
	write_res(index_header,individuals,fout)		
	fout.close()
	return

def main():
	input_name,output_name = input_output_name(argv[1])
	var_check(input_name,output_name)


if __name__=='__main__':
	main()
