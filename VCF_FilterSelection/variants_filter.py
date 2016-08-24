#!/usr/local/bin/python
# Leave out the variants that not match our requirements
# python variants_filter.py input.vcf output.vcf

from sys import argv


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

# genotype is not 0/0 or ./.
def has_allele(GT):
	if GT == '0/0' or GT == './.':
		return False
	else:
		return True
"""
# Modity the genotype of individual
def update_indi_genotype(item,items,indi_has_the_var,min_fre,DP):
	x,y = items[0].split('/')
	if x == y: #homozygous
		if items[2]== '0' or items[2] == '.' :
			#print "No AD at this site"
			return './.:.:.:.:.',indi_has_the_var
		return ':'.join(items),indi_has_the_var+1
	else: #heterozygous
		if items[1] == '.':
			return './.:.:.:.:.',indi_has_the_var
		AD = items[1].split(',')
		#print item,AD,int(x),int(y)
		dx,dy = map(float,[AD[int(x)],AD[int(y)]])
		
		if min(dx,dy) >= DP: # minor allele also larger than DP, trust this variant
			indi_has_the_var += 1
		elif dx + dy <= DP: # allele depth less than DP, not treat as a variant.
			#print "remove the site for DP <=",DP
			items[0] = '0/0'
		elif dy/(dx+dy) < min_fre:
			#print "modify the site",items,dy/(dx+dy)
			items[0] = x + '/' + x
		elif dx/(dx+dy) < min_fre:
			#print "modify the site",items,dx/(dx+dy)
			items[0] = y + '/' + y
		else:
			indi_has_the_var += 1
		return ':'.join(items),indi_has_the_var
"""
def update_indi_genotype(item,items,indi_has_the_var,DP,Allele_ratio):
	try:
		x,y = items[0].split('/')
		if x==y: #homozygous
			if items[2] == '0' or items[2] == '.':
				return './.:.:.:.:.',indi_has_the_var
			elif int(items[2]) < DP:
				return './.:.:.:.:.',indi_has_the_var
			return ':'.join(items),indi_has_the_var + 1
		else: #heterozygous
			if items[1] == '.':
				return './.:.:.:.:.',indi_has_the_var
			AD = items[1].split(',')
			dx,dy = map(float,[AD[int(x)],AD[int(y)]])
			if min(dx,dy) >= DP: # minor allele also larger than DP, trust this variant
				indi_has_the_var += 1
				return ':'.join(items),indi_has_the_var + 1
			elif dx+dy == 0: #no called
				return './.:.:.:.:.',indi_has_the_var
			if dy/(dx+dy) < Allele_ratio:
				if dx < DP:
					return './.:.:.:.:.',indi_has_the_var
				else:
					items[0] = x + '/' + x
					return ':'.join(items),indi_has_the_var + 1
			elif dx/(dx+dy) < Allele_ratio:
				x = y
				if dy < DP:
					return './.:.:.:.:.',indi_has_the_var
				else:
					items[0] = y + '/' + y
					return ':'.join(items),indi_has_the_var + 1
			else:
				if dx + dy <= DP: # allele depth less then DP, not treat as a variant.
					return './.:.:.:.:.',indi_has_the_var
				else:
					indi_has_the_var += 1
					return ':'.join(items),indi_has_the_var
	except ValueError:
		print "Error Occur!!"
		print items
		return './.:.:.:.:.',indi_has_the_var

def pass_allele_depth_check(l,min_fre,DP):
	data = l.strip().split()
	INFO = make_info(data[7])
	
	indi_has_the_var = 0 # how many individuals has this variant. if None of them pass the check, this variant will be discard.
	# Check on each individual
	for i in xrange(9,len(data)):
		item = data[i] # item: current individual's genotype
		# items[GT,AD,DP,GQ,PL] 
		items = item.split(':') # split according to format
		flag_site_filt = False # whether filter out this genotype or not
		if has_allele(items[0]):
			data[i],indi_has_the_var = update_indi_genotype(item,items,indi_has_the_var,DP,min_fre)

	if indi_has_the_var == 0:
		return False
	else:
		data[0] = data[0].strip('chr')
		new_l = '\t'.join(data)+"\n"
		return new_l

def pass_cadd_check(flag_cadd_check,INFO, cadd_cutoff):
	if flag_cadd_check == False:
		return True
	if INFO['CADD_phred'] == '.':
		return True
	elif float(INFO['CADD_phred']) < cadd_cutoff:
		return False
	else:
		return True

def my_filter(vcf_file, out_vcf, min_fre, DP, ref_func, exonic_func, maf, flag_cadd_check=False, cadd_cutoff=0):
	fout = open(out_vcf,'wb')
	fin = open(vcf_file,'rb')

	pass_ref_func = False
	pass_exonic_func = False
	pass_cadd = not flag_cadd_check
	pass_maf = False
	
	counter = 0
	for l in fin:
		
		pass_ref_func = False
		pass_exonic_func = False
		pass_cadd = not flag_cadd_check
		pass_maf = False
		
		counter += 1
		if counter % 1000 == 0:
			print "process",counter
		
		if l[0:2] == "##": #meta info
			fout.write(l)
		elif l[0] == "#": #header
			fout.write(l)
			header = l.strip().split('\t')
		else: #data
			data = l.strip().split('\t')
			if data[6] != "PASS":
				continue
			INFO = make_info(data[7])
			# First filter on Func.refGene
			if INFO['Func.refGene'] in ref_func:
				pass_ref_func = True
			# Second filter on ExonicFunc.refGene
			if INFO['Func.refGene'] == 'exonic' and INFO['ExonicFunc.refGene'] in exonic_func:
				pass_exonic_func = False
			else:
				pass_exonic_func = True
			# Third filter on flag_cadd_check
			pass_cadd = pass_cadd_check(flag_cadd_check, INFO, cadd_cutoff)
			# Fourth filter on maf
			if INFO[maf[0]] == '.':
				pass_maf = True
			elif float(INFO[maf[0]]) < maf[1]:
				pass_maf = True
			# If pass previous filter
			if pass_ref_func and pass_exonic_func and pass_cadd and pass_maf:
				# Finally filter on allele depth
				res = pass_allele_depth_check(l,min_fre,DP)
				if res != False:
					fout.write(res)
	fout.close()
	fin.close()

def main():
	#take_argv(argv)
	vcf_file = argv[1]
	out_vcf = argv[2]
	
	######### Filter Parameters
	allele_ratio = 0.3 
	DP = 5

	flag_cadd_check = False
	cadd_cutoff = 10
	
	ref_func = ['exonic',r'splicing',r'exonic\x3bsplicing'] ## ref.func in this types
	exonic_func = ['synonymous_SNV']
	#exonic_func = ['unknown','synonymous_SNV','nonframeshift_insertion','nonframeshift_deletion'] ## exnoic func not in these types
	
	maf = ['ExAC_NFE',0.01]
	
	print "MIN_DEPTH =",allele_ratio
	print "DP =",DP
	print "CADD_Phred check:",flag_cadd_check, "CADD cutoff",cadd_cutoff
	print "MAF type",maf[0],"MAF cutoff",maf[1]
	print "ref_func",ref_func
	print "exonic_func",exonic_func

	######### Run Filter
	my_filter(vcf_file, out_vcf, allele_ratio, DP, ref_func, exonic_func, maf, flag_cadd_check, cadd_cutoff)

if __name__=='__main__':
	main()
