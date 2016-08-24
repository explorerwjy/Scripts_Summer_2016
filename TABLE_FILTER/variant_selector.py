#!/usr/local/apps/spark-submit
#Author: Jiayao Wang <jw1509@georgetown.edu>
#"Select variants record match the certain contidion (gene, maf, individual) from big ANNOVAR text table"
import csv
from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext
import os
from pyspark.sql.functions import *
from sys import argv
import sys
import re

def usage():
	print "Select variants record match the certain contidion (gene, maf, individual)"
	print "Usages: spark-submit variants_selector.py -i big_annovar_table -h header -o filtered_table -g gene_list -c config_file -ind individual_list"
	exit()	

#return a dictionary of header item and it's index for query of data
def make_header(header_file):
	res = {}
	header = open(header_file,'rb')
	r = csv.reader(header,delimiter='\t')
	header = r.next()
	for i,item in enumerate(header):
		res[item.strip()] = i
	#print header,res
	return res,header

#return dictionary of header:index and rdd of variants
def load_data(sc,input_file,header_file):
	data_file_large = input_file
	header_file = header_file
	print "header file:",header_file
	print "big table file:",data_file_large
	header_dic,header_list = make_header(header_file)
	
	#header = clean_header(header) # make header clear
	print "Creating RDD from",data_file_large,"..."
	#if large_data_file has header in it, use this line
	#data = sc.textFile(data_file_large).zipWithIndex().filter(lambda x:x[1] > 0).map(lambda x:x[0]).map(lambda x: x.strip().split("\t"))
	data = sc.textFile(data_file_large).map(lambda x: x.strip().split("\t"))
	print "Done"
	return header_dic,header_list,data

# Read a file contains a list of genes as criteria and return a list of genes, if file not provided, return None
# Return type list/None
def read_gene_list(gene_list_file):
	if gene_list_file == "":
		return None
	gene_file = open(gene_list_file,'rb')
	gene_list = []
	for l in gene_file:
		gene_list.append(l.strip())
	gene_file.close()
	return gene_list

def read_maf_list(maf_file,header_dic): #Not in use at this version
	if maf_file == "":
		return None,None
	maf_file = open(maf_file,'rb')
	try:
		maf_list = maf_file.readline().strip().split()
		maf_cutoff = float(maf_file.readline().strip())
	except:
		print "Invalid format of maf file"
		print "first line should contains maf name seperated by space or tab, second line should be a float number as threshold."
		exit()
	for i in maf_list:
		if i not in header_dic:
			print "maf name",i,"not in header, please check input"
			print i,"will be removed from maf list"
			maf_list.remove(i)
	return maf_list,maf_cutoff

def read_individual(individual_file,header_dic):
	if individual_file == "":
		return None
	individual_file = open(individual_file,'rb')
	try:
		individuals = [indi.strip() for indi in individual_file.readlines()]
		for i in xrange(0,len(individuals)):
			individuals[i] = get_individual_name(individuals[i],header_dic)
	except:
		print "Invalid format of individual file"
		print "Each line of the file should be id of an individual"
		exit()
	for i in individuals:
		if i not in header_dic:
			print "individual id",i,"not in header, please check input"
			print i,"will be remove from individual list"
			individuals.remove(i)
	return individuals 

# Return True if ref.Gene in gene list or don't filter on gene list.
def select_on_gene(record,header_dic,gene_list):
	if gene_list == None:
		return True
	else:
		if record[6] in gene_list:
			return True
		else:
			return False

def select_on_maf(record,header_dic,maf_list):
	if maf_list == None:
		return True
	else:
		for maf,cutoff in maf_list:
			tmp = record[header_dic[maf]]
			#print maf,header_div[maf],record[header_dic[maf]]
			if tmp != '.':
				if float(tmp) > float(cutoff):
					return False
		return True

def has_allele(GT):
	if GT == '0/0' or GT == './.':
		return False
	else:
		return True

def select_on_ind(record,header_dic,individuals):
	if individuals == None:
		return True
	else:
		for ind in individuals:
			tmp = record[header_dic[ind]]
			GT = tmp.split(':')[0]
			if has_allele(GT):	
				return True
		return False

def select_on_type(record,header_dic,FUNC_TYPE,EXONIC_TYPE):
	#if FUNC_TYPE == None:
	#	return True
	#if record[header_dic['Func.refGene']] == 'exonic' and record[header_dic['ExonicFunc.refGene']] == 'synonymous SNV':
	#	return False
	#if record[header_dic['Func.refGene']] in FUNC_TYPETYPE:
	#	return True
	#elif 
	#return False

	if FUNC_TYPE == None:
		return True
	if (record[header_dic['Func.refGene']] == 'exonic') and ('exoinc' in FUNC_TYPE):
		if record[header_dic['ExonicFunc.refGene']] in EXONIC_TYPE:
			return True
		else:
			return False
	elif record[header_dic['Func.refGene']] in FUNC_TYPE:
		return True
	else:
		return False

def select_on_CADD(record,header_dic,CADD):
	if record[header_dic['CADD_phred']] == '.':
		return True
	elif record[header_dic['ExonicFunc.refGene']] == 'nonsynonymous SNV' and  float(record[header_dic['CADD_phred']]) < CADD:
		return False
	else:
		return True

# Main filter contains 
# Filter on gene list
# Filter on maf&cutoff
# Filter on individual of interested
def my_filter(record,header_list,header_dic,gene_list,individuals,FUNC_TYPE,EXONIC_TYPE,CADD,AD,Allele_ratio,maf_list):
	#try:
	if record == None:
		return False
	if record[0] == 'Chr':
		return False
	if not select_on_type(record,header_dic,FUNC_TYPE,EXONIC_TYPE):
		return False
	if not select_on_gene(record,header_dic,gene_list):
		return False
	if not select_on_maf(record,header_dic,maf_list):
		return False
	if not select_on_ind(record,header_dic,individuals):
		return False
	#if not select_on_AD(record,header_dic,AD,Allele_ratio):
	#	return False
	if not select_on_CADD(record,header_dic,CADD):
		return False
	#except:
	#	print "Unexpected error:", sys.exc_info()[0]
	return True

def trim_on_header(header_list,header_dic,individuals):
	if individuals == None:
		return header_list
	res = header_list[:63]
	for ind in individuals:
		res.append(ind)
	return res

def get_individual_name(indi,header_dic):
	
	for name in header_dic.keys():
		if indi in name:
			return name
	return indi

# remove the individuals that not in individual list
def trim_on_ind(x,header_list,header_dic,individuals):
	#split by general information and individuals
	if individuals == None:
		return x
	res = x[:63]
	samples = []
	for ind in individuals:
		samples.append(x[header_dic[ind]])
	res.extend(samples)
	return res

def my_min(a,b):
	if a >= b:
		return b
	else:
		return a
def my_max(a,b):
	if a >= b:
		return a
	else:
		return b

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

			if my_min(dx,dy) >= DP: # minor allele also larger than DP, trust this variant
				indi_has_the_var += 1
				return ':'.join(items),indi_has_the_var + 1
			elif dx+dy == 0: #no called
				return './.:.:.:.:.',indi_has_the_var
			if dy/(dx+dy) < Allele_ratio: #alts take less than Allele_ratio consider it as homo of wild (0/0)
				if dx < DP: # wild reads less than DP, consider it no called
					return './.:.:.:.:.',indi_has_the_var
				else: 
					items[0] = x + '/' + x
					return ':'.join(items),indi_has_the_var + 1
			elif dx/(dx+dy) < Allele_ratio: #wild take less than Allele_ratio consider it as homo of alt (1/1)
				#x = y
				if dy < DP: # alt reads less than DP, consider it no called
					return './.:.:.:.:.',indi_has_the_var
				else:
					items[0] = y + '/' + y
					return ':'.join(items),indi_has_the_var + 1
			else:
				if dx + dy <= DP: # allele depth less then DP, not treat as a variant.
					return './.:.:.:.:.',indi_has_the_var
					
				else: # Trust this genotype
					indi_has_the_var += 1
					return ':'.join(items),indi_has_the_var
	except ValueError:
		print items
		return './.:.:.:.:.',indi_has_the_var

def revise_genotype(record,header_list,header_dic,AD,Allele_ratio):
	indi_has_the_var = 0
	for i in xrange(63,len(record)):
		item = record[i] # item: current individual's genotype info
		items = item.split(':') # split according to format
		# items[GT,AD,DP,GQ,PL] 
		flag_site_filt = False # whether filter out this genotype or not
		if has_allele(items[0]):
			#Choose whether filter on allele ratio
			#continue
			#print record[i]
			record[i],indi_has_the_var = update_indi_genotype(item,items,indi_has_the_var,AD,Allele_ratio)
	if indi_has_the_var == 0:
		return None
	else:
		return record

def revise_variants(record):
	if record == None:
		return False
	else:
		indi_has_the_var = 0
		for i in xrange(63,len(record)):
			item = record[i] # item: current individual's genotype info
			items = item.split(':') # split according to format
			# items[GT,AD,DP,GQ,PL] 
			if has_allele(items[0]):
				indi_has_the_var += 1
		if indi_has_the_var == 0:
			return False
		else:
			return True

# Find variants records match the condition provided
# Return list of variatns (each variants also is a list)
def find_variant(df,header_list,header_dic,gene_list,individuals,FUNC_TYPE,EXONIC_TYPE,CADD,AD,Allele_ratio,maf_list):
	df = df.filter(lambda x : my_filter(x,header_list,header_dic,gene_list,individuals,FUNC_TYPE,EXONIC_TYPE,CADD,AD,Allele_ratio,maf_list)) # Filter by gene list, maf list, individuals
	df = df.map(lambda x:trim_on_ind(x,header_list,header_dic,individuals)) # Trimming the record according to individual list
	df = df.map(lambda x:revise_genotype(x,header_list,header_dic,AD,Allele_ratio))
	header_list = trim_on_header(header_list,header_dic,individuals)
	df = df.filter(lambda x: revise_variants(x))
	"""
	header_dic = {}
	for i,item in enumerate(header_list):
		header_dic[item] = i
	df = df.filter(lambda x : my_filter(x,header_list,header_dic,gene_list,individuals,FUNC_TYPE,EXONIC_TYPE,CADD,AD,Allele_ratio,maf_list))
	"""
	return df.collect(),header_list

# Write the result on disk according out_file
def write_out(df,out_file,header_list):
	fout = open(out_file,'wb')
	fout.write("\t".join(header_list) + "\n")
	#print df[0:10]
	for item in df:
		if item != None:
			fout.write('\t'.join(item) + "\n")
	fout.close()

def take_argv(argv):
	input_file = ""
	output_file = ""
	gene_list_file = ""
	maf_file = ""
	individual_file = ""
	header_file = ""
	i = 1
	while i < len(argv):
		if argv[i] == '-i':
			input_file = argv[i+1]
		if argv[i] == '-h':
			header_file = argv[i+1]
		if argv[i] == '-o':
			output_file = argv[i+1]
		if argv[i] == '-c':
			maf_file = argv[i+1]
		if argv[i] == '-ind':
			individual_file = argv[i+1]
		if argv[i] == '-g':
			gene_list_file = argv[i+1]
		i += 2
	flag = True
	if input_file == "":
		print "ARGV missing: Big annovar table file must be provided"
		flag = False
	if output_file == "":
		print "ARGV missing: File after filtered must be provided"
		flag = False
	if header_file == "":
		print "ARGV missing: Header file contains headers of table must provided"
		flag = False
	if (maf_file == "") and (individual_file == "") and (gene_list_file == ""):
		print "ARGV missing: At least one criteria must be provided"
		flag = False
	if flag == False:
		usage()
	else:
			return input_file,header_file,output_file,gene_list_file,maf_file,individual_file

def more_criteria(config_file,header_dic):
	FUNC_TYPE,EXONIC_TYPE,CADD,ALLELE_DEPTH,ALLELE_RATIO,MAF = None,None,None,None,None,None
	def read_maf(value):
		res = []
		value = value.split(';')
		for item in value:
			maf,cutoff = item.split(':')
			if maf not in header_dic:
				print "maf name",maf,"not in header, please check input"
				print maf,"will not included in maf list"
			else:
				res.append((maf,cutoff))
		return res
	if config_file == "":
		return FUNC_TYPE,EXONIC_TYPE,CADD,ALLELE_DEPTH,ALLELE_RATIO,MAF
	hand = open(config_file,'rb')
	try:
		for l in hand:
			if l.startswith('#'):
				continue
			Key,Value = l.strip().split('=')
			if Key == 'FUNC_TYPE':
				FUNC_TYPE = Value.split(';')
			elif Key == 'EXONIC_TYPE':
				EXONIC_TYPE = Value.split(';')
			elif Key == 'CADD':
				CADD = float(Value)
			elif Key == 'ALLELE_DEPTH':
				ALLELE_DEPTH = int(Value)
			elif Key == 'ALLELE_RATIO':
				ALLELE_RATIO = float(Value)
			elif Key == 'MAF':
				MAF = read_maf(Value)
			else:
				print 'Unknow Paramter:%s'%Key
	except ValueError:
		print "Invalid format of config file"
		print "Please follow the format of config file like this:\nTYPE=splicing;exonic;exonic,splicing\nCADD=20\nALLELE_DEPTH=5\nALLELE_RATIO=0.3\nMAF=ExAC_ALL:0.01\n"
		exit()
	for l in hand:
		if i not in header_dic:
			print "maf name",i,"not in header, please check input"
			print i,"will be removed from maf list"
			maf_list.remove(i)
	print 'FUNC_TYPE: ',FUNC_TYPE
	print 'EXONIC_TPYE: ',EXONIC_TYPE
	print 'CADD: ',CADD
	print 'ALLELE_DEPTH: ',ALLELE_DEPTH
	print 'ALLELE_RATIO: ',ALLELE_RATIO
	print 'MAF: ',MAF
	return FUNC_TYPE,EXONIC_TYPE,CADD,ALLELE_DEPTH,ALLELE_RATIO,MAF

def main():
	conf = (SparkConf().set("spark.driver.maxResultSize", "4g").set("appName","VARIANTS_FILTER").set("spark.executor.memory","2g").set("spark.driver.memory","2g"))
	sc = SparkContext(conf = conf)
	#sqlContext = SQLContext(sc)
	
	input_file,header_file,output_file,gene_list_file,config_file,individual_file = take_argv(argv)
	header_dic,header_list,df = load_data(sc,input_file,header_file) #load data as spark dataframes
	gene_list = read_gene_list(gene_list_file) #load gene list from file
	#maf_list,maf_cutoff = read_maf_list(maf_file,header_dic) #load maf and cutoff from file
	individuals = read_individual(individual_file,header_dic) # load individual list from file
	FUNC_TYPE,EXONIC_TYPE,CADD,ALLELE_DEPTH,ALLELE_RATIO,MAF = more_criteria(config_file,header_dic)
	df,header_list = find_variant(df,header_list,header_dic,gene_list,individuals,FUNC_TYPE,EXONIC_TYPE,CADD,ALLELE_DEPTH,ALLELE_RATIO,MAF) #select variants according to gene/maf/individuals/other_criteria 
	#df.filter(lambda x: x!=None)
	write_out(df,output_file,header_list)

if __name__=="__main__":
	main()

