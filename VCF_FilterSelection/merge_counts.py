#!/usr/local/bin/python
#Merge All the count result from each chromosome into one
from sys import argv
from os import path



#get all the variant files
def get_files():
	res  = []
	for i in range(1,23):
		prefix = str(i)
		if len(prefix) == 1:
			prefix = "0" + prefix
		res.append("chr" + prefix + "_exomefiltered_recalibrated_SNP_INDEL_hg19_multianno_count.txt")
	print res
	return res

def get_init(_file):

	individuals = {}
	indi_sequence = []
	with open(_file,'rb') as fin:
		for l in fin:
			meta = l.strip().split('\t')
			indi = meta[0]
			if indi not in individuals:
				individuals[indi] = {}
				indi_sequence.append(indi)
			for v_type in meta[1:]:
				v_type = v_type.split(":")
				if v_type[0] not in individuals[indi]:
					individuals[indi][v_type[0]] = int(v_type[1])

				else:
					individuals[indi][v_type[0]] += int(v_type[1])
	fin.close()
	return individuals,indi_sequence

def update(individuals,_file):

	with open(_file,'rb') as fin:
		for l in fin:
			meta = l.strip().split('\t')
			indi = meta[0]
			for v_type in meta[1:]:
				v_type = v_type.split(":")
				if v_type[0] not in individuals[indi]:
					individuals[indi][v_type[0]] = int(v_type[1])
				else:
					individuals[indi][v_type[0]] += int(v_type[1])
	fin.close()

def merge_write(individuals,indi_sequence):
	header = ['individual_name','synonymous_SNV','nonsynonymous_SNV','stopgain','stoploss','frameshift_insertion','frameshift_deletion','nonframeshift_insertion','nonframeshift_deletion','unknown']
	fout = open("individual_variant_count.txt",'wb')
	fout.write("\t".join(header)+"\n")
	
	for indi in indi_sequence:
		record = [indi]
		for v_type in header[1:]:
			record.append(str(individuals[indi][v_type]))
		record = '\t'.join(record) + '\n'
		
		fout.write(record)


#Read each file, load results into dictionary, merge with previous one
def merge(files):
	individuals,indi_sequence = get_init(files[0])
	for _file in files[1:]:
		update(individuals,_file)
	merge_write(individuals,indi_sequence)

	return

def main():
	files = get_files()
	merge(files)

if __name__=="__main__":
	main()
