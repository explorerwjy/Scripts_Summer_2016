#Select gene set in test result according to a gene name list file
#Usage: python filter_by_retnet.py skat.epacts.txt retnet.txt

import os
from sys import argv

def get_filter(filter_set):
	hand = open(filter_set,'r')
	res = {}
	for l in hand:
		gene = l.strip()
		res[gene] = 1
	hand.close()
	return res

def getfilter(aimfile,filter_set):
	hand = open(aimfile,'r')
	outname = '.'.join(aimfile.split(".")[:-1]) + "_filter_by_retnet.txt"
	fout = open(outname,'w')
	fout.write(hand.readline())

	for l in hand:
		gene = l.split("\t")[3].split("_")[-1]
		if gene in filter_set:
			fout.write(l)
	hand.close()
	fout.close()


def main():
	aimfile = argv[1]
	filter_set = argv[2]
	filter_set = get_filter(filter_set)
	getfilter(aimfile,filter_set)

if __name__=="__main__":
	main()
