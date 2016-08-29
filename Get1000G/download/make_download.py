import sys
from sys import argv
import os

def usage():
	print "python make_download.py --s supur_population --p population"
	print "supur_population is required"
	print "Two file must place under same folder with it."
	print "1000_GenomeProjectPhase3_SampleDetails.txt with information of individual and their population"
	print "Filterd_file.index with information of individual and which file to download and other informations"
	exit()

def take_argv():
	
	super_pop = ""
	pop = ""
	i = 1 
	while i < len(argv):
		if argv[i] == "--s":
			super_pop = argv[i+1]
		if argv[i] == "--p":
			pop = argv[i+1]

def make_ind_pop(hand):
	hand.seek(0)
	hand.readline()
	ind_pop = {}
	for l in hand:
		A = l.split("\t")
		ind_pop[A[0]] = (A[1],A[2])
	return ind_pop

def make_ind_file(hand):
	hand.seek(0)
	hand.readline()
	ind_file = {}
	for l in hand:
		A = l.split("\t")
		if A[2] in ind_file:
			continue
		else:
			ind_file[A[2]] = (A[0],A[1])
	
	return ind_file



def downlist(hand,pop):
	hand.seek(0)
	hand.readline()
	dic = {}
	for l in hand:
		split = l.split()
		if split[2] == pop:
			if split[1] in dic:
				dic[split[1]].append(split[0])
			else:
				dic[split[1]] = [split[0]]
	"""
	for k in dic.keys():
		print k
		print dic[k]
	"""
	return dic

def make_download_sh(downloads,ind_file):
	for k in downloads.keys():
		fout = open(str(k)+"_download.sh",'w')
		for indi in downloads[k]:
			if indi in ind_file:
				addr1 = ind_file[indi][0]
				addr2 = ind_file[indi][1]
				fout.write("wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/" + addr1 + " -P /data/wangj36/1000genomes_phase3_exome/" + argv[1] + "/" + k + "\n")
				fout.write("wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/" + addr2 + " -P /data/wangj36/1000genomes_phase3_exome/" + argv[1] + "/" + k + "\n")
			else:
				print "ERROR!",indi,"not in",k,"in address file"
		fout.close()
		#os.system('mv '+str(k)+'_download.sh '+'/data/wangj36/1000genomes_phase3_exome/' + argv[1] + '/' + k)


def main():
	#take_argv()
	hand1 = open("1000_GenomeProjectPhase3_SampleDetails.txt",'r')
	hand2 = open("Filterd_file.index",'r')
	downloads = downlist(hand1,argv[1])
	
	
	#ind_pop	= make_ind_pop(hand1)
	ind_file = make_ind_file(hand2)
	make_download_sh(downloads,ind_file)


if __name__=="__main__":
	main()
