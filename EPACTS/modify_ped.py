from sys import argv
import csv

def usage():
	print "This Script is used to modify the ped file for association test."
	print "python modify_ped.py ped_file_from_vcftools meta_sheet"
	print "Two file needed as input. "
	print "Filtered metadata file contains individual file name as cases"
	print "Seems not necessary need 1000G controls list is needed"
	print "outliers.txt will contains individuals that we have the bam file but not have disease we interested. Further used to filter individual in .vcf file"

def add_line_remove_genotype(head,infile):
	hand = open(argv[1],'rb')
	print hand
	out_name = infile.split(".")[0]+"_with_head.ped"
	fout = open(out_name,'wb')
	print fout
	fout.write(head)
	for l in hand:
		fout.write("\t".join(l.split()[:6])+"\n")
	hand.close()
	fout.close()

def make_pair_fam():
	meta_sheet = open("/data/wangj36/common_used_files/AMR_metasheet.csv",'rb')
	data = csv.reader(meta_sheet, delimiter=',')
	# modified.ID_2:[Affection, FamID, FatherID, MotherID, SKAT] 
	#meta_sheet.readline()
	pair = {}
	for l in data:
		#record = l.strip().split(",")
		record = l
		print record[25],record[26]
		if ((record[25] == 'AMR') and (record[26] == 'YES') and record[27] == 'IRD'):
			if record[13] == '0': record[13] = record[4]
			if record[14] == '':record[14] = '0'
			if record[15] == '':record[15] = '0'
			pair[record[4]] = [record[8],record[13],record[14],record[15]]
			print [record[8],record[13],record[14],record[15]]
			if pair[record[4]][1] == '0' or pair[record[4]][1] == '1': # If FamID eq 0, FamID eq individual_ID (just this patient)
				pair[record[4]][1] == record[4]
	#print pair
	meta_sheet.close()
	return pair


def make_pair_ind(meta_sheet):
	#meta_sheet = open("/data/wangj36/common_used_files/AMR_metasheet.csv",'rb')
	data = csv.reader(meta_sheet, delimiter=',')
	# modified.ID_2:[Affection, FamID, FatherID, MotherID, SKAT] 
	pair = {}
	for l in data:
		record = l
		if ((record[25] == 'AMR') and (record[26] == 'YES') and record[27] == 'IRD'):
			record[13] = record[4]
			record[14] = '0'
			record[15] = '0'
			pair[record[4]] = [record[8],record[13],record[14],record[15]]
			print [record[8],record[13],record[14],record[15]]
			if pair[record[4]][1] == '0' or pair[record[4]][1] == '1': # If FamID eq 0, FamID eq individual_ID (just this patient)
				pair[record[4]][1] == record[4]
	#print pair
	meta_sheet.close()
	return pair

def control_list():
	control = open('1000_GenomeProjectPhase3_SampleDetails.txt','rb')
	control.readline()
	res = {}
	for l in control:
		res[l.strip().split()[0]] = 1
	control.close()
	return res

def modify_ped(head,infile,case,control):
	
	hand = open(argv[1],'rb')
	out_name = infile.split(".")[0]+"_with_head.ped"
	fout = open(out_name,'wb')
	#This file contains individuals that we have the bam file but not have disease we interested.
	fout2 = open("outliers.txt",'wb')
	#fout2.write("#This file contains individuals that we have the bam file but not have disease we interested.")
	fout.write(head)
	for l in hand:
		record = l.split()[:6]
		if record[1] in case: #Case
			record[0] = case[record[1]][1] #FamID
			record[1] = record[1] #Individual ID
			record[2] = case[record[1]][2] #FatherID
			record[3] = case[record[1]][3] #MotherID
			record[5] = case[record[1]][0] #Affection
			fout.write("\t".join(record)+"\n")
		#elif record[1] in control:
		elif not record[1].startswith('eG'): #Control	
			record[0] = record[1]
			record[5] = "1"
			fout.write("\t".join(record)+"\n")
		else: #individual not in case and control
			fout2.write(record[1]+'\n')
			
	hand.close()
	fout.close()
	fout2.close()

def main():
	head = "#FAM_ID	IND_ID	FAT_ID	MOT_ID	SEX	AFFECTION\n"
	meta_sheet = argv[2]
	case = make_pair_ind(meta_sheet)
	#control = control_list()
	control = None
	modify_ped(head,argv[1],case,control)

if __name__=="__main__":
	main()
