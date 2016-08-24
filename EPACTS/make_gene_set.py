#!/usr/local/bin/python
# make gene set from a list of genes and randomly selected genes exclude genes from that set

from sys import argv
import re
from random import shuffle

def usage():
	print "python make_geneset.py -i group_file.grp -l gene_list.txt -e exclude_list.txt -o geneset.grp"
	exit()

#take argv needed in this script
def take_arg(argv):
	gene_grp = "" # variants grouped by genes
	gene_list = "" # A list of genes will grouped into interested set
	gene_exclude = "" # genes will exclude from the gene_list
	geneset = "" # output gene set file name

	i = 1
	while i<len(argv):
		if argv[i] == '-i':
			gene_grp = argv[i+1]
		if argv[i] == '-l':
			gene_list = argv[i+1]
		if argv[i] == '-o':
			geneset = argv[i+1]
		i += 2
	if gene_grp == "" or gene_list == "" or geneset == "":
		usage()
	return gene_grp,gene_list,gene_exclude,geneset

def read_grp(gene_grp):
	hand = open(gene_grp,'rb')
	res = {}
	for l in hand:
		data = l.strip().split('\t')
		gene = data[0]
		variants = data[1:]
		res[gene] = variants
	return res

# Pick genes in retnet set form a new gene: variants dic, and remove those from gene_grp
# INPUT: big group
def make_aim(gene_grp,gene_list,exclude):
	aim_set = {}
	exclude_set = {}
	#Read the exclude gene as dictionary
	if exclude != "":
		hand = open(exclude,'rb')
		for l in hand:
			exclude_set[l.strip()] = 1
	#Read the gene list exclude the exclude_set as dictionary of gene:variants
	hand = open(gene_list,'rb')
	for l in hand:
		gene = l.strip()
		if (gene in gene_grp) and (gene not in exclude_set):
			aim_set[gene] = gene_grp[gene]
			gene_grp.pop(gene, None)
	return gene_grp,aim_set

def sort_by_g(res):
	# res: [setid,[var1,var2,var3]]
	# vari: chrNum:pos_Base1/Base2
	# sort first on chrNum and then on pos
	var_list = res[1]
	setid = res[0]
	tmp = []
	for item in var_list:
		chrNum,other = item.split(':')[0],item.split(':')[1]
		pos,other = other.split('_')[0],other.split('_')[1]
		tmp.append([chrNum,pos,other])
	var_list = sorted(tmp,key=lambda x:[int(x[0]),int(x[1])])
	# reconstruct var format
	res = []
	for var in var_list:
		res.append(var[0] + ':' + var[1] + '_' + var[2])
	return (setid,res)
		

def genes2set(setid,tmp,fout):
	#tmp:[ [gene1,variants], [gene2,variants], [gene3,variants] ... ,[genen,variants]]
	res = [setid,[]]
	# add vairants from each gene into set
	for item in tmp:
		fout.write(item[0]+'\t')
		res[1].extend(item[1])
	# sort the variants according to increasing order of genomic coordinate
	#print res
	res = sort_by_g(res)
	# res: setid,[var1,var2,...,varN]
	fout.write('\n')
	return res

def make_gene_sets(gene_grp,set_aim,geneset):
	window = len(set_aim.keys())
	fout = open('Set2Gene.txt','wb') #write a file of randomly selected genes
	fout.write('Aim Genes'+"\n")
	geneset = open(geneset,'wb')
	tmp = genes2set('RetnetGenes',set_aim.items(),fout)
	geneset.write(tmp[0]+'\t'+'\t'.join(tmp[1])+'\n')	
	gene_grp_tuple = gene_grp.items() #Gene_grp from dic to tuple list
	shuffle(gene_grp_tuple) # In place shuffle the list
	set_random = []
	for i in xrange(0,10):
		tmp = gene_grp_tuple[i*window:i*window+window]
		#print tmp
		setid = "RandomSet_"+str(i)
		fout.write(setid+"\n")
		tmp = genes2set(setid,tmp,fout) # setid,[var1...varN]
		set_random.append(tmp)
		#print tmp
		geneset.write(tmp[0]+'\t'+'\t'.join(tmp[1])+'\n')
	geneset.close()
	fout.close()
		

def main():
	gene_grp,gene_list,exclude,geneset = take_arg(argv)
	gene_grp = read_grp(gene_grp)
	gene_grp, set_aim = make_aim(gene_grp,gene_list,exclude)
	make_gene_sets(gene_grp,set_aim,geneset)	
	
if __name__=="__main__":
	main()
