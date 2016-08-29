#make the meta_data_meta.csv file according to FASTQ FOLDER
import os

def walk(some_dir,level=1):
	some_dir = some_dir.rstrip(os.path.sep)
	assert os.path.isdir(some_dir)
	num_sep = some_dir.count(os.path.sep)
	for root, dirs, files in os.walk(some_dir):
		#print root, dirs, files
		num_sep_this = root.count(os.path.sep)
		#print files
		pair = {}
		for _file in files:
			if _file.split('_')[0] not in pair:
				pair[_file.split('_')[0]] = [_file]
			else:
				pair[_file.split('_')[0]].append(_file)
		if num_sep + level <= num_sep_this:
			del dirs[:]
	return pair

def make_meta(pair):
	fout = open("meta_data_meta.csv",'w')
	fout.write("R1,R2\n")
	for k,v in pair.items():
		if len(v) == 2:
			v[0],v[1] = sorted([v[0],v[1]])
			fout.write(v[0]+","+v[1]+"\n")
		elif len(v) == 1:
			fout.write(v[0]+","+"NA"+"\n")
		else:
			print "Unexpected Error pair"

def main():
	some_dir="./FASTQ"
	pairs = walk(some_dir,1)
	make_meta(pairs)	

if __name__=="__main__":
	main()
