import sys

#get all Exome samples and thier runs
def get_sample(hand):
	hand.readline()
	samples = {}
	for l in hand:
		cells = l.strip().split("\t")
		_filename = cells[0]
		_pairfile = cells[19]
		_sample = cells[9]
		_run = cells[2]
		_reads = cells[23]
		_type = cells[25]
		if _type != "exome":
			continue
		else:
			record = [_run,_filename,_pairfile,_reads,_type]
			if _sample not in samples:
				samples[_sample] = [record]
			else:
				samples[_sample].append(record)
	return samples

def check_data_in_dic(samples):
	for k,v in samples.items()[0:100]:
		print k
		print v

#get best run from all runs from same sample
def get_best_run(records):
	runs = {}
	#group records by run
	#record:[run, filename, pairfile, reads, type]
	for record in records:
		if record[0] not in runs:
			runs[record[0]] = [record[1:]]
		else:
			runs[record[0]].append(record[1:])

	#leave out the "fragments where the other end failed qc"
	for k,v in runs.items():
		if len(v) <2:
			print "Unecpected Situation, maybe a single end file"
			del runs[k]
			continue
		files = []
		for _file in v:
			name = _file[0].split("/")[-1]
			if "_" in name:
				files.append(_file)
		#print files
		files.append(int(files[0][2])+int(files[1][2])) #This is the sum of reads number of two pair-end sequence file
		runs[k] = files

	best_run = None
	best_reads = 0
	for k,v in runs.items():
		print v
		if v[2] >= best_reads:
			best_run = [k,v[0],v[1]]
			best_reads = v[2]
	return best_run
		

def sample2run(samples):
	runs = {}
	for k,v in samples.items():
		#print v
		runs[k] = get_best_run(v)
	runs = runs.items()
	res = []
	for run in runs:
		#print run
		if run[1] == None:
			continue
			#print run[0]
		(_sample,[_run,[_filename1,_pairfile1,_reads,_type],[_filename2,_pairfile2,_reads,_type]]) = run
		res.append([_filename1,_pairfile1,_sample,_run,_reads,_type])
		res.append([_filename2,_pairfile2,_sample,_run,_reads,_type])
	return res

def step_1():
	indexFile = "20130502.phase3.analysis.sequence.index"
	hand = open(indexFile,'r')
	samples = get_sample(hand)
	#check_data_in_dic(samples)
	runs = sample2run(samples)
	#check_data_in_runs(runs)
	print "Writing new index ..."
	#for i in runs[0:100]:
		#print i
	with open("./Filterd_file.index",'w') as fout:
		for _filename,_pairfile,_sample,_run,_reads,_type in runs:
			fout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(_filename,_pairfile,_sample,_run,_reads,_type))
	fout.close()
	hand.close()

#leave out samples that don't include in "1000_GenomeProjectPhase3_SampleDetails.txt"
#Join file with groups.
def step2():
	hand = open("1000_GenomeProjectPhase3_SampleDetails.txt",'r')
	hand.readline()
	names = {}
	for l in hand:
		info = l.strip().split("\t")
		names[l.split("\t")[0]] = [info[1],info[2],info[3]]
	hand.close()

	fout = open("Filterd_with_POP.index",'w')
	fin = open("Filterd_file.index",'r')
	for l in fin:
		key = l.split("\t")[2]
		if key in names:
			new_line = l.strip() + "\t" + names[key][0] + "\t"+ names[key][1] + "\t" + names[key][2] + "\n" 
			fout.write(new_line)
	fin.close()
	fout.close()


def main():
	#step1()
	step2()

if __name__=="__main__":
	main()
