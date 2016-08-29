from sys import argv
import commands
import os

def count_size():
	fin=open("filelist.txt",'r')
	fin.readline()

	count = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	cutoff = [0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000]
	for l in fin:
		size = l.split()[4]
		_size = size[:-1]
		_size = int(_size)
		if _size < cutoff[1]:
			count[0] += 1
		elif _size > cutoff[1] and _size < cutoff[2]:
			count[1] += 1
		elif _size > cutoff[2] and _size < cutoff[3]:
			count[2] += 1
		elif _size > cutoff[3] and _size < cutoff[4]:
			count[3] += 1
		elif _size > cutoff[4] and _size < cutoff[5]:
			count[4] += 1
		elif _size > cutoff[5] and _size < cutoff[6]:
			count[5] += 1
		elif _size > cutoff[6] and _size < cutoff[7]:
			count[6] += 1
		elif _size > cutoff[7] and _size < cutoff[8]:
			count[7] += 1
		elif _size > cutoff[8] and _size < cutoff[9]:
			count[8] += 1
		elif _size > cutoff[9] and _size < cutoff[10]:
			count[9] += 1
		elif _size > cutoff[10] and _size < cutoff[11]:
			count[10] += 1
		elif _size > cutoff[11] and _size < cutoff[12]:
			count[11] += 1
		elif _size > cutoff[12] and _size < cutoff[13]:
			count[12] += 1
		elif _size > cutoff[13] and _size < cutoff[14]:
			count[13] += 1
		elif _size > cutoff[14]:
			count[14] += 1
	for i in range(0,len(count[0:-1])):
		print "number of files less than",cutoff[i+1],":",count[i]
	print "number of files more than 15000",count[-1]

def filter_file():
	os.system('mkdir left_out')
	fin=open("filelist.txt",'r')
	fin.readline()
	for l in fin:
		info = l.strip().split()
		size = int(info[4][:-1])
		file_name = info[8]
		if size < 1000 or size > 8000:
			os.system('mv FASTQ/' + file_name + ' left_out') 


def main():
	count_size()
	filter_file()

if __name__=="__main__":
	main()
