#!/local/usr/bin/python
import pandas as pd 
from sys import argv

# Calculate total coverage from coverage file from Bedtools
def calculate_bam_coverage(coverage_file):
	df = pd.read_csv(coverage_file,sep='\t',header=None,names=['CHR','START_POS','END_POS','MARKER','NUM_READS','REGION_COVERED','REGION_TOTAL','COVERED_RATIO'])
	Total = list(df['REGION_TOTAL'])
	Cover = list(df['REGION_COVERED'])
	Cover,Total = sum(Cover),sum(Total)
	Total_Ratio = float(Cover)/Total
	Reads_Num = sum(list(df['NUM_READS']))
	Coverage_Cover = float(Reads_Num)*125/Cover
	Coverage_Total = float(Reads_Num)*125/Total
	out_name = coverage_file + '.tmp'
	res = '\t'.join(map(str,[coverage_file.strip('.coverage')+'.bam',Cover,Total,Total_Ratio,Reads_Num,Coverage_Cover,Coverage_Total]))
	#open(out_name,'wb').write('Covered_Region\tTotal_Region\tRegion_covered_Ratio\tTotal_Reads\tCovered_Depth\tTotal_Depth\n'+res+'\n')
	open(out_name,'wb').write(res+'\n')


def main():
	coverage_file = argv[1]
	calculate_bam_coverage(coverage_file)

if __name__=='__main__':
	main()
