#!/usr/local/bin/python
# Author: Jiayao Wang
# Call Variants From Exome-seq data

from sys import argv
import subprocess as sub
import os
import pandas as pd

def usage():
	print "Create directory and Generate scripts to call variants from Whole Exome sequencing"
	print "Put this script under a 'root' folder to run whole pipeline, for example ../path/to/ExomePipeline/"
	print "./ExomePipeline -h/-help Print Hele Message"
	print "./ExomePipeline -step 0  Create folders to store files needed"
	print "./ExomePipeline -step 1  Run the first part Data pre-processing"
	print "./ExomePipeline -step 2  Run the second part Variant Discovery"

def get_root():
	p = sub.Popen(['pwd'],stdout=sub.PIPE,stderr=sub.PIPE)	
	output, errors = p.communicate()
	return output.strip()

# Make dirctories for Pipeline
def step_0():
	print "Prepare directory for pipeline"
	print "Root directory of Exome pipeline is: %s"%get_root()
	if not os.path.exists('FASTQ'):
		os.makedirs('FASTQ')
		print "FASTQ directory not exists, create it"
	else:
		print "Faseq directory already exists."
	if not os.path.exists('BAM'):
		os.makedirs('BAM')
		print "BAM directory not exists, create it"
	else:
		print "BAM directory already exists."
	if not os.path.exists('Variant'):
		os.makedirs('Variant')
		print "Variant directory not exists, create it"
	else:
		print "Variant directory already exists."
	if not os.path.exists('Script'):
		os.makedirs('Script')
		print "Script directory not exists, create it"
	else:
		print "Script directory already exists."
	
	print "Please move all the fastq file or reads into ./FASTQ folder"
	return

def get_fastq():
 	some_dir="./FASTQ"
	level=1
	some_dir = some_dir.rstrip(os.path.sep)
	assert os.path.isdir(some_dir)
	num_sep = some_dir.count(os.path.sep)
	pair = {}
	for root, dirs, files in os.walk(some_dir):
		num_sep_this = root.count(os.path.sep)
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

def create_bam(root,f1,f2,f_bam):
	f_bam.write('cd '+root+'; \\\n')
	f_bam.write('module load trimmomatic bwa samtools; \\\n')
	
	if f2 != 'NA': #paired end reads
		#Trimming
		f_bam.write('java -jar $TRIMMOJAR PE  -threads 16 \\\n')
		f_bam.write('./FASTQ/'+f1+' \\\n')
		f_bam.write('./FASTQ/'+f2+' \\\n')
		tmp1,tmp2 = f1.split('.')[0] + '_1_Paired.fq',f2.split('.')[0] + '_2_Paired.fq'
		f_bam.write('/lscratch/${SLURM_JOBID}/'+ tmp1 + ' \\\n')
		f_bam.write('/lscratch/${SLURM_JOBID}/'+ f1.split('.')[0]+ '_1_UnPaired.fq \\\n')
		f_bam.write('/lscratch/${SLURM_JOBID}/'+ tmp2 + ' \\\n')
		f_bam.write('/lscratch/${SLURM_JOBID}/'+ f2.split('.')[0]+ '_2_UnPaired.fq \\\n')
		f_bam.write('ILLUMINACLIP:/data/wangj36/genomes/NexteraPE-PE.fa:2:30:10:8:TRUE SLIDINGWINDOW:4:15 MINLEN:36; \\\n')
		#Mapping
		f_bam.write('bwa mem /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa \\\n')
		f_bam.write('-R "@RG\\tID:pairend\\tSM:'+ tmp1.strip('_Paired.fq') + '\\tPL:ILLUMINA\"' ' \\\n')
		f_bam.write('-t 24 /lscratch/${SLURM_JOBID}/' + tmp1 + ' \\\n')
		f_bam.write('/lscratch/${SLURM_JOBID}/' + tmp2 + ' \\\n')
		f_bam.write(' > /lscratch/${SLURM_JOBID}/' + tmp1.strip('_Paired.fq') + tmp2.strip('_Paired.fq') + '.sam; \\\n')
		f_bam.write('samtools view -b -S /lscratch/${SLURM_JOBID}/'+ tmp1.strip('_Paired.fq') + tmp2.strip('_Paired.fq') + '.sam \\\n')
		f_bam.write(' > ./BAM/' + tmp1.strip('_Paired.fq') + tmp2.strip('_Paired.fq') + '_trimmed.bam;\n')
		return tmp1.strip('_Paired.fq') + tmp2.strip('_Paired.fq') + '_trimmed.bam'

	else: # Single end reads
		f_bam.write('./FASTQ/'+f1+' \\\n')
		tmp1,tmp2 = f1.split('.')[0] + '_trimmed.fq','NA'
		f_bam.write('/lscratch/${SLURM_JOBID}/'+ tmp1 + ' \\\n')
		f_bam.write('SLIDINGWINDOW:4:15 MINLEN:36; \\\n')
		#Mapping
		f_bam.write('bwa mem /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa')
		f_bam.write('-R "@RG\\tID:pairend\\tSM:'+f1.split('.')[0]+'\\tPL:ILLUMINA\" \\\n')
		f_bam.write('-t 24 ','/lscratch/${SLURM_JOBID}/'+f1.split('.')[0]+'_trimmed.fq \\\n')
		f_bam.write(' > ','/lscratch/${SLURM_JOBID}/'+f1.split('.')[0]+'.sam; \\\n')
		f_bam.write('samtools view -b -S /lscratch/${SLURM_JOBID}/'+f1.split('.')[0]+'.sam \\\n')
		f_bam.write(' > ','./BAM/'+tmp1+';\n')
		return f1.split('.')[0] + '_trimmed.bam'

def create_dedup(root,bam,f_dedup):
	f_dedup.write('cd ' + root + '/BAM; \\\n')
	f_dedup.write('module load samtools/0.1.19; \\\n')
	sort_bam = bam.strip('_trimmed.bam') + '_sorted.bam'
	f_dedup.write('samtools sort ./'+bam+' /lscratch/${SLURM_JOBID}/' + sort_bam + '; \\\n')
	dedup_bam = sort_bam.strip('_sorted.bam') + '_DeDup.bam'
	f_dedup.write('samtools rmdup -S /lscratch/${SLURM_JOBID}/' + sort_bam + ' ./' + dedup_bam + '; \\\n' )
	f_dedup.write('samtools index ./' + dedup_bam + '; \\\n')
	f_dedup.write('rm ./' + bam + ';\n')
	
	return dedup_bam

def create_recal(root,dedup,f_recal):
	# Recal
	f_recal.write('cd '+root+'/BAM; \\\nmodule load GATK; \\\n')
	f_recal.write('java -Xmx24g -Djava.io.tmpdir=/lscratch/${SLURM_JOBID} -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
	f_recal.write('-R /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \\\n')
	f_recal.write('-T RealignerTargetCreator -nt 24 -I ./' + dedup + ' \\\n')
	list_name = dedup.strip('_DeDup.bam') + '_DeDup_target_intervals.list' # interval list
	f_recal.write('-o /lscratch/${SLURM_JOBID}/'+ list_name +'; \\\n')
	
	f_recal.write('java -Xmx24g -Djava.io.tmpdir=/lscratch/${SLURM_JOBID} -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
	f_recal.write('-R /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \\\n')
	f_recal.write('-T IndelRealigner -targetIntervals /lscratch/${SLURM_JOBID}/' + list_name + ' \\\n')
	f_recal.write('-I ./' + dedup + ' \\\n')
	realign = dedup.strip('_DeDup.bam') + '_DeDup_realigned.bam'
	f_recal.write('-o /lscratch/${SLURM_JOBID}/' + realign + '; \\\n')

	f_recal.write('java -Xmx24g -Djava.io.tmpdir=/lscratch/${SLURM_JOBID} -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
	f_recal.write('-R /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \\\n')
	f_recal.write('-T BaseRecalibrator -nct 8 -knownSites /fdb/GATK_resource_bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz \\\n')
	f_recal.write('-I /lscratch/${SLURM_JOBID}/'+ realign + ' \\\n')
	recal_grp = realign.strip('_realigned.bam') + '_recal.grp'
	f_recal.write('-o /lscratch/${SLURM_JOBID}/' + recal_grp + '; \\\n')

	f_recal.write('java -Xmx24g -Djava.io.tmpdir=/lscratch/${SLURM_JOBID} -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
	f_recal.write('-R /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \\\n')
	f_recal.write('-T PrintReads -nct 8 \\\n')
	f_recal.write('-I /lscratch/${SLURM_JOBID}/' + realign + ' \\\n')
	f_recal.write('-BQSR /lscratch/${SLURM_JOBID}/' + recal_grp +' \\\n')
	recal = realign.strip('_realigned.bam') + '_recal.bam'
	f_recal.write('-o ./' + recal + '; \\\n') 
	
	f_recal.write('rm ./'+dedup+'; \\\n')
	f_recal.write('rm ./'+dedup+'.bai;\n')

def step_1():
	print "Generate data preprocessing scripts.\nPlease make sure fastq file are putted under FASTQ directory."
	fastqs = get_fastq() # get all the name of fastq files as tuple
	f_bam = open('./Script/Bam.swarm','wb')
	f_dedup = open('./Script/Dedup.swarm','wb')
	f_recal = open('./Script/Recal.swarm','wb')
	make_meta(fastqs)
	fastqs = pd.read_csv('meta_data_meta.csv',sep=',')
	root = get_root() 
	for row in fastqs.itertuples():
		f1,f2 = row[1],row[2]
		bam = create_bam(root,f1,f2,f_bam) #Create bam file
		dedup = create_dedup(root,bam,f_dedup)
		recal = create_recal(root,dedup,f_recal)
	f_part1 = open('./Script/ExomePipe_part1.sh','wb')
	f_part1.write('#!/bin/bash\n')
	f_part1.write('cd '+root+'/Script/;\n')
	f_part1.write('module load trimmomatic bwa samtools GATK; module list 2> toolVersions.txt;\n')
	f_part1.write('BamID=`swarm -f Bam.swarm -g 60 -t 24 --gres=lscratch:256`\n')
	f_part1.write('DedupID=`swarm -f Dedup.swarm -g 40 --gres=lscratch:128 --dependency afterany:$BamID`\n')
	f_part1.write('RecalID=`swarm -f Recal.swarm -g 40 -t 24 --gres=lscratch:128 --dependency afterany:$RmDupID`\n')
	f_bam.close()
	f_dedup.close()
	f_recal.close()
	f_part1.close()
	return

def get_bam(root):
	some_dir="./BAM"
	level=1
	some_dir = some_dir.rstrip(os.path.sep)
	assert os.path.isdir(some_dir)
	num_sep = some_dir.count(os.path.sep)
	#bamfiles = []
	"""
	for root, dirs, files in os.walk(some_dir):
		#yield root, dirs, files
		num_sep_this = root.count(os.path.sep)
		bamfiles = []
		count = 0
		print len(files)
		for my_file in files:
			if my_file.endswith('.bam'):
				print my_file
				count += 1
				bamfiles.append(str(my_file))
		if num_sep + level <= num_sep_this:
			del dirs[:]
	"""
	bamfiles = []
	count = 0
	files = [f for f in os.listdir(root+'/BAM') if os.path.isfile(os.path.join(root+'/BAM', f))]
	for my_file in files:
		if my_file.endswith('.bam'):
			count += 1
			bamfiles.append(my_file)
	print len(bamfiles),"BAM FILES"
	#bamfiles.sort(reverse=False)
	#print bamfiles

	return bamfiles

def var_call(root,bams,f_varcall):
	#Create chunks to split task of variants call
	chr_names = ['ch'+ str(i) for i in range(1,23)]
	chr_names.extend(['chrX','chrY'])
	#print chr_names
	chr_sizes = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566]
	chr_breaks = [10, 10, 8, 8, 8, 7, 7, 6, 6, 6, 6, 6, 5, 5, 4, 4, 4, 3, 3, 2, 2, 2, 7, 3]
	chunks = [None]*sum(chr_breaks)
	current_chunk = 0
	for i in range(0,len(chr_sizes)):
		chunk_length = int(round(chr_sizes[i]/chr_breaks[i]))
		num_chr_chunks = chr_breaks[i]
		chunk_start = 1
		for j in range(1,num_chr_chunks+1):
			if j < num_chr_chunks:
				chunk_end = chunk_start + chunk_length - 1
				chunks[current_chunk] = chr_names[i]+':'+str(chunk_start)+'-'+str(chunk_end)
				chunk_start = chunk_end + 1
				current_chunk += 1
			else:
				chunk_end = chr_sizes[i]
				chunks[current_chunk] = chr_names[i]+':'+str(chunk_start)+'-'+str(chunk_end)
				current_chunk += 1
	#print chunks
	for i,chunk in enumerate(chunks):
		f_varcall.write('cd '+root+'/BAM ;\\\n')
		f_varcall.write('module load GATK; \\\n')
		f_varcall.write('java -Xmx64g -Djava.io.tmpdir=/lscratch/${SLURM_JOBID} -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
		f_varcall.write('-R /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \\\n')
		f_varcall.write('-T HaplotypeCaller -nct 64 -A Coverage -A QualByDepth \\\n')
		f_varcall.write('-A FisherStrand -A MappingQualityRankSumTest -A ReadPosRankSumTest -minPruning 3 \\\n')
		for bam in bams:
			f_varcall.write('-I ./' + bam + ' \\\n')
		f_varcall.write('-L '+chunk+' \\\n')
		f_varcall.write('-o ../Variant/Variant_Chunk'+str(i+1)+'.vcf;\n')
	return chunks

def merge_var(root,raw_var,f_mergevar):
	f_mergevar.write('cd '+root+'/Variant ;\\\n')
	f_mergevar.write('module load GATK; \\\n')
	f_mergevar.write('java -Xmx56g -cp $GATK_HOME/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \\\n')
	f_mergevar.write('-R /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \\\n')
	for i in range(len(raw_var)):
		f_mergevar.write('-V ./Variant_Chunk'+str(i+1)+'.vcf \\\n')
	f_mergevar.write('-out GATK_All_Variants.vcf \\\n')
	f_mergevar.write('-assumeSorted;\n')

def SNP_recal(root,var_whole,f_snp):
	f_snp.write('cd ' + root + '\Variant; \\\n')
	f_snp.write('module load GATK; \\\n')
	f_snp.write('java -Xmx64g -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
	f_snp.write('-T VariantRecalibrator \\\n')
	f_snp.write('-R /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \\\n')
	f_snp.write('-input ./GATK_All_Variants.vcf \\\n')
	## Set Training need consideration
	f_snp.write('-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /fdb/GATK_resource_bundle/hg19/hapmap_3.3.hg19.vcf.gz \\\n')
	f_snp.write('-resource:omni,known=false,training=true,truth=true,prior=12.0 /fdb/GATK_resource_bundle/hg19/1000G_omni2.5.hg19.vcf.gz \\\n')
	f_snp.write('-resource:1000G,known=false,training=true,truth=false,prior=10.0 /fdb/GATK_resource_bundle/hg19/1000G_phase1.snps.high_confidence.hg19.vcf.gz \\\n')
	f_snp.write('-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /fdb/GATK_resource_bundle/hg19/dbsnp_138.hg19.vcf.gz \\\n')
	f_snp.write('-an QD \\\n-an MQ \\\n-an MQRankSum \\\n-an ReadPosRankSum \\\n-an FS \\\n-an SOR \\\n-an DP \\\n')
	f_snp.write('-mode SNP \\\n-tranche 100.0 \\\n-tranche 99.9 \\\n-tranche 99.0 \\\n-tranche 90.0 \\\n')
	f_snp.write('-recalFile ./recalibrate_SNP.recal \\\n')
	f_snp.write('-tranchesFile ./recalibrate_SNP.traches \\\n')
	f_snp.write('-rscriptFile ./recalibrate_SNP_plots.R; \\\n')

	f_snp.write('java -Xmx24g -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
	f_snp.write('-T ApplyRecalibration \\\n')
	f_snp.write('-R /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \\\n')
	f_snp.write('-input ./GATK_All_Variants.vcf  \\\n')
	f_snp.write('-mode SNP \\\n')
	f_snp.write('-ts_filter_level 99.5 \\\n')
	f_snp.write('-recalFile ./recalibrate_SNP.recal \\\n')
	f_snp.write('-tranchesFile ./recalibrate_SNP.traches \\\n')
	f_snp.write('-o recalibrated_SNPS.vcf;\n')

	return "recalibrated_SNPS.vcf"


def INDEL_recal(root,var_SNP,f_indel):
	f_indel.write('cd ' + root + '\Variant; \\\n')
	f_indel.write('module load GATK; \\\n')
	f_indel.write('java -Xmx64g -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
	f_indel.write('-T VariantRecalibrator \\\n')
	f_indel.write('-R /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \\\n')
	f_indel.write('-input ./recalibrated_SNPS.vcf \\\n')
	## Set Training need consideration
	f_indel.write('-resource:mills,known=true,training=true,truth=true,prior=12.0 /fdb/GATK_resource_bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz \\\n')
	f_indel.write('-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /fdb/GATK_resource_bundle/hg19/dbsnp_138.hg19.vcf.gz \\\n')
	f_indel.write('-an QD \\\n-an MQ \\\n-an MQRankSum \\\n-an ReadPosRankSum \\\n-an FS \\\n-an SOR \\\n-an DP \\\n')
	f_indel.write('-mode INDEL \\\n-tranche 100.0 \\\n-tranche 99.9 \\\n-tranche 99.0 \\\n-tranche 90.0 \\\n--maxGaussians 4 \\\n')
	f_indel.write('-recalFile ./recalibrate_SNP_INDEL.recal \\\n')
	f_indel.write('-tranchesFile ./recalibrate_SNP_INDEL.traches \\\n')
	f_indel.write('-rscriptFile ./recalibrate_SNP_INDEL_plots.R; \n')

	f_indel.write('java -Xmx24g -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
	f_indel.write('-T ApplyRecalibration \\\n')
	f_indel.write('-R /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \\\n')
	f_indel.write('-input ./recalibrated_SNPS.vcf  \\\n')
	f_indel.write('-mode INDEL \\\n')
	f_indel.write('-ts_filter_level 99.0 \\\n')
	f_indel.write('-recalFile ./recalibrate_SNP_INDEL.recal \\\n')
	f_indel.write('-tranchesFile ./recalibrate_SNP_INDEL.traches \\\n')
	f_indel.write('-o recalibrated_SNP_INDEL.vcf;\n')

	return "recalibrated_SNP_INDEL.vcf"

def Coverage(root,bams,f_coverage,bed_file):
	f_coverage.write('#!/bin/bash\n')
	f_coverage.write('cd '+root+'/Script\n')
	f_coverage.write('CoverID=`swarm -f coverage.swarm -g 60 -t 8 --partition quick`\n')
	f_coverage.write('ConcatID=`swarm -f concat_coverage.swarm -g 60 -t 8 --partition quick --dependency afterany:$CoverID`')
	f_coverID = open('./Script/coverage.swarm','wb')
	f_coverID.write('# Calculate Coverage and Depth of Bedfile')
	for bam in bams:
		f_coverID.write('module load bedtools; \\\ncd '+root+'; \\\n')
		f_coverID.write('bedtools bamtobed -i ./BAM/'+bam+' | bedtools coverage -a '+bed_file+' -b - > ./BAM/'+bam.strip('.bam')+'.coverage; \\\n')
		f_coverID.write('cd '+root+'/Script; \\\n')
		f_coverID.write('python Calculate_Bam_coverage.py ../BAM/'+bam.strip('.bam')+'.coverage;\n' )
	f_concatID = open('./Script/concat_coverage.swarm','wb')
	f_concatID.write('cd '+root+'/BAM; \\\n')
	f_concatID.write('echo \'Bam_File\tCovered_Region\tTotal_Region\tRegion_covered_Ratio\tTotal_Reads\tCovered_Depth\tTotal_Depth\' | cat - \\\n')
	for bam in bams:
		f_concatID.write(bam.strip('.bam')+'.coverage.tmp \\\n')
	#f_concatID.write('*.coverage.tmp')
	f_concatID.write('> Coverage.tsv')


def step_2():
	print "Generate variants call scripts.\nPlease make sure bam file are putted under BAM directory."
	root = get_root()
	bams = get_bam(root) # get all the name of bam files
	
	f_coverage = open('./Script/Coverage.sh','wb')
	f_varcall = open('./Script/VarCall.swarm','wb')
	f_mergevar = open('./Script/VarMerge.swarm','wb')
	f_snp = open('./Script/RunVsqrSNP.swarm','wb')
	f_indel = open('./Script/RunVsqrINDEL.swarm','wb')
	bed_file = '/data/wangj36/common_used_files/OverlappingExome.bed'
	Coverage(root,bams,f_coverage,bed_file)
	raw_var = var_call(root,bams,f_varcall)
	var_whole = merge_var(root,raw_var,f_mergevar)
	var_SNP = SNP_recal(root,var_whole,f_snp)
	var_INDEL = INDEL_recal(root,var_SNP,f_indel)

	
	f_part2 = open('./Script/ExomePipe_part2.sh','wb')
	f_part2.write('#!/bin/bash\n')
	f_part2.write('cd '+root+'/Script/;\n')
	f_part2.write('VarCallID=`swarm -f VarCall.swarm -g 72 -t 24 --gres=lscratch:256 --time 24:00:00`\n')
	f_part2.write('MergeVariantsID=`swarm -f MergeVariants.swarm -g 72 --dependency afterany:$VarCallID`\n')
	f_part2.write('RunVsqrSNPID=`swarm -f RunVsqrSNP.swarm -g 72 -t 24 --dependency afterany:$MergeVariantsID`\n')
	f_part2.write('RunVsqrINDELID=`swarm -f RunVsqrINDEL.swarm -g 72 -t 24 --dependency afterany:$RunVsqrSNPID`')
	f_part2.close()

	f_varcall.close()
	f_mergevar.close()
	f_snp.close()
	f_indel.close()
	return

def take_argv():
	step = 0
	i = 1
	while i < len(argv):
		if argv[i] == '-s':
			step = argv[i+1]
		i += 2
	return step

def main():
	step = take_argv()
	if step == '0':
		step_0()
	if step == '1':
		step_1()
	if step == '2':
		step_2()
	
if __name__=="__main__":
	main()

