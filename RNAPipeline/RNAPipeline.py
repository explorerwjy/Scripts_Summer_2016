#!/usr/local/bin/python
# Author: Jiayao Wang
# Call Variants From RNA-seq data

from sys import argv
import subprocess as sub
import os
import pandas as pd

def usage():
	print "Create directory and Generate scripts to call variants from rna-seq"
	print "Put this script under a 'root' folder to run whole pipeline, for example ../path/to/RNAPipeline/"
	print "./RNAPipeline -h/-help Print Hele Message"
	print "./RNAPipeline -step 0  Create folders to store files needed"
	print "./RNAPipeline -step 1  Run the first part Data pre-processing"
	print "./RNAPipeline -step 2  Run the second part Variant Discovery"

def get_root():
	p = sub.Popen(['pwd'],stdout=sub.PIPE,stderr=sub.PIPE)	
	output, errors = p.communicate()
	return output.strip()

# Make dirctories for Pipeline
def step_0():
	print "Prepare directory for pipeline"
	print "Root directory of rna pipeline is: %s"%get_root()
	if not os.path.exists('Fastq'):
		os.makedirs('Fastq')
		print "Fastq directory not exists, create it"
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
	
	print "Please move all the fastq file or reads into ./Fastq folder"
	return

def get_fastq():
 	some_dir="./Fastq"
	level=1
	some_dir = some_dir.rstrip(os.path.sep)
	assert os.path.isdir(some_dir)
	num_sep = some_dir.count(os.path.sep)
	for root, dirs, files in os.walk(some_dir):
		num_sep_this = root.count(os.path.sep)
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

def create_bam(root,f1,f2,f_bam):
	f_bam.write('cd '+root+'; \\\n')
	f_bam.write('module load trimmomatic STAR samtools; \\\n')
	#Trimming
	f_bam.write('java -jar $TRIMMOJAR PE  -threads 16 \\\n')
	if f2 != 'NA': #paired end reads
		f_bam.write('./Fastq/'+f1+' \\\n')
		f_bam.write('./Fastq/'+f2+' \\\n')
		tmp1,tmp2 = f1.split('.')[0] + '_1.Paired.fq',f2.split('.')[0] + '_2.Paired.fq'
		f_bam.write('/lscratch/${SLURM_JOBID}/'+ tmp1 + ' \\\n')
		f_bam.write('/lscratch/${SLURM_JOBID}/'+ f1.split('.')[0]+ '_1.UnPaired.fq \\\n')
		f_bam.write('/lscratch/${SLURM_JOBID}/'+ tmp2 + ' \\\n')
		f_bam.write('/lscratch/${SLURM_JOBID}/'+ f2.split('.')[0]+ '_2.UnPaired.fq \\\n')
		f_bam.write('ILLUMINACLIP:/data/wangj36/genomes/TruSeq3-PE-2.fa:2:30:10:1:TRUE SLIDINGWINDOW:4:5 TRAILING:5 MINLEN:25; \\\n')
	else: # Single end reads
		f_bam.write('./Fastq/'+f1+' \\\n')
		tmp1,tmp2 = f1.split('.')[0] + '_1.Paired.fq','NA'
		f_bam.write('/lscratch/${SLURM_JOBID}/'+ tmp1 + ' \\\n')
		f_bam.write('/lscratch/${SLURM_JOBID}/'+ f1.split('.')[0]+ '_1.UnPaired.fq \\\n')
		f_bam.write('ILLUMINACLIP:/data/wangj36/genomes/TruSeq3-PE-2.fa:2:30:10:8:TRUE SLIDINGWINDOW:4:15 MINLEN:36; \\\n')
	#Mapping
	#f_bam.write('find ' + 'BAM' + ' -type d -exec chmod 777 {} \; \\\n')
	f_bam.write('STAR --genomeDir /data/wangj36/genomes/star_indices/star_ensembl_GRCh38.82/ \\\n')
	if tmp2 != 'NA': # Paired end reads
		f_bam.write('--readFilesIn /lscratch/${SLURM_JOBID}/' + tmp1 +' /lscratch/${SLURM_JOBID}/' + tmp2 + ' \\\n')
	else: # Single end reads
		f_bam.write('--readFilesIn /lscratch/${SLURM_JOBID}/' + tmp1 + ' \\\n')
	f_bam.write('--readFilesCommand - --runThreadN 16 --outSAMtype BAM Unsorted --outReadsUnmapped Fastx \\\n')
	f_bam.write('--twopassMode Basic --genomeSAindexNbases 11 --outFilterType BySJout \\\n')
	f_bam.write('--outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \\\n')
	f_bam.write('--outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 \\\n')
	f_bam.write('--alignMatesGapMax 1000000 --quantMode TranscriptomeSAM \\\n')
	f = '_'.join(f1.split('_')[:-1]) + '_'
	f_bam.write('--outFileNamePrefix ./BAM/' + f + ';\n\n')
	return f + 'Aligned.out.bam'

def create_dedup(root,bam,f_dedup):
	f_dedup.write('cd ' + root + '; \\\n')
	f_dedup.write('module load picard; \\\n')
	f_dedup.write('java -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups \\\n')
	f_dedup.write('I=./BAM/'+ bam +' \\\n')
	tmp = bam.strip('.bam') + '.sorted.bam'
	f_dedup.write('O=/lscratch/${SLURM_JOBID}/' + tmp + ' \\\n')
	f_dedup.write('SO=coordinate RGID=id RGLB=library RGPL=illumina RGPU=machine RGSM=sample; \\\n')

	f_dedup.write('java -jar $PICARDJARPATH/picard.jar MarkDuplicates \\\n')
	f_dedup.write('I=/lscratch/${SLURM_JOBID}/' + tmp + ' \\\n')
	tmp = tmp.strip('t.sroted.bam') + 't.dedup.bam'
	f_dedup.write('O=./BAM/' + tmp + ' \\\n')
	f_dedup.write('CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics;\n\n')
	
	return tmp

def create_recal(root,dedup,f_recal):
	# SplitNcigars
	f_recal.write('cd ' + root + '/BAM/; \\\n')
	f_recal.write('module load GATK; \\\n')
	f_recal.write('java -Xmx24g -Djava.io.tmpdir=/lscratch/${SLURM_JOBID} -jar $GATK_JAR \\\n')
	f_recal.write('-T SplitNCigarReads -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS \\\n')
	f_recal.write('-R /data/wangj36/genomes/ensembl/ensembl_GRCh38.82/ensembl_homo_sapiens_GRCh38.82_dna_primary_assembly.fa \\\n')
	f_recal.write('-I ' + dedup + ' \\\n')
	tmp = dedup.strip('bam') + 'splitN.bam' # Filename after SplitNcigars
	f_recal.write('-o /lscratch/${SLURM_JOBID}/'+ tmp + '; \\\n')

	# Recal
	f_recal.write('java -Xmx24g -Djava.io.tmpdir=/lscratch/${SLURM_JOBID} -jar $GATK_JAR \\\n')
	f_recal.write('-T RealignerTargetCreator -nt 24 \\\n')
	f_recal.write('-R /data/wangj36/genomes/ensembl/ensembl_GRCh38.82/ensembl_homo_sapiens_GRCh38.82_dna_primary_assembly.fa \\\n')
	f_recal.write('-I /lscratch/${SLURM_JOBID}/' + tmp + ' \\\n')
	list_name = tmp.strip('.bam') + '_intervals.list' # interval list
	f_recal.write('-o /lscratch/${SLURM_JOBID}/'+ list_name +'; \\\n')
	
	f_recal.write('java -Xmx24g -Djava.io.tmpdir=/lscratch/${SLURM_JOBID} -jar $GATK_JAR \\\n')
	f_recal.write('-R /data/wangj36/genomes/ensembl/ensembl_GRCh38.82/ensembl_homo_sapiens_GRCh38.82_dna_primary_assembly.fa \\\n')
	f_recal.write('-T IndelRealigner -targetIntervals /lscratch/${SLURM_JOBID}/' + list_name + ' \\\n')
	f_recal.write('-I /lscratch/${SLURM_JOBID}/' + tmp + ' \\\n')
	realign = tmp.strip('.bam') + '.realigned.bam'
	f_recal.write('-o /lscratch/${SLURM_JOBID}/' + realign + '; \\\n')

	f_recal.write('java -Xmx24g -Djava.io.tmpdir=/lscratch/${SLURM_JOBID} -jar $GATK_JAR \\\n')
	f_recal.write('-R /data/wangj36/genomes/ensembl/ensembl_GRCh38.82/ensembl_homo_sapiens_GRCh38.82_dna_primary_assembly.fa \\\n')
	f_recal.write('-T BaseRecalibrator -nct 8 \\\n')
	f_recal.write('-knownSites /fdb/GATK_resource_bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz \\\n') # KnowSites
	f_recal.write('-I /lscratch/${SLURM_JOBID}/'+ realign + ' \\\n')
	recal_grp = tmp.strip('.bam') + '.recal.grp'
	f_recal.write('-o /lscratch/${SLURM_JOBID}/' + recal_grp + '; \\\n')

	f_recal.write('java -Xmx24g -Djava.io.tmpdir=/lscratch/${SLURM_JOBID} -jar $GATK_JAR \\\n')
	f_recal.write('-R /data/wangj36/genomes/ensembl/ensembl_GRCh38.82/ensembl_homo_sapiens_GRCh38.82_dna_primary_assembly.fa \\\n')
	f_recal.write('-T PrintReads -nct 8 \\\n')
	f_recal.write('-I /lscratch/${SLURM_JOBID}/' + realign + ' \\\n')
	f_recal.write('-BQSR /lscratch/${SLURM_JOBID}/' + recal_grp +' \\\n')
	recal = tmp.strip('.bam') + '.recal.bam'
	f_recal.write('-o ./' + recal + '; \n\n') 
	# Add rm cmd
	#f_recal.write('rm ./SRR077397_2_1SRR077397_1_2_DeDup.bam; \\\n')
	#f_recal.write('rm ./SRR077397_2_1SRR077397_1_2_DeDup.bam.bai; \n\n')

def step_1():
	print "Generate data preprocessing scripts.\nPlease make sure fastq file are putted under Fastq directory."
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
	f_part1 = open('./Script/RNAPipe_part1.sh','wb')
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
	bamfiles = []
	files = [f for f in os.listdir(root+'/BAM') if os.path.isfile(os.path.join(root+'/BAM', f))]
	for my_file in files:
		if my_file.endswith('_Aligned.out.dedup.splitN.recal.bam'):
			bamfiles.append(my_file)
	for bam in bamfiles:
		print bam
	print len(bamfiles),"BAM FILES"
	return bamfiles

def var_call(root,bam,f_varcall):
	tmp = bam
	f_varcall.write('cd ' + root + '; \\\n')
	f_varcall.write('module load GATK; \\\n')
	f_varcall.write('java -Xmx64g -Djava.io.tmpdir=/lscratch/${SLURM_JOBID} -jar $GATK_JAR \\\n')
	f_varcall.write('-T HaplotypeCaller -nct 32 \\\n')
	f_varcall.write('-R /data/wangj36/genomes/ensembl/ensembl_GRCh38.82/ensembl_homo_sapiens_GRCh38.82_dna_primary_assembly.fa \\\n')
	f_varcall.write('-I ./BAM/' + bam + ' \\\n')
	tmp = bam.strip('.bam') + '.vcf'
	f_varcall.write('-o ./Variant/' + tmp +' \\\n')
	f_varcall.write('-A Coverage -A QualByDepth -A FisherStrand -A MappingQualityRankSumTest -A ReadPosRankSumTest \\\n')
	f_varcall.write('-minPruning 3 -dontUseSoftClippedBases \\\n')
	f_varcall.write('-stand_call_conf 20.0 -stand_emit_conf 20.0; \\\n')
	f_varcall.write('cd ./Variant; \\\n')
	f_varcall.write('bgzip ' + tmp + '; \\\n')
	f_varcall.write('tabix -p vcf ' + tmp + '.gz;\n')
	return tmp+'.gz'

def merge_var(root,raw_var,f_mergevar):
	vcf_files = ' '.join(raw_var)
	f_mergevar.write('cd ' + root + '\Variant; \\\n')
	f_mergevar.write('module load vcftools; \\\n')
	f_mergevar.write('vcf-merge ' + vcf_files + ' > GATK_All_Variants.vcf')
	return 'GATK_All_Variants.vcf'

def SNP_recal(root,var_whole,f_snp):
	f_snp.write('cd ' + root + '\Variant; \\\n')
	f_snp.write('module load GATK; \\\n')
	f_snp.write('java -Xmx64g -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
	f_snp.write('-T VariantRecalibrator \\\n')
	f_snp.write('-R /data/wangj36/genomes/ensembl/ensembl_GRCh38.82/ensembl_homo_sapiens_GRCh38.82_dna_primary_assembly.fa \\\n')
	f_snp.write('-input ./' + var_whole + '\\\n')
	## Set Training need consideration
	f_snp.write('-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ')
	f_snp.write('-resource:omni,known=false,training=true,truth=true,prior=12.0 ')
	f_snp.write('-resource:1000G,known=false,training=true,truth=false,prior=10.0 ')
	f_snp.write('-resource:dbsnp ')
	f_snp.write('-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \\\n')
	f_snp.write('-mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\\n')
	f_snp.write('-recalFile ./recalibrate_SNP.recal \\\n')
	f_snp.write('-tranchesFile ./recalibrate_SNP.traches \\\n')
	f_snp.write('-rscriptFile ./recalibrate_SNP_plots.R; \n')

	f_snp.write('java -Xmx24g -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
	f_snp.write('-T ApplyRecalibration \\\n')
	f_snp.write('-R /data/wangj36/genomes/ensembl/ensembl_GRCh38.82/ensembl_homo_sapiens_GRCh38.82_dna_primary_assembly.fa \\\n')
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
	f_indel.write('-R /data/wangj36/genomes/ensembl/ensembl_GRCh38.82/ensembl_homo_sapiens_GRCh38.82_dna_primary_assembly.fa \\\n')
	f_indel.write('-input ./' + var_SNP + '\\\n')
	## Set Training need consideration
	f_indel.write('-resource:mills,known=true,training=true,truth=true,prior=12.0 ')
	#f_indel.write('-resource:omni,known=false,training=true,truth=true,prior=12.0 ')
	#f_indel.write('-resource:1000G,known=false,training=true,truth=false,prior=10.0 ')
	f_indel.write('-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ')
	f_indel.write('-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \\\n')
	f_indel.write('-mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 \\\n')
	f_indel.write('-recalFile ./recalibrate_SNP_INDEL.recal \\\n')
	f_indel.write('-tranchesFile ./recalibrate_SNP_INDEL.traches \\\n')
	f_indel.write('-rscriptFile ./recalibrate_SNP_INDEL_plots.R; \n')

	f_indel.write('java -Xmx24g -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
	f_indel.write('-T ApplyRecalibration \\\n')
	f_indel.write('-R /data/wangj36/genomes/ensembl/ensembl_GRCh38.82/ensembl_homo_sapiens_GRCh38.82_dna_primary_assembly.fa \\\n')
	f_indel.write('-input ./recalibrated_SNPS.vcf  \\\n')
	f_indel.write('-mode INDEL \\\n')
	f_indel.write('-ts_filter_level 99.0 \\\n')
	f_indel.write('-recalFile ./recalibrate_SNP_INDEL.recal \\\n')
	f_indel.write('-tranchesFile ./recalibrate_SNP_INDEL.traches \\\n')
	f_indel.write('-o recalibrated_SNP_INDEL.vcf;\n')

	return "recalibrated_SNP_INDEL.vcf"

def Filter(root,var_whole,f_filter):
	f_filter.write('cd ' + root + '\Variant; \\\n')
	f_filter.write('module load GATK; \\\n')
	f_filter.write('gunzip GATK_All_Variants.vcf.gz; \\\n')
	f_filter.write('java -Xmx64g -jar $GATK_HOME/GenomeAnalysisTK.jar \\\n')
	f_filter.write('-T VariantFiltration \\\n')
	f_filter.write('-R /data/wangj36/genomes/ensembl/ensembl_GRCh38.82/ensembl_homo_sapiens_GRCh38.82_dna_primary_assembly.fa \\\n')
	f_filter.write('-V GATK_All_Variants.vcf \\\n')
	f_filter.write('-window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" \\\n')
	f_filter.write('-o GATK_All_Variants_filtered.vcf')

def step_2():
	print "Generate variants call scripts.\nPlease make sure bam file are putted under BAM directory."
	root = get_root()
	bams = get_bam(root) # get all the name of bam files
	#Genereating Scripts
	f_varcall = open('./Script/VarCall.swarm','wb')
	f_mergevar = open('./Script/VarMerge.swarm','wb')
	#f_snp = open('./Script/RunVsqrSNP.swarm','wb')
	#f_indel = open('./Script/RunVsqrINDEL.swarm','wb')
	f_filter = open('./Script/RunFilter.swarm','wb')
	
	raw_var = []
	for bam in bams:
		raw_var.append(var_call(root,bam,f_varcall))
	var_whole = merge_var(root,raw_var,f_mergevar)
	#var_SNP = SNP_recal(root,var_whole,f_snp)
	#var_INDEL = INDEL_recal(root,var_SNP,f_indel)
	Filter(root,var_whole,f_filter)
	
	f_part2 = open('./Script/RNAPipe_part2.sh','wb')
	f_part2.write('#!/bin/bash \n')
	f_part2.write('cd '+root+'/Script; \n')
	f_part2.write('VarCallID=`swarm -f VarCall.swarm -g 60 -t 32 --gres=lscratch:256`\n')
	f_part2.write('MergeID=`swarm -f VarMerge.swarm -g 60 --dependency afterany:$VarCallID --partition quick`\n')
	f_part2.write('FilterID=`swarm -f RunFilter.swarm -g 60 --gres=lscratch:128 --dependency afterany:$MergeID --partition quick`\n')
	f_part2.close()
	
	f_varcall.close()
	f_mergevar.close()
	#f_snp.close()
	#f_indel.close()
	f_filter.close()
	return

def take_argv():
	step = None
	i = 1
	while i < len(argv):
		if argv[i] == '-step':
			step = argv[i+1]
		i += 2
	if not step:
		usage()
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

