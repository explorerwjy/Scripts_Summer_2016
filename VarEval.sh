#!/bin/bash
module load GATK;
java -jar $GATK_JAR \
	-T VariantEval \
	-R /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \
	-eval:myCalls recalibrated_SNP_INDEL.vcf -L chr1\
	-eval:1000G_SNP /fdb/GATK_resource_bundle/hg19/1000G_phase1.snps.high_confidence.hg19.vcf.gz \
	-eval:1000G_INDEL /fdb/GATK_resource_bundle/hg19/1000G_phase1.indels.hg19.vcf.gz \
	-EV TiTvVariantEvaluator -EV CompOverlap \
	-o Chr1.eval
	#-o recalibrated_SNP_INDEL.eval

