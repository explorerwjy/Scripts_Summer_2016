echo "Generate swarm filtering"
echo "#!/bin/bash
for entry in \`ls -a chr*.vcf\`
do
	echo \"python variants_filter.py \$entry Exonic_filt_\$entry;\"
done" > TMP.TiTvRaio.Gen.sh;
bash TMP.TiTvRaio.Gen.sh > variants_filter.swarm;
rm TMP.TiTvRaio.Gen.sh;
echo "Submit swarm filtering"
FILTER_ID=`swarm -f variants_filter.swarm -g 5 --partition quick`
echo -e '#!/bin/bash\ntouch DONE' \
	| sbatch --time=5 --partition=quick --dependency=afterany:$FILTER_ID
while true; do
	[[ -f DONE ]] && break
	echo 'waiting swarm tasks finish ...';
	sleep 1m
done
rm DONE
rm variants_filter.swarm

module load vcftools;
echo "Concat vcf file together..."
vcf-concat Exonic_filt_*.vcf | gzip -c > All_filtered.multianno.recode.vcf.gz;
gunzip All_filtered.multianno.recode.vcf.gz;
module load GATK;
java -jar $GATK_JAR \
	-T VariantEval \
	-R /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \
	-eval:myCalls All_filtered.multianno.recode.vcf\
	-EV TiTvVariantEvaluator -EV CompOverlap \
	-o All_filtered.multianno.recode.eval
grep 'TiTvVariantEvaluator' All_filtered.multianno.recode.eval > All_filtered.multianno.recode.titv

