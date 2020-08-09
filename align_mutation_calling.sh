#/bin/sh
#MINIMAP2=/home/xiaofei/theSoftware/miniconda3/bin/minimap2
BWA=/home/xiaofei/theSoftware/bwa/bwa
SAMTOOLS=/home/xiaofei/theSoftware/miniconda3/bin/samtools
BCFTOOLS=/home/xiaofei/theSoftware/miniconda3/bin/bcftools

REFDIR=/home/xiaofei/data/2019_ncov/variation
REF=$REFDIR/ref_genome.fasta

#$MINIMAP2 -ax asm5 -DP ref_genome.fasta all_2020_0305.complete.fasta | $SAMTOOLS view -b > all_SARS-COV-2.bam
FADIR=/home/xiaofei/data/2019_ncov/variation/theVirusFa
OUTDIR=/home/xiaofei/data/2019_ncov/variation/theVariance
$BWA index $REF

cd $OUTDIR
for i in `ls $FADIR`; do
	echo $i
	$BWA mem $REF $FADIR/$i | samtools view -b -o $i.bam
	$SAMTOOLS mpileup -uf $REF $i.bam | bcftools call -mv > $i.var.raw.vcf
	$BCFTOOLS filter -s LowQual -e '%QUAL<20 || DP>100' $i.var.raw.vcf  > $i.var.flt.vcf
done

