#Firstly we must create an index of the deduplicate file
samtools index align.dedup.bam

#Create our variant calling file
freebayes -f Homo_sapiens.GRCh38.dna.primary_assembly.fa --min-coverage 1 align.dedup.bam > variants.vcf

#Obtain variant calling stats
rtg vcfstats variants.vcf > variants.vcfstats

#Prepare two files to do a separate analysis: Indels and SNPs. Indels are more likely to contain errors, so the min quality required is lower than SNPs
vcftools --vcf variants.vcf --keep-only-indels --minQ 30 --recode --recode-INFO-all --out variants_indels.vcf
vcftools --vcf variants.vcf --remove-indels --minQ 20 --recode --recode-INFO-all --out variants_snvs.vcf

#Zip the files so we can introduce them in VEP website and get our final variant calling stats of the project
bgzip variants_snps.vcf.recode.vcf
bgzip variants_indels.vcf.recode.vcf