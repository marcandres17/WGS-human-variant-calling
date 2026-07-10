
#Now we have to transform our SAM file into a more-readable BAM file using samtools:
samtools view -bS align.sam > align.bam

#Then we sort the coordenates
samtools sort align.bam align.sorted.bam

#We have to index the BAM file in order to get the stats of the mapping
samtools index align.sorted.bam

#Obtain the stats
samtools flagstat align.sorted.bam