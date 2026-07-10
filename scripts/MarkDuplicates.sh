
#The next step is getting rid of the duplicates and try to obtain the most curated sequence we can get and build a proper file so variant calling tools can read it. To do so, firstly, we have to create the read groups with information about our reads. We can obtain some info where we extracted our reads (NCBI) and we can also explore the first lines of the file:
zcat out1_clean.fq.gz | head -n 1

#These first lines tell us:
#HWI-EAS110_103327062-Sequenciator ID(--RGID/--RGPU)
#6-Sequentiatior chip lane (--RGID/--RGPU)
#NCBI WEB
#28.456.850- Real name/label (--RGLB)
#Illumina Genome Analyzer II-Illumina platform (--RGLP)

#Now we cann add our groups with PicardTools:
picard AddOrReplaceReadGroups --INPUT align.sorted.bam --OUTPUT align.sorted.rg.bam --RGID HWI-EAS110_103327062.6 --RGLB 2845856850 --RGPL illumina --RGPU HWI-EAS110_103327062 --RGSM SRR062634

#Mark duplicates
picard MarkDuplicates --INPUT align.sorted.rg.bam --OUTPUT align.dedup.bam --METRICS_FILE markDuplicatesMetrics.txt --ASSUME_SORTED True