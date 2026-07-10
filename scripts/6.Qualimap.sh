
#We use Qualimap to observe the analysis of the quality of our mapping:
qualimap bamqc -bam align.sorted.bam -nt 4 --java-mem-size=12G

#We move onto the new folder we created:
cd align.sorted_stats

#Open the .html file with the core information of our mapping
firefox *.HTML
