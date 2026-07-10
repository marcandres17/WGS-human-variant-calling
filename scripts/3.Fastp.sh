
#Trimm your reads using FastP using the next command:
fastp -i SRR062634_1.fastq.gz -I SRR062634_2.fastq.gz -o out_1.clean.fq.gz -O out_2.cleanfq.gz --detect_adapter_for_pe --trim_poly_x --trim_poly_g --cut_front 20 --cut_tail 20 --cut_mean_quality 20 -h out_FastP.html
