
#To create my WGS_variant_calling environment, just go to Terminal and introduce:
conda create env -f environment.yml

#Or create your environment manually introducing:
conda create -n wgs_variant_calling -c bioconda -c conda-forge fastqc fastp bwa samtools qualimap picard freebayes vcftools igv

#Download the reads
-wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SR062/SRR062634/SRR062634_1.fastq.gz
-wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SR062/SRR062634/SRR062634_2.fastq.gz

#Download the Homo sapiens GRCh38 reference genome
wget -c ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GCRh38.dna.primary_assembly.fa.gz