
#Before the mapping step, we need the indexed file of the reference genome:
bwa index Homo_sapiens.GCRh38.dna.primary_assembly.fa

#Now, we can proceed with the mapping
bwa mem -a Homo_sapiens.GCRh38.dna.primary_assembly.fa out1.clean.fq.gz out2.clean.fq.gz -o align.sam 2> stderror.out
