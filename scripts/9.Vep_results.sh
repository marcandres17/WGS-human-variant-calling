#First, we download the .vcf files in VEP website. Those files were based on our VCF files previously introduced
#Then, we create a script to obtain the HIGH value varieties of both Indels and SNPs files

# Indels script
awk -F '\t' '!/^#/{
	if($8 ~ /CSQ=/) {
#Isolate the CSQ block
		split($8, t1, "CSQ=");
		split(t1[2], cb, ";");

#Separate the different transcripts (comma-separated)
		n= split(cb[1], tr, ",");
#Save already printed combinations to avoid duplicates at the same position
		p = "";
		for(i=1; i<=n; i++) {
#Separate the internal VEP fields (pipe-separated). Standard VEP order is Allele|Consequence|Impact|Gene...
			split(tr[i], f, "|");
			cons = f[2];
			imp = f[3];
			gen = f[4];

			if (gen == "") gen = "Unknown";
#Filter consequences of interest
			if (imp == "HIGH" && cons ~ /frameshift_variant|stop_gained|stop_lost|missense_variant/) {
#Create a unique key per position-consequence-gene to avoid cluttering the file
				key = $1 "_" $2 "_" cons "_" gen;

				if (p !~ key) {
					print $1 "\t" $2 "\t" cons "\t" gen;
					p = p " " key
				}
			}
		}
	}
}' results_indels.vcf | sort -V -u > results_indels.txt

#SNPs script
awk -F '\t' '!/^#/{
	if($8 ~ /CSQ=/) {
#Isolate the CSQ block
		split($8, t1, "CSQ=");
		split(t1[2], cb, ";");

#Separate the different transcripts (comma-separated)
		n= split(cb[1], tr, ",");
#Save already printed combinations to avoid duplicates at the same position
		p = "";
		for(i=1; i<=n; i++) {
#Separate the internal VEP fields (pipe-separated). Standard VEP order is Allele|Consequence|Impact|Gene...
			f_n = split(tr[i], f, "|");
			cons = f[2];
			imp = f[3];
			gen = f[4];

			if (gen == "") gen = "Unknown";
#Filter importance and consequences of interest
			if (imp == "HIGH" && cons ~ /stop_gained|stop_lost/) {
#Create a unique key per position-consequence-gene to avoid cluttering the file
				key = $1 "_" $2 "_" cons "_" gen;

				if (p !~ key) {
					print $1 "\t" $2 "\t" cons "\t" imp "\t" gen;
					p = p " " key
				}
			}
		}
	}
}' results_snps.vcf | sort -V -u > results_snps.txt

# Once both files are created, the following command allows us to see which genes overlap between them. This helps us to identify potencial affected genes who are not false positives.
comm -12 <(awk -F'\t' '{print $4}' results_snps.txt | sort -u) <(awk -F'\t' '{print $4}' results_indels.txt | sort -u)

#We have obtained that CYP4B1 and ZNF717 are in both files. Now we search for their coordenates to see where are the mutations so we can search them with IGV
cat results_indels.txt | grep -E "CYP4B1|ZNF717"
cat results_snps.txt | grep -E "CYP4B1|ZNF717"

#Now we start IGV with 10GB RAM and search for our coordenates
_JAVA_OPTIONS='-Xmx10g' igv

#Load the following files in IGV
variants_indels.vcf.recode.vcf.gz
variants_snps.vcf.recode.vcf.gz
align.sorted.bam
align.sorted.bam.bai
