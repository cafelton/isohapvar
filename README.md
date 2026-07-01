# isohapvar
Take a haplotyped vcf and isoform bed and create a variant-aware isoform fasta

Longshot and longcallR variant calling from long-read RNA-seq have been tested, but this should work with variants called from WGS as well.

This has been tested using FLAIR isoforms, but should work with any bed12 file.
```
python3 add_vars_to_isos_070126.py --vcf SAMPLE.longshot.vcf --isoform_bed /SAMPLE.isoforms.bed --genome GENOME.fa --output TEST
```
To visualize the resulting variant-aware isoforms, I recommend aligning them to the genome.

```
minimap2 -ax splice --secondary=no genome.fa test.isoswithvars.fa | samtools view -hb - | samtools sort - > test.isoswithvars.bam; samtools index test.isoswithvars.bam
```
