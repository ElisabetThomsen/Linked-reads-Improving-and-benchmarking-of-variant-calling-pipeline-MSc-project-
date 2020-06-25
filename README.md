# Linked-reads: Improving and benchmarking of variant calling pipeline
## MSc project

This repository serves the purpose of being a supplement to my hand-in of my MSc thesis.
Here the majority of commands and scripts used undir the MSc project can be found.

Please see the MSc thesis in the doc folder for more details of the projects purpose, methods and results.

The content of the bash scripts and other unix commands are found beneath here. The python and R scripts are found in the script folder. The used nextflow.config-files and timelines of the pipelines are in the nextflow folder.

## Rancher
To get access to the server rancher was used. Following commands were used:

To list available containers:
```
rancher kubectl get pods -n <namespace>
```
Go into a container:
```
rancher kubectl exec -it <container_name> -n <namespace> -- /bin/bash
```

Transfering files between the author's laptop and the server:
```
rancher kubectl cp <copy_this_file> <to_this_file>
```
For files from the server "<name_space>/<container_name>:/home/path/to/file" should be used.


## Conda environments

To install software a .yml file containing the software was made and the conda environment was made by running the command
```
conda env create -f <env_name.yml>
```
The software were found by searching on https://anaconda.org/


## Downloading data and making subsets

10X exome data of NA12878
```
wget http://cf.10xgenomics.com/samples/genome/2.1.4/NA12878_WES_v2/NA12878_WES_v2_fastqs.tar
tar -xvf NA12878_WES_v2_fastqs.tar
```
Make subset (1000 reads) of data for testing
```
less /path/to/NA12878_WES_v2_fastqs/NA12878_WES_v2_S1_L004_R1_001.fastq.gz | head -1000 > new_name_L004_R1.fastq
```

Make subset with only barcode contaminated reads
```
mkdir NA12878_bctrimmed
mkdir NA12878_withbc

# For these R2 files ..
for file in $(ls /path/to/10xdata/NA12878_WES_v2/*R2*fastq.gz)
do

# This is the belonging R1
r1=$(echo $file | sed 's/R2/R1/')

# Change name, so they will be seperated into bctrimmed and withbc folders
out2=$(echo $file | sed 's/.*_v2\//NA12878_bctrimmed\//' | sed 's/.gz//')
out_org=$(echo $file | sed 's/.*_v2\//NA12878_withbc\//' | sed 's/.gz//')

# Trim them, and save only reads that had barcode contamination (all contaminated reads are printed before and after trim, into seperate files)
echo Trimming $file ...
time ./trimR2bc_whitelist.py $r1 $file /path/to/whitelist/4M-with-alts-february-2016.txt $out2 $out_org 1>>NA12878_bctrim_stats.txt

#### Filter the R1s that belong to the trimmed R2s ####

mkdir tmp

# Name all the files, so we can use them in repair.sh
r1_bctrim_sync=$(echo $out2 | sed 's/R2/R1/')
r2_bctrim_sync=$(echo $out2 | sed 's/NA12878_bctrimmed/tmp/')
s_bctrim_sync=$(echo $out2 | sed 's/R2/s/' | sed 's/NA12878_bctrimmed/tmp/')

# Syncronize the R1 and R2 files
echo Synchronizing $r1 and $out2 ...
time repair.sh in1=$r1 in2=$out2 out1=$r1_bctrim_sync out2=$r2_bctrim_sync outs=$s_bctrim_sync repair 2>>repair.log

# Gzip them
echo gzipping $r1_bctrim_sync and $out2 ...
gzip $r1_bctrim_sync
gzip $out2

# Do the same for the untrimmed file
r1_withbc_sync=$(echo $out_org | sed 's/R2/R1/')
r2_withbc_sync=$(echo $out_org | sed 's/NA12878_withbc/tmp/')
s_withbc_sync=$(echo $out_org | sed 's/R2/s/' | sed 's/NA12878_withbc/tmp/')
echo Synchronizing $r1 and $out_org ...
time repair.sh in1=$r1 in2=$out_org out1=$r1_withbc_sync out2=$r2_withbc_sync outs=$s_withbc_sync repair 2>>repair.log
echo gzipping $r1_withbc_sync and $out_org ...
gzip $r1_withbc_sync
gzip $out_org

# Clean up
rm -rf tmp

done
```

## Downloading references

GRCh38 were downloaded with the gatk_bundle.sh provided in linkseq repository (https://github.com/olavurmortensen/linkseq/tree/master/reference).

b37
```
# Reference genome
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.fai.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.dict.gz

# Known variants
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase3_v4_20130502.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase3_v4_20130502.sites.vcf.gz.tbi

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_omni2.5.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_omni2.5.b37.vcf.idx.gz

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.idx.gz

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.idx.gz

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.idx.gz

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/hapmap_3.3.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/hapmap_3.3.b37.vcf.idx.gz
```
All files except 1000G_phase3_v4_20130502.sites.vcf.gz and 1000G_phase3_v4_20130502.sites.vcf.gz.tbi were then gunzipped with
```
gunzip <filename>
```
index b37 reference
```
bwa index -a bwtsw human_g1k_v37.fasta
```

### Modifying bedfiles
The targets BED-files were downloaded from Agilent’s webpage ( https://earray.chem.agilent.com/suredesign/ ) and then modified.

<b>GRCh38</b> - adding 0 and . in the 5th and 6th column
```
awk 'BEGIN{OFS="\t"}{ if(NR <= 2) { print "#"$0 } else { print $1,$2,$3,$4,0,"." } }' ../sureselect_human_all_exon_v6_utr_grch38/S07604624_Padded.bed > S07604624_Padded_UTR_modified.bed
```
<b>hg19</b> - adding 0 and . and changing chrX to X
```
awk 'BEGIN{OFS="\t"}{ if(NR <= 2) { gsub(/chr/, ""); print "#"$0 } else { print substr($1,4),$2,$3,$4,0,"." } }' ../SureSelect_Human_All_Exon_V6_r2_hg19/S07604514_Padded.bed > S07604514_Padded_noUTR_hg19_modified.bed
```

## GitHub
First make a fork. Then cloning LinkSeq
```
git clone https://github.com/ElisabetThomsen/linkseq.git
```
Make a new branch and go onto it
```
git checkout SNP-INDEL-fix
```
Make changes to a script
```
vim joint_genotyping.nf
```
Push changes to GitHub
```
git add joint_genotyping.nf
git commit -m "<the_changes_make_and_why>"
git push
```
After this go to GitHub - make New pull request


## Inserting Read Groups

Convert BAM to SAM
```
in_bam=<path/to/bam>
time samtools view -h -o noRG.sam $in_bam
```

Add lane information to SAM
```
time ./addRG.py noRG.sam RG.sam :4,:5 1,4
```

Clean up and convert SAM to BAM
```
rm noRG.sam
time samtools view -S -b RG.sam > sortedRG.bam
rm RG.sam
```

## Plot hard filter parameters

Annotate the variants with if they are TP or FP.<br />
The -E uses "giab".<expression>, from the --resource:"giab", so they must be the same
```
gatk VariantAnnotator \
    -V /path/to/indel_or_snp.vcf \
    --resource:giab /home/elisabetht/resources/GiaB/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf \
    -E giab.callsets \
    -O giab_ann.vcf
```
Make a new tab-seperated file, that is more easy for R to handle
```
gatk VariantsToTable \
    -V giab_ann.vcf \
    -F CHROM -F POS -F QUAL \
    -F BaseQRankSum -F MQRankSum -F ReadPosRankSum \
    -F DP -F FS -F MQ -F QD -F SOR \
    -F giab.callsets \
    -GF GQ \
    -O giab_ann.txt
 ``` 
 Seperate TP and FP into seperate files
```
./seperateTPFP.py giab_ann.txt
```
Make the plots. They will be saved to "TPFPplots.pdf".<br />
Note: If you change the name "giab_ann.txt" and "giab.callsets", you need to change it inside TPFPplots.R also.
```
Rscript TPFPhisto.R
```

## Running LinkSeq
Run the main script
```
time nextflow run /path/to/linkseq/main.nf -c nextflow.config --sample "NA12878" --fastq_r1 "/path/to/fastqs/*R1*" --fastq_r2 "/path/to/fastqs/*R2*" -with-report NA12878_report_main.html -with-timeline NA12878_timeline_main.html -resume
```

Run the joint_genotyping script
```
time nextflow run /path/to/linkseq/joint_genotyping.nf -c nextflow.config --gvcf_path "/path/to/gvcf_folder" -with-report NA12878_report_joint.html -with-timeline NA12878_timeline_joint.html -resume
```

## Running Long Ranger
```
time ./../software/longranger-2.2.2/longranger targeted \
    --id=NA12878 \
    --reference=/path/to/refdata-GRCh38-2.1.0 \
    --fastqs=/path/to/10xdata/NA12878_WES_v2 \
    --targets=/path/to/SureSelect_Human_All_Exon_V6_r2_hg38/S07604514_Padded.bed \
    --vcmode=freebayes
```

## Variant correctness benchmark

hap.py
```
export HGREF=/home/elisabetht/resources/gatk_bundle/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta

mkdir NA12878_happy
time hap.py \
    /path/to/GiaB/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf \
    /path/to/phased.vcf.gz \
    -f /path/to/GiaB/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed \
    -o NA12878_happy/NA12878 \
    -r /path/to/Homo_sapiens_assembly38.fasta \
    -T /path/to/modified_BEDs/S07604514_Covered_noUTR_modified.bed
```

Split the multi-sample called vcf, before benchmark
```
gatk SelectVariants \
    -R /path/to/human_g1k_v37.fasta \
    -V /path/to/variants.vcf.gz \
    --sample-name <Samplename> \
    -O <samplename>.vcf.gz
```
GATK's Concordance
```
gatk Concordance \
    -R /path/to/human_g1k_v37.fasta \
    -eval /path/to/<samplename>.vcf.gz \
    --truth /path/to/Old.vcf.gz \
    --summary concordance/concordance_<Samplename>_LR.tsv \
    -L /path/to/modified_BEDs/S07604514_Padded_noUTR_hg19_modified.bed
```

## Phase evaluation

error_rates.py
```
python /path/to/error_rates.py /path/to/GiaB/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf /path/to/blockfile > phase_errors.txt
```

Whatshap
```
whatshap stats --tsv=whatshap_out/stats.txt --block-list=whatshap_out/blocklist.txt --chr-lengths=chr_lengths_tab.txt $vcf
```

## Venn diagram
Note: this script was made in last moment, so the coding is not so good. It could have been improved by making loop and a tsv file to read into R.

First get the VCF-files
```
# LinkSeq VCF
LS="/path/to/phased.vcf.gz"
# Long Ranger VCF
LR="/path/to/phased_variants.vcf.gz"
# GiaB truth VCF
GB="/path/to/GiaB/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf"
# Targets
T="/path/to/SureSelect_Human_All_Exon_V6_r2_hg38/S07604514_Padded.bed"
```

Filter to keep only variants from target regions
```
bedtools intersect -header -a $GB -b $T -wa > GB_T.vcf
GB="GB_T.vcf"
bedtools intersect -header -a $LR -b $T -wa > LR_T.vcf
LR="LR_T.vcf"
bedtools intersect -header -a $LS -b $T -wa > LS_T.vcf
LS="LS_T.vcf"
```

Isolate INDELs (or SNPs)
```
gatk SelectVariants \
    -V $GB \
    -select-type INDEL \
    -O GBindel.vcf

gatk SelectVariants \
    -V $LR \
    -select-type INDEL \
    -O LRindel.vcf

gatk SelectVariants \
    -V $LS \
    -select-type INDEL \
    -O LSindel.vcf
    
GB="GBindel.vcf"
LR="LRindel.vcf"
LS="LSindel.vcf"
```

Print number of variants needed for the venn diagram R-script
```
echo Variants in LS
less $LS | grep -v "^#" | uniq | wc -l

echo Variants in LR
less $LR | grep -v "^#" | uniq | wc -l

echo Variants in GB
less $GB | grep -v "^#" | uniq | wc -l

echo LS intersect with GB
bedtools intersect -a $GB -b $LS -wa | uniq | wc -l

echo LR intersect with GB
bedtools intersect -a $GB -b $LR -wa | uniq | wc -l

echo LS intersect with LR
bedtools intersect -a $LS -b $LR -wa | uniq | wc -l

echo LS intersect with LR and GB
bedtools intersect -header -a $GB -b $LS -wa | bedtools intersect -a stdin -b $LR -wa | uniq | wc -l
```

To upload:
addRG.py
seperateTPFP.py
TPFPhisto.R
Modified demux
newflow.config (main)
modified trimR2bc.py
venndiagram.R
