# Linked-reads: Improving and benchmarking of variant calling pipeline
## MSc project

This repository serves the purpose of being a supplement to my hand-in of my MSc thesis.
Here the majority of commands and scripts used undir the MSc project can be found.

Please see the MSc thesis in the doc folder for more details of the projects purpose, methods and results.

Short scripts are found benight here and the larger scripts are found in the script folder.

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


## Downloading data

Make subset (1000 reads) of data for testing 

## Downloading references

### Modifying bedfiles
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
vim joint_genotyping.nf (ger broytingar tú vilt)
```
Push changes to GitHub
```
git add joint_genotyping.nf
git commit -m "<the_changes_make_and_why>"
git push
```
After this go to GitHub - make New pull request


## Inserting Read Groups


To upload:
addRG.py
