# Run riboviz on Eddie

This page describes how you can run riboviz on Eddie

---
### Log in
Connect to the cluster using ssh from a terminal window (Linux and Mac OS) or use a client such as MobaXterm (Windows)

`$ ssh -X <YOUR UUN>@eddie.ecdf.ed.ac.uk`

### Ask to work inside a node interactively

`$qlogin -l h_vmem=8G`

This means that I ask for 8GB RAM

Much more info available from: https://www.wiki.ed.ac.uk/display/ResearchServices/Interactive+Sessions

### Change to group directory

`cd /exports/csce/eddie/biology/groups/wallace_rna`

### Load necessary modules on node

```
$module load igmm/apps/BEDTools 
$module load igmm/apps/bowtie
$module load igmm/apps/hdf5 
$module load igmm/apps/HISAT2
$module load igmm/apps/pigz
$module load R/3.5.3
$module load anaconda
```
### Activate environment

`$source activate riboviz`

### Run the workflow

Remember to change to the `riboviz` directory

To run the Python workflow:

```console
$ python -m riboviz.tools.prep_riboviz -c vignette/vignette_config.yaml
Running under Python: 3.7.6 | packaged by conda-forge | (default, Jun  1 2020, 18:57:50)
[GCC 7.5.0]
Created by: RiboViz Date: 2020-06-06 02:09:31.484844 Command-line tool: /exports/csce/eddie/biology/groups/wallace_rna/riboviz/riboviz/tools/prep_riboviz.py File: /exports/csce/eddie/biology/groups/wallace_rna/riboviz/riboviz/tools/prep_riboviz.py Version: commit 0fe0191f585d55763d00283d129a12e4c3c1e5c1 date 2020-06-04 00:43:45-07:00
Configuration file: vignette/vignette_config.yaml
Command file: run_riboviz_vignette.sh
Number of processes: 1
Build indices for alignment, if necessary/requested
Build indices for alignment (vignette/input/yeast_rRNA_R64-1-1.fa). Log: vignette/logs/20200606-020931/hisat2_build_r_rna.log
Build indices for alignment (vignette/input/yeast_YAL_CDS_w_250utrs.fa). Log: vignette/logs/20200606-020931/hisat2_build_orf.log
Processing samples
Processing sample: WTnone
Processing file: vignette/input/SRR1042855_s1mi.fastq.gz
Cut out sequencing library adapters. Log: vignette/logs/20200606-020931/WTnone/01_cutadapt.log
Remove rRNA or other contaminating reads by alignment to rRNA index files. Log: vignette/logs/20200606-020931/WTnone/02_hisat2_rrna.log
Align remaining reads to ORFs index files using hisat2. Log: vignette/logs/20200606-020931/WTnone/03_hisat2_orf.log
Trim 5' mismatches from reads and remove reads with more than 2 mismatches. Log: vignette/logs/20200606-020931/WTnone/04_trim_5p_mismatch.log
Convert SAM to BAM and sort on genome. Log: vignette/logs/20200606-020931/WTnone/05_samtools_view_sort.log
Index BAM file. Log: vignette/logs/20200606-020931/WTnone/06_samtools_index.log
Calculate transcriptome coverage for + strand and save as a bedgraph. Log: vignette/logs/20200606-020931/WTnone/07_bedtools_genome_cov_plus.log
Calculate transcriptome coverage for - strand and save as a bedgraph. Log: vignette/logs/20200606-020931/WTnone/08_bedtools_genome_cov_minus.log
Make length-sensitive alignments in H5 format. Log: vignette/logs/20200606-020931/WTnone/09_bam_to_h5.log
Create summary statistics, and analyses and QC plots for both RPF and mRNA datasets. Log: vignette/logs/20200606-020931/WTnone/10_generate_stats_figs.log
Finished processing sample: vignette/input/SRR1042855_s1mi.fastq.gz
Processing sample: WT3AT
Processing file: vignette/input/SRR1042864_s1mi.fastq.gz
Cut out sequencing library adapters. Log: vignette/logs/20200606-020931/WT3AT/01_cutadapt.log
Remove rRNA or other contaminating reads by alignment to rRNA index files. Log: vignette/logs/20200606-020931/WT3AT/02_hisat2_rrna.log
Align remaining reads to ORFs index files using hisat2. Log: vignette/logs/20200606-020931/WT3AT/03_hisat2_orf.log
Trim 5' mismatches from reads and remove reads with more than 2 mismatches. Log: vignette/logs/20200606-020931/WT3AT/04_trim_5p_mismatch.log
Convert SAM to BAM and sort on genome. Log: vignette/logs/20200606-020931/WT3AT/05_samtools_view_sort.log
Index BAM file. Log: vignette/logs/20200606-020931/WT3AT/06_samtools_index.log
Calculate transcriptome coverage for + strand and save as a bedgraph. Log: vignette/logs/20200606-020931/WT3AT/07_bedtools_genome_cov_plus.log
Calculate transcriptome coverage for - strand and save as a bedgraph. Log: vignette/logs/20200606-020931/WT3AT/08_bedtools_genome_cov_minus.log
Make length-sensitive alignments in H5 format. Log: vignette/logs/20200606-020931/WT3AT/09_bam_to_h5.log
Create summary statistics, and analyses and QC plots for both RPF and mRNA datasets. Log: vignette/logs/20200606-020931/WT3AT/10_generate_stats_figs.log
Finished processing sample: vignette/input/SRR1042864_s1mi.fastq.gz
File not found: vignette/input/example_missing_file.fastq.gz
Finished processing 3 samples, 1 failed
Collate TPMs across sample results. Log: vignette/logs/20200606-020931/collate_tpms.log
Count reads. Log: vignette/logs/20200606-020931/count_reads.log
Completed
```

To run the Nextflow workflow:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.04.1
Launching `prep_riboviz.nf` [big_shirley] - revision: 281e0a4d55
No such sample file (NotHere): example_missing_file.fastq.gz
[eb/da58ed] Submitted process > cutAdapters (WT3AT)
[40/456217] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[c7/83b32a] Submitted process > cutAdapters (WTnone)
[1c/9c76ba] Submitted process > buildIndicesrRNA (yeast_rRNA)
[75/c7ac3b] Submitted process > hisat2rRNA (WT3AT)
[8a/52ead5] Submitted process > hisat2rRNA (WTnone)
[7c/be343a] Submitted process > hisat2ORF (WT3AT)
[57/d2ac14] Submitted process > hisat2ORF (WTnone)
[db/6b0cd9] Submitted process > trim5pMismatches (WT3AT)
[0e/f163f9] Submitted process > trim5pMismatches (WTnone)
[db/80a57c] Submitted process > samViewSort (WT3AT)
[6b/4e6432] Submitted process > samViewSort (WTnone)
[17/0557d6] Submitted process > outputBams (WT3AT)
[5e/236c46] Submitted process > outputBams (WTnone)
[b0/be630a] Submitted process > makeBedgraphs (WT3AT)
[4a/7dc7de] Submitted process > bamToH5 (WT3AT)
[50/cb98ce] Submitted process > makeBedgraphs (WTnone)
[0d/ab3d58] Submitted process > bamToH5 (WTnone)
[f1/7ba436] Submitted process > generateStatsFigs (WT3AT)
[af/08787f] Submitted process > generateStatsFigs (WTnone)
Finished processing sample: WT3AT
[68/53d480] Submitted process > renameTpms (WT3AT)
Finished processing sample: WTnone
[4d/815c11] Submitted process > renameTpms (WTnone)
[77/8cebac] Submitted process > collateTpms (WT3AT, WTnone)
[d8/a8fcb5] Submitted process > countReads
Workflow finished! (OK)
```
