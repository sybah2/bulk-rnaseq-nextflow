# <p align=center>Bulk RNA-Seq analysis</p>

**<p align=left>What the workflow does</p>**
This nextflow workflow is for the QC, mapping and read counting of bulk RNA-Seq. The workflow accept and analyze both single and paired end RNA-Seq data.  
The quality of the FastQ file are determine using FastQC and trimming is done used trimmomatic with these parameters `LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20` taking into account the quality score of the reads.
<br />
<br />
The are mapped to the reference genome using HISAT. To enable fasting mapping the FastQ files are split to smaller chuck which are mapped individually and the respective bam files merge into one and sorted by coordinate. Then the mapping quality of the bam files is generated using Samtools. 
<br />
<br />
HTSeq is then used to count reads generating four outputs (count files) for stranded libraries: `genes.htseq-union.firststrand.counts`, `genes.htseq-union.secondstrand.count`, `genes.htseq-union.firststrand.nonunique.counts` and `genes.htseq-union.secondstrand.nonunique.counts`. And two count files for un-stranded library: `genes.htseq-union.unstranded.counts` and `genes.htseq-union.unstranded.nonunique.counts` representing unique and non-unique respectively. 
<br />
<br />
The workflow accept already downloaded FastQ files and also SRA accession number of FastQ files which will be automatically downloaded and analyze.
<br />

**<p align=left>Get Started</p>**
To run the work the following dependencies need to be install
* Docker
> `https://docs.docker.com/engine/install/`
* Nextflow
> `curl https://get.nextflow.io | bash`

* The pull the git hub repo using the following command
> `git pull https://github.com/VEuPathDB/bulk-rnaseq-nextflow.git`

* Alternatively the workflow can be run directly using nextflow which pull down the repo. 
> `nextflow run VEuPathDB/bulk-rnaseq-nextflow -with-trace -c  <config_file> -r main`


<br />
<br />

**<p align=left>Input Data</p>**
Example data can be found in the data directory. 
* nextflow configA
> The nextflow config file (`nextflow.config`) contain the configuration for analysis. File paths to where the sequence reads, reference genome and output directory are specify in the config file. See example in workflow parent directory. 
* Input Fastq files
> If the Fastq files are already downloaded they should be store under the data folder. If the data is not downloaded the SRA accession number need to be specified in a csv file and its path set in the config file. 
* Reference genome
> The reference genome should be added in a folder under data and its path specify in config file (See example in the data folder).

**<p align=left>Output Data</p>**
Depending on weather the library is stranded or not the following output are generated and found in the results folder

* Stranded
> `genes.htseq-union.firststrand.counts`, `genes.htseq-union.secondstrand.count`, `genes.htseq-union.firststrand.nonunique.counts` and `genes.htseq-union.secondstrand.nonunique.counts`

* Un-stranded 
> `genes.htseq-union.unstranded.counts` and `genes.htseq-union.unstranded.nonunique.counts`
***<p align=center>Nextflow workflow diagram</p>*** 
```mermaid
flowchart TD
    p0((Channel.fromFilePairs))
    p1(( ))
    p2(( ))
    p3[rna_seq:createIndex]
    p4[rna_seq:qualityControl]
    p5[rna_seq:fastqcCheck]
    p6([first])
    p7[rna_seq:pairedEndTrimming]
    p10([splitFastq])
    p11(( ))
    p12[rna_seq:hisatMappingPairedEnd]
    p13[rna_seq:sortSam]
    p14([groupTuple])
    p15[rna_seq:mergeSams]
    p17[rna_seq:sortBams]
    p18(( ))
    p19(( ))
    p20(( ))
    p21[rna_seq:htseqCounting]
    p22(( ))
    p23[rna_seq:bedBamStats]
    p24(( ))
    p25[rna_seq:spliceCrossingReads]
    p26(( ))
    p0 -->|reads_ch| p4
    p1 -->|organismAbbv| p3
    p2 -->|reference| p3
    p3 --> p12
    p3 --> p12
    p4 --> p5
    p4 --> p7
    p5 --> p6
    p6 -->|check_fastq| p7
    p7 --> p10
    p10 -->|reads| p12
    p6 -->|check_fastq| p12
    p11 -->|intronLength| p12
    p12 --> p13
    p13 --> p14
    p14 -->|samSet| p15
    p15 --> p17
    p17 --> p21
    p18 -->|annotation| p21
    p19 -->|isCds| p21
    p20 -->|isStranded| p21
    p21 --> p22
    p15 --> p23
    p23 --> p24
    p15 --> p25
    p25 --> p26
```
