# <p align=center>Bulk RNA-Seq analysis</p>

**<p align=left>What the workflow does</p>**
This nextflow workflow is for the QC, mapping and read counting of bulk RNA-Seq. The workflow accept and analyze both single and paired end RNA-Seq data.  
The quality of the FastQ file are determine using FastQC and trimming is done used trimmomatic with these parameters `LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20` taking into account the quality score of the reads.
<br />
The are mapped to the reference genome using HISAT. To enable fasting mapping the FastQ files are split to smaller chuck which are mapped individually and the respective bam files merge into one and sorted by coordinate. Then the mapping quality of the bam files is generated using Samtools. 
<br />
HTSeq is then used to count reads generating four outputs (count files) for stranded libraries: `genes.htseq-union.firststrand.counts`, `genes.htseq-union.secondstrand.count`, `genes.htseq-union.firststrand.nonunique.counts` and `genes.htseq-union.secondstrand.nonunique.counts`. And two count files for un-stranded library: `genes.htseq-union.unstranded.counts` and `genes.htseq-union.unstranded.nonunique.counts` representing unique and non-unique respectively. 


**<p align=left>Running the workflow</p>**
To run the work the following dependencies need to be install
* Docker
* Nextflow
> url https://get.nextflow.io | bash



<br />
<br />
<br />

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
