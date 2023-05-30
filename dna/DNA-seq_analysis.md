Invaded inbred DNA-seq
================
Matthew Beaumont
2023-05-27

We were provided with the following BAM files containing the multiplexed
DNA sequencing reads for all of our samples (plus the i-seq probe).

``` bash

cd /Volumes/Data/Projects/invaded_inbred_lines/dna/raw
ls -h *.bam
```

    ## HY5WYDRX2_1_20230419B_20230420.bam
    ## HY5WYDRX2_2_20230419B_20230420.bam
    ## i-seq_BSB09408-2634_1_20230329B_20230330.bam
    ## iil_merged.bam

These two files were merged using the built in samtools command.

``` bash

samtools merge -o iil_merged.bam HY5WYDRX2_1_20230419B_20230420.bam HY5WYDRX2_2_20230419B_20230420.bam
```

Then, we ran the merged BAM file through the following two scripts,
respectively. This was in order to demultiplex the files and extract out
each read by their index primer sequence.

``` bash

samtools view --threads 10 /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/iil_merged.bam | paste -d "|" - - | tee >(grep BC:Z:ATCACG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dmel_1.bam) >(grep BC:Z:CGATGT |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dmel_2.bam) >(grep BC:Z:TTAGGC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dmel_3.bam) >(grep BC:Z:TGACCA |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dmel_4.bam) >(grep BC:Z:ACAGTG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dmel_5.bam) >(grep BC:Z:GCCAAT |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dmel_6.bam) >(grep BC:Z:CAGATC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dmel_7.bam) >(grep BC:Z:ACTTGA |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dmel_8.bam) >(grep BC:Z:GATCAG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dmel_9.bam) >(grep BC:Z:TAGCTT |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dmel_10.bam) >(grep BC:Z:GGCTAC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dsim_11.bam) >(grep BC:Z:CTTGTA |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dsim_12.bam) >(grep BC:Z:AGTCAA |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dsim_13.bam) >(grep BC:Z:AGTTCC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dsim_14.bam) >(grep BC:Z:ATGTCA |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dsim_15.bam) >(grep BC:Z:CCGTCC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dsim_16.bam) >(grep BC:Z:GTAGAG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dsim_17.bam) >(grep BC:Z:GTCCGC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dsim_18.bam) >(grep BC:Z:GTGAAA |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dsim_19.bam) >(grep BC:Z:GTGGCC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dsim_20.bam) >(grep BC:Z:GTTTCG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dyak_21.bam) >(grep BC:Z:CGTACG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dyak_22.bam) >(grep BC:Z:GAGTGG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dyak_23.bam) >(grep BC:Z:GGTAGC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dyak_24.bam) >(grep BC:Z:ACTGAT |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dyak_25.bam) >(grep BC:Z:ATGAGC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dyak_26.bam) >(grep BC:Z:ATTCCT |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dyak_27.bam) >(grep BC:Z:CAAAAG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dyak_28.bam) >(grep BC:Z:CAACTA |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dyak_29.bam) >(grep BC:Z:CACCGG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dyak_30.bam) >(grep BC:Z:CACGAT |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dmel_N1.bam) >(grep BC:Z:CACTCA |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dsim_N2.bam) >(grep BC:Z:CAGGCG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dyak_N3.bam) >(grep BC:Z:CATTTT |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/Dere_A2.bam) >(grep BC:Z:GGGGGG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/GGGGGG.bam) | grep BC:Z:NNNNNN |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/demultiplexed/NNNNNN.bam

samtools view --threads 10 /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/HY5WYDRX2_1_20230419B_20230420.bam | paste -d "|" - - | tee   >(grep BC:Z:CATGGC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer34.bam) >(grep BC:Z:CCAACA |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer36.bam) >(grep BC:Z:TAATCG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer42.bam) >(grep BC:Z:CGGAAT |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer37.bam) >(grep BC:Z:CTAGCT |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer38.bam) >(grep BC:Z:CTATAC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer39.bam) >(grep BC:Z:GTGATC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer40.bam) >(grep BC:Z:GACGAC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer41.bam) >(grep BC:Z:TACAGC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer43.bam) >(grep BC:Z:TATAAT |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer44.bam) >(grep BC:Z:TCATTC |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer45.bam) >(grep BC:Z:TCCCGA |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer46.bam) >(grep BC:Z:TCGAAG |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer47.bam) >(grep BC:Z:TCGGCA |tr "|" "\n" | samtools view --threads 2 -b > /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/primer48.bam)
```

Bam to fq.gz

``` bash

for i in D*bam;
do n=${i%.bam};
samtools view $i| tee >(awk 'NR%2 == 1'| awk '{print "@" $1"/1"; print $10; print "+" $1"/1";
print $11}'| gzip -c > ${n}_1.fq.gz)| awk 'NR%2 == 0'| awk '{print "@" $1"/2"; print $10; print  
"+" $1"/2"; print $11}'| gzip -c > ${n}_2.fq.gz;   
done
```

# Quality control

We then ran teh FastQC tool to reassess the quality of the fq.gz files
for each individual sample, granting the following results.

``` bash
```

# DeviaTE analysis

Fastq-miner

``` bash

nohup zsh fastq-miner.sh iil /Volumes/Data/Projects/invaded_inbred_lines/dna/raw/output/FastQ > /Volumes/Data/Projects/invaded_inbred_lines/logs/IIlines.log &
```

``` bash

nohup zsh deviate-family.sh iil PPI251 > /Volumes/Data/Projects/invaded_inbred_lines/logs/IIlines.log &
```

# Attempt 2

De-multiplexing as before …

Converting directory of BAM files to gzipped FASTQ files and splitting
PE into separate forward and reverse reads.

``` bash

for file in /path/to/directory/*.bam; do
    base=$(basename "$file" .bam)
    samtools fastq -1 "${base}_1.fq.gz" -2 "${base}_2.fq.gz" -0 /dev/null -s /dev/null -n "$file"
done
```

\#PopoolationTE2

First, we needed to construct the required ppileup file.

To begin, we concatenated the D. mel reference genome (r6.51) with the
P-element to generate the following FASTA file -

``` bash

cd /Volumes/Data/Tools/RefGenomes/dmel
ls -h dmelall*.fasta
```

    ## dmelallchrom-r6.51_PPI251.fasta

Then created an index and /map directory for later steps.

``` bash

bwa index walkthrough-refgenome/2R-603-consensusTE.fasta
mkdir map
```

Then we ran the following PPileupGen.sh script on the pairs of forward
and reverse reads, generating the required PPileup file.

``` bash

nohup zsh PopoolationTE2/PPileupGen.sh > logs/PPileupGen.log
```

``` bash

#!/bin/bash

ref_genome="RefGenomes/dmel/dmelallchrom-r6.51_PPI251.fasta"
popte2_jar="PopoolationTE2/popte2-v1.10.03.jar"
hier_file="pelement.hier"
samples=("Dmel_1" "Dmel_2" "Dmel_3" "Dmel_4" "Dmel_5" "Dmel_6" "Dmel_7" "Dmel_8" "Dmel_9" "Dmel_10" "Dmel_N1")

for sample in "${samples[@]}"; do
    bwa bwasw -t 10 "$ref_genome" "/Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/fastq/${sample}_1.fq.gz" > "RefGenomes/map/${sample}_1.sam" &
    bwa bwasw -t 10 "$ref_genome" "/Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/fastq/${sample}_2.fq.gz" > "RefGenomes/map/${sample}_2.sam" &
done

wait

for sample in "${samples[@]}"; do
    java -jar "$popte2_jar" se2pe --fastq1 "/Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/fastq/${sample}_1.fq.gz" \
        --fastq2 "/Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/fastq/${sample}_2.fq.gz" \
        --bam1 "RefGenomes/map/${sample}_1.sam" \
        --bam2 "RefGenomes/map/${sample}_2.sam" \
        --sort --output "RefGenomes/${sample}.sort.bam" &
done

wait

java -jar "$popte2_jar" ppileup --bam RefGenomes/Dmel_1.sort.bam --bam RefGenomes/Dmel_2.sort.bam \
    --bam RefGenomes/Dmel_3.sort.bam --bam RefGenomes/Dmel_4.sort.bam --bam RefGenomes/Dmel_5.sort.bam \
    --bam RefGenomes/Dmel_6.sort.bam --bam RefGenomes/Dmel_7.sort.bam --bam RefGenomes/Dmel_8.sort.bam \
    --bam RefGenomes/Dmel_9.sort.bam --bam RefGenomes/Dmel_10.sort.bam --bam RefGenomes/Dmel_N1.sort.bam \
    --map-qual 15 --hier "$hier_file" --output RefGenomes/Dmel.ppileup.gz
```