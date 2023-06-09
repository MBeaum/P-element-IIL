Invaded-inbred lines DNA-seq analysis
================
Matthew Beaumont
2023-06-12

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

# BAM to FASTQ

Converting directory of BAM files to gzipped FASTQ files and splitting
PE into separate forward and reverse reads.

``` bash
for file in /path/to/directory/*.bam; do
    base=$(basename "$file" .bam)
    samtools fastq -1 "demultiplexed/${base}_1.fq.gz" -2 "demultiplexed/${base}_2.fq.gz" -0 /dev/null -s /dev/null -n "$file"
done
```

# Quality control

We then ran the FastQC tool to reassess the quality of the fq.gz files
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

# Renaming FastQ files

We then realised that the identifier lines of our PE fastq.gz read files
were not annotated with a /1 or /2, for forward and reverse reads
respectively, which would cause some errors down the line. So we ran the
following script in order to fix this.

``` bash
# Specify the output directory
output_directory="/Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/renamed_fastq2"

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Process forward reads
for file in /Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/fastq2/*_1.fq.gz; do
    output_file="$output_directory/renamed_${file##*/}"
    gzip -cd "$file" | paste - - - - | awk '{print $1"/1"; print $2; print $3; print $4}' | tr '\t' '\n' | gzip -c > "$output_file"
done

# Process reverse reads
for file in /Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/fastq2/*_2.fq.gz; do
    output_file="$output_directory/renamed_${file##*/}"
    gzip -cd "$file" | paste - - - - | awk '{print $1"/2"; print $2; print $3; print $4}' | tr '\t' '\n' | gzip -c > "$output_file"
done
```

Perhaps unnecessary.

# PopoolationTE2

First, we needed to construct the required ppileup file.

To begin, we concatenated the D. mel reference genome (r6.51) with the
P-element to generate the following FASTA file -

``` bash
cd /Volumes/Data/Tools/RefGenomes/dmel/dna/dmel_PPI251/
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

ref_genome="/Volumes/Data/Tools/RefGenomes/dyak/dyak_prin_Tai18E2_2.1_PPI251.fasta"
popte2_jar="/Volumes/Data/Tools/PopoolationTE2/popte2-v1.10.03.jar"
hier_file="/Volumes/Data/Tools/RefGenomes/hierarchies/pelement.hier"
samples=("Dyak_21" "Dyak_22" "Dyak_23" "Dyak_24" "Dyak_25" "Dyak_26" "Dyak_27" "Dyak_28" "Dyak_29" "Dyak_30" "Dyak_N3")

for sample in "${samples[@]}"; do
   bwa mem -M -t 2 "$ref_genome" "/Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/fastq/${sample}_1.fq.gz" > "/Volumes/Data/Projects/invaded_inbred_lines/dna/map/unpaired/${sample}_1.sam" &
   bwa mem -M -t 2 "$ref_genome" "/Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/fastq/${sample}_2.fq.gz" > "/Volumes/Data/Projects/invaded_inbred_lines/dna/map/unpaired/${sample}_2.sam" &
done
wait

for sample in "${samples[@]}"; do
    java -Duser.country=US -Duser.language=en -jar "$popte2_jar" se2pe --fastq1 "/Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/fastq/${sample}_1.fq.gz" \
        --fastq2 "/Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/fastq/${sample}_2.fq.gz" \
        --bam1 "/Volumes/Data/Projects/invaded_inbred_lines/dna/map/unpaired/${sample}_1.sam" \
        --bam2 "/Volumes/Data/Projects/invaded_inbred_lines/dna/map/unpaired/${sample}_2.sam" \
        --sort --output "/Volumes/Data/Projects/invaded_inbred_lines/dna/map/paired/${sample}.sort.bam" &
done
wait

java -Duser.country=US -Duser.language=en -jar "$popte2_jar" ppileup \
    --bam /Volumes/Data/Projects/invaded_inbred_lines/dna/map/paired/Dyak_21.sort.bam \
    --bam /Volumes/Data/Projects/invaded_inbred_lines/dna/map/paired/Dyak_22.sort.bam \
    --bam /Volumes/Data/Projects/invaded_inbred_lines/dna/map/paired/Dyak_23.sort.bam \
    --bam /Volumes/Data/Projects/invaded_inbred_lines/dna/map/paired/Dyak_24.sort.bam \
    --bam /Volumes/Data/Projects/invaded_inbred_lines/dna/map/paired/Dyak_25.sort.bam \
    --bam /Volumes/Data/Projects/invaded_inbred_lines/dna/map/paired/Dyak_26.sort.bam \
    --bam /Volumes/Data/Projects/invaded_inbred_lines/dna/map/paired/Dyak_27.sort.bam \
    --bam /Volumes/Data/Projects/invaded_inbred_lines/dna/map/paired/Dyak_28.sort.bam \
    --bam /Volumes/Data/Projects/invaded_inbred_lines/dna/map/paired/Dyak_29.sort.bam \
    --bam /Volumes/Data/Projects/invaded_inbred_lines/dna/map/paired/Dyak_30.sort.bam \
    --bam /Volumes/Data/Projects/invaded_inbred_lines/dna/map/paired/Dyak_N3.sort.bam \
    --map-qual 15 --hier "$hier_file" --output /Volumes/Data/Projects/invaded_inbred_lines/dna/ppileup/Dyak.ppileup.gz
```

Then we ran the basic PoPoolationTE2 pipeline to assess P-element
insertion location and abundance in the different populations.

``` bash
# We need to subsample the ppileup file to a chosen coverage depth to generate am unbiased comparison of P-element abundance.
# java -Duser.country=US -Duser.language=en -jar popte2.jar subsampleppileup --ppileup dmel.ppileup.gz --target-coverage 100 --output output.ss100.ppileup.gz

# First we generate a file of the found P-element insertions and their locations.
java -Duser.country=US -Duser.language=en -jar /Volumes/Data/Tools/PopoolationTE2/popte2-v1.10.03.jar identifySignatures --ppileup /Volumes/Data/Projects/invaded_inbred_lines/dna/ppileup/Dmel.ppileup.gz --mode separate --output /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel.signatures --min-count 3

java -Duser.country=US -Duser.language=en -jar /Volumes/Data/Tools/PopoolationTE2/popte2-v1.10.03.jar identifySignatures --ppileup /Volumes/Data/Projects/invaded_inbred_lines/dna/ppileup/Dsim.ppileup.gz --mode separate --output /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim.signatures --min-count 3

java -Duser.country=US -Duser.language=en -jar /Volumes/Data/Tools/PopoolationTE2/popte2-v1.10.03.jar identifySignatures --ppileup /Volumes/Data/Projects/invaded_inbred_lines/dna/ppileup/Dyak.ppileup.gz --mode separate --output /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dyak.signatures --min-count 3

# Then we look at the frequency of the found signatures.
java -Duser.country=US -Duser.language=en -jar /Volumes/Data/Tools/PopoolationTE2/popte2-v1.10.03.jar frequency --ppileup /Volumes/Data/Projects/invaded_inbred_lines/dna/ppileup/Dmel.ppileup.gz --signature /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel.signatures --output /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel.freqsig

java -Duser.country=US -Duser.language=en -jar /Volumes/Data/Tools/PopoolationTE2/popte2-v1.10.03.jar frequency --ppileup /Volumes/Data/Projects/invaded_inbred_lines/dna/ppileup/Dsim.ppileup.gz --signature /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim.signatures --output /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim.freqsig

java -Duser.country=US -Duser.language=en -jar /Volumes/Data/Tools/PopoolationTE2/popte2-v1.10.03.jar frequency --ppileup /Volumes/Data/Projects/invaded_inbred_lines/dna/ppileup/Dyak.ppileup.gz --signature /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dyak.signatures --output /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dyak.freqsig

# Finally, we combine it all together. 
java -Duser.country=US -Duser.language=en -jar /Volumes/Data/Tools/PopoolationTE2/popte2-v1.10.03.jar pairupSignatures --signature /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel.freqsig --ref-genome /Volumes/Data/Tools/RefGenomes/dmel/dna/dmel_PPI251/dmelallchrom-r6.51_PPI251.fasta --hier /Volumes/Data/Tools/RefGenomes/hierarchies/pelement.hier --min-distance -200 --max-distance 300 --output /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel.teinsertions

java -Duser.country=US -Duser.language=en -jar /Volumes/Data/Tools/PopoolationTE2/popte2-v1.10.03.jar pairupSignatures --signature /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim.freqsig --ref-genome /Volumes/Data/Tools/RefGenomes/dsim/dsim_ASM75419v3_PPI251.fasta --hier /Volumes/Data/Tools/RefGenomes/hierarchies/pelement.hier --min-distance -200 --max-distance 300 --output /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim.teinsertions

java -Duser.country=US -Duser.language=en -jar /Volumes/Data/Tools/PopoolationTE2/popte2-v1.10.03.jar pairupSignatures --signature /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dyak.freqsig --ref-genome /Volumes/Data/Tools/RefGenomes/dyak/dyak_prin_Tai18E2_2.1_PPI251.fasta --hier /Volumes/Data/Tools/RefGenomes/hierarchies/pelement.hier --min-distance -200 --max-distance 300 --output /Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dyak.teinsertions
```

This provides us with a list of P-element insertions found in all 11
samples for each species, their location and population frequency.

# Visualisation

# Manhattan plots

We then visualised the .teinsertion files in Manhattan plots for all
samples in each species.

``` r
library(viridisLite)
library(ggplot2)
library(gridExtra)
theme_set(theme_bw())

dm <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel/dmel.teinsertions")
names(dm) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

ds <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim/dsim.teinsertions")
names(ds) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

#dy <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dyak.teinsertions")
#names(dy) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

dm <- subset(dm, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
ds <- subset(ds, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
#dy <- subset(dy, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")

dm$Chromosome <- factor(dm$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4"))
lim <- c(0.0, 0.51)
ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
ds$Chromosome <- factor(ds$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4"))
lim <- c(0.0, 0.51)
ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
#dy$Chromosome <- factor(dy$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4"))
#lim <- c(0.0, 0.51)
#ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

sample_colors <- viridisLite::turbo(11)

dm$Sample <- as.factor(dm$Sample)
ds$Sample <- as.factor(ds$Sample)
#dy$Sample <- as.factor(dy$Sample)

dmp <- ggplot(dm, aes(x = Position, y = Frequency, color = Sample)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +
  scale_color_manual(values = sample_colors) +  
  theme(legend.position = "right", 
          plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(title = "P-element insertion frequencies in D. melanogaster populations")


dsp <- ggplot(ds, aes(x = Position, y = Frequency, color = Sample)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(ds$Frequency), max(ds$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +
  scale_color_manual(values = sample_colors) +  
  theme(legend.position = "right", 
          plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(title = "P-element insertion frequencies in D. simulans populations")


#dyp <- ggplot(dy, aes(x = Position, y = Frequency, color = Sample)) +
  #geom_point(size = 0.5) +
  #facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  #scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     #labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  #scale_y_continuous(name = "Population frequency", limits = c(min(dy$Frequency), max(dy$Frequency)),
                     #breaks = seq(0, 1, by = 0.1)) +
  #scale_color_manual(values = sample_colors) +  
  #theme(legend.position = "right")

plot(dmp)
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
plot(dsp)
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

``` r
#dyp
```

Then we looked at each sample individually, starting with D. mel.

``` r
library(ggplot2)
library(gridExtra)
theme_set(theme_bw())

dm1 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel/replicates/1.teinsertions")
names(dm1) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

dm2 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel/replicates/2.teinsertions")
names(dm2) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

dm3 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel/replicates/3.teinsertions")
names(dm3) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

dm4 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel/replicates/4.teinsertions")
names(dm4) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

dm5 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel/replicates/5.teinsertions")
names(dm5) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

dm6 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel/replicates/6.teinsertions")
names(dm6) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

dm7 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel/replicates/7.teinsertions")
names(dm7) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

dm8 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel/replicates/8.teinsertions")
names(dm8) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

dm9 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel/replicates/9.teinsertions")
names(dm9) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

dm10 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dmel/replicates/10.teinsertions")
names(dm10) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

dm1 <- subset(dm1, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
dm2 <- subset(dm2, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
dm3 <- subset(dm3, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
dm4 <- subset(dm4, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
dm5 <- subset(dm5, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
dm6 <- subset(dm6, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
dm7 <- subset(dm7, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
dm8 <- subset(dm8, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
dm9 <- subset(dm9, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
dm10 <- subset(dm10, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")

dm1$Chromosome <- factor(dm1$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
dm2$Chromosome <- factor(dm2$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
dm3$Chromosome <- factor(dm3$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
dm4$Chromosome <- factor(dm4$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
dm5$Chromosome <- factor(dm5$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
dm6$Chromosome <- factor(dm6$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
dm7$Chromosome <- factor(dm7$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
dm8$Chromosome <- factor(dm8$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
dm9$Chromosome <- factor(dm9$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
dm10$Chromosome <- factor(dm10$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

dmp1 <- ggplot(dm1, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 1")

dmp2 <- ggplot(dm2, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 2")

dmp3 <- ggplot(dm3, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 3")

dmp3 <- ggplot(dm3, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 3")

dmp3 <- ggplot(dm3, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 3")

dmp4 <- ggplot(dm4, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 4")

dmp5 <- ggplot(dm5, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 5")

dmp6 <- ggplot(dm6, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 6")

dmp7 <- ggplot(dm7, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 7")

dmp8 <- ggplot(dm8, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 8")

dmp9 <- ggplot(dm9, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 9")

dmp10 <- ggplot(dm10, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 10")

dmp1
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
dmp2
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
dmp3
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->

``` r
dmp4
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

``` r
dmp5
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-15-5.png)<!-- -->

``` r
dmp6
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-15-6.png)<!-- -->

``` r
dmp7
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-15-7.png)<!-- -->

``` r
dmp8
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-15-8.png)<!-- -->

``` r
dmp9
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-15-9.png)<!-- -->

``` r
dmp10
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-15-10.png)<!-- -->

``` r
library(ggplot2)
library(gridExtra)
theme_set(theme_bw())

ds1 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim/replicates/1.teinsertions")
names(ds1) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

ds2 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim/replicates/2.teinsertions")
names(ds2) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

ds3 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim/replicates/3.teinsertions")
names(ds3) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

ds4 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim/replicates/4.teinsertions")
names(ds4) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

ds5 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim/replicates/5.teinsertions")
names(ds5) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

ds6 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim/replicates/6.teinsertions")
names(ds6) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

ds7 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim/replicates/7.teinsertions")
names(ds7) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

ds8 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim/replicates/8.teinsertions")
names(ds8) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

ds9 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim/replicates/9.teinsertions")
names(ds9) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

ds10 <- read.table("/Volumes/Data/Projects/invaded_inbred_lines/dna/popTE2/dsim/replicates/10.teinsertions")
names(ds10) <- c("Sample", "Chromosome", "Position", "Strand", "TE", "Order", "FR", "Comment", "Frequency")

ds1 <- subset(ds1, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
ds2 <- subset(ds2, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
ds3 <- subset(ds3, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
ds4 <- subset(ds4, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
ds5 <- subset(ds5, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
ds6 <- subset(ds6, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
ds7 <- subset(ds7, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
ds8 <- subset(ds8, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
ds9 <- subset(ds9, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")
ds10 <- subset(ds10, Chromosome == "X" | Chromosome == "2L" | Chromosome == "2R" | Chromosome == "3L" | Chromosome == "3R" | Chromosome == "4")

ds1$Chromosome <- factor(ds1$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
ds2$Chromosome <- factor(ds2$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
ds3$Chromosome <- factor(ds3$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
ds4$Chromosome <- factor(ds4$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
ds5$Chromosome <- factor(ds5$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
ds6$Chromosome <- factor(ds6$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
ds7$Chromosome <- factor(ds7$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
ds8$Chromosome <- factor(ds8$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
ds9$Chromosome <- factor(ds9$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
ds10$Chromosome <- factor(ds10$Chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4")) 
    lim <- c(0.0, 0.51) 
    ybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

dsp1 <- ggplot(ds1, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. simulans - replicate 1")

dsp2 <- ggplot(ds2, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. simulans - replicate 2")

dsp3 <- ggplot(ds3, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. simulans - replicate 3")

dsp3 <- ggplot(ds3, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. simulans - replicate 3")

dsp3 <- ggplot(ds3, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. simulans - replicate 3")

dsp4 <- ggplot(ds4, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. simulans - replicate 4")

dsp5 <- ggplot(ds5, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. simulans - replicate 5")

dsp6 <- ggplot(ds6, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. simulans - replicate 6")

dsp7 <- ggplot(ds7, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. simulans - replicate 7")

dsp8 <- ggplot(ds8, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. melanogaster - replicate 8")

dsp9 <- ggplot(ds9, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. simulans - replicate 9")

dsp10 <- ggplot(ds10, aes(x = Position, y = Frequency, color = Frequency)) +
  geom_point(size = 0.5) +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0", "5m", "10m", "15m", "20m", "25m")) +
  scale_y_continuous(name = "Population frequency", limits = c(min(dm$Frequency), max(dm$Frequency)),
                     breaks = seq(0, 1, by = 0.1)) +  scale_color_gradient(low = "blue", high = "red") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) + labs(title = "P-element insertion frequencies in D. simulans - replicate 10")

dsp1
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
dsp2
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

``` r
dsp3
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->

``` r
dsp4
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-16-4.png)<!-- -->

``` r
dsp5
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-16-5.png)<!-- -->

``` r
dsp6
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-16-6.png)<!-- -->

``` r
dsp7
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-16-7.png)<!-- -->

``` r
dsp8
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-16-8.png)<!-- -->

``` r
dsp9
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-16-9.png)<!-- -->

``` r
dsp10
```

![](DNA-seq_analysis_files/figure-gfm/unnamed-chunk-16-10.png)<!-- -->

# IGV

We wanted to assess whether the results from PopoolationTE2 were
expected based on the coverage of input files relative to a reference.
To do this, we simply gave a quick visual assessment on IGV.

First, we needed to map the files against the indexed ISO1 reference
with bwa, convert them to BAM files and sort them.

``` bash
#!/bin/bash

ref="/Volumes/Data/Tools/RefGenomes/dmel/dna/dmel-all-chromosome-r6.51.fasta"
if="/Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/fastq"
of="/Volumes/Data/Projects/invaded_inbred_lines/dna/demultiplexed/map-bwamem"

bwa mem -t 12 $ref $if/Dmel_1_1.fq.gz $if/Dmel_1_2.fq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_1.sort.bam
bwa mem -t 12 $ref $if/Dmel_2_1.fq.gz $if/Dmel_2_2.fq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_2.sort.bam
bwa mem -t 12 $ref $if/Dmel_3_1.fq.gz $if/Dmel_3_2.fq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_3.sort.bam
bwa mem -t 12 $ref $if/Dmel_4_1.fq.gz $if/Dmel_4_2.fq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_4.sort.bam
bwa mem -t 12 $ref $if/Dmel_5_1.fq.gz $if/Dmel_5_2.fq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_5.sort.bam
bwa mem -t 12 $ref $if/Dmel_6_1.fq.gz $if/Dmel_6_2.fq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_6.sort.bam
bwa mem -t 12 $ref $if/Dmel_7_1.fq.gz $if/Dmel_7_2.fq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_7.sort.bam
bwa mem -t 12 $ref $if/Dmel_8_1.fq.gz $if/Dmel_8_2.fq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_8.sort.bam
bwa mem -t 12 $ref $if/Dmel_9_1.fq.gz $if/Dmel_9_2.fq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_9.sort.bam
bwa mem -t 2 $ref $if/Dmel_10_1.fq.gz $if/Dmel_10_2.fq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_10.sort.bam
bwa mem -t 2 $ref $if/Dmel_N1_1.fq.gz $if/Dmel_N1_2.fq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_N1.sort.bam
```

We then open the BAM files in IGV and compare to the reference.

(.png of coverage)
