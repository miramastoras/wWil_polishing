# Wwil nanopore assembly polishing

This document contains methods used for polishing the wWil nanopore genome assembly

Useful paper comparing polishing tools for microbial genome assembly: https://www.nature.com/articles/s41598-021-00178-w#Sec2


Spreadsheet of results: https://docs.google.com/spreadsheets/d/1Rbb5gen7m2lebTOjihjBFObD-WStY5hlO6VrOfL8lvU/edit#gid=0

Slides: https://docs.google.com/presentation/d/1uRFTi6YXE-Jekr8PCdXKE4QVoxRLXadJIZi0G_wNA-U/edit#slide=id.g2cb56d0e0f1_0_90

## 1. Small variant polishing (SNP and INDELs <50bp)

### Step 1: Establish assembly QC metrics on unpolished assembly

Downsample illumina bam from 760x to 40x for QV calculation
```
#!/bin/bash
#SBATCH --job-name=wWil_downsample
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=3:00:00

samtools view -s 0.055 -b -h -@ 32 /private/groups/russelllab/cade/wwil_polishing/JW18wWil0703A_wWil-only.bam > /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only.downsample.40x.bam
```

#### Merqury

Estimate the best k size for this genome:
```
# Get number of bases in unpolished assembly
seqkit stats /private/groups/patenlab/mira/wWil_polishing/data/wWil_Nanopore_assembly.fasta
file                          format  type  num_seqs    sum_len    min_len    avg_len    max_len
wWil_Nanopore_assembly.fasta  FASTA   DNA          1  1,266,954  1,266,954  1,266,954  1,266,954

# calculate merqury estimate of best k size
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups juklucas/hpp_merqury:latest /opt/merqury/best_k.sh 1266954

genome: 1266954
tolerable collision rate: 0.001
15.1186
```

Build a kmer db for merqury: Also testing k21 to see how it goes
```
# convert bam to fastq
samtools fastq -@16 /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only.downsample.40x.bam > /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only.downsample.40x.fastq

# meryl count
cd /private/groups/patenlab/mira/wWil_polishing/data/illumina
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups juklucas/hpp_merqury:latest meryl k=15 count /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only.downsample.40x.fastq output /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k15.meryl

cd /private/groups/patenlab/mira/wWil_polishing/data/illumina
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups juklucas/hpp_merqury:latest meryl k=21 count /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only.downsample.40x.fastq output /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k21.meryl
```

Run merqury on unpolished assembly
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/wWil_unpolished:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k15.meryl /private/groups/patenlab/mira/wWil_polishing/data/wWil_Nanopore_assembly.fasta

# wWil_Nanopore_assembly	65437	1266940	24.5233	0.00352918

docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/wWil_unpolished_k21:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k21.meryl /private/groups/patenlab/mira/wWil_polishing/data/wWil_Nanopore_assembly.fasta wWil_unpolished_merqury_k21

# wWil_Nanopore_assembly	78792	1266934	25.1529	0.0030529
```

kmer completeness:
```
cat /private/groups/patenlab/mira/wWil_polishing/merqury/wWil_unpolished/wWil_unpolished_merqury_k15.completeness.stats
wWil_Nanopore_assembly  all     1099602 1114097 98.6989
```
QV estimate is 24.523, which means an error every 283 bases

#### Busco gene completeness:

from paper: BUSCO17 v5.1.1 was used with enterobacterales_odb10 database

```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco/unpolished:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/data/wWil_Nanopore_assembly.fasta \
    -c 16 \
    -l bacteria_odb10
```
Results:
```
***** Results: *****

	C:42.7%[S:42.7%,D:0.0%],F:37.1%,M:20.2%,n:124
	53	Complete BUSCOs (C)
	53	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	46	Fragmented BUSCOs (F)
	25	Missing BUSCOs (M)
	124	Total BUSCO groups searched

Assembly Statistics:
	1	Number of scaffolds
	1	Number of contigs
	1266954	Total length
	0.000%	Percent gaps
	1 MB	Scaffold N50
	1 MB	Contigs N50
```

### Step 2: Run multiple polishing tools:

Align nanopore reads to unpolished assembly

Reads from Jodie: `/private/groups/russelllab/jodie/wWil_genome/Wwil_fastq.tar.gz`
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 -ax map-ont -t 14 \
    /private/groups/patenlab/mira/wWil_polishing/data/wWil_Nanopore_assembly.fasta \
    /private/groups/patenlab/mira/wWil_polishing/data/nanopore/Wwil_nanopore.all.fastq.gz \
    > /private/groups/patenlab/mira/wWil_polishing/data/nanopore/Wwil.R10.mm2.wWil_unpolished.sam
```

Racon polishing: Using params suggested for Racon + medaka combo
https://www.nature.com/articles/s41598-021-00178-w#Sec2
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    staphb/racon:1.4.20 racon \
    -m 8 -x -6 -g -8 -w 500 \
    /private/groups/patenlab/mira/wWil_polishing/data/nanopore/Wwil_nanopore.all.fastq.gz \
    /private/groups/patenlab/mira/wWil_polishing/data/nanopore/Wwil.R10.mm2.wWil_unpolished.sam \
    /private/groups/patenlab/mira/wWil_polishing/data/wWil_Nanopore_assembly.fasta \
    > /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.fasta
```
Medaka following Racon:
```
docker run -it --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    staphb/medaka:1.2.0 medaka_consensus \
    -t 14 -m r103_min_high_g345 \
    -i /private/groups/patenlab/mira/wWil_polishing/data/nanopore/Wwil_nanopore.all.fastq.gz \
    -d /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.fasta \
    -o /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/medaka/

```
Racon + Homopolish polishing

```
conda activate homopolish

python3 /private/home/mmastora/progs/homopolish/homopolish.py polish \
    -a /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.fasta \
    -s /private/home/mmastora/progs/homopolish/bacteria.msh \
    -m /private/home/mmastora/progs/homopolish/R10.3.pkl \
    -o /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.homopolish.fasta
```

Pilon only short read polishing

```
java -Xmx16G -jar /private/home/mmastora/progs/pilon-1.24.jar \
    --genome /private/groups/patenlab/mira/wWil_polishing/data/wWil_Nanopore_assembly.fasta \
    --bam /private/groups/russelllab/cade/wwil_polishing/JW18wWil0703A_wWil-only.bam \
    --output wWil_Nanopore_assembly.pilon.polished \
    --outdir /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/
```

Racon+Medaka+Pilon
```
# bwa index
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa index -p \
    /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/medaka/consensus.fasta \
    /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/medaka/consensus.fasta

# align illumina reads to racon+medaka polished assembly
docker run -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa mem -t32 \
    /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/medaka/consensus.fasta \
    /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only_R1.fastq.gz \
    /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only_R2.fastq.gz \
    | samtools view -b -h \
    > /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.Racon_Medaka_polished.bam

# run pilon
java -Xmx16G -jar /private/home/mmastora/progs/pilon-1.24.jar \
    --genome /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/medaka/consensus.fasta \
    --bam /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.Racon_Medaka_polished.srt.bam\
    --output wWil_Nanopore_assembly.racon.medaka.pilon.polished \
    --outdir /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/
```


Racon+homopolish+Pilon
```
# bwa index
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa index -p \
    /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.homopolish.fasta/wWil_Nanopore_assembly_homopolished.fasta \
    /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.homopolish.fasta/wWil_Nanopore_assembly_homopolished.fasta

# align illumina reads to racon+homopolish polished assembly
docker run -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa mem -t32 \
    /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.homopolish.fasta/wWil_Nanopore_assembly_homopolished.fasta \
    /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only_R1.fastq.gz \
    /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only_R2.fastq.gz \
    | samtools view -b -h \
    > /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.Racon_homopolish_polished.bam

# run pilon
java -Xmx16G -jar /private/home/mmastora/progs/pilon-1.24.jar \
    --genome /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.homopolish.fasta/wWil_Nanopore_assembly_homopolished.fasta \
    --bam /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.Racon_homopolish_polished.srt.bam \
    --output wWil_Nanopore_assembly.racon.homopolish.pilon.polished \
    --outdir /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/
```


Pilon x 2
```
# bwa index
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa index -p \
    /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.pilon.polished.fasta \
    /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.pilon.polished.fasta

# align illumina reads to pilon polished assembly
docker run -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa mem -t32 \
    /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.pilon.polished.fasta \
    /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only_R1.fastq.gz \
    /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only_R2.fastq.gz \
    | samtools view -b -h \
    > /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.pilon_polished.bam

samtools sort /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.pilon_polished.bam > /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.pilon_polished.srt.bam

samtools index /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.pilon_polished.srt.bam

# run pilon
java -Xmx16G -jar /private/home/mmastora/progs/pilon-1.24.jar \
    --genome /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.pilon.polished.fasta \
    --bam /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.pilon_polished.srt.bam \
    --output wWil_Nanopore_assembly.pilonx2.polished \
    --outdir /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/
```

### Step 3: Compare QC metrics across different polishing tools

#### QV with merqury

Racon only: QV improves by .2
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Racon_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.fasta wWil_Racon_merqury_k15

wWil_Nanopore_assembly.Racon    61574   1264923 24.7869 0.00332132
```

Racon + medaka: another .6 QV improvement
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Racon_Medaka_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/medaka/consensus.fasta wWil_Racon_Medaka_merqury_k15

wWil_Wolbachia_chromosome       54950   1268440 25.3045 0.00294813
```
Racon + homopolish:
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Racon_homopolish_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.homopolish.fasta/wWil_Nanopore_assembly_homopolished.fasta wWil_Racon_homopolish_merqury_k15

cat wWil_Racon_homopolish_merqury_k15.qv
wWil_Nanopore_assembly_homopolished     55086   1265680 25.2839 0.00296216
```

Pilon only
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Pilon_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.pilon.polished.fasta wWil_Pilon_merqury_k15

wWil_Nanopore_assembly.pilon.polished   52071   1268400 25.5429 0.0027907
```

Pilon x 2
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Pilonx2_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.pilonx2.polished.fasta wWil_Pilonx2_merqury_k15
```

Racon+Medaka+Pilon
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Racon_Medaka_Pilon_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.racon.medaka.pilon.polished.fasta wWil_racon_medaka_Pilon_merqury_k15

wWil_Nanopore_assembly.racon.medaka.pilon.polished      51928   1268472 25.5553 0.00278272
```

Racon+homopolish+pilon
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Racon_homopolish_Pilon_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.racon.homopolish.pilon.polished.fasta wWil_racon_homopolish_Pilon_merqury_k15

```
#### Busco

Racon only
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco/Racon:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/busco/Racon/wWil_Nanopore_assembly.Racon.fasta \
    -c 16 \
    -l bacteria_odb10

    ---------------------------------------------------
    |Results from dataset bacteria_odb10               |
    ---------------------------------------------------
    |C:60.5%[S:60.5%,D:0.0%],F:25.8%,M:13.7%,n:124     |
    |75    Complete BUSCOs (C)                         |
    |75    Complete and single-copy BUSCOs (S)         |
    |0    Complete and duplicated BUSCOs (D)           |
    |32    Fragmented BUSCOs (F)                       |
    |17    Missing BUSCOs (M)                          |
    |124    Total BUSCO groups searched                |
    ---------------------------------------------------
```

Racon + Medaka
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco/Racon_Medaka:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/medaka/consensus.fasta \
    -c 16 -f \
    -l bacteria_odb10

---------------------------------------------------
|Results from dataset bacteria_odb10               |
---------------------------------------------------
|C:77.4%[S:77.4%,D:0.0%],F:10.5%,M:12.1%,n:124     |
|96    Complete BUSCOs (C)                         |
|96    Complete and single-copy BUSCOs (S)         |
|0    Complete and duplicated BUSCOs (D)           |
|13    Fragmented BUSCOs (F)                       |
|15    Missing BUSCOs (M)                          |
|124    Total BUSCO groups searched                |
---------------------------------------------------
```

Racon + homopolish
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco/Racon_homopolish:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.homopolish.fasta/wWil_Nanopore_assembly_homopolished.fasta \
    -c 16 -f \
    -l bacteria_odb10

---------------------------------------------------
|Results from dataset bacteria_odb10               |
---------------------------------------------------
|C:83.1%[S:83.1%,D:0.0%],F:5.6%,M:11.3%,n:124      |
|103    Complete BUSCOs (C)                        |
|103    Complete and single-copy BUSCOs (S)        |
|0    Complete and duplicated BUSCOs (D)           |
|7    Fragmented BUSCOs (F)                        |
|14    Missing BUSCOs (M)                          |
|124    Total BUSCO groups searched                |
---------------------------------------------------
```

Pilon only
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco/Pilon:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.pilon.polished.fasta \
    -c 16 -f \
    -l bacteria_odb10


    ---------------------------------------------------
    |Results from dataset bacteria_odb10               |
    ---------------------------------------------------
    |C:83.9%[S:83.9%,D:0.0%],F:4.8%,M:11.3%,n:124      |
    |104    Complete BUSCOs (C)                        |
    |104    Complete and single-copy BUSCOs (S)        |
    |0    Complete and duplicated BUSCOs (D)           |
    |6    Fragmented BUSCOs (F)                        |
    |14    Missing BUSCOs (M)                          |
    |124    Total BUSCO groups searched                |
    ---------------------------------------------------
```

Pilon x 2
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco/Pilonx2:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.pilonx2.polished.fasta \
    -c 16 -f \
    -l bacteria_odb10

```

Racon+Medaka+Pilon
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco/Racon_Medaka_Pilon:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.racon.medaka.pilon.polished.fasta \
    -c 16 -f \
    -l bacteria_odb10

```

```
---------------------------------------------------
  |Results from dataset bacteria_odb10               |
  ---------------------------------------------------
  |C:82.3%[S:82.3%,D:0.0%],F:6.5%,M:11.2%,n:124      |
  |102    Complete BUSCOs (C)                        |
  |102    Complete and single-copy BUSCOs (S)        |
  |0    Complete and duplicated BUSCOs (D)           |
  |8    Fragmented BUSCOs (F)                        |
  |14    Missing BUSCOs (M)                          |
  |124    Total BUSCO groups searched                |
  ---------------------------------------------------
```

Racon+homopolish+Pilon
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco/Racon_homopolish_Pilon:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.racon.homopolish.pilon.polished.fasta \
    -c 16 -f \
    -l bacteria_odb10
```

```
---------------------------------------------------
|Results from dataset bacteria_odb10               |
---------------------------------------------------
|C:83.1%[S:83.1%,D:0.0%],F:5.6%,M:11.3%,n:124      |
|103    Complete BUSCOs (C)                        |
|103    Complete and single-copy BUSCOs (S)        |
|0    Complete and duplicated BUSCOs (D)           |
|7    Fragmented BUSCOs (F)                        |
|14    Missing BUSCOs (M)                          |
|124    Total BUSCO groups searched                |
---------------------------------------------------
```


Align nanopore reads to pilonx1 polished assembly to look in IGV

```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 -ax map-ont -t 14 \
    /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.pilon.polished.fasta \
    /private/groups/patenlab/mira/wWil_polishing/data/nanopore/Wwil_nanopore.all.fastq.gz \
    | samtools view -bh \
    > /private/groups/patenlab/mira/wWil_polishing/data/nanopore/Wwil.R10.mm2.wWil_pilon_polished.bam
```
