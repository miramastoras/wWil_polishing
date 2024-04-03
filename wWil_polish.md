# Wwil nanopore assembly polishing

This document contains analysis used for polishing wWil nanopore genome assembly with illumina reads

Useful paper comparing polishing tools for microbial genome assembly: https://www.nature.com/articles/s41598-021-00178-w#Sec2
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
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Racon_homopolish_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.homopolish.fasta/wWil_Nanopore_assembly.Racon.homopolish.fasta wWil_Racon_homopolish_merqury_k15

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
## 2. Check for any SVs that need polishing
