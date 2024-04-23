# Wwil nanopore assembly polishing

This document contains methods used for polishing the wWil nanopore genome assembly

Useful paper comparing polishing tools for microbial genome assembly: https://www.nature.com/articles/s41598-021-00178-w#Sec2


Spreadsheet of results: https://docs.google.com/spreadsheets/d/1Rbb5gen7m2lebTOjihjBFObD-WStY5hlO6VrOfL8lvU/edit#gid=0

Slides: https://docs.google.com/presentation/d/1uRFTi6YXE-Jekr8PCdXKE4QVoxRLXadJIZi0G_wNA-U/edit#slide=id.g2cb56d0e0f1_0_90

### Step 0: Establish filtering for illumina bam files

Align reads to assembly with minimap2
```
samtools fastq -@ 16 /private/groups/russelllab/cade/wwil_polishing/no_filtering/JW18wWil0703A_wWil-only-no-filtering.bam > /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only-no-filtering.fastq

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 -ax sr --cs --eqx -t 14 \
    /private/groups/patenlab/mira/wWil_polishing/data/wWil_nanopore_assembly_decontaminated.fasta \
    /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only-no-filtering.fastq \
    | samtools view -bh - > /private/groups/patenlab/mira/wWil_polishing/data/illumina/Wwil.ilm.mm2.wWil_nanopore_assembly_decontaminated.bam

```

Plot read divergences (de)
```
python3 /private/home/mmastora/progs/element_polishing/scripts/plot_de_distribution.py -b /private/groups/patenlab/mira/wWil_polishing/data/illumina/Wwil.ilm.mm2.wWil_nanopore_assembly_decontaminated.srt.bam -p /private/groups/patenlab/mira/wWil_polishing/data/illumina/de_plot_ilm
```

Remove reads with > .04 divergence

```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl
```

```
{
  "runCorrectBam.correctBam.dockerImage": "mobinasri/secphase:dev-v0.2.0-hom",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/wWil_polishing/data/illumina/Wwil.ilm.mm2.wWil_nanopore_assembly_decontaminated.srt.bam",
  "runCorrectBam.correctBam.options": "--maxDiv 0.04 --minReadLen 0 --minAlignmentLen 0",
  "runCorrectBam.correctBam.suffix": "maxDiv.04"
}
```

```
cd /private/groups/patenlab/mira/wWil_polishing/data/illumina/correct_bam

time java -jar /private/home/mmastora/progs/cromwell-85.jar \
    run /private/groups/hprc/polishing/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl \
    --inputs /private/groups/patenlab/mira/wWil_polishing/data/illumina/correct_bam/correct_bam_inputs.json
```

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

samtools view -s 0.025 -b -h -@ 32 /private/groups/patenlab/mira/wWil_polishing/data/illumina/correct_bam/Wwil.ilm.mm2.wWil_nanopore_assembly_decontaminated.srt.maxDiv.04.bam > /private/groups/patenlab/mira/wWil_polishing/data/illumina/correct_bam/Wwil.R10.mm2.wWil_nanopore_assembly_decontaminated.srt.maxDiv.04.40x.bam
```

#### Merqury

Estimate the best k size for this genome:
```
# Get number of bases in unpolished assembly
seqkit stats /private/groups/patenlab/mira/wWil_polishing/data/wWil_nanopore_assembly_decontaminated.fasta
file                                                                                           format  type  num_seqs    sum_len    min_len    avg_len    max_len
/private/groups/patenlab/mira/wWil_polishing/data/wWil_nanopore_assembly_decontaminated.fasta  FASTA   DNA          1  1,268,501  1,268,501  1,268,501  1,268,501

# calculate merqury estimate of best k size
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups juklucas/hpp_merqury:latest /opt/merqury/best_k.sh 1268501

genome: 1268501
tolerable collision rate: 0.001
15.1195
```

Build a kmer db for merqury: Also testing k21 to see how it goes
```
# convert bam to fastq
samtools fastq -@16 /private/groups/patenlab/mira/wWil_polishing/data/illumina/correct_bam/Wwil.R10.mm2.wWil_nanopore_assembly_decontaminated.srt.maxDiv.04.40x.bam > /private/groups/patenlab/mira/wWil_polishing/data/illumina/correct_bam/Wwil.R10.mm2.wWil_nanopore_assembly_decontaminated.srt.maxDiv.04.40x.fastq

# meryl count
cd /private/groups/patenlab/mira/wWil_polishing/data/illumina
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups juklucas/hpp_merqury:latest meryl k=15 count /private/groups/patenlab/mira/wWil_polishing/data/illumina/correct_bam/Wwil.R10.mm2.wWil_nanopore_assembly_decontaminated.srt.maxDiv.04.40x.fastq output /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.de.04.k15.meryl
```

Run merqury on unpolished assembly
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/wWil_unpolished:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.de.04.k15.meryl /private/groups/patenlab/mira/wWil_polishing/data/wWil_nanopore_assembly_decontaminated.fasta wWil_nanopore_assembly_decontaminated.k15.merqury
```

#### Busco gene completeness:

```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco_rickettsiales_odb10/unpolished:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/data/wWil_nanopore_assembly_decontaminated.fasta \
    -c 16 \
    -f \
    -l rickettsiales_odb10

    ---------------------------------------------------
     |Results from dataset rickettsiales_odb10          |
     ---------------------------------------------------
     |C:98.6%[S:98.6%,D:0.0%],F:0.8%,M:0.6%,n:364       |
     |359    Complete BUSCOs (C)                        |
     |359    Complete and single-copy BUSCOs (S)        |
     |0    Complete and duplicated BUSCOs (D)           |
     |3    Fragmented BUSCOs (F)                        |
     |2    Missing BUSCOs (M)                           |
     |364    Total BUSCO groups searched                |
     ---------------------------------------------------
```

### Step 2: Run multiple polishing tools:

Align nanopore reads to unpolished assembly

Reads from Jodie: `/private/groups/russelllab/jodie/wWil_genome/Wwil_fastq.tar.gz`
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 -ax map-ont -t 14 \
    /private/groups/patenlab/mira/wWil_polishing/data/wWil_nanopore_assembly_decontaminated.fasta \
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
    /private/groups/patenlab/mira/wWil_polishing/data/wWil_nanopore_assembly_decontaminated.fasta \
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
    --genome /private/groups/patenlab/mira/wWil_polishing/data/wWil_nanopore_assembly_decontaminated.fasta \
    --bam /private/groups/patenlab/mira/wWil_polishing/data/illumina/correct_bam/Wwil.ilm.mm2.wWil_nanopore_assembly_decontaminated.srt.maxDiv.04.bam \
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

# align illumina reads to racon+homopolish polished assembly
docker run -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa mem -t32 \
    /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/medaka/consensus.fasta \
    /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only_R1.fastq.gz \
    /private/groups/patenlab/mira/wWil_polishing/data/illumina/JW18wWil0703A_wWil-only_R2.fastq.gz | samtools view -b -h > /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.Racon_medaka_polished.bam

samtools sort /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.Racon_medaka_polished.bam > /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.Racon_medaka_polished.srt.bam
samtools index /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.Racon_medaka_polished.srt.bam

# run pilon
java -Xmx16G -jar /private/home/mmastora/progs/pilon-1.24.jar \
    --genome /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/medaka/consensus.fasta \
    --bam /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil_ilm.bwa.Racon_medaka_polished.srt.bam \
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


### Step 3: Compare QC metrics across different polishing tools

#### QV with merqury

Racon only: QV improves by .2
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Racon_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.de.04.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.fasta wWil_Racon_merqury_k15

wWil_Nanopore_assembly.Racon    61574   1264923 24.7869 0.00332132
```

Racon + medaka: another .6 QV improvement
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Racon_Medaka_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.de.04.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/medaka/consensus.fasta wWil_Racon_Medaka_merqury_k15

```
Racon + homopolish:
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Racon_homopolish_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.de.04.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.homopolish.fasta/wWil_Nanopore_assembly_homopolished.fasta wWil_Racon_homopolish_merqury_k15

```

Pilon only
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Pilon_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.de.04.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.pilon.polished.fasta wWil_Pilon_merqury_k15

```

Racon+Medaka+Pilon
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Racon_Medaka_Pilon_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.de.04.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.racon.medaka.pilon.polished.fasta wWil_racon_medaka_Pilon_merqury_k15

```

Racon+homopolish+pilon
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/Racon_homopolish_Pilon_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.de.04.k15.meryl /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.racon.homopolish.pilon.polished.fasta wWil_racon_homopolish_Pilon_merqury_k15

```
#### Busco

Racon only
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco_rickettsiales_odb10/Racon:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.fasta \
    -f \
    -c 16 \
    -l rickettsiales_odb10

    ---------------------------------------------------
        |Results from dataset rickettsiales_odb10          |
        ---------------------------------------------------
        |C:80.2%[S:80.2%,D:0.0%],F:13.7%,M:6.1%,n:364      |
        |292    Complete BUSCOs (C)                        |
        |292    Complete and single-copy BUSCOs (S)        |
        |0    Complete and duplicated BUSCOs (D)           |
        |50    Fragmented BUSCOs (F)                       |
        |22    Missing BUSCOs (M)                          |
        |364    Total BUSCO groups searched                |
        ---------------------------------------------------
```

Racon + Medaka
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco_rickettsiales_odb10/Racon_Medaka:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/medaka/consensus.fasta \
    -c 16 -f \
    -l rickettsiales_odb10

    ---------------------------------------------------
       |Results from dataset rickettsiales_odb10          |
       ---------------------------------------------------
       |C:94.0%[S:94.0%,D:0.0%],F:3.3%,M:2.7%,n:364       |
       |342    Complete BUSCOs (C)                        |
       |342    Complete and single-copy BUSCOs (S)        |
       |0    Complete and duplicated BUSCOs (D)           |
       |12    Fragmented BUSCOs (F)                       |
       |10    Missing BUSCOs (M)                          |
       |364    Total BUSCO groups searched                |
       ---------------------------------------------------
```

Racon + homopolish
```

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco_rickettsiales_odb10/Racon_homopolish:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.Racon.homopolish.fasta/wWil_Nanopore_assembly_homopolished.fasta \
    -c 16 -f \
    -l rickettsiales_odb10

    ---------------------------------------------------
    |Results from dataset rickettsiales_odb10          |
    ---------------------------------------------------
    |C:98.9%[S:98.9%,D:0.0%],F:0.3%,M:0.8%,n:364       |
    |360    Complete BUSCOs (C)                        |
    |360    Complete and single-copy BUSCOs (S)        |
    |0    Complete and duplicated BUSCOs (D)           |
    |1    Fragmented BUSCOs (F)                        |
    |3    Missing BUSCOs (M)                           |
    |364    Total BUSCO groups searched                |
    ---------------------------------------------------
```

Pilon only
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco_rickettsiales_odb10/Pilon:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.pilon.polished.fasta \
    -c 16 -f \
    -l rickettsiales_odb10

---------------------------------------------------
|Results from dataset rickettsiales_odb10          |
---------------------------------------------------
|C:99.5%[S:99.5%,D:0.0%],F:0.0%,M:0.5%,n:364       |
|362    Complete BUSCOs (C)                        |
|362    Complete and single-copy BUSCOs (S)        |
|0    Complete and duplicated BUSCOs (D)           |
|0    Fragmented BUSCOs (F)                        |
|2    Missing BUSCOs (M)                           |
|364    Total BUSCO groups searched                |
---------------------------------------------------
```

Racon+Medaka+Pilon
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco_rickettsiales_odb10/Racon_Medaka_Pilon:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.racon.medaka.pilon.polished.fasta \
    -c 16 -f \
    -l rickettsiales_odb10
---------------------------------------------------
 |Results from dataset rickettsiales_odb10          |
 ---------------------------------------------------
 |C:99.5%[S:99.5%,D:0.0%],F:0.0%,M:0.5%,n:364       |
 |362    Complete BUSCOs (C)                        |
 |362    Complete and single-copy BUSCOs (S)        |
 |0    Complete and duplicated BUSCOs (D)           |
 |0    Fragmented BUSCOs (F)                        |
 |2    Missing BUSCOs (M)                           |
 |364    Total BUSCO groups searched                |
 ---------------------------------------------------
```

Racon+homopolish+Pilon
```

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco_rickettsiales_odb10/Racon_homopolish_Pilon:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/polished_assemblies/wWil_Nanopore_assembly.racon.homopolish.pilon.polished.fasta \
    -c 16 -f \
    -l rickettsiales_odb10

---------------------------------------------------
|Results from dataset rickettsiales_odb10          |
---------------------------------------------------
|C:99.5%[S:99.5%,D:0.0%],F:0.0%,M:0.5%,n:364       |
|362    Complete BUSCOs (C)                        |
|362    Complete and single-copy BUSCOs (S)        |
|0    Complete and duplicated BUSCOs (D)           |
|0    Fragmented BUSCOs (F)                        |
|2    Missing BUSCOs (M)                           |
|364    Total BUSCO groups searched                |
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

### Debugging busco discrepancy with Anne's assembly

full flye output from https://drive.google.com/drive/folders/1Zndl7JpjZ3AuX-rAJZcuI4gZjG-4XlfU
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco/assembly.fastas:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/data/assembly.fasta \
    -c 16 \
    -l rickettsiales_odb10

```

Anne's assembly
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/wWil_polishing/busco_rickettsiales_odb10/anne_assembly.fasta:/busco_wd \
    ezlabgva/busco:v5.7.0_cv1 busco \
    --mode genome \
    -i /private/groups/patenlab/mira/wWil_polishing/data/anne_assembly.fasta \
    -c 16 \
    -l rickettsiales_odb10
```

```
---------------------------------------------------
|Results from dataset rickettsiales_odb10          |
---------------------------------------------------
|C:98.6%[S:98.6%,D:0.0%],F:0.8%,M:0.6%,n:364       |
|359    Complete BUSCOs (C)                        |
|359    Complete and single-copy BUSCOs (S)        |
|0    Complete and duplicated BUSCOs (D)           |
|3    Fragmented BUSCOs (F)                        |
|2    Missing BUSCOs (M)                           |
|364    Total BUSCO groups searched                |
---------------------------------------------------
```

Align Anne's assembly to unpolished assembly
```
docker run --rm -u `id -u`:`id -g` -v /private/groups:/private/groups mobinasri/long_read_aligner:v0.3.3 minimap2 -a --eqx -x asm5 -t128 -c --cs /private/groups/patenlab/mira/wWil_polishing/data/wWil_Nanopore_assembly.fasta /private/groups/patenlab/mira/wWil_polishing/data/anne_assembly.fasta -o /private/groups/patenlab/mira/wWil_polishing/data/anne_asm_to_wWil_Nanopore_assembly.sam
```

Anne's assembly merqury
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups -v /private/groups/patenlab/mira/wWil_polishing/merqury/annes_asm_k15:/data juklucas/hpp_merqury:latest merqury.sh /private/groups/patenlab/mira/wWil_polishing/data/illumina/wWil.ilm_40x.k15.meryl /private/groups/patenlab/mira/wWil_polishing/data/anne_assembly.fasta annes_asm_merqury_k15
```
