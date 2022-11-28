# Introduction
# Phylogenetic analysis
## GTDB-tk (0.3.2)
Create the GTDB-Tk environment and install
```
conda create -n gtdbtk -c conda-forge -c bioconda gtdbtk
```
download database:
```
download-db.sh
```
Verifying your installation:
```
gtdbtk check_install
```
output be like:
```
[2020-11-04 09:35:16] INFO: GTDB-Tk v1.4.0
[2020-11-04 09:35:16] INFO: gtdbtk check_install
[2020-11-04 09:35:16] INFO: Using GTDB-Tk reference data version r95: /release95
[2020-11-04 09:35:16] INFO: Running install verification
[2020-11-04 09:35:16] INFO: Checking that all third-party software are on the system path:
[2020-11-04 09:35:16] INFO:          |-- FastTree         OK
[2020-11-04 09:35:16] INFO:          |-- FastTreeMP       OK
[2020-11-04 09:35:16] INFO:          |-- fastANI          OK
[2020-11-04 09:35:16] INFO:          |-- guppy            OK
[2020-11-04 09:35:16] INFO:          |-- hmmalign         OK
[2020-11-04 09:35:16] INFO:          |-- hmmsearch        OK
[2020-11-04 09:35:16] INFO:          |-- mash             OK
[2020-11-04 09:35:16] INFO:          |-- pplacer          OK
[2020-11-04 09:35:16] INFO:          |-- prodigal         OK
[2020-11-04 09:35:16] INFO: Checking /release95
[2020-11-04 09:35:16] INFO:          |-- pplacer          OK
[2020-11-04 09:35:16] INFO:          |-- masks            OK
[2020-11-04 09:35:17] INFO:          |-- markers          OK
[2020-11-04 09:35:17] INFO:          |-- radii            OK
[2020-11-04 09:35:20] INFO:          |-- msa              OK
[2020-11-04 09:35:20] INFO:          |-- metadata         OK
[2020-11-04 09:35:20] INFO:          |-- taxonomy         OK
[2020-11-04 09:47:36] INFO:          |-- fastani          OK
[2020-11-04 09:47:36] INFO:          |-- mrca_red         OK
[2020-11-04 09:47:36] INFO: Done.
```
check test file:
```
gtdbtk test --out_dir /tmp/test --cpus 8
```
output be like:
```
[2020-04-13 09:50:58] INFO: GTDB-Tk v1.1.0
[2020-04-13 09:50:58] INFO: gtdbtk test --out_dir /tmp/test --cpus 3
[2020-04-13 09:50:58] INFO: Using GTDB-Tk reference data version r89: /release89
[2020-04-13 09:50:58] INFO: Command: gtdbtk classify_wf --genome_dir /tmp/test/genomes --out_dir /tmp/test/output --cpus 3
[2020-04-13 09:52:35] INFO: Test has successfully finished.
```
Output files: [prefix].warnings.log, output/, test_execution.log
What’s in output/: gtdbtk.ar122.classify.tree, gtdbtk.ar122.summary.tsv, gtdbtk.ar122.markers_summary.tsv, gtdbtk.ar122.msa.fasta, gtdbtk.ar122.user_msa.fasta

Species annotation:
```
gtdbtk classify_wf --genome_dir <my_genomes> --out_dir <output_dir> --cpus 64 -x fna --force
```
Upload \*.unrooted.tree to itol (https://itol.embl.de/) for visualization.
## Gene trees of sadABC
Using blastp to align coding sequences (CDS) of each bin against amino acid sequences of enzymes SadA, SadB and SadC.
```
makeblastdb -in sad.faa -dbtype prot -out sad
blastp -query bin.1.faa -out SDB1_sad.txt -db sad -outfmt 6 -evalue 1e-5 -num_threads 8
```
CDS that identity > 30% were chosen as functional genes. Phylogenetic and molecular evolutionary analyses of sadABC genes were conducted using MEGA X. Tree files (sadA.nwk, sadB.nwk and sadC.nwk) were put in files/ and uploaded to iTOL for visualization.

# RNA-seq
## Install required software
It is strongly recommended to create a new conda environment for metatranscriptome before installation. The software required, the version used and installation commands are listed:  
fastp, version 0.20.1  
Bowtie2, version 2.4.4  
samtools, version 1.14  
featureCounts, version 2.0.1  
```
conda create -n RNA-seq
conda activate RNA-seq
conda install -c bioconda fastp
conda install -c bioconda bowtie2
conda install -c bioconda samtools
```
enter the direction to install featureCounts
```
wget https://jaist.dl.sourceforge.net/project/subread/subread-2.0.1/subread-2.0.1-Linux-x86_64.tar.gz
tar -zxvf subread-2.0.1-Linux-x86_64.tar.gz
```
# count reads
Raw data was deposited in NCBI under the accession number of PRJNA799876. Four triplicate samples were collected during SMX degradation at 0, 4, 8 and 12 h.
Quality control:
```
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
```
Build sam database
```
bowtie2-build ASSEMBLY/megahit/final.contigs.fa smx
bowtie2 -x smx -1 example_1.1.fastq -2 example_1.2.fastq -S example.sam -p 32
```
sam to bam:
```
samtools view -bS smx.sam > smx.bam
samtools sort -n smx.bam -o smx.sorted.bam
samtools view -h smx.sorted.bam|less -N
```
Count mapped reads using the sorted bam file
```
featureCounts -p -t gene -g gene_id -a total.gtf -o counts.txt smx.sorted.bam
```
Sort 12 counts files into 1 file manually, see detailed in files/total_counts.txt.
Select degradation genes shown in Figure 5 and calculate their average counts, see detailed in files/degradation_genes_expression.txt.

## Heatmap (based on R)
Install and library the package:
```
install.packages("pheatmap")
library(pheatmap)
```
read data and transform to log2:
```
data <- read.table ("degradation_genes_expression.txt", header = T, sep ="\t", row.names = 1)
datalg = log2(data+1)
```
draw heatmap:
```
pheatmap (datalg, cluster_cols = F, legend = T, angle_col = 0, clustering_method = "complete", annotation_legend = F, >
show_rownames = T, cellwidth = 25, cellheight = 25, border_color = "black")
```
# Plasmid
## Plasflow (1.1.0)
Install:
```
conda create -n plasflow
conda install -c bioconda plasflow
```
filter contigs less than 1000 bp:
```
filter_sequences_by_length.pl -input SDB1.fa -output SDB1filter.fasta -thresh 1000
```
run plasflow:
```
PlasFlow.py --input SDB1filter.fasta --output SDB1.tsv --threshold 0.7
```
The results contain four files: SDB1.tsv, SDB1.tsv_chromosomes.fasta, SDB1.tsv_plasmids.fasta, SDB1.tsv_unclassified.fasta

## Plascad
Install:
```
conda create -n plascad
conda install -c pianpianyouche plascad
```
take fasta file of plasmid sequences as input:
```
Plascad -i SDB1.tsv_plasmids.fasta
```
Output
\*\_classification_sum.txt: name, plasmid type, ARGs (type_subtype)
\*\_conj_plasimids_loc_sum.txt: name, marker genes, c-value & e-value (c-value for hmmsearch, e-value for blastp), genetic location, strand
A plasmid visualization component using AngularJS is integrated into our pipeline, all the plasmid maps are in HTML formats. In order to view the map locally, you need to download the js folder in addtion to the HTML files.

## SPAdes
Require a 64-bit Linux system and Python (supported versions are Python 2.7, and Python3: 3.2 and higher) on it. Go to the directory in which you wish SPAdes to be installed and run:

```
wget http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Linux.tar.gz
tar -xzf SPAdes-3.15.5-Linux.tar.gz
cd SPAdes-3.15.5-Linux/bin/
```
In case of successful installation the following files will be placed in the bin directory:
```
spades.py (main executable script)
metaspades.py (main executable script for metaSPAdes)
plasmidspades.py (main executable script for plasmidSPAdes)
metaplasmidspades.py (main executable script for metaplasmidSPAdes)
metaviralspades.py (main executable script for metaviralSPAdes)
rnaspades.py (main executable script for rnaSPAdes)
truspades.py (main executable script for truSPAdes, DEPRECATED)
rnaviralspades.py (main executable script for rnaviralSPAdes)
coronaspades.py (wrapper script for coronaSPAdes mode)
spades-core (assembly module)
spades-gbuilder (standalone graph builder application)
spades-gmapper (standalone long read to graph aligner)
spades-kmercount (standalone k-mer counting application)
spades-hammer (read error correcting module for Illumina reads)
spades-ionhammer (read error correcting module for IonTorrent reads)
spades-bwa (BWA alignment module which is required for mismatch correction)
spades-corrector-core (mismatch correction module)
spades-truseq-scfcorrection (executable used in truSPAdes pipeline)
```
adding SPAdes installation directory to the PATH variable.
```
export PATH=/data/home/a204a/SPAdes/SPAdes-3.15.5-Linux/bin:$PATH
```
Verifying your installation
```
spades.py --test
```
If the installation is successful, you will find the following information at the end of the log:
```
===== Assembling finished. Used k-mer sizes: 21, 33, 55

 * Corrected reads are in spades_test/corrected/
 * Assembled contigs are in spades_test/contigs.fasta
 * Assembled scaffolds are in spades_test/scaffolds.fasta
 * Assembly graph is in spades_test/assembly_graph.fastg
 * Assembly graph in GFA format is in spades_test/assembly_graph_with_scaffolds.gfa
 * Paths in the assembly graph corresponding to the contigs are in spades_test/contigs.paths
 * Paths in the assembly graph corresponding to the scaffolds are in spades_test/scaffolds.paths

======= SPAdes pipeline finished.

========= TEST PASSED CORRECTLY.

SPAdes log can be found here: spades_test/spades.log

Thank you for using SPAdes!
```
search for plasmids from metagenomic clean reads:
```
spades.py --metaplasmid -o /path/to/output -1 /CLEAN_READS/ALL_READS_1.fastq -2 /CLEAN_READS/ALL_READS_2.fastq --threads 32
```
Output:
Contigs and scaffolds format (scaffolds.fasta is recomended)
Contigs/scaffolds names in SPAdes output FASTA files have the following format: >NODE_3_length_237403_cov_243.207 Here 3 is the number of the contig/scaffold, 237403 is the sequence length in nucleotides and 243.207 is the k-mer coverage for the last (largest) k value used. Note that the k-mer coverage is always lower than the read (per-base) coverage.
In order to distinguish contigs with putative plasmids detected at different cutoff levels, contig name in FASTA file was added with cutoff value used for this particular contig (in format \_cutoff_N). For metaplasmid mode only circular putative plasmids was output.

# Calculate coverage of contigs
```
conda activate RNA-seq
cd /path/to/assemble/contigs
bowtie2-build --thread 40 final.contigs.fa all_contig_build
bowtie2 -x all_contig_build -1 ALL_READS_1.fastq -2 ALL_READS_2.fastq -S all_contig.sam --threads 1 
samtools view -bS all_contig.sam > all_contig.bam
samtools sort  all_contig.bam -o all_contig_sorted.bam
samtools index all_contig_sorted.bam
```
put three files: fasta, sorted.bam, sorted.bam.bai into one folder checkM_all
```
checkm coverage checkM_all/ all_coverage.out checkM_all/all_contig_sorted.bam -x fasta -m 20
```
