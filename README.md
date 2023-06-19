# Contamination of Genomes using BlobToolKit
## Mapping network drive for the database server
In order to access the server with the databases stored on it and to move the finished files to your server you need to map the network drive.
This can be done by clicking This PC > Computer > Map network drive then chose the drive you want. The folder needed is \\biostore01\bioblob. Then click finish.

## Programs needed for using BlobToolKit and how to install them:
### Diamond 
Diamond is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data (https://github.com/bbuchfink/diamond). It can be installed using the following commands:
```
conda create -n diamond
conda activate diamond
conda install -c bioconda diamond
```
### 
### Blast
### MiniMap2
MiniMap2 is a sequence alignment tool (https://github.com/lh3/minimap2). It can be installed using the following commands:
```
conda create -n minimap2
conda activate minimap2
conda install -c bioconda minimap2
```

### BUSCO
BUSCO is a tool used for quality control, gene prediction and phylogenomics (https://busco.ezlab.org/). It can be installed using the following commands:
```
conda create -n busco
conda activate busco
 conda install -c bioconda busco
```

## Accessing the Databases
1. First change directories to where the databases are stored. This can be done using the following command:
```
cd /nfs/bioblob/btk
```
2. Create a new BlobToolKit directory using mkdir and make sure to include your name so you know which one is yours e.g. ```blobdir_sgriffin```:
```
mkdir blobdir_<your_name>
```
Do all of the following steps below in the btk directory then move to your own blobdir directory once finished. When you have finished everything move your data out to your own server.

## Example of How to Create a BlobToolKit dataset
### 1. Create a metadata file:
Open a text editor, e.g. ```nano```, and create a file called ```<your_species_name.yaml```. For example ```M_fallax.yaml```. Write the corresponding information for your species in your yaml file, e.g.:
```
assembly:
  alias: M_fallax
  record_type: contig
taxon:
  name: Meloidogyne fallax
 ```
Make sure to change the content to match your species information.

### 2. Create database input files:
#### 2.1. Create a Diamond database
Using your assembly.fasta or contigs.fasta file from your genome assembly of your study species run the following command:
```
diamond blastx --query contigs.fasta --db ./uniprot/reference_proteomes.dmnd --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 30 > diamond.blastx.out
```
Be sure to change the name of the contig.fasta file in the command to match your file name.
#### 2.2. Create a blastn database
Using your assembly.fasta or contigs.fasta file from your genome assembly of your study species run the following command:
```
diamond blastx --query contigs.fasta --db ./uniprot/reference_proteomes.dmnd --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 30 > diamond.blastx.out
```
Be sure to change the name of the contig.fasta file in the command to match your file name.
### 3. Create a coverage file:
Using minimap2 create a coverage file using your assembly.fasta/contigs.fasta file from your genome assembly and your two Illumina reads. Be sure to change the file names in the line below to match your file names. An example is included for *Meloidogyne fallax*:
```
minimap2 -ax sr -t 30 contigs.fasta \
MF1_S12_R1_001_val_1.fq.gz MF1_S12_R2_001_val_2.fq.gz \
| samtools sort -@30 -O BAM -o coverage.bam -
```
### Create a BUSCO Summary File:
Use the correct lineage for your species. This can be found using the following link: https://busco.ezlab.org/list_of_lineages.html 
An example for *Meloidogyne fallax* is below:
```
busco -i contigs.fasta -l nematoda_odb10 -o Meloidogyne -m genome --cpu 30
```
## Add Files to BlobToolKit Dataset
