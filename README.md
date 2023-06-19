# Contamination of Genomes using BlobToolKit
## Mapping network drive for the database server
In order to access the server with the databases stored on it and to move the finished files to your server you need to map the network drive.
This can be done by clicking This PC > Computer > Map network drive then chose the drive you want. The folder needed is \\biostore01\bioblob. Then click finish.

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
