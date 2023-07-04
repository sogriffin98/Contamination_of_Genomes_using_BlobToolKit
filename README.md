# Contamination of Genomes using BlobToolKit
I used BlobToolKit to look at contamination in the genome assembly of *Meloidogyne fallax*. I used the a Future Learn tutorial and edited the commands to make it work for my species (https://github.com/blobtoolkit/tutorials/tree/main/futurelearn). I have written instructions using my study species as an example for others to try.

## Mapping network drive for the database server
In order to access the server with the databases stored on it and to move the finished files to your server you need to map the network drive.
This can be done by clicking This PC > Computer > Map network drive then chose the drive you want. The folder needed is \\biostore01\bioblob. Then click finish.

## Programs needed for using BlobToolKit and how to install them:
### 1. BlobToolKit
BlobToolKit is a modular command-line solution for visualisation, quality control and taxonomic partitioning of genome datasets (https://blobtoolkit.genomehubs.org/install/). It can be installed using the following commands:
```
conda create -y -n btk python=3.9
conda activate btk
pip install blobtoolkit[full]
conda deactivate
``` 
### 2. Diamond 
Diamond is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data (https://github.com/bbuchfink/diamond). It can be installed using the following commands:
```
conda create -n diamond
conda activate diamond
conda install -c bioconda diamond
conda deactivate
``` 
### 3. Blast
Blast is a tool used for finding regions of local similarity between sequences (https://github.com/ncbi/blast_plus_docs). It can be installed using the following commands:
```
conda create -n blast
conda activate blast
conda install -c bioconda blast
conda deactivate
``` 
### 4. MiniMap2
MiniMap2 is a sequence alignment tool (https://github.com/lh3/minimap2). It can be installed using the following commands:
```
conda create -n minimap2
conda activate minimap2
conda install -c bioconda minimap2
conda deactivate
```
### 5. Samtools
Samtools is a tool used for handling SAM, BAM and CRAM files (https://github.com/samtools/samtools). It can be installed using the following commands:
```
conda create -n samtools
conda activate samtools
 conda install -c bioconda samtools
conda deactivate
```
### 6. BUSCO
BUSCO is a tool used for quality control, gene prediction and phylogenomics (https://busco.ezlab.org/). It can be installed using the following commands:
```
conda create -n busco
conda activate busco
conda install -c bioconda busco
conda deactivate
```
## Downloading Databases
In order to blast search your genome assembly you need to download the diamond, blast and uniprot databases.
### 1. First make a directory for all of your databases:
```
mkdir btk
```
### 2. Fetch the NCBI Taxdump:
```
mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;
```
### 3. Fetch the nt database:
```
mkdir -p nt
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz" -P nt/ && \
        for file in nt/*.tar.gz; \
            do tar xf $file -C nt && rm $file; \
        done
```
### 4. Fetch and format the UniProt database:
```
mkdir -p uniprot
wget -q -O uniprot/reference_proteomes.tar.gz \
 ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
     -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
     awk '/tar.gz/ {print $9}')
cd uniprot
tar xf reference_proteomes.tar.gz

touch reference_proteomes.fasta.gz
find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz

echo -e "accession\taccession.version\ttaxid\tgi" > reference_proteomes.taxid_map
zcat */*/*.idmapping.gz | grep "NCBI_TaxID" | awk '{print $1 "\t" $1 "\t" $3 "\t" 0}' >> reference_proteomes.taxid_map

diamond makedb -p 16 --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map --taxonnodes ../taxdump/nodes.dmp -d reference_proteomes.dmnd
cd -
```

Once you have all of the programs and databases installed follow the instructions below:
## Accessing the Databases Needed for BlobToolKIt
1. Change directories to where the databases are stored. This can be done using the following command:
```
cd btk
```
2. Create a new BlobToolKit directory using mkdir and make sure to include your name so you know which one is yours e.g. ```blobdir_sgriffin```:
```
mkdir blobdir_<your_name>
```
Do all of the following steps below in the btk directory then move to your own blobdir directory once finished. When you have finished everything move your data out to your own directory.

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
Using your assembly.fasta or contigs.fasta file from your genome assembly of your study species run the following commands:
```
conda activate diamond
diamond blastx --query contigs.fasta --db ./uniprot/reference_proteomes.dmnd --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 30 > diamond.blastx.out
conda deactivate
```
Be sure to change the name of the contig.fasta file in the command to match your file name.
#### 2.2. Create a blastn database
Using your assembly.fasta or contigs.fasta file from your genome assembly of your study species run the following command:
```
conda activate blast
blastn -db ./nt/nt -query contigs.fasta -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 30 -out blastn.out
conda deactivate
```
Be sure to change the name of the contig.fasta file in the command to match your file name.
### 3. Create a coverage file:
Using minimap2 create a coverage file using your assembly.fasta/contigs.fasta file from your genome assembly and your two Illumina reads. Be sure to change the file names in the line below to match your file names. An example is included for *Meloidogyne fallax*:
```
conda activate minimap2
minimap2 -ax sr -t 30 contigs.fasta MF1_S12_R1_001_val_1.fq.gz MF1_S12_R2_001_val_2.fq.gz
conda deactivate
conda activate samtools
samtools sort -@30 -O BAM -o coverage.bam -
conda deactivate
```
### Create a BUSCO Summary File:
Use the correct lineage for your species (https://busco.ezlab.org/list_of_lineages.html). An example for *Meloidogyne fallax* is below:
```
conda activate busco
busco -i contigs.fasta -l nematoda_odb10 -o Meloidogyne -m genome --cpu 30
conda deactivate
```
## Add Files to BlobToolKit Dataset
### 1. Add assembly and a metadata files
This example uses the M_fallax.yaml file. Ensure to change it to match the name of your .yaml file. It also uses a contigs.fasta file which might be different for you. Your genome assembly file might be a assembly.fasta file instead or it might have a different name.
```
blobtools create \
--fasta contigs.fasta \
--meta M_fallax.yaml \
--taxid 71801 \
--taxdump ./taxdump/ \
./btk
```
### 2. Add hits
Again as before these files are the ones created for *Meloidogyne fallax*, make sure you have the correct names of your diamond and blast output files.
```
blobtools add \
--hits blastn.out \
--hits diamond.blastx.out \
--taxrule bestsumorder \
--taxdump ./taxdump/ \
./btk
```
### 3. Add Coverage
Again as before this coverage file was created for *Meloidogyne fallax*, make sure you have the correct name for your coverage.bam file.
```
blobtools add \
--cov coverage.bam \
--threads 30 \
./btk
```
### 4. Add BUSCO
Again as before this BUSCO file was created for *Meloidogyne fallax*, make sure you have the correct name for your BUSCO file.
```
blobtools add \
--busco ./Meloidogyne/run_nematoda_odb10/full_table.tsv \
./btk
```
## View the dataset on BlobToolKit Viewer
### 1. Run the following command to initialize the viewer
```
blobtools view \
--remote \
./btk
```
### 2. SHH Instructions
1. Open PuTTY Configuration Panel
   
![MicrosoftTeams-image (1)](https://github.com/sogriffin98/BlobToolKitMfallax/assets/117275745/8daa0f45-b714-408b-b922-6ff241885249)

2. In Host Name type the name of the server you are using. In this example the server name is 'bio-prod02'
3. In the Saved Sessions type a name that makes sense for you. In this example I have used 'bioprod02 with blob tunnel'.
4. Click Save
5. Click Connection > SSH > Tunnels on the left side

![MicrosoftTeams-image (2)](https://github.com/sogriffin98/BlobToolKitMfallax/assets/117275745/33502f2e-47ee-4049-ab55-9d7b6247f8bd)

5. In the Add new forwarded port section put the numbers that were outputted on the terminal window with step 1. An example is below:

   5.1. Source port: 8001, Destination: 127.0.0.1:8001
   
   5.2. Source port: 8002, Destination: 127.0.0.1:8002

### 3. Open your browser (i.e. Google Chrome, Internet Explorer, Firefox) and go to the URL indicated in the terminal window (step 1)

## Filtering blobtoolkit dataset
Datasets can be filtered based on the values in any variable (e.g. GC proportion and length) or category field (e.g. assigned phylum), or by using a list of identifiers (sequence IDs). Filters may be applied to a complete dataset to allow for use of a reduced dataset without repeating analyses or applied to assembly FASTA and read FASTQ files to allow for reassembly and reanalysis. Filter parameters are all shared between BlobTools and the BlobToolKit Viewer, allowing interactive sessions to be reproduced on the command line.

Configuration options for blobtools filter can be considered in three groups:
### 1. Setting filter parameters:
```--param``` – String of type param=value to specify individual filter parameters.
```--query-string``` – List of param=value pairs (separated by &) pairs from url query string.
```--json``` – JSON format list file as generated by BlobtoolKit Viewer.
```--list``` – Text file containing a space or newline separated list of identifiers.
```--invert```– Flag to invert the filter (exclude matching records).

### 2. Specifying files to filter:
```--output``` – Path to directory to save a filtered copy of the BlobDir.
```--fasta``` – FASTA format assembly file to be filtered.
```--fastq``` – FASTQ format read file to be filtered (requires --cov).
```--cov``` – BAM/SAM/CRAM read alignment file.
```--text``` – Generic text file to be filtered.
```--text-delimiter``` – Text file delimiter. [Default: whitespace]
```--text-id-column``` – Index of column containing identifiers (1-based). [Default: 1]
```--text-header``` – Flag to indicate first row of text file contains field names. [Default: False]
```--suffix``` -STRING String to be added to filtered filename. [Default: filtered]

### 3. Generating summary data:
```--summary``` – Filename for a JSON-format summary of the filtered dataset.
```--summary-rank``` – Taxonomic level for summary. [Default: phylum]
```--taxrule``` – Taxrule used when processing hits. [Default: bestsumorder]
```--table``` – Filename for a tabular output of filtered dataset.
```--table-fields``` – Comma separated list of field IDs to include in the table output. Use ‘plot’ to include all plot axes. [Default: plot]

## Other Filtering Options
### 1. Filter an assembly FASTA file to exclude sequences shorter than 1000 bp, to also filter the assembly to remove sequences with no hit in the reference databases and print a set of summary statistics for the filtered assembly to the terminal:
```
blobtools filter \
     --param length--Min=1000 \
     --param bestsumorder_phylum--Keys=no-hit \
     --fasta assembly.fasta \
     --summary STDOUT \
     ./btk/

