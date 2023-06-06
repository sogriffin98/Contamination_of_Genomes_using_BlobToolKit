# Contamination of Genomes using BlobToolKit
## Mapping network drive for the database server
In order to access the server with the databases stored on it and to move the finished files to your server you need to map the network drive.
This can be done by clicking This PC > Computer > Map network drive then chose the drive you want. The folder needed is \\biostore01\bioblob. Then click finish.

## Accessing the Databases
First change directories to where the databases are stored. This can be done using the following command:
```
cd /nfs/bioblob/btk
```
Do all of the following steps in the btk directory then move to your own directory once finished

## Example of How to Create a BlobToolKit dataset
1. Create a new BlobToolKit directory using mkdir:
```
mkdir blobdir
```
2. Create a metadata file by opening a text editor and create a file called B_cinera.yaml. Then write the following in your yaml file:
```
assembly:
  alias: B_cinera_112
  record_type: contig
taxon:
  name: Botrytis cinerea
 ```
Make sure to change the content to match your species information
2. Create hit files:
