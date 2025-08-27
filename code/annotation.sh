#!/bin/bash

## Building a database

#Step 1: Configure a new genome in SnpEff's config file snpEff.config.
# In order to tell SnpEff that there is a new genome available, you must update SnpEff's configuration file snpEff.config.
# You must add a new genome entry to snpEff.config.

#vi snpEff.config
# Homo sapiens (hg19) custom
#"hg19_custom.genome : Homo_sapiens (custom hg19)
#"	hg19_custom.reference : data/hg19_custom/sequences.fa
#""	hg19_custom.M.codonTable : Vertebrate_Mitochondrial
#""	hg19_custom.MT.codonTable : Vertebrate_Mitochondrial
#""	hg19_custom.coordinates : GRCh37



#Step 2: Build using gene annotations and reference sequences
#Option 1: Building a database from GTF files (recommended for large genomes)
#Option 2: Building a database from GenBank files (recommended for small genomes)
#Option 3: Building a database from GFF files
#Option 4: Building a database from RefSeq table from UCSC

## Dowloading Genome Annotations and Reference Sequences
SNPEFF_DATA="/home/prospecmol/anaconda3/envs/DNASeq/share/snpeff-5.2-1/data/hg19_custom"
mkdir -p $SNPEFF_DATA
cd $SNPEFF_DATA

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip -f hg19.fa.gz
mv hg19.fa sequences.fa

wget -c ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip -f Homo_sapiens.GRCh37.75.gtf.gz
mv Homo_sapiens.GRCh37.75.gtf genes.gtf

wget -c ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.75.cds.all.fa.gz
gunzip -f Homo_sapiens.GRCh37.75.cds.all.fa.gz
mv Homo_sapiens.GRCh37.75.cds.all.fa cds.fa

wget -c ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.75.pep.all.fa.gz
gunzip -f Homo_sapiens.GRCh37.75.pep.all.fa.gz
mv Homo_sapiens.GRCh37.75.pep.all.fa protein.fa

#Step 3: Checking the database: SnpEff will check the database by comparing predicted protein sequences and CDS sequences with ones provided by the user.
cd /home/prospecmol/anaconda3/envs/DNASeq/share/snpeff-5.2-1
java -jar snpEff.jar build -gtf22 -v hg19_custom

echo " Dtabase hg19_custom completed"
#Checking CDS sequences
#Checking Protein sequences

