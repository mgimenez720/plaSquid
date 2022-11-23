# plaSquid: **Pla**smid **S**e**q**uences **Id**entification in metagenomic assemblies.  

### Description

- **plaSquid** is a Nextflow pipeline for plasmid detection and classification in genomic and metagenomic data. This pipeline accepts either genomic or metagenomic *assemblies* as input (.fasta). It uses two different approaches to detect plasmids sequences: alignment with minimap2 against a plasmidic database (minidist) and HMM dependent search of plasmid specific genes (repsearch).

- **plaSquid** also classifies plasmids into replicon types and MOB groups by comparing RIPs or Relaxases against custom HMMs. 

- **plaSquid** can extract plasmids RIP or MOB sequences in order to further analyze these proteins.   

- **plaSquid** summarises the information gathered by the two complementary approaches in a single output table and allows further analysis as it outputs plasmidic contigs in a single multifasta file (*"Result.fasta"*)

![Pipeline overview](./img/plaSquid_pipeline.png)

### Installation

    git clone https://github.com/mgimenez720/plaSquid/
    cd plaSquid/
    
You need nextflow installed in order to run plaSquid. Documentation is available [here](https://www.nextflow.io/docs/latest/getstarted.html)     
    
PlaSquid can be ran using docker or conda.  

If you want to generate a permanent docker image you can try:
    
    docker pull mgimenez720/plasquid:latest
    
If you want to generate a permanent conda environment you can try:

    conda env create -f environments/plaSquid.yml


### Dependencies

All dependencies are provided within the containers available. Manual installation is discouraged. 

[hmmer 3.3.1](http://hmmer.org/download.html),
[infernal 1.1.3](http://eddylab.org/infernal/),
[minimap2 2.17](https://github.com/lh3/minimap2),
[prodigal 2.6.3](https://github.com/hyattpd/Prodigal),
R packages:
dplyr 1.0.4,
tidyverse 1.3.0,
seqinr 4.2.5,
biostrings 2.58.0.


### Usage 


    nextflow run main.nf --contigs {testdata/test.fasta} --outdir {plaSquid_result} -profile docker

    arguments:

    --contigs       Path to input assemblies.
    --mmi           Path to Minimap2 indexed (.mmi) or fasta (.fasta/.fna) plsdb database.
    --outdir        Path to output directory where results are written.
    --help          Print help message and exit

    subworkflows:

    --minidist      Run mapping of contigs against plsdb database.
    --repsearch     Run search and classification of RIP and MOB (Rel) genes.
    --ripextract    Extract replication initiator proteins sequences.
    --mobextract    Extract relaxases sequences.
    
    profiles:

    -profile conda  Installs dependencies using a conda environment
    -profile docker Installs dependencies within a docker image
    -profile server runs using 15 cpus and 50 Gb
    -profile test   tests dependencies and normal functioning


    Authors:

    Matías Giménez
    Ignacio Ferrés
    Gregorio Iraola


    Microbial Genomics Laboratory
    Institut Pasteur Montevideo (Uruguay)

### Output

- **plaSquid** outputs consist of a fasta file "Result.fasta" with plasmids contigs detected along with a table "Result.tsv" with the following fields: 

>**"Contig"**: contig id for plaSquid   
**"name"**: contig name in the assembly file  
**"Sim-dist"**: S value obtained by Minidist workflow  
**"plsdb_match"**: plasmid matched at plsdb database  
**"Match_length"**: length of the plasmid matched at plsdb  
**"RIP_domain"**: RIP-domain found in that contig   
**"MOB_group"**: MOB group classification of relaxase found in that contig   
**"Rep_type"**: Rep-type classifiation of the contig detected  
**"Contig_length"**: size of the contig detected  

###Citation 

Available preprint 
https://doi.org/10.1101/2022.08.04.502827

### Note

This is a beta version, please report bugs or misfunctions detected.
